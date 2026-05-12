seq Analysis Pipeline (Plant Species)

|> [!TIP]
> Complete workflow for single-cell RNA-seq analysis of non-model plant species using cross-species cell type annotation with Arabidopsis thaliana reference and LLM-assisted functional annotation.

## Overview

This pipeline enables cell type annotation for non-model plant species by:
1. Quality control with automated organelle gene identification (ATC/ATM prefixes)
2. Standard Seurat clustering workflow across multiple resolutions
3. Cross-species homology mapping to Arabidopsis markers (Ath assay)


### Pipeline Summary

| Step | Description | Key Output |
|------|-------------|------------|
| **Step 1** | Data loading & Prefix-based QC | `Step1_Species_raw.rds` |
| **Step 1.5** | QC filtering | `Step1_5_Species_filtered.rds` |
| **Step 2** | Normalization & Multi-res Clustering | `Step2_Species_res[0.5/1/1.5/2].rds` |
| **Step 3** | Cross-species  | `Step5_Species_annotated.rds` |


---

## Step 1 - Step 2:  Code Template

> [!IMPORTANT]
> The following is a production-ready template that can be adapted by changing only the **configuration variables** at the top.

```r
#!/usr/bin/env Rscript
# Cross-Species scRNA-seq Annotation Workflow for Non-Model Plants
# ============================================================================

# ============================================================================
# CONFIGURATION - Modify these for your project
# ============================================================================
workdir <- "/path/to/your/project"
species_name <- "YourSpecies"  # Used for output file naming

# Input paths
matrix_dir <- file.path(workdir, "step4_cellranger_count/YOUR_SAMPLE/outs/filtered_feature_bc_matrix")
gene_mapping_file <- file.path(workdir, "step5_homolog_mapping/Species_Ath_best_hit.txt")
pcmdb_file <- file.path(workdir, "PCMDB_Ath.csv")  # Plant Cell Marker Database
func_anno_file <- file.path(workdir, "gene_annotation.txt")  # Functional annotations

# Output directory (all results go here)
output_dir <- file.path(workdir, "step5_clustering_annotation")

# ============================================================================
# Environment Setup
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(RColorBrewer)
  library(Matrix)
})

setwd(workdir)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

cat("=", rep("=", 70), "\n", sep="")
cat("SCRNA-seq CROSS-SPECIES ANNOTATION WORKFLOW\n")
cat("Species:", species_name, "\n")
cat("Output directory:", output_dir, "\n")
cat("=", rep("=", 70), "\n\n", sep="")

# ============================================================================
# Step 1: Data Loading & QC Metrics
# ============================================================================
cat("STEP 1: Loading 10X Data & Calculating QC Metrics\n")
cat("-", rep("-", 70), "\n", sep="")

data <- Read10X(data.dir = matrix_dir)
seurat_obj <- CreateSeuratObject(counts = data, project = species_name,
                                  min.cells = 3, min.features = 200)
cat("Initial Seurat object:\n")
cat("  - Cells:", ncol(seurat_obj), "\n")
cat("  - Genes:", nrow(seurat_obj), "\n\n")

# Load gene mapping for organelle identification
gene_map <- read.table(gene_mapping_file, header = FALSE, sep = "\t", 
                       stringsAsFactors = FALSE)
colnames(gene_map) <- c("TargetGene", "AthGene")

# Identify organelle genes by Arabidopsis ID prefixes
# ATC = Chloroplast, ATM = Mitochondrion
target_chloroplast <- unique(gene_map$TargetGene[grepl("^ATC", gene_map$AthGene)])
target_mito <- unique(gene_map$TargetGene[grepl("^ATM", gene_map$AthGene)])
target_chloroplast <- target_chloroplast[target_chloroplast %in% rownames(seurat_obj)]
target_mito <- target_mito[target_mito %in% rownames(seurat_obj)]

# Calculate QC metrics
seurat_obj[["percent.chloro"]] <- PercentageFeatureSet(seurat_obj, features = target_chloroplast)
seurat_obj[["percent.mito"]] <- PercentageFeatureSet(seurat_obj, features = target_mito)

# Generate QC violin plot
pdf(file.path(output_dir, "Step1_QC_ViolinPlot.pdf"), width = 12, height = 6)
print(VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", 
                                        "percent.chloro", "percent.mito"), ncol = 4, pt.size = 0.1))
dev.off()

saveRDS(seurat_obj, file.path(output_dir, paste0("Step1_", species_name, "_raw.rds")))
cat("  Saved: Step1_", species_name, "_raw.rds\n\n", sep="")

# ============================================================================
# Step 1.5: QC Filtering
# ============================================================================
cat("STEP 1.5: QC Filtering\n")
cat("-", rep("-", 70), "\n", sep="")

cells_before <- ncol(seurat_obj)
seurat_filtered <- subset(seurat_obj,
                          subset = nFeature_RNA > 300 & nFeature_RNA < 6000 &
                                   percent.mito < 5 & percent.chloro < 5)
cells_after <- ncol(seurat_filtered)
retention <- round(cells_after / cells_before * 100, 1)

cat("Filtering results:\n")
cat("  - Before:", cells_before, "cells\n")
cat("  - After:", cells_after, "cells\n")
cat("  - Retention rate:", retention, "%\n")

saveRDS(seurat_filtered, file.path(output_dir, paste0("Step1_5_", species_name, "_filtered.rds")))
cat("  Saved: Step1_5_", species_name, "_filtered.rds\n\n", sep="")

# ============================================================================
# Step 2: Normalization & Multi-Resolution Clustering
# ============================================================================
cat("STEP 2: Normalization & Multi-Resolution Clustering\n")
cat("-", rep("-", 70), "\n", sep="")

seurat_obj <- NormalizeData(seurat_filtered, normalization.method = "LogNormalize", scale.factor = 1e4)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = 'vst', nfeatures = 3000)
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("percent.chloro", "percent.mito"))
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:25, verbose = FALSE)

resolutions <- c(0.5, 1.0, 1.5, 2.0)

for (res in resolutions) {
  cat(sprintf("Processing resolution: %.1f\n", res))
  seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:25, verbose = FALSE)
  
  Idents(seurat_obj) <- paste0("RNA_snn_res.", res)
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, 
                            logfc.threshold = 0.58, verbose = FALSE)
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  
  write.csv(markers, file.path(output_dir, sprintf("Step2_All_Markers_res%.1f.csv", res)), row.names = FALSE)
  
  pdf(file.path(output_dir, sprintf("Step2_UMAP_res%.1f.pdf", res)), width = 8, height = 6)
  print(DimPlot(seurat_obj, reduction = "umap", label = TRUE) + ggtitle(paste0("UMAP Res: ", res)))
  dev.off()
  
  pdf(file.path(output_dir, sprintf("Step2_Heatmap_res%.1f.pdf", res)), width = 12, height = 8)
  print(DoHeatmap(seurat_obj, features = top10$gene) + NoLegend())
  dev.off()
  
  saveRDS(seurat_obj, file.path(output_dir, sprintf("Step2_%s_res%.1f.rds", species_name, res)))
  cat("  - Completed resolution", res, "\n")
}

# ============================================================================
# Step 3: Cell Type Annotation (Based on res=0.5)
# ============================================================================
cat("\nSTEP 3: Cell Type Annotation (Resolution 0.5)\n")
cat("-", rep("-", 70), "\n", sep="")

# Load res=0.5 object from output_dir
seurat_obj <- readRDS(file.path(output_dir, "Step3_species_name_res0.5.rds"))
Idents(seurat_obj) <- "RNA_snn_res.0.5"

# --- Method 1: Ath Assay ---
gene_map_filtered <- gene_map[gene_map$TargetGene %in% rownames(seurat_obj), ]
raw_counts <- seurat_obj[["RNA"]]$counts
sub_counts <- raw_counts[gene_map_filtered$TargetGene, ]
mapped_matrix <- rowsum(as.matrix(sub_counts), group = gene_map_filtered$AthGene)
mapped_matrix <- as(mapped_matrix, "sparseMatrix")

seurat_obj[["Ath"]] <- CreateAssayObject(counts = mapped_matrix)
seurat_obj <- NormalizeData(seurat_obj, assay = "Ath")
ath_markers <- FindAllMarkers(seurat_obj, assay = "Ath", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58, verbose = FALSE)

# Match with PCMDB (Plant Cell Marker Database)
pcmdb <- read.csv(pcmdb_file, stringsAsFactors = FALSE)

# add details 
cluster_annotations <- data.frame(cluster=character(), cell_type=character(), 
                                   n_matched_types=integer(), details=character(), markers=character())

for (clust in unique(ath_markers$cluster)) {
  cluster_markers <- ath_markers %>% filter(cluster == clust) %>% pull(gene)
  
  #  Cluster to PCMDB 
  overlap_data <- pcmdb %>% filter(Gene_id %in% cluster_markers)
  
  if (nrow(overlap_data) > 0) {
    
    # find Cell-Marker gene and Tissue enriched genes
    t1_data <- overlap_data %>% 
      filter(Note_gene_type %in% c("Cell-Marker gene", "Tissue enriched genes"))
    
    if (nrow(t1_data) > 0) {
      #---
      type_counts <- t1_data %>%
        group_by(Cell_type) %>%
        summarise(
          n_cell_marker = sum(Note_gene_type == "Cell-Marker gene", na.rm = TRUE),
          n_tissue = sum(Note_gene_type == "Tissue enriched genes", na.rm = TRUE),
          total_markers = n_distinct(Gene_id),
          .groups = 'drop'
        ) %>%
        #  sort as Cell-Marker and Tissue 
        arrange(desc(n_cell_marker), desc(n_tissue)) %>%
        # =================  =================
        head(3) 
      
      # 
      final_type <- paste(type_counts$Cell_type, collapse = "-")
      
      # detail as "root cortex(M:2, T:1) | root endodermis(M:0, T:3)"
      details_str <- paste(
        sprintf("%s(Marker:%d, Tissue:%d)", type_counts$Cell_type, type_counts$n_cell_marker, type_counts$n_tissue),
        collapse = " | "
      )
      
      all_matched_genes <- unique(t1_data$Gene_id)
      marker_str <- paste(all_matched_genes[1:min(10, length(all_matched_genes))], collapse = ", ")
      n_types <- nrow(type_counts)
      
    } else {
      # --- DEG as planB ---
      t2_data <- overlap_data %>% filter(Note_gene_type == "Differentially expressed gene")
      
      if (nrow(t2_data) > 0) {
        type_counts <- t2_data %>%
          group_by(Cell_type) %>%
          summarise(n_deg = n_distinct(Gene_id), .groups = 'drop') %>%
          arrange(desc(n_deg)) %>%
          # =================  =================
          head(3)
        
        final_type <- paste(type_counts$Cell_type, collapse = "-")
        details_str <- paste(sprintf("%s(DEG:%d)", type_counts$Cell_type, type_counts$n_deg), collapse = " | ")
        
        all_matched_genes <- unique(t2_data$Gene_id)
        marker_str <- paste(all_matched_genes[1:min(10, length(all_matched_genes))], collapse = ", ")
        n_types <- nrow(type_counts)
      } else {
        # 
        final_type <- "Unknown"; n_types <- 0; details_str <- "None"; marker_str <- ""
      }
    }
  } else {
    # 
    final_type <- "Unknown"; n_types <- 0; details_str <- "None"; marker_str <- ""
  }
  
  cluster_annotations <- rbind(cluster_annotations, data.frame(
    cluster = clust, 
    cell_type = final_type, 
    details = details_str,
    markers = marker_str
  ))
}

#  back to Seurat 
seurat_obj$cell_type <- cluster_annotations$cell_type[match(as.character(Idents(seurat_obj)), 
                                                             as.character(cluster_annotations$cluster))]
write.csv(cluster_annotations, file.path(output_dir, "Step3_Cluster_Annotations.csv"), row.names = FALSE)

# ============================================================================
#  PCMDB annotation Visualization 
# ============================================================================

# 1. 构建超大容量的 CNS 级高辨识度色库 (共 36 种独立颜色)
cns_mega_palette <- c(
  # Nature Publishing Group (NPG)
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
  # Science / D3.js (高饱和离散色)
  "#1F77B4", "#FF7F0E", "#2CA02C", "#9467BD", "#8C564B", "#E377C2", "#BCBD22", "#17BECF",
  # Cell / IGV (温和对照色)
  "#CE3D32", "#749B58", "#F0E685", "#466983", "#BA6338", "#5DB1DD", "#802268", "#6BD76B", "#D595A7", "#924822",
  # 补充高对比度色
  "#FF9896", "#C5B0D5", "#C49C94", "#F7B6D2", "#DBDB8D", "#9EDAE5", "#393B79", "#5254A3"
)

# 动态取色函数：优先使用原生高辨识度颜色，超载时才启动平滑插值
get_distinct_colors <- function(n, palette = cns_mega_palette) {
  if (n <= length(palette)) {
    return(palette[1:n])
  } else {
    message(sprintf("Warning: %d colors requested, exceeding the distinct palette limit (%d). Fallback to interpolation applied.", n, length(palette)))
    return(colorRampPalette(palette)(n))
  }
}

# 提取 Cluster 和 Cell Type 的数量
num_clusters <- length(unique(seurat_obj$RNA_snn_res.0.5))
num_celltypes <- length(unique(seurat_obj$cell_type))

# 安全获取配色
colors_clusters <- get_distinct_colors(num_clusters)
colors_celltypes <- get_distinct_colors(num_celltypes)

# 将 "Unknown" 细胞群强制设为高级灰 (从颜色池中安全覆盖)
names(colors_celltypes) <- unique(seurat_obj$cell_type)
if ("Unknown" %in% names(colors_celltypes)) {
  colors_celltypes["Unknown"] <- "#CCCCCC" 
}

# 2. 构建独立图层
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "RNA_snn_res.0.5", 
              label = TRUE, label.size = 5, cols = colors_clusters) + 
      ggtitle("Clusters (res=0.5)") +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "none"
      )

p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type", 
              label = FALSE, cols = colors_celltypes) + 
      ggtitle("Annotated Cell Types") +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right",
        legend.text = element_text(size = 9),
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm")
      ) +
      guides(color = guide_legend(override.aes = list(size = 4), ncol = 1)) # 强制图例单列显示，避免长文字错位

# 3. 组合并输出 (根据数量动态微调画布宽度)
# 如果细胞类型超过 15 种，右侧图例会很长，适当增加整体图片宽度以保证 UMAP 不被挤压
final_width <- ifelse(num_celltypes > 15, 20, 18)

pdf(file.path(output_dir, "Step3_Final_UMAP_Combined_CNS.pdf"), width = final_width, height = 7)
print(cowplot::plot_grid(p1, p2, ncol = 2, rel_widths = c(1, 1.6)))
dev.off()

# ... 后续的 saveRDS 和 summary 代码保持不变 ...
saveRDS(seurat_obj, file.path(output_dir, "Step3_species_name_PCMDB.rds"))
cat("  Saved: Step6 PCMDB files \n\n")

# Summary Statistics
cat("\n", rep("=", 70), "\n", sep="")
cat("ANALYSIS COMPLETE\n")
cat("Summary:\n")
cat("  Final cells:", ncol(seurat_obj), "\n")
cat("  Cell types identified:", length(unique(seurat_obj$cell_type)), "\n")
cat("  Annotation rate:", round(sum(seurat_obj$cell_type != "Unknown") / ncol(seurat_obj) * 100, 1), "%\n")
cat("  Output directory:", output_dir, "\n")
cat("\nCell Type Distribution:\n")
print(table(seurat_obj$cell_type))

```



---

## Key Parameters Reference

| Parameter | Value | Description |
|-----------|-------|-------------|
| `min.cells` | 3 | Minimum cells for gene retention |
| `min.features` | 200 | Minimum genes per cell (lower bound) |
| `nFeature_RNA` | 300-6000 | QC filtering range |
| `percent.mito/chloro` | < 5% | Organelle contamination threshold |
| `nfeatures` | 3000 | Variable features count |
| `dims` | 1:25 | PCs for clustering |
| `resolutions` | 0.5, 1.0, 1.5, 2.0 | Multi-resolution testing |
| `min.pct` | 0.25 | Marker detection threshold |
| `logfc.threshold` | 0.58 | Log2FC cutoff for markers |


## Input File Formats

### Gene Homology Mapping
```
TargetGene1    AT1G01010
TargetGene2    AT2G32910
...
```

### PCMDB Database Format
CSV with columns: `Cell_type`, `Gene_id`, plus additional metadata columns.

### Functional Annotation Format
Tab-separated with gene ID as first column, followed by annotation columns (NR, SwissProt, GO, KEGG, etc.)

---

## Expected Outputs

| File | Description |
|------|-------------|
| `Step1_Species_raw.rds` | Raw Seurat object with QC metrics |
| `Step1_5_Species_filtered.rds` | Filtered Seurat object |
| `Step2_Species_res[0.5-2.0].rds` | Clustering results per resolution |
| `Step2_All_Markers_res[0.5-2.0].csv` | Marker genes per resolution |
| `Step5_Cluster_Annotations.csv` | Cluster → CellType mapping (PCMDB) |



## Troubleshooting

| Issue | Solution |
|-------|----------|
| Low annotation rate (<50%) | Check gene mapping quality; consider using multiple reference species |
| Too many Unknown clusters | Lower threshold to ≥2 markers or check PCMDB format |
| Empty Ath assay | Verify gene mapping file format and TargetGene IDs match Seurat object |
| Memory issues with large datasets | Process markers in batches; use `verbose=FALSE` |


---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 2.0 | 2024-04-06 | Complete code template |
| 1.3 | 2024-04-05 | Prefix-based QC & LLM updates |
| 1.2 | 2024-04-01 | Added multi-resolution logic |
| 1.0 | 2024-03-20 | Initial documentation |
