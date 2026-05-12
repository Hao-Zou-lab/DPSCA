# LLM annotation Pipeline (Plant Species)

|> [!TIP]
> workflow for single-cell RNA-seq analysis of non-model plant species using LLM-assisted functional annotation.

## Overview

This pipeline enables cell type annotation for non-model plant species by:
1. identified marker gene for LLM 
2. Get gene function of these marker genes
3. Cell Type Inference by LLM



---

## step 6 LLM-Based Cell Type Annotation Workflow


```r
# ---  LLM-Based Cell Type Annotation ---
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
output_dir <- file.path(workdir, "step6_LLM")
last_dir  <- file.path(workdir, "step5_clusting_annotation")
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
# ============================================================================
# ---  LLM Input Preparation ---
# ============================================================================
func_anno <- read.table(func_anno_file, header = TRUE, sep = "\t", 
                        stringsAsFactors = FALSE, comment.char = "", 
                        quote = "", fill = TRUE)
rna_markers <- read.csv(file.path(last_dir, "Step2_All_Markers_res0.5.csv"), stringsAsFactors = FALSE)
top20_rna <- rna_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

# Auto-detect gene ID column
gene_col <- intersect(colnames(func_anno), c("query-id", "GeneID", "gene_id", "Gene", "gene"))[1]
llm_input <- left_join(top20_rna, func_anno, by = setNames(gene_col, "gene"))
write.csv(llm_input, file.path(output_dir, "Step6_LLM_Input_Markers.csv"), row.names = FALSE)




# LLM-based annotation uses functional descriptions (NR/SwissProt) of marker genes to infer cell types based on biological knowledge, providing an alternative to homology-based annotation.
# Read LLM input markers
llm_markers <- read.csv("Step6_LLM_Input_Markers.csv", stringsAsFactors = FALSE)

# Preview structure
cat("Columns:", colnames(llm_markers), "\n")
cat("Total rows:", nrow(llm_markers), "\n")
cat("Clusters:", length(unique(llm_markers$cluster)), "\n")

# Extract key functional annotation columns
# Adjust column names based on your annotation file
nr_col <- "NR_Annotation"
swissprot_col <- "SwissProt_Annotation"

# View markers for specific cluster
cluster_to_check <- 0

write.csv(llm_markers %>% 
  filter(cluster == cluster_to_check) %>% 
  select(cluster, gene, avg_log2FC, all_of(nr_col, swissprot_col)), file.path(output_dir, "Step6_LLM_Markers_Function.csv"), row.names = FALSE)
cat("  Saved: LLM Input files\n")
```

# ============================================================================
### Step:  Cell Type Inference
# ============================================================================
> [!IMPORTANT]
> Not in R environment
### Input Requirements
- `Step6_LLM_Markers_Function.csv` (Functional annotation file with NR/SwissProt columns)


For each cluster, analyze:
1. **Functional themes**: Group markers by biological function
2. **Tissue-specificity**: Identify markers known to be expressed in specific cell types
3. **Pathways**: Check for pathway enrichment (e.g., lignin biosynthesis for xylem)

## LLM Prompt Template

Use this prompt with AI assistants (Claude, GPT-4, etc.) to automate cell type inference:

```markdown
You are an expert in plant single-cell transcriptomics and cell biology.

I will provide you with Top 20 Marker genes for each cluster from a non-model plant species, along with their functional annotations from NR/SwissProt databases.

**Task**: Infer the likely cell type for each cluster based on:
1. Functional descriptions of marker genes
2. Known tissue-specific expression patterns in plants
3. Biological pathways and processes represented

**Output Format**: Provide a Markdown table with:
- Cluster number
- Inferred Cell Type
- Confidence level (High/Medium/Low)
- Key Markers (top 5-10 genes)
- Justification (biological reasoning)

**Context**: This is root tissue data. Consider common root cell types:
- Trichoblast (root hair cells)
- Root cortex/epidermis/endodermis
- Root cap/columella
- Phloem and Xylem (vascular tissues)
- Root apical meristem (dividing cells)
- Quiescent center (stem cell niche)

Focus on identifying unique functional signatures that distinguish cell types.

```
## Interpretation Guidelines

### High-Confidence Indicators by Cell Type

| Cell Type | Key Indicator Genes | Interpretation |
|-----------|---------------------|----------------|
| **Trichoblast** | Aquaporins (TIP, PIP), Nitrate transporters (NRT), Amino acid transporters | High membrane transport activity |
| **Root cortex** | Peroxidases (PER), Chitinases, Cytochrome P450 | Defense response, secondary metabolism |
| **Root endodermis** | Auxin transporters (LAX, PIN), Thionins, NIP aquaporins | Barrier function, hormone transport |
| **Xylem** | Laccases (LAC), Berberine bridge enzyme, Dirigent proteins, PAL | Lignin polymerization, secondary cell wall |
| **Phloem** | β-1,3-glucosidases, NO VEIN protein, Dof TFs, Tetraspanins | Sugar transport, vascular development |
| **RAM** | Histones (H2A, H2B, H3, H4), RNR, DNA replication | Active cell division |
| **Quiescent center** | CLAVATA1, CML proteins, Highly expressed unknowns | Stem cell niche |
| **Root cap** | ALMT, Glycine-rich proteins, RBOH | Metal transport, stress response |
| **Columella** | Proline-rich proteins, Pectinesterases, FLA proteins | Gravity sensing |
| **Root stele** | NPF transporters, GASA proteins | Vascular transport |

### Interpretation Workflow

```
For each cluster:
1. Extract top 20 markers with functional annotations
2. Identify dominant functional themes:
   - Transport activity (aquaporins, transporters)
   - Cell wall proteins (pectinesterases, extensins, TBL)
   - Defense proteins (peroxidases, chitinases, thionins)
   - Cell division (histones, RNR)
   - Secondary metabolism (P450, laccases, PAL)
   - Hormone signaling (auxin-responsive, gibberellin-regulated)
3. Match to known cell type signatures
4. Assign confidence level:
   - High: Multiple specific markers with consistent theme
   - Medium: Some specific markers but mixed signals
   - Low: Few informative markers or conflicting signals
5. Document justification with biological reasoning
```

### Common Pitfalls

> [!CAUTION]
> **Common interpretation mistakes**:
>
> 1. **Multiple cell types in one cluster**: If markers show mixed functions, the cluster may contain multiple cell types
> 2. **Generic housekeeping genes**: Ribosomal proteins, mitochondrial genes are NOT cell-type-specific
> 3. **Unknown annotations**: If most markers lack functional info, treat as "Unknown" with low confidence
> 4. **Cross-species homology**: Arabidopsis annotations may not perfectly map to your species
> 5. **Validation**: Always compare with UMAP visualization to check spatial consistency


### Output File Format

After LLM annotation, generate `Step6_LLM_CellType_Annotation.md` with:

```markdown
# LLM-Based Cell Type Annotation Report

**Method**: Functional annotation analysis of Top 20 marker genes per cluster
**Input**: Step5_LLM_Input_Markers.csv
**Date**: YYYY-MM-DD

---

## Annotation Results Summary

| Cluster | Inferred Cell Type | Confidence | Key Evidence |
|---------|-------------------|------------|--------------|
| 0 | Trichoblast | High | Aquaporins, transporters - high transport activity |
| ... | ... | ... | ... |

---

## Comparison with PCMDB Annotation

| Cluster | PCMDB Result | LLM Result | Agreement |
|--------|-------------|-----------|-----------|
| 0 | trichoblast | Trichoblast | ✅ |
| ... | ... | ... | ... |

---

## Final Recommendations

1. **High-confidence annotations** (both methods agree): Use these directly
2. **Discrepancies**: Manually review marker genes and UMAP spatial distribution
3. **Unknown**: Requires additional validation or marker gene analysis
```
### Output
- `Step6_LLM_CellType_Annotation.md` - LLM annotation results with justifications
---

### Step: Create Annotation Table
> [!IMPORTANT]
> based on `Step6_LLM_CellType_Annotation.md`
create annotation table based on marker analysis:

```r
# Create annotation data frame
llm_annotations <- data.frame(
  cluster = character(),
  cell_type = character(),
  confidence = character(),
  key_markers = character(),
  justification = character(),
  stringsAsFactors = FALSE
)

# Add annotations manually based on analysis
# Example:
llm_annotations <- rbind(llm_annotations, data.frame(
  cluster = "0",
  cell_type = "Trichoblast",
  confidence = "High",
  key_markers = "Aquaporin TIP2-2, NRT3.1, Amino acid transporter",
  justification = "High transport activity characteristic of root hair cells"
))

# Continue for all clusters...

# Save annotations
write.csv(llm_annotations, file.path(output_dir, "Step6_LLM_CellType_Annotation.csv"), row.names = FALSE)
```





# ============================================================================
# Integration with Seurat Object (PCMDB Top3 + LLM)
# ============================================================================
> [!IMPORTANT]
> After generating `Step6_LLM_CellType_Annotation.csv`, combined LLM-based cell type annotation and PCMDB results.
```r
# ---  LLM-Based Cell Type Annotation ---
# ============================================================================
# CONFIGURATION - Modify these for your project
# ============================================================================
workdir <- "/path/to/your/project"
species_name <- "YourSpecies"  # Used for output file naming

# Input paths
matrix_dir <- file.path(workdir, "step3_cellranger_count/YOUR_SAMPLE/outs/filtered_feature_bc_matrix")
gene_mapping_file <- file.path(workdir, "step4_homolog_mapping/Species_Ath_best_hit.txt")
pcmdb_file <- file.path(workdir, "PCMDB_Ath.csv")  # Plant Cell Marker Database
func_anno_file <- file.path(workdir, "gene_annotation.txt")  # Functional annotations

# Output directory (all results go here)
output_dir <- file.path(workdir, "step6_LLM")
last_dir  <- file.path(workdir, "step5_clusting_annotation")
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

# 1. 加载 LLM 注释结果
llm_annotations <- read.csv(file.path(output_dir, "Step6_LLM_CellType_Annotation.csv"), stringsAsFactors = FALSE)
seurat_obj<-readRDS(file.path(last_dir, "Step3_Platiflora_PCMDB.rds"))
# 2. 综合拼接：PCMDB Top3 / LLM
clusters <- as.character(unique(Idents(seurat_obj)))
combined_map <- sapply(clusters, function(cl) {
  
  # 获取当前 cluster 的标记基因
  cl_markers <- ath_markers$gene[ath_markers$cluster == cl]
  
  # 提取 PCMDB 前 3 名匹配结果 (基于 Step 5 的重叠逻辑)
  # 筛选重叠基因数 >= 3 的类型并按重叠数降序排列
  pcmdb_top3_df <- pcmdb %>% 
    filter(Gene_id %in% cl_markers) %>% 
    count(Cell_type, sort = TRUE) %>% 
    filter(n >= 3) %>% 
    head(3)
  
  # 拼接 PCMDB 结果
  pcmdb_str <- if(nrow(pcmdb_top3_df) > 0) {
    paste(pcmdb_top3_df$Cell_type, collapse = "/")
  } else {
    "Unknown"
  }
  
  # 提取 LLM 注释结果
  # 匹配 cluster 列，获取 cell_type 列
  llm_val <- llm_annotations$cell_type[match(cl, llm_annotations$cluster)]
  llm_str <- if(!is.na(llm_val) && llm_val != "") llm_val else "Unknown"
  
  # 最终格式: PCMDB1/PCMDB2/PCMDB3/LLM
  return(paste0(pcmdb_str, "/", llm_str))
})

# 3. 将最终注释映射回 Seurat 对象
seurat_obj$cell_type_final <- combined_map[as.character(Idents(seurat_obj))]

# 4. 基于最终结果进行可视化绘图
pdf("Step6_Final_Combined_Annotation_UMAP.pdf", width = 16, height = 9)
# 使用 repel = TRUE 处理较长的标签
print(DimPlot(seurat_obj, group.by = "cell_type_final", label = TRUE, repel = TRUE) + 
      ggtitle("Final Combined Annotation (PCMDB Top3 / LLM)"))
dev.off()

# 5. 导出最终注释结果表格
final_anno_table <- data.frame(
  cluster = names(combined_map),
  combined_annotation = as.character(combined_map)
)
write.csv(final_anno_table, "Step6_Final_Annotation_Summary.csv", row.names = FALSE)

# 设置为默认标识符方便后续分析
Idents(seurat_obj) <- "cell_type_final"
```
saveRDS(seurat_obj, file.path(output_dir, "combine_annotation.rds"))
| **annotation threshold** | **≥3** | Overlapping markers with PCMDB |
| **LLM markers** | **Top 20** | Markers per cluster for LLM analysis |

---

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
|------|--------------|
| `Step5_LLM_Input_Markers.csv` | Top 20 markers with functional annotations |
| `Step5_LLM_CellType_Annotation.md` | LLM annotation results with justifications |
| `Step6_Final_UMAP_Combined.pdf` | Publication-ready dual-panel UMAP |

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Low annotation rate (<50%) | Check gene mapping quality; consider using multiple reference species |
| Too many Unknown clusters | Lower threshold to ≥2 markers or check PCMDB format |
| Empty Ath assay | Verify gene mapping file format and TargetGene IDs match Seurat object |
| Memory issues with large datasets | Process markers in batches; use `verbose=FALSE` |
| LLM annotations disagree with PCMDB | Manually review marker genes; check UMAP spatial distribution |

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 2.1 | 2024-04-06 | **Added LLM annotation workflow**: Detailed functional annotation-based cell type inference workflow with manual guide, LLM prompt template, and interpretation guidelines |
| 2.0 | 2024-04-06 | Complete code template |
| 1.3 | 2024-04-05 | Prefix-based QC & LLM updates |
| 1.2 | 2024-04-01 | Added multi-resolution logic |
| 1.0 | 2024-03-20 | Initial documentation |

