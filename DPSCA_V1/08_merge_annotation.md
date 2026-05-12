# SKILL: Combine PCMDB and LLM Cell Types & Generate UMAP

> [!TIP]
> This SKILL script is used to merge the cell annotations based on PCMDB database matching and the LLM-generated cell annotations exactly as they are into a new column, `merge_celltype`. It then plots a UMAP with corner axes using a custom advanced light color palette (`light_palette`) and `scRNAtoolVis`.

## Dependency Preparation
Please ensure the following R packages are installed:
```R
install.packages("grDevices")
# If scRNAtoolVis is not installed, install it via GitHub: devtools::install_github("junjunlab/scRNAtoolVis")
```

---

## Core Integration and Plotting Code

```r
# ============================================================================
# Step 1: Load necessary packages
# ============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(grDevices)
  library(scRNAtoolVis)
})

# ============================================================================
# Step 2: Load data
# ============================================================================
# Load Seurat object (please adjust according to the actual path)
seurat_path <- '/mnt/97d9f81e-6799-4cb1-a909-27a485ae1cf6/Agent_files/Platiflora_NGS_SCRNA/step5_clusting_annotation/Step6_Platiflora_merge_celltype.rds'
seurat <- readRDS(seurat_path)

# Load LLM annotation results (please adjust according to the actual path)
# Assuming the LLM output file is Step6_LLM_CellType_Annotation.csv
llm_annotations <- read.csv("step6_LLM/Step6_LLM_CellType_Annotation.csv", stringsAsFactors = FALSE)

# ============================================================================
# Step 3: Synthesize merge_celltype exactly as is
# ============================================================================
# Assuming PCMDB annotations are stored in seurat$cell_type, and clustering info in seurat$seurat_clusters
clusters <- as.character(unique(seurat$seurat_clusters))

# Extract existing PCMDB annotation mapping from the Seurat object
pcmdb_map <- seurat@meta.data %>% 
  select(seurat_clusters, cell_type) %>% 
  distinct()

# Build mapping dictionary
merge_map <- sapply(as.character(seurat$seurat_clusters), function(cl) {
  # Get PCMDB annotation
  pcmdb_val <- pcmdb_map$cell_type[pcmdb_map$seurat_clusters == cl]
  pcmdb_str <- if(length(pcmdb_val) > 0 && !is.na(pcmdb_val[1])) pcmdb_val[1] else "Unknown"
  
  # Get LLM annotation
  llm_val <- llm_annotations$cell_type[llm_annotations$cluster == cl]
  llm_str <- if(length(llm_val) > 0 && !is.na(llm_val[1]) && llm_val[1] != "") llm_val[1] else "Unknown"
  
  # Concatenate exactly as is (Format: PCMDB / LLM)
  return(paste0(pcmdb_str, " / ", llm_str))
})

# Add the synthesized annotation results back to the Seurat object
seurat$merge_celltype <- merge_map

# ============================================================================
# Step 4: Palette definition and UMAP plotting
# ============================================================================

# Generate n light colors: evenly distributed hues, saturation S=0.3, value V=0.95
light_palette <- function(n, s = 0.3, v = 0.95) {
  hues <- seq(0, 1, length.out = n + 1)[1:n]
  hsv(h = hues, s = s, v = v)
}

# Generate an exclusive color palette for merge_celltype
unique_merge_types <- unique(seurat$merge_celltype)
n_merge_celltype <- length(unique_merge_types)
merge_celltype_colors <- light_palette(n_merge_celltype)

# Plot UMAP
p_merge <- clusterCornerAxes(object = seurat, reduction = 'umap',
                       clusterCol = "merge_celltype",  # Group by the newly generated merge_celltype
                       pSize = 0.1,
                       arrowType = 'open',      # Axis arrow type
                       lineTextcol = 'grey50',  # Corner line and label color
                       cornerTextSize = 3.5,
                       keySize = 1,             # Legend size
                       show.legend = TRUE,      # Show legend
                       cellLabel = FALSE,       # Do not display cell labels directly on the plot (to avoid overlap)
                       cellLabelSize = 5, 
                       noSplit = TRUE,
                       addCircle = FALSE,       # Enable subpopulation background cloud circles
                       cicAlpha = 0.1,          # Circle transparency
                       cicDelta = 0.5,
                       nbin = 200) + 
  scale_color_manual(values = merge_celltype_colors) +
  guides(color = guide_legend(ncol = 1,         # Single-column legend to avoid truncation of long text
                              title = "Merge Celltype\n(PCMDB / LLM)", 
                              override.aes = list(size = 2.5))) + 
  theme(legend.text = element_text(size = 10),  # Legend label font size
        legend.title = element_text(size = 12), # Legend title font size
        legend.position = "right")              # Legend position

# Display plot
print(p_merge)

# ============================================================================
# Step 5: Save results
# ============================================================================
# Since the synthesized names might be long, it's recommended to increase width to fully display the legend
output_dir <- '/home/a9/Plot_agent'
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

ggsave(plot = p_merge, 
       filename = file.path(output_dir, 'umap_merge_celltype.png'),
       width = 12, height = 6) 

# (Optional) Save the updated Seurat object
# saveRDS(seurat, file.path(output_dir, "Step6_Platiflora_merge_celltype_updated.rds"))
```

> [!CAUTION]
> **Plotting Fine-tuning Tips**:
> If the generated `merge_celltype` label text is too long, causing the right legend to overflow or squeeze the UMAP plot area, please appropriately increase the `width` parameter in the `ggsave` function (e.g., to 14 or 16), or adjust the font size in `theme(legend.text = element_text(size = ...))`.