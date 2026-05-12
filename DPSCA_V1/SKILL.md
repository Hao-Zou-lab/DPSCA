---
id: DPSCA_index
name: De novo Plant Single-Cell Annotator Workflow
description: |
  Complete 8-step workflow for non-model plant single-cell RNA-seq analysis, 
  from transcriptome reference preparation to cell type annotation using 
  Arabidopsis markers and LLM. Covers transcriptome assembly/correction, 
  ORF finding, CellRanger quantification, Seurat clustering, and dual annotation.
tags: [plant, scrna-seq, cross-species, annotation, arabidopsis, cellranger, seurat, non-model-species, de-novo, llm]
---

# De novo Plant Single-Cell Annotator Workflow

Complete 8-step workflow for non-model plant single-cell RNA-seq analysis, from transcriptome reference preparation to cell type annotation using Arabidopsis markers and Large Language Models (LLM).

## Workflow Overview

```text
[A: hybrid correction]
PacBio Long Reads + NGS Short Reads ──→ LORDEC Correction ───────┐
                                                                     │
[B: NGS denovo assembly]                                             │
NGS Short Reads ──────────────────────────→ Trinity De Novo Assembly ┼─→ transcriptome reference
                                                                     │             ↓
[C:PacBio consensusFLNC]                                             │      TransDecoder ORF
PacBio Long Reads ────────────────────────→ Iso-Seq Pipeline ────────┘             ↓
                                                                          CellRanger Reference
                                                                                   ↓
                                                                           CellRanger Count
                                                                                   ↓
                                                                      Finding Best Hit Homologs
                                                                                   ↓
                                                                           Seurat Clustering
                                                                                   ↓
                                                                 label transfer cell-type Annotation
                                                                                   ↓
                                                                        LLM cell-type Annotation
                                                                                   ↓
                                                                    Merge cell-type Annotation Results
```

## When to use

- Non-model plant species without reference genome
- Plant cell type annotation

---

## Skills Index

### Step 1: Transcriptome Reference Preparation

[01_transcriptome_reference.md](./DPSCA/01_transcriptome_reference.md)

**Description**: Create a baseline transcriptome reference for non-model plants using one of three approaches: NGS short-read *de novo* assembly (e.g., Trinity), PacBio long-read consensus (Iso-Seq), or hybrid error correction (LORDEC) combining both.

**Keywords**: transcriptome, de novo assembly, Trinity, Iso-Seq, LORDEC, PacBio, NGS

---

### Step 2: ORF Finding

[02_Transdecoder.md](./DPSCA/02_Transdecoder.md)

**Description**: Open Reading Frame (ORF) prediction using TransDecoder on the transcriptome reference. Produces high-quality CDS and protein sequences for each transcript.

**Keywords**: TransDecoder, ORF-prediction, CDS, protein-sequences

---

### Step 3: CellRanger Reference Preparation

[03_prepare_cellranger_reference.md](./DPSCA/03_prepare_cellranger_reference.md)

**Description**: Convert TransDecoder ORF predictions to CellRanger-compatible reference files. Selects the best ORF per transcript, generating reference FASTA and GTF files for single-cell analysis.

**Keywords**: CellRanger, reference-generation, GTF, FASTA

---

### Step 4: CellRanger Count

[04_cellranger_count.md](./DPSCA/04_cellranger_count.md)

**Description**: Build the reference database and quantify gene expression using CellRanger count for the non-model plant species.

**Keywords**: CellRanger, count, 10X, quantification, single-cell

---

### Step 5: Best Hit Homolog Finding

[05_best_hit_homolog.md](./DPSCA/05_best_hit_homolog.md)

**Description**: Identify the best orthologs between the target species and a reference species (e.g., Arabidopsis) using BLAST/DIAMOND and OrthoFinder.

**Keywords**: homolog, ortholog, BLAST, OrthoFinder, cross-species, Arabidopsis

---

### Step 6: scRNA Clustering and Label Transfer Cell-Type Annotation

[06_clusting_transfer.md](./DPSCA/06_clusting_transfer.md)

**Description**: Complete Seurat workflow with plant-specific QC (chloroplast/mitochondrion filtering) and multi-resolution clustering. Includes cross-species homology-based annotation (label transfer) using datasets like PlantscRNAdb (Ath assay).

**Keywords**: Seurat, clustering, label-transfer, cross-species, plant, QC

---

### Step 7: LLM Cell-Type Annotation

[07_LLM_annotation.md](./DPSCA/07_LLM_annotation.md)

**Description**: Utilize Large Language Models (LLM) for cell-type annotation based on gene functional descriptions (NR/SwissProt) to provide a secondary, independent layer of biological evidence (RNA assay).

**Keywords**: LLM, annotation, functional-description, cell-type

---

### Step 8: Merge Cell-Type Annotation Results

[08_merge_annotation.md](./DPSCA/08_merge_annotation.md)

**Description**: Merge results from the label transfer (Step 6) and LLM annotation (Step 7) to finalize robust cell type identification. Generates publication-quality visualizations.

**Keywords**: merge, consensus-annotation, visualization

---

## ⚠️ Critical Best Practices

> [!CAUTION]
> **Key lessons learned from real-world analysis**

### 1. Dual Annotation Verification
Always ensure the merged annotation (Step 8) cross-references both the Ath assay (homology) and RNA assay (LLM). Combined annotations should use an underscore format (e.g., "Root cap_Vascular tissue") when resolving discrepancies.

### 2. Use Merged Clusters for Downstream Analysis
After cluster merging, always use the `merged_cluster` column for all subsequent analyses (annotation, visualization, statistics). Never use original `seurat_clusters` after merging.

### 3. Verify Data Statistics at Every Step
Always verify cell counts directly from RDS files using `ncol(obj)` before writing reports. Never use cached or assumed values. Common errors include:
- Reporting wrong cell counts (e.g., 5,405 vs actual 3,958)
- Mixing raw vs filtered vs merged statistics
- Using marker counts instead of cell counts