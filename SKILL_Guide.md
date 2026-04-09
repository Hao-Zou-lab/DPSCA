## De novo Plant Single-Cell Annotator Workflow

Complete 6-step workflow for non-model plant single-cell RNA-seq analysis, from
 transcriptome reference preparation to cell type annotation using Arabidopsis markers.

**Skill index**: [DPSCA/SKILL.md](./DPSCA/SKILL.md)



**Workflow Overview**:
```
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

**When to use**:
- Non-model plant species without reference genome
- plant cell type annotation



**Individual Steps**:

| Step | Skill File | Description |
|------|------------|-------------|
| 1 | [01_transcriptome_reference.md](./DPSCA/01_transcriptome_reference.md) | NGS assembly or consensused of PacBio or hybrid correction(NGS+ Pacbio)|
| 2 | [02_Transdecoder.md](./DPSCA/02_Transdecoder.md) | ORF finding |
| 3 | [03_prepare_cellranger_reference.md](./DPSCA/03_prepare_cellranger_reference.md) | CellRanger reference preparation |
| 4 | [04_cellranger_count.md](./DPSCA/04_cellranger_count.md) | CellRanger count |
| 5 | [05_best_hit_homolog.md](./DPSCA/05_best_hit_homolog.md) | Best hit homolog finding |
| 6 | [06_clusting_transfer.md](./DPSCA/06_clusting_transfer.md) | scRNA clusting and label transfer cell-type Annotation |
| 7 | [07_LLM_annotation.md](./DPSCA/07_LLM_annotation.md) | LLM cell-type Annotation |
| 8 | [08_merge_annotation.md](./DPSCA/08_merge_annotation.md) | Merge cell-type Annotation Results |

---
