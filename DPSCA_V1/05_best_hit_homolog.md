# Best-Hit Homolog Identification Pipeline for Plant Species

## Overview

This skill provides a comprehensive workflow for identifying best-hit homologous genes between a target plant species and reference species using OrthoFinder, BLASTP, and phyloMapper.py.

**Use Case**: Identify high-confidence orthologs between a newly sequenced plant species and well-annotated model species (e.g., Arabidopsis, Rice, Maize).

**Workflow Duration**: 4-8 hours (depending on dataset size)

---

## Prerequisites

### Software Requirements
- **Conda environment**: `orthofinder` with OrthoFinder 2.5.4+
- **Python packages**: pandas, biopython
- **Tools**: BLAST+, OrthoFinder

### Input Files
1. **Target species protein sequences**: FASTA format (e.g., `target_species.fa`)
2. **Reference species protein sequences**: Provided in `Ref_data/` directory
   - `Ath.fasta` (Arabidopsis thaliana)
   - `Osativa.fasta` (Oryza sativa)
   - `Zmays.fasta` (Zea mays)
3. **phyloMapper.py script**: Provided in `Ref_data/` directory

---

## Workflow Steps

### Step 1: Environment Setup

```bash
# Create working directory
WORKDIR="workdir_best_homolog_<project_name>"
mkdir -p $WORKDIR
cd $WORKDIR

# Activate conda environment
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate orthofinder

# Verify installation
orthofinder -h
makeblastdb -version
```

### Step 2: Prepare Input Files

#### 2.1 Create target gene list
```bash
# Extract gene IDs from target species protein file
grep '^>' target_species.fa | tr -d '>' | cut -d' ' -f1 > Target_gene.list

# Verify gene count
wc -l Target_gene.list
```

#### 2.2 Merge reference species proteins
```bash
# Concatenate all reference species protein files
cat ref_species1.fa ref_species2.fa ref_species3.fa > Ref_pep.fasta

# Verify merged file
ls -lh Ref_pep.fasta
grep -c '^>' Ref_pep.fasta
```

#### 2.3 Build BLAST database
```bash
# Create BLAST protein database
makeblastdb -in Ref_pep.fasta -dbtype prot -out Ref_pep.fasta

# Verify database creation (should create .phr, .pin, .psq files)
ls -lh Ref_pep.fasta.*
```

#### 2.4 Organize OrthoFinder input directory
```bash
# Create orthofinder input directory
mkdir -p orthofinder_input

# Copy protein files with clear naming
cp target_species.fa orthofinder_input/Target_species.fa
cp ref_species1.fa orthofinder_input/Ath.fa      # Arabidopsis
cp ref_species2.fa orthofinder_input/Osativa.fa  # Rice
cp ref_species3.fa orthofinder_input/Zmays.fa    # Maize

# Verify all files
ls -lh orthofinder_input/
```

---

### Step 3: Run OrthoFinder (Manual Execution Required)

> [!CAUTION]
> **Known Issue**: OrthoFinder may encounter bugs (TypeError, memory errors). Manual execution is recommended.

```bash
# Run OrthoFinder with optimal parameters
cd $WORKDIR

orthofinder -f orthofinder_input \
  -t 16 \
  -a 8 \
  -S diamond

# Parameters:
# -f : Input directory
# -t : Number of BLAST threads (adjust based on CPU)
# -a : Number of analysis threads
# -S : Sequence search program (diamond is faster)
```

**Expected Output:**
- Results directory: `orthofinder_input/OrthoFinder/Results_<date>/`
- Orthogroups: `Orthogroups/Orthogroups.tsv`
- Gene trees: `Gene_Trees/` directory
- Orthologues: `Orthologues/` directory

**Monitoring:**
```bash
# Check progress
ls -lh orthofinder_input/OrthoFinder/

# Verify orthogroups
wc -l orthofinder_input/OrthoFinder/Results_*/Orthogroups/Orthogroups.tsv

# Count gene trees
ls orthofinder_input/OrthoFinder/Results_*/Gene_Trees/ | wc -l
```

---

### Step 4: Run BLASTP Alignment (Manual Execution Required)

> [!WARNING]
> **Long Runtime**: BLASTP alignment takes 2-4 hours for ~50K queries against ~100K references.

```bash
cd $WORKDIR

# Run BLASTP search
blastp -query orthofinder_input/Target_species.fa \
  -db Ref_pep.fasta \
  -out blast_results.tsv \
  -evalue 1e-5 \
  -num_threads 16 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"

# Monitor progress
ls -lh blast_results.tsv
wc -l blast_results.tsv
```

**Expected Output:**
- File: `blast_results.tsv` (~30-40MB for 50K queries)
- Format: Tab-separated (12 columns)
- Alignment count: ~450K-500K hits

---

### Step 5: Identify Best Homologs with phyloMapper.py

```bash
cd $WORKDIR

# Run phyloMapper.py to integrate evidence
python /path/to/phyloMapper.py \
  --tree_dir orthofinder_input/OrthoFinder/Results_Mar13/Gene_Trees \
  --og_file orthofinder_input/OrthoFinder/Results_Mar13/Orthogroups/Orthogroups.tsv \
  --blast_file blast_results.tsv \
  --target_ids Target_gene.list \
  --target_species Target_species \
  --ref_species Ath,Osativa,Zmays \
  --output Query_Refspeices_homolog.tsv

# Verify output
ls -lh Query_Refspeices_homolog.tsv
wc -l Query_Refspeices_homolog.tsv
head -5 Query_Refspeices_homolog.tsv
```

**Scoring System:**
- OG membership: 0.2 points
- Phylogenetic distance (closest in tree): 0.4 points
- Highest BLAST bitscore: 0.35 points
- Single-copy ortholog bonus: 0.05 points
- **Maximum score**: 1.0

**Output Format:**
```
Target_Gene  OG_ID  Ath_Best_Match  Ath_Score  Ath_Other_Candidates  Osativa_Best_Match  ...
```

---

### Step 6: Extract High-Confidence Homologs

```bash
cd $WORKDIR

# Extract gene pairs with score >= 0.6 for each reference species
awk -F'\t' 'NR>1 && $4!="NA" && $4>=0.6 {print $1"\t"$3}' Query_Refspeices_homolog.tsv > Target_species_Ath_best_hit.txt

awk -F'\t' 'NR>1 && $7!="NA" && $7>=0.6 {print $1"\t"$6}' Query_Refspeices_homolog.tsv > Target_species_Osativa_best_hit.txt

awk -F'\t' 'NR>1 && $10!="NA" && $10>=0.6 {print $1"\t"$9}' Query_Refspeices_homolog.tsv > Target_species_Zmays_best_hit.txt

# Verify outputs
ls -lh Target_species_*_best_hit.txt
wc -l Target_species_*_best_hit.txt
```

**Output Format:**
```
query_gene    ref_gene
transcript1000    AT5G04930
transcript10000    LOC_Os03g21530
```

---

## Expected Results

### Summary Statistics

| Metric | Typical Value |
|--------|--------------|
| **Total target genes** | ~50,000 |
| **Orthogroups identified** | ~23,000 |
| **Gene trees generated** | ~14,000 |
| **BLAST alignments** | ~450,000 |
| **High-confidence matches (score >= 0.6)** | |
| - Arabidopsis | ~26,000 (52%) |
| - Rice | ~24,000 (48%) |
| - Maize | ~22,000 (44%) |

### Output Files

1. **Target_gene.list**: List of target species gene IDs
2. **Ref_pep.fasta**: Merged reference species protein sequences
3. **Ref_pep.fasta.***: BLAST database files
4. **orthofinder_input/**: Organized input directory
5. **orthofinder_input/OrthoFinder/Results_*/**: OrthoFinder results
6. **blast_results.tsv**: BLASTP alignment results
7. **Query_Refspeices_homolog.tsv**: Integrated homolog scores
8. **Target_species_*_best_hit.txt**: High-confidence gene pairs for each reference species

---

## Troubleshooting

### Common Issues

#### 1. OrthoFinder Memory Errors
**Symptom**: Process terminates during gene tree construction
**Solution**:
```bash
# Reduce thread count
orthofinder -f orthofinder_input -t 8 -a 4

# Or use resume mode if interrupted
orthofinder -f orthofinder_input -b previous_results/
```

#### 2. Empty Gene Tree Files
**Symptom**: Gene tree files exist but are empty (0 bytes)
**Solution**:
- This is expected for some orthogroups
- phyloMapper.py handles missing trees gracefully
- BLAST and OG evidence still provide robust scoring

#### 3. Low Homolog Coverage
**Symptom**: <70% of target genes have matches
**Possible causes**:
- Target species is evolutionarily distant from references
- Protein sequences are incomplete or fragmented
- E-value threshold too stringent

**Solution**:
```bash
# Relax BLAST e-value threshold
blastp ... -evalue 1e-3 ...

# Include more reference species
```

#### 4. BLASTP Hanging
**Symptom**: No output file growth
**Monitoring**:
```bash
# Check process
top -u $USER

# Monitor file growth
watch -n 60 'ls -lh blast_results.tsv'
```

**Solution**:
- Reduce threads: `-num_threads 8`
- Split query file into smaller batches
- Use DIAMOND instead of BLAST: `-S diamond` in OrthoFinder

---

## Quality Control Checklist

- [ ] All input files exist and are properly formatted
- [ ] Gene count matches expectations
- [ ] BLAST database created successfully
- [ ] OrthoFinder completed without errors
- [ ] Gene trees generated for >50% of orthogroups
- [ ] BLASTP completed with expected alignment count
- [ ] phyloMapper.py output has correct format
- [ ] High-confidence matches extracted for all reference species
- [ ] Output files contain expected number of gene pairs

---

## Advanced Options

### Custom Score Thresholds

```bash
# Extract with different thresholds
# Score >= 0.8 (very high confidence)
awk -F'\t' 'NR>1 && $4>=0.8 {print $1"\t"$3}' Query_Refspeices_homolog.tsv > high_confidence_pairs.txt

# Score >= 0.4 (moderate confidence)
awk -F'\t' 'NR>1 && $4>=0.4 {print $1"\t"$3}' Query_Refspeices_homolog.tsv > moderate_confidence_pairs.txt
```

### Filter by Orthogroup Type

```bash
# Extract single-copy orthologs only
awk -F'\t' 'NR>1 && $4>=0.6 && $2 ~ /OG/ {print $1"\t"$3}' Query_Refspeices_homolog.tsv
```

### Add Additional Reference Species

```bash
# Simply add more protein files to orthofinder_input/
cp new_species.fa orthofinder_input/NewSpecies.fa

# Re-run OrthoFinder (will analyze all species together)
orthofinder -f orthofinder_input -t 16 -a 8
```

---

## Performance Optimization

### For Large Datasets (>100K genes)

```bash
# Use DIAMOND for faster search
orthofinder -f orthofinder_input -S diamond_ultra_sens

# Reduce memory usage
orthofinder -f orthofinder_input -t 8 -a 4 -M 32000
```

### For Quick Testing

```bash
# Run on subset of data
head -1000 target_species.fa > test_subset.fa
# ... repeat workflow with test_subset.fa
```

---

## References

- **OrthoFinder**: Emms, D.M., Kelly, S. (2019). Genome Biology 20:238
- **BLAST+**: Camacho et al. (2009). BMC Bioinformatics 10:421
- **phyloMapper.py**: Custom script integrating phylogeny, orthology, and BLAST evidence

---

## Contact & Support

For issues or questions:
1. Check the troubleshooting section above
2. Verify input file formats and paths
3. Review log files for error messages
4. Consult OrthoFinder documentation: https://github.com/davidemms/OrthoFinder

---

## Version History

- **v1.0** (2024-03): Initial workflow documentation
- Tested with: OrthoFinder 2.5.4, BLAST 2.15.0, Python 3.8+
- Validated on: Plant species (target: 50K genes, references: Arabidopsis, Rice, Maize)
