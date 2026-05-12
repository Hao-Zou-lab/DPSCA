# Reference Data for Best-Hit Homolog Identification and Cell Type Annotation 

This directory contains essential reference data for the best-hit homolog identification pipeline.

## Contents


### Analysis Script

- **`phyloMapper.py`**: Python script for integrated homolog identification
  - Combines orthogroup, phylogenetic, and BLAST evidence
  - Generates comprehensive scoring (max score: 1.0)
  - Output: TSV format with best matches for each reference species

## Usage

### Quick Start

```bash
# Navigate to skill directory
cd /path/to/.pantheon/skills/omics/

# Copy reference data to your working directory
cp -r Ref_data/ /path/to/your/workdir/

# Use with the workflow
cd /path/to/your/workdir/
# ... follow steps in best_hit_homolog.md
```

### Integration with Workflow

1. **Step 2.2**: Use these files directly for merging reference proteins
   ```bash
   cat Ref_data/Ath.fasta Ref_data/Osativa.fasta Ref_data/Zmays.fasta > Ref_pep.fasta
   ```

2. **Step 5**: Use phyloMapper.py for integrated analysis
   ```bash
   python Ref_data/phyloMapper.py --tree_dir ... --og_file ... --blast_file ...
   ```

## Adding Custom Reference Species

To add your own reference species:

1. Obtain protein sequences in FASTA format
2. Name file as `<SpeciesAbbrev>.fasta` (e.g., `Gmax.fasta` for soybean)
3. Place in `Ref_data/` directory
4. Update workflow command to include new species

```bash
# Example: Adding soybean (Glycine max)
cp /path/to/soybean_proteins.fasta Ref_data/Gmax.fasta

# Merge with existing references
cat Ref_data/*.fasta > Ref_pep.fasta
```

## Data Sources

- **Arabidopsis thaliana**: TAIR database (https://www.arabidopsis.org/)
- **Oryza sativa**: Ensembl Plants (http://plants.ensembl.org/)
- **Zea mays**: MaizeGDB (https://www.maizegdb.org/)

## Version Information

- **Created**: 2024-03
- **OrthoFinder Version**: 2.5.4
- **Python Version**: 3.8+
- **Last Updated**: 2024-03

## Requirements for phyloMapper.py

```bash
# Python packages
pip install pandas biopython

# System requirements
- Python 3.8+
- Memory: 8GB+ RAM recommended
- Disk: 100MB+ free space
```

## Troubleshooting

### Issue: Memory errors with phyloMapper.py
**Solution**: Process smaller batches or increase available RAM

### Issue: Missing dependencies
**Solution**: Install required packages
```bash
pip install pandas biopython
```

### Issue: Incorrect file paths
**Solution**: Use absolute paths when calling phyloMapper.py

## Citation

If using this pipeline, please cite:
- OrthoFinder: Emms & Kelly (2019) Genome Biology 20:238
- BLAST: Camacho et al. (2009) BMC Bioinformatics 10:421
- phyloMapper.py: Custom script (included)
