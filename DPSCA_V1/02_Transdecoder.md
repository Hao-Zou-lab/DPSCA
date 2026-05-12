---
id: Transdecoder
name: TransDecoder ORF Prediction
description: |
   ORF prediction using TransDecoder.
  Produces high-quality CDS and protein sequences for each transcript.
tags: [ORF-prediction, transcriptome]
---

# TransDecoder ORF Prediction
 ORF prediction.

---

## Overview

This workflow performs:
1. **TransDecoder ORF prediction** - Identify coding regions in corrected transcripts
2. **Best ORF selection** - Select the longest, highest-scoring ORF per transcript
3. **Format Conversion** - Convert GFF3 to GTF format using gffread
4. **ID Simplification** - Use transcript IDs instead of ORF IDs (e.g., `transcript0` instead of `transcript0.p1`)

**Input**:
- transcripts (FASTA)


**Output**:
- Predicted CDS sequences (FASTA)
- Predicted protein sequences (FASTA)
- Reference FASTA (one sequence per transcript)
- Reference GTF (gene annotations in GTF format)
- Statistics file
---

## Environment Setup

### Required Software

```bash
# Check and install required tools
# Using Conda/Mamba (recommended)
conda install -c bioconda seqtk seqkit
conda install -c bioconda transdecoder gffread

# Or using system package manager
# LORDEC: http://www.atgc-montpellier.fr/downloads/software/Lordec/
# TransDecoder: https://github.com/TransDecoder/TransDecoder
```

### Verify Installation

```bash
# Check versions
seqtk version                      # seqtk 1.3+
TransDecoder.LongOrfs --version   # TransDecoder 5.7+
gffread --version  # gffread 0.12.7+
python --version  # Python 3.7+
```

### Environment File Template

Create `environment.yml`:

```yaml
name: Transdecoder
channels:
  - bioconda
  - conda-forge
dependencies:
  - python=3.11
  - lordec=0.6
  - seqtk=1.3
  - seqkit=2.13
  - transdecoder=5.7.1
  - gzip
  - gffread=0.12.7
```

Activate:
```bash
conda env create -f environment.yml
conda activate Transdecoder
```

---

## Workflow Steps

### Step 1: TransDecoder ORF Prediction

Identify candidate coding regions within transcripts.

#### 1.1 Predict Long ORFs

```bash
#!/bin/bash
WORKDIR="/path/to/workdir"
MIN_PROTEIN=100  # Minimum ORF length (amino acids)

# Create output directory
mkdir -p ${WORKDIR}/02_Transdecoder
cd ${WORKDIR}/02_Transdecoder

# Step 1: Identify long ORFs
TransDecoder.LongOrfs \
    -t ${WORKDIR}/01_transcriptome_reference/transcriptome.fasta \
    -m ${MIN_PROTEIN} \
    2>&1 | tee transdecoder_longorfs.log

# Count predicted ORFs
ORF_COUNT=$(find . -name "*.transdecoder_dir" -exec cat {}/*.faa \; | grep -c "^>")
echo "Initial ORFs predicted: ${ORF_COUNT}"
```

#### 1.2 (Optional) Homology Search for Improved Prediction

```bash
# Optional: Use BLAST or Pfam to improve ORF prediction
# This step increases accuracy but requires more time

# Option A: BLAST against protein database
blastp -query longest_orfs.pep \
    -db nr \
    -max_target_seqs 1 \
    -outfmt 6 \
    -evalue 1e-5 \
    -num_threads ${THREADS} \
    > blastp.outfmt6

# Option B: Pfam domain search
hmmscan --cpu ${THREADS} \
    --domtblout pfam.domtblout \
    /path/to/Pfam-A.hmm \
    longest_orfs.pep

# Then run TransDecoder.Predict with homology results
TransDecoder.Predict \
    -t ${WORKDIR}/01_transcriptome_reference/transcriptome.fasta \
    --retain_blastp_hits blastp.outfmt6 \
    --retain_pfam_hits pfam.domtblout
```

#### 1.3 Predict Final ORFs

```bash
# Without homology search (faster)
TransDecoder.Predict \
    -t ${WORKDIR}/01_transcriptome_reference/transcriptome.fasta \
    2>&1 | tee transdecoder_predict.log

# Output files:
# - transcriptome.fasta.transdecoder.pep  (protein sequences)
# - transcriptome.fasta.transdecoder.cds  (CDS sequences)
# - transcriptome.fasta.transdecoder.gff3 (coordinates)
# - transcriptome.fasta.transdecoder.bed  (BED format)

# Count final predictions
PEP_COUNT=$(grep -c "^>" transcriptome.fasta.transdecoder.pep)
CDS_COUNT=$(grep -c "^>" transcriptome.fasta.transdecoder.cds)
echo "Final ORFs: ${PEP_COUNT} proteins, ${CDS_COUNT} CDS"
```

---

### Step 2: Select Best ORF per Transcript

Each transcript may have multiple ORFs. Select the best one based on length and score.

#### 2.1 Selection Script

Create `select_best_orf.py`:

```python
#!/usr/bin/env python3
"""
Select the best ORF for each transcript based on length and score.
Best ORF = longest with highest score.
"""

import sys
from collections import defaultdict
from Bio import SeqIO

def parse_gff(gff_file):
    """Parse GFF to get ORF information."""
    orf_info = defaultdict(list)
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            transcript_id = parts[0]
            attributes = parts[8]
            
            # Parse attributes
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key] = value
            
            orf_id = attr_dict.get('ID', '')
            score = float(parts[5]) if parts[5] != '.' else 0.0
            length = int(parts[4]) - int(parts[3]) + 1
            
            orf_info[transcript_id].append({
                'orf_id': orf_id,
                'score': score,
                'length': length
            })
    
    return orf_info

def select_best_orf(orf_list):
    """Select the longest ORF with highest score."""
    if not orf_list:
        return None
    
    # Sort by length (descending), then by score (descending)
    sorted_orfs = sorted(orf_list, 
                        key=lambda x: (x['length'], x['score']), 
                        reverse=True)
    return sorted_orfs[0]['orf_id']

def filter_fasta(input_fa, output_fa, selected_ids):
    """Filter FASTA file to keep only selected ORFs."""
    with open(output_fa, 'w') as out:
        for record in SeqIO.parse(input_fa, 'fasta'):
            if record.id in selected_ids:
                SeqIO.write(record, out, 'fasta')

def main():
    if len(sys.argv) != 4:
        print("Usage: python select_best_orf.py <input.pep> <input.cds> <gff>")
        sys.exit(1)
    
    pep_file = sys.argv[1]
    cds_file = sys.argv[2]
    gff_file = sys.argv[3]
    
    # Parse GFF
    print("Parsing GFF file...")
    orf_info = parse_gff(gff_file)
    
    # Select best ORF for each transcript
    print("Selecting best ORFs...")
    selected_orfs = set()
    for transcript, orfs in orf_info.items():
        best_orf = select_best_orf(orfs)
        if best_orf:
            selected_orfs.add(best_orf)
    
    print(f"Selected {len(selected_orfs)} ORFs")
    
    # Filter PEP file
    print("Filtering PEP file...")
    filter_fasta(pep_file, 'final_pep.fa', selected_orfs)
    
    # Filter CDS file
    print("Filtering CDS file...")
    filter_fasta(cds_file, 'final_cds.fa', selected_orfs)
    
    print("Done! Output files:")
    print("  - final_pep.fa")
    print("  - final_cds.fa")

if __name__ == '__main__':
    main()
```

#### 2.2 Run Selection

```bash
#!/bin/bash
WORKDIR="/path/to/workdir"


# Copy TransDecoder outputs
cp ${WORKDIR}/02_Transdecoder/transcriptome.fasta.transdecoder.* .

# Run selection script
python select_best_orf.py \
    transcriptome.fasta.transdecoder.pep \
    transcriptome.fasta.transdecoder.cds \
    transcriptome.fasta.transdecoder.gff3

# Generate statistics
PEP_SEQS=$(grep -c "^>" final_pep.fa)
CDS_SEQS=$(grep -c "^>" final_cds.fa)
INPUT_SEQS=$(grep -c "^>" ${WORKDIR}/01_transcriptome_reference/transcriptome.fasta)

cat > orf_statistics.txt << EOF
Total transcripts with ORF: ${PEP_SEQS}
Selected best ORFs: ${PEP_SEQS}
Final PEP sequences: ${PEP_SEQS}
Final CDS sequences: ${CDS_SEQS}
ORF length range: $(seqkit stats -a final_cds.fa | awk 'NR==2 {print $6" - "$7}')
ORF average length: $(seqkit stats -a final_cds.fa | awk 'NR==2 {print $8}')
EOF

cat orf_statistics.txt
```

---

## Output Files

### Directory Structure

```
workdir/
── 02_transdecoder/
   ├── transcriptome.fasta.transdecoder.pep
   ├── transcriptome.fasta.transdecoder.cds
   ├── transcriptome.fasta.transdecoder.gff3
   ├── transcriptome.fasta.transdecoder.bed
   ├── transdecoder_longorfs.log
   ├── transdecoder_predict.log
   ├── final_cds.fa ⭐
   ├── final_pep.fa ⭐
   ├── orf_statistics.txt
   └── select_best_orf.py
```

### Final Output Files

| File | Format | Description |
|------|--------|-------------|
| `final_cds.fa` | FASTA | CDS nucleotide sequences (one per transcript) |
| `final_pep.fa` | FASTA | Protein amino acid sequences (one per transcript) |
| `orf_statistics.txt` | Text | Summary statistics of ORF prediction |

---

## Quality Control

### Check Data Quality

```bash
echo "Final ORFs:"
grep -c "^>" ${WORKDIR}/02_Transdecoder/final_pep.fa

# Check sequence lengths
seqkit stats ${WORKDIR}/02_Transdecoder/final_cds.fa
seqkit stats ${WORKDIR}/02_Transdecoder/final_pep.fa

# Verify no duplicate IDs
seqkit rmdup -s ${WORKDIR}/02_Transdecoder/final_pep.fa | \
    awk '{print "Duplicates found: " $2}'
```

### Expected Quality Metrics

| Metric | Expected Value |
|--------|---------------|
| LORDEC correction rate | >95% |
| ORF identification rate | 85-95% |
| Average ORF length | 500-2000 bp |
| Duplicate sequences | 0 |

---

## Troubleshooting

### Common Issues


#### 1. TransDecoder Slow Performance

**Problem**: ORF prediction takes too long

**Solutions**:
- Use homology search only for important transcripts
- Increase thread count
- Use a compute cluster

#### 2. Inconsistent CDS/PEP Counts

**Problem**: CDS and PEP files have different numbers of sequences

**Solutions**:
- Check for sequence duplicates
- Verify selection script correctness
- Re-run from TransDecoder step

---

## Advanced Usage

### Parameter Tuning


#### TransDecoder Parameters

```bash
# Stricter minimum ORF length (200 aa)
TransDecoder.LongOrfs -m 200 ...

# More permissive minimum (50 aa)
TransDecoder.LongOrfs -m 50 ...
```

### Batch Processing

```bash
# Process multiple transcriptome assemblies
for ASSEMBLY in assembly_*.fa; do
    SAMPLE=$(basename ${ASSEMBLY} .fa)
    
    lordec-correct -i ${ASSEMBLY} ... -o ${SAMPLE}_corrected.fa
    TransDecoder.LongOrfs -t ${SAMPLE}_corrected.fa
    TransDecoder.Predict -t ${SAMPLE}_corrected.fa
done
```

---

## References

1. **TransDecoder**: Haas BJ et al. (2013) TransDecoder: Identifying coding regions within transcripts
2. **seqtk**: https://github.com/lh3/seqtk

---

## Version History

- v1.0 (2026-03-11): Initial version with complete workflow
