---
id: prepare_cellranger_from_transdecoder
name: Prepare CellRanger Reference from TransDecoder
description: |
  Convert TransDecoder ORF predictions to CellRanger-compatible reference files.
  Selects best ORF per transcript, generates reference FASTA and GTF files for single-cell analysis.
tags: [cellranger, transdecoder, reference-generation, single-cell, transcriptome]
---

# Prepare CellRanger Reference from TransDecoder

Convert TransDecoder ORF predictions to CellRanger-compatible reference files for single-cell RNA-seq analysis.

---

## Overview

This workflow converts TransDecoder output (from de novo transcriptome assembly) into CellRanger-compatible reference files:
1. **Best ORF Selection** - Select the highest-quality ORF for each transcript
2. **Sequence Extraction** - Extract CDS sequences for selected ORFs
3. **Format Conversion** - Convert GFF3 to GTF format using gffread
4. **ID Simplification** - Use transcript IDs instead of ORF IDs (e.g., `transcript0` instead of `transcript0.p1`)

**Input**:
- TransDecoder CDS file (`.transdecoder.cds`)
- TransDecoder GFF3 file (`.transdecoder.gff3`)
- (Optional) TransDecoder PEP file (`.transdecoder.pep`)

**Output**:
- Reference FASTA (one sequence per transcript)
- Reference GTF (gene annotations in GTF format)
- Statistics file

---

## Environment Setup

### Required Software

```bash
# Core tools
conda install -c bioconda gffread  # For GFF3 to GTF conversion

# Python packages (usually pre-installed)
# - Python 3.7+
# - Biopython (optional, for advanced filtering)
```

### Verify Installation

```bash
# Check gffread
gffread --version  # gffread 0.12.7+

# Check Python
python --version  # Python 3.7+
```

### Environment File Template

Create `environment.yml`:

```yaml
name: cellranger_ref_prep
channels:
  - bioconda
  - conda-forge
dependencies:
  - python=3.11
  - gffread=0.12.7
  - biopython  # Optional
```

Activate:
```bash
conda env create -f environment.yml
conda activate cellranger_ref_prep
```

---

## Workflow Steps

### Step 1: Filter Best ORFs from TransDecoder Output

Each transcript may have multiple ORFs. Select the best one based on:
1. **ORF Type Priority**: complete > 5prime_partial > 3prime_partial > internal
2. **Score**: Higher score = better ORF quality
3. **One representative per transcript**: Each transcript gets only one ORF

#### 1.1 Create Filtering Script

Create `filter_best_orf.py`:

```python
#!/usr/bin/env python3
"""
Filter TransDecoder GFF3 to keep only the best ORF per transcript.
Strategy:
1. Prioritize complete ORFs over partial
2. Among same type, select highest score
3. Generate filtered GFF3 and extract corresponding sequences
"""

import re
from collections import defaultdict

def parse_orf_info(attributes):
    """Extract ORF type and score from GFF3 attributes."""
    type_match = re.search(r'ORF type:([^\s,]+)', attributes)
    orf_type = type_match.group(1) if type_match else "unknown"
    
    score_match = re.search(r'score[=:](\d+\.?\d*)', attributes)
    score = float(score_match.group(1)) if score_match else 0.0
    
    return orf_type, score

def orf_priority(orf_type):
    """Assign priority to ORF types. Lower = higher priority."""
    if "complete" in orf_type:
        return 1
    elif "5prime_partial" in orf_type or "3prime_partial" in orf_type:
        return 2
    elif "internal" in orf_type:
        return 3
    return 4

def main():
    import sys
    
    if len(sys.argv) != 5:
        print("Usage: python filter_best_orf.py <input.gff3> <input.cds> <input.pep> <output_dir>")
        sys.exit(1)
    
    gff3_file = sys.argv[1]
    cds_file = sys.argv[2]
    pep_file = sys.argv[3]
    output_dir = sys.argv[4]
    
    from pathlib import Path
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    filtered_gff3 = output_dir / "filtered_orfs.gff3"
    filtered_cds = output_dir / "genes.fasta"
    filtered_pep = output_dir / "proteins.fasta"
    stats_file = output_dir / "filtering_stats.txt"
    
    print("Step 1: Parsing GFF3 file...")
    transcripts = defaultdict(list)
    
    with open(gff3_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'gene':
                continue
            
            seqid = parts[0]
            attributes = parts[8]
            
            orf_id_match = re.search(r'ID=GENE\.([^;]+)', attributes)
            if not orf_id_match:
                continue
            
            orf_id = orf_id_match.group(1)
            orf_type_str, score_val = parse_orf_info(attributes)
            score_val = float(score_val)
            
            transcripts[seqid].append({
                'line': line,
                'orf_id': orf_id,
                'orf_type': orf_type_str,
                'score': score_val
            })
    
    print(f"  Found {len(transcripts)} transcripts")
    
    print("\nStep 2: Selecting best ORF per transcript...")
    selected_orfs = set()
    orf_type_counts = defaultdict(int)
    
    for transcript_id, orf_list in transcripts.items():
        # Sort by priority then by score (descending)
        for orf in orf_list:
            orf['score'] = float(orf['score'])
        
        best_orf = sorted(orf_list, 
                         key=lambda x: (orf_priority(x['orf_type']), -float(x['score'])))[0]
        selected_orfs.add(best_orf['orf_id'])
        orf_type_counts[best_orf['orf_type']] += 1
    
    print(f"  Selected {len(selected_orfs)} best ORFs")
    print(f"  ORF type distribution:")
    for orf_type, count in sorted(orf_type_counts.items(), key=lambda x: -x[1]):
        print(f"    {orf_type}: {count}")
    
    print("\nStep 3: Extracting filtered GFF3 entries...")
    output_lines = []
    keep_current = False
    
    with open(gff3_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            feature_type = parts[2]
            attributes = parts[8]
            
            if feature_type == 'gene':
                orf_id_match = re.search(r'ID=GENE\.([^;]+)', attributes)
                if orf_id_match:
                    current_orf_id = orf_id_match.group(1)
                    keep_current = current_orf_id in selected_orfs
                    if keep_current:
                        output_lines.append(line)
            else:
                if keep_current:
                    output_lines.append(line)
    
    with open(filtered_gff3, 'w') as f:
        f.writelines(output_lines)
    print(f"  Written {len(output_lines)} lines to {filtered_gff3}")
    
    print("\nStep 4: Filtering CDS sequences...")
    with open(cds_file, 'r') as f_in, open(filtered_cds, 'w') as f_out:
        keep = False
        for line in f_in:
            if line.startswith('>'):
                orf_id = line[1:].split()[0]
                keep = orf_id in selected_orfs
                if keep:
                    f_out.write(line)
            elif keep:
                f_out.write(line)
    
    print("\nStep 5: Filtering protein sequences...")
    with open(pep_file, 'r') as f_in, open(filtered_pep, 'w') as f_out:
        keep = False
        for line in f_in:
            if line.startswith('>'):
                orf_id = line[1:].split()[0]
                keep = orf_id in selected_orfs
                if keep:
                    f_out.write(line)
            elif keep:
                f_out.write(line)
    
    # Write statistics
    with open(stats_file, 'w') as f:
        f.write("TransDecoder Best ORF Filtering Statistics\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total transcripts: {len(transcripts)}\n")
        f.write(f"Selected ORFs: {len(selected_orfs)}\n\n")
        f.write("ORF type distribution:\n")
        for orf_type, count in sorted(orf_type_counts.items(), key=lambda x: -x[1]):
            f.write(f"  {orf_type}: {count}\n")
    
    print(f"\n✅ Filtering complete!")
    print(f"\nOutput files:")
    print(f"  GFF3: {filtered_gff3}")
    print(f"  CDS:  {filtered_cds}")
    print(f"  PEP:  {filtered_pep}")
    print(f"  Stats: {stats_file}")

if __name__ == "__main__":
    main()
