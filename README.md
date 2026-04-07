# <div style="display: flex; align-items: center;font-family: 'Arial';">Antisense-Discovery Scripts </div>

This repository provides downstream scripts for **antisense candidate discovery** in *Escherichia coli* from modification-aware site tables. In my full workflow, these tables are generated from Nanopore DRS data after dorado basecalling and Modkit pileup processing.
The complete upstream pipeline is still under development, but the scripts in this repository are designed to work on the resulting site tables and can be adapted to other file structures.

## Environment Requirements 
```plaintext
pandas
numpy
biopython
```
Clone the repository and install dependencies:
```plaintext
git clone https://github.com/ChelseaSeqs/Antisense-Discovery.git
cd Antisense-Discovery
pip install -r requirements.txt
```

## 1. **`build_antisense_candidates.py`** 
Builds antisense candidate regions around modified positions and extracts candidate sequences in transcript 5'→3' orientation.

### Inputs:
- a reference FASTA file
- a TSV file containing modified positions
 

The modification table must contain, at minimum, columns for:
- chromosome / reference name
- modified position
- strand
- modification type

Because different workflows may place these fields in different columns, the script lets you specify column indices through command-line arguments.

### Example usage

    python scripts/build_antisense_candidates.py \
      --ref data/fasta.fa \
      --input results/modkit_sites.tsv \
      --output results/antisense_candidates.tsv \
      --chrom-col 0 \
      --pos-col 2 \
      --strand-col 5 \
      --modtype-col 3 \
      --window-size 150 \
      --merge-gap 50

## 2. **`scan_candidate_architecture.py`**  
   Scans each candidate sequence for:
   - promoter-like motifs at the 5' end
   - intrinsic terminator-like features at the 3' end

### Example usage

    python scripts/scan_candidate_architecture.py \
      --input results/antisense_candidates.tsv \
      --output results/antisense_candidate_architecture.tsv


## Outputs
**`build_antisense_candidates.py`** produces a TSV of merged antisense candidate regions with extracted sequences.

Example columns:

candidate_id | chrom | start | end | strand | length | n_modified_sites | modified_positions | mod_type | sequence

**Important notes:**
- candidate regions are built on the opposite strand of the modified sites
- nearby windows are merged into larger candidate intervals
- sequences are returned in transcript 5'→3' orientation for both plus and minus strand candidates

**`scan_candidate_architecture.py`** takes the candidate table from script A and produces a TSV with motif-based architecture annotations, including:
1. Promoter detection results
2. Promoter motif coordinates within the candidate
3. Estimated TSS position in sequence space
4. Estimated TSS genomic coordinate
5. Intrinsic terminator detection results
6. Hairpin and T-run features
7. Estimated TTS position in sequence space
8. Estimated TTS genomic coordinate


## Author

Chelsea Uju Amajirionwu

MSc Molecular & Computational Biology
Focus: RNA-seq, transcriptomics, reproducible workflows, and RNA modification-aware analysis


























