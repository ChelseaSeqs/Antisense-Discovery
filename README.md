# <div style="display: flex; align-items: center;font-family: 'Arial';">Antisense-Discovery Scripts </div>

This repository provides downstream scripts for antisense candidate discovery from modification-aware site tables. In my full workflow, these tables are generated from Nanopore DRS data after dorado basecalling and Modkit pileup processing.

## Environment Requirements 
```plaintext
pandas
numpy
biopython
```

## 1. Build antisense candidates

Example usage

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

## 2. Scan candidates

Example usage

    python scripts/scan_candidate_architecture.py \
      --input results/antisense_candidates.tsv \
      --output results/antisense_candidate_architecture.tsv
