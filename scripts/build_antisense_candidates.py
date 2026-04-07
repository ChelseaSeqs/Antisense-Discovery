import argparse
from pathlib import Path
import pandas as pd
from Bio import SeqIO

# -----------------------------------------------------------------------------------
# This script is part of an ongoing project but can be modified for stand-alone use
# Author: Chelsea Uju Amajirionwu
# Created Date: 04/2026
# ------------------------------------------------------------------------------------


def load_genome(fasta_file):
    """Load reference sequences into a dictionary and store chromosome lengths."""
    genome = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    chrom_lengths = {chrom: len(record.seq) for chrom, record in genome.items()}
    return genome, chrom_lengths


def load_modifications(tsv_file, chrom_col, pos_col, strand_col, modtype_col):
    """Load modification table and keep required columns."""
    df = pd.read_csv(tsv_file, sep="\t", header=None)

    mods = df[[chrom_col, pos_col, strand_col, modtype_col]].copy()
    mods.columns = ["chrom", "mod", "mod_strand", "mod_type"]

    mods["mod"] = mods["mod"].astype(int)
    mods["mod_strand"] = mods["mod_strand"].astype(str).str.strip()
    mods["mod_type"] = mods["mod_type"].astype(str).str.strip()

    valid_strands = {"+", "-"}
    mods = mods[mods["mod_strand"].isin(valid_strands)].copy()

    return mods


def build_candidate_intervals(mods, chrom_lengths, window_size, merge_gap):
    """
    Build antisense candidate regions around modified sites and merge nearby intervals.
    """
    mods = mods.copy()

    mods["antisense_strand"] = mods["mod_strand"].map({"+": "-", "-": "+"})
    mods["start"] = mods["mod"] - window_size
    mods["end"] = mods["mod"] + window_size

    mods["start"] = mods.apply(lambda x: max(1, x["start"]), axis=1)
    mods["end"] = mods.apply(
        lambda x: min(chrom_lengths[x["chrom"]], x["end"]), axis=1
    )

    merged_candidates = []

    for (chrom, strand), sub in mods.groupby(["chrom", "antisense_strand"]):
        sub = sub.sort_values("start").reset_index(drop=True)
        if sub.empty:
            continue

        current_start = sub.loc[0, "start"]
        current_end = sub.loc[0, "end"]
        current_sites = [sub.loc[0, "mod"]]
        current_mods = [sub.loc[0, "mod_type"]]

        for i in range(1, len(sub)):
            row = sub.loc[i]

            if row["start"] <= current_end + merge_gap:
                current_end = max(current_end, row["end"])
                current_sites.append(row["mod"])
                current_mods.append(row["mod_type"])
            else:
                merged_candidates.append(
                    {
                        "chrom": chrom,
                        "strand": strand,
                        "start": current_start,
                        "end": current_end,
                        "n_modified_sites": len(current_sites),
                        "modified_positions": ",".join(
                            map(str, sorted(set(current_sites)))
                        ),
                        "mod_type": ",".join(map(str, sorted(set(current_mods)))),
                    }
                )

                current_start = row["start"]
                current_end = row["end"]
                current_sites = [row["mod"]]
                current_mods = [row["mod_type"]]

        merged_candidates.append(
            {
                "chrom": chrom,
                "strand": strand,
                "start": current_start,
                "end": current_end,
                "n_modified_sites": len(current_sites),
                "modified_positions": ",".join(map(str, sorted(set(current_sites)))),
                "mod_type": ",".join(map(str, sorted(set(current_mods)))),
            }
        )

    return pd.DataFrame(merged_candidates)


def extract_transcript_sequence(genome, chrom, start, end, strand):
    """
    Extract sequence in transcript 5'->3' orientation.
    Coordinates are assumed to be 1-based inclusive.
    """
    seq = genome[chrom].seq[start - 1 : end]
    if strand == "-":
        seq = seq.reverse_complement()
    return str(seq)


def add_sequences(candidates, genome):
    """Add sequence, length, and candidate IDs."""
    candidates = candidates.copy()

    candidates["length"] = candidates["end"] - candidates["start"] + 1
    candidates["sequence"] = candidates.apply(
        lambda x: extract_transcript_sequence(
            genome, x["chrom"], x["start"], x["end"], x["strand"]
        ),
        axis=1,
    )

    candidates = candidates.reset_index(drop=True)
    candidates["candidate_id"] = [
        f"antisense_candidate_{i+1}" for i in range(len(candidates))
    ]

    candidates = candidates[
        [
            "candidate_id",
            "chrom",
            "start",
            "end",
            "strand",
            "length",
            "n_modified_sites",
            "modified_positions",
            "mod_type",
            "sequence",
        ]
    ]

    return candidates


def parse_args():
    parser = argparse.ArgumentParser(description="Build antisense candidate regions around modified positions.")
    parser.add_argument('-r', "--ref", required=True, help="Reference genome/transcriptome FASTA file.",)
    parser.add_argument('-i', "-input", required=True, help="Input TSV file containing modified positions.",)
    parser.add_argument('-o', "--output", required=True, help="Output TSV file for antisense candidates.",)
    parser.add_argument("--chrom-col", type=int, default=0, help="Column index for chromosome/reference name in the input TSV.",)
    parser.add_argument("--pos-col", type=int, default=2, help="Column index for modification position in the input TSV.",)
    parser.add_argument("--strand-col", type=int, default=5, help="Column index for strand in the input TSV.",)
    parser.add_argument("--modtype-col", type=int, default=3, help="Column index for modification type in the input TSV.",)
    parser.add_argument("--window-size", type=int, default=150, help="Number of nt to extend upstream and downstream of each modified site.",)
    parser.add_argument("--merge-gap", type=int, default=50, help="Merge candidate intervals if they overlap or are within this many nt.",)
    return parser.parse_args()


def main():
    args = parse_args()

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print("Loading genome...")
    genome, chrom_lengths = load_genome(args.fasta)

    print("Loading modification table...")
    mods = load_modifications(
        args.input,
        chrom_col=args.chrom_col,
        pos_col=args.pos_col,
        strand_col=args.strand_col,
        modtype_col=args.modtype_col,
    )

    if mods.empty:
        raise ValueError("No valid modification records found after input parsing.")

    missing_chroms = sorted(set(mods["chrom"]) - set(chrom_lengths.keys()))
    if missing_chroms:
        raise ValueError(
            "The following references were found in the TSV but not in the FASTA: "
            + ", ".join(missing_chroms[:10])
        )

    print("Building merged antisense candidate intervals...")
    candidates = build_candidate_intervals(
        mods,
        chrom_lengths=chrom_lengths,
        window_size=args.window_size,
        merge_gap=args.merge_gap,
    )

    if candidates.empty:
        print("No antisense candidates were generated.")
        candidates.to_csv(output_path, sep="\t", index=False)
        return

    print("Extracting sequences...")
    candidates = add_sequences(candidates, genome)

    candidates.to_csv(output_path, sep="\t", index=False)
    print(f"Saved {len(candidates)} antisense transcript candidates to {output_path}")
    print(candidates.head())


if __name__ == "__main__":
    main()
