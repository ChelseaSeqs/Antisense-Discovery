import argparse
from pathlib import Path
import pandas as pd

# -----------------------------------------------------------------------------------
# This script is part of an ongoing project but can be modified for stand-alone use
# Author: Chelsea Uju Amajirionwu
# Created Date: 04/2026
# ------------------------------------------------------------------------------------


def hamming_distance(s1, s2):
    """Count mismatches between known promoter motif and candidate promoter"""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def scan_promoter(
    seq,
    promoter_window=120,
    minus10_consensus="TATAAT",
    minus35_consensus="TTGACA",
    max_minus10_mismatches=2,
    max_minus35_mismatches=2,
    spacing_range=(15, 19),
):
    """
    Scan the 5' end of the candidate sequence for a sigma70-like promoter.
    Sequence must already be in transcript 5'->3' orientation.
    Returns the best promoter hit, or a null-result dictionary.
    """
    scan_seq = seq[:promoter_window].upper()
    best_hit = None

    minus10_hits = []
    for i in range(len(scan_seq) - len(minus10_consensus) + 1):
        kmer = scan_seq[i : i + len(minus10_consensus)]
        mm = hamming_distance(kmer, minus10_consensus)

        if mm <= max_minus10_mismatches:
            minus10_hits.append(
                {
                    "minus10_seq": kmer,
                    "minus10_start": i,
                    "minus10_end": i + len(minus10_consensus) - 1,
                    "minus10_mismatches": mm,
                }
            )

    for m10 in minus10_hits:
        m10_start = m10["minus10_start"]

        for spacing in range(spacing_range[0], spacing_range[1] + 1):
            m35_start = m10_start - spacing - len(minus35_consensus)
            m35_end = m35_start + len(minus35_consensus)

            if m35_start < 0:
                continue

            kmer35 = scan_seq[m35_start:m35_end]
            mm35 = hamming_distance(kmer35, minus35_consensus)

            if mm35 <= max_minus35_mismatches:
                estimated_tss = m10["minus10_end"] + 1 + 7
                score = (
                    (len(minus10_consensus) - m10["minus10_mismatches"])
                    + (len(minus35_consensus) - mm35)
                )

                hit = {
                    "promoter_found": True,
                    "promoter_score": score,
                    "minus35_seq": kmer35,
                    "minus35_start": m35_start,
                    "minus35_end": m35_end - 1,
                    "minus35_mismatches": mm35,
                    "minus10_seq": m10["minus10_seq"],
                    "minus10_start": m10["minus10_start"],
                    "minus10_end": m10["minus10_end"],
                    "minus10_mismatches": m10["minus10_mismatches"],
                    "spacing": spacing,
                    "estimated_tss_pos_in_seq": estimated_tss,
                }

                if best_hit is None or hit["promoter_score"] > best_hit["promoter_score"]:
                    best_hit = hit

    if best_hit is None:
        return {
            "promoter_found": False,
            "promoter_score": None,
            "minus35_seq": None,
            "minus35_start": None,
            "minus35_end": None,
            "minus35_mismatches": None,
            "minus10_seq": None,
            "minus10_start": None,
            "minus10_end": None,
            "minus10_mismatches": None,
            "spacing": None,
            "estimated_tss_pos_in_seq": None,
        }

    return best_hit


def can_pair(a, b):
    """Simple Watson-Crick DNA pairing."""
    pairs = {("A", "T"), ("T", "A"), ("G", "C"), ("C", "G")}
    return (a, b) in pairs


def stem_matches(left, right_reversed):
    """Count pairable positions between left stem and reverse-complemented right stem."""
    return sum(can_pair(a, b) for a, b in zip(left, right_reversed))


def scan_intrinsic_terminator(
    seq,
    terminator_window=120,
    min_t_run=4,
    max_t_run_search_len=10,
    stem_len_range=(5, 9),
    loop_len_range=(3, 8),
    upstream_search_before_t=50,
):
    """
    Scan the 3' end of the candidate sequence for a simple intrinsic terminator:
    hairpin-like inverted repeat + T-run.
    Returns the best terminator hit, or a null-result dictionary.
    """
    scan_start = max(0, len(seq) - terminator_window)
    scan_seq = seq[scan_start:].upper()
    best_hit = None

    for i in range(len(scan_seq)):
        for run_len in range(max_t_run_search_len, min_t_run - 1, -1):
            if i + run_len > len(scan_seq):
                continue

            fragment = scan_seq[i : i + run_len]
            if set(fragment) != {"T"}:
                continue

            t_run_start = i
            t_run_end = i + run_len - 1

            region_end = t_run_start
            region_start = max(0, region_end - upstream_search_before_t)
            upstream_region = scan_seq[region_start:region_end]

            for stem_len in range(stem_len_range[0], stem_len_range[1] + 1):
                for loop_len in range(loop_len_range[0], loop_len_range[1] + 1):
                    total_len = stem_len + loop_len + stem_len
                    if total_len > len(upstream_region):
                        continue

                    for j in range(len(upstream_region) - total_len + 1):
                        left = upstream_region[j : j + stem_len]
                        loop = upstream_region[j + stem_len : j + stem_len + loop_len]
                        right = upstream_region[j + stem_len + loop_len : j + total_len]

                        pair_score = stem_matches(left, right[::-1])

                        if pair_score >= stem_len - 1:
                            score = pair_score + run_len

                            hit = {
                                "terminator_found": True,
                                "terminator_score": score,
                                "hairpin_left_seq": left,
                                "hairpin_loop_seq": loop,
                                "hairpin_right_seq": right,
                                "hairpin_start": scan_start + region_start + j,
                                "hairpin_end": scan_start + region_start + j + total_len - 1,
                                "stem_len": stem_len,
                                "loop_len": loop_len,
                                "stem_pair_score": pair_score,
                                "t_run_seq": fragment,
                                "t_run_start": scan_start + t_run_start,
                                "t_run_end": scan_start + t_run_end,
                                "estimated_tts_pos_in_seq": scan_start + t_run_end,
                            }

                            if best_hit is None or hit["terminator_score"] > best_hit["terminator_score"]:
                                best_hit = hit

            break

    if best_hit is None:
        return {
            "terminator_found": False,
            "terminator_score": None,
            "hairpin_left_seq": None,
            "hairpin_loop_seq": None,
            "hairpin_right_seq": None,
            "hairpin_start": None,
            "hairpin_end": None,
            "stem_len": None,
            "loop_len": None,
            "stem_pair_score": None,
            "t_run_seq": None,
            "t_run_start": None,
            "t_run_end": None,
            "estimated_tts_pos_in_seq": None,
        }

    return best_hit


def convert_seq_to_genomic(seq_pos, row_start, row_end, strand):
    """
    Convert a 0-based sequence-relative position to genomic coordinate.
    Candidate sequence is assumed to be in transcript 5'->3' orientation.
    """
    if seq_pos is None:
        return None

    if strand == "+":
        return row_start + seq_pos
    if strand == "-":
        return row_end - seq_pos

    return None


def scan_candidate_architecture(
    candidates_df,
    promoter_window=120,
    minus10_consensus="TATAAT",
    minus35_consensus="TTGACA",
    max_minus10_mismatches=2,
    max_minus35_mismatches=2,
    spacing_range=(15, 19),
    terminator_window=120,
    min_t_run=4,
    max_t_run_search_len=10,
    stem_len_range=(5, 9),
    loop_len_range=(3, 8),
    upstream_search_before_t=50,
):
    """
    Scan candidate sequences for promoter-like motifs and intrinsic terminators.
    Input dataframe must contain:
    candidate_id, chrom, start, end, strand, sequence
    """
    required_cols = {"candidate_id", "chrom", "start", "end", "strand", "sequence"}
    missing = required_cols - set(candidates_df.columns)
    if missing:
        raise ValueError(
            "Input candidate table is missing required columns: "
            + ", ".join(sorted(missing))
        )

    results = []

    for _, row in candidates_df.iterrows():
        seq = str(row["sequence"]).upper()

        promoter = scan_promoter(
            seq=seq,
            promoter_window=promoter_window,
            minus10_consensus=minus10_consensus,
            minus35_consensus=minus35_consensus,
            max_minus10_mismatches=max_minus10_mismatches,
            max_minus35_mismatches=max_minus35_mismatches,
            spacing_range=spacing_range,
        )

        terminator = scan_intrinsic_terminator(
            seq=seq,
            terminator_window=terminator_window,
            min_t_run=min_t_run,
            max_t_run_search_len=max_t_run_search_len,
            stem_len_range=stem_len_range,
            loop_len_range=loop_len_range,
            upstream_search_before_t=upstream_search_before_t,
        )

        out = {
            "candidate_id": row["candidate_id"],
            "chrom": row["chrom"],
            "start": row["start"],
            "end": row["end"],
            "strand": row["strand"],
            "length": len(seq),
        }

        if "n_modified_sites" in row.index:
            out["n_modified_sites"] = row["n_modified_sites"]
        if "modified_positions" in row.index:
            out["modified_positions"] = row["modified_positions"]
        if "mod_type" in row.index:
            out["mod_type"] = row["mod_type"]

        out.update(promoter)
        out.update(terminator)

        out["estimated_tss_genomic"] = convert_seq_to_genomic(
            promoter["estimated_tss_pos_in_seq"],
            row["start"],
            row["end"],
            row["strand"],
        )

        out["estimated_tts_genomic"] = convert_seq_to_genomic(
            terminator["estimated_tts_pos_in_seq"],
            row["start"],
            row["end"],
            row["strand"],
        )

        results.append(out)

    return pd.DataFrame(results)


def parse_args():
    parser = argparse.ArgumentParser(description="Scan antisense candidate sequences for promoter-like motifs and intrinsic terminators.")
    parser.add_argument('-i', "-input", required=True, help="Input TSV file from build_antisense_candidates.py",)
    parser.add_argument('-o', "--output", required=True, help="Output TSV file with promoter/terminator scan results.",)
    parser.add_argument("--promoter-window", type=int, default=120, help="Number of nt from the 5' end to scan for promoter motifs.",)
    parser.add_argument("--minus10-consensus", default="TATAAT", help="Consensus sequence for the -10 box.",)
    parser.add_argument("--minus35-consensus", default="TTGACA", help="Consensus sequence for the -35 box.",)
    parser.add_argument("--max-minus10-mismatches", type=int, default=2, help="Maximum mismatches allowed for -10 box matches.",)
    parser.add_argument("--max-minus35-mismatches", type=int, default=2, help="Maximum mismatches allowed for -35 box matches.",)
    parser.add_argument("--spacing-min", type=int, default=15, help="Minimum allowed spacing between -35 and -10 boxes.",)
    parser.add_argument("--spacing-max", type=int, default=19, help="Maximum allowed spacing between -35 and -10 boxes.",)
    parser.add_argument("--terminator-window", type=int, default=120, help="Number of nt from the 3' end to scan for intrinsic terminators.",)
    parser.add_argument("--min-t-run", type=int, default=4, help="Minimum poly-T run length for terminator detection.",)
    parser.add_argument("--max-t-run-search-len", type=int, default=10, help="Maximum T-run length to consider.",)
    parser.add_argument("--stem-min", type=int, default=5, help="Minimum stem length for hairpin detection.",)
    parser.add_argument("--stem-max", type=int, default=9, help="Maximum stem length for hairpin detection.",)
    parser.add_argument("--loop-min", type=int, default=3, help="Minimum loop length for hairpin detection.",)
    parser.add_argument("--loop-max", type=int, default=8, help="Maximum loop length for hairpin detection.",)
    parser.add_argument("--upstream-search-before-t", type=int, default=50, help="How far upstream of the T-run to search for a hairpin.",)
    return parser.parse_args()


def main():
    args = parse_args()

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print("Loading antisense candidate table...")
    candidates = pd.read_csv(args.input, sep="\t")

    print("Scanning candidate architecture...")
    scan_results = scan_candidate_architecture(
        candidates_df=candidates,
        promoter_window=args.promoter_window,
        minus10_consensus=args.minus10_consensus,
        minus35_consensus=args.minus35_consensus,
        max_minus10_mismatches=args.max_minus10_mismatches,
        max_minus35_mismatches=args.max_minus35_mismatches,
        spacing_range=(args.spacing_min, args.spacing_max),
        terminator_window=args.terminator_window,
        min_t_run=args.min_t_run,
        max_t_run_search_len=args.max_t_run_search_len,
        stem_len_range=(args.stem_min, args.stem_max),
        loop_len_range=(args.loop_min, args.loop_max),
        upstream_search_before_t=args.upstream_search_before_t,
    )

    scan_results.to_csv(output_path, sep="\t", index=False)
    print(f"Saved scan results for {len(scan_results)} candidates to {output_path}")


if __name__ == "__main__":
    main()
