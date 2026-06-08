"""
Comparative K-mer Frequency Regulatory Analyzer

Computes and compares normalized k-mer frequencies in the upstream regulatory
regions of target and regulator genes to identify potential transcription
factor binding sites (TFBS).

This tool extracts upstream sequences for two genes, computes k-mer frequency
distributions, normalizes them by sequence length (CPK - Counts per Kilobase),
and identifies motifs enriched in one region versus the other.

Author: Jan Ephraim R. Vallente (ephrvallente@gmail.com)
Date: 2026-06-05
License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).

Example Usage:
    $ python3 kmer_analyzer.py -i genome.gbk -t ctg1_50 -r ctg1_74 -o analysis.tsv
"""

__version__ = "1.0.1"

import sys
import argparse
from pathlib import Path
from collections import Counter
from utils import base_parser, extract_upstream_sequence


def get_kmer_counts(sequence: str, k: int) -> Counter:
    """
    Returns the frequency count of k-mers in a sequence.

    Args:
        sequence: A nucleotide or amino acid string.
        k: The length of each k-mer. Must be >= 1.

    Returns:
        A Counter object mapping k-mer strings to their integer counts.

    Raises:
        ValueError: If k < 1 or if the sequence is shorter than k.
    """
    if k < 1:
        raise ValueError(f"k must be at least 1, got {k}.")
    if len(sequence) < k:
        raise ValueError(
            f"Sequence length ({len(sequence)} bp) is shorter than k ({k}). "
            "Reduce k or increase the upstream window."
        )

    seq = sequence.upper()
    kmers = [seq[i : i + k] for i in range(len(seq) - k + 1)]
    return Counter(kmers)


def get_args() -> argparse.Namespace:
    """Configures the CLI and returns parsed arguments."""
    parser = base_parser("Comparative K-mer Frequency Analyzer for Regulatory Regions")
    parser.add_argument(
        "-t", "--target", required=True, help="Locus tag for target gene"
    )
    parser.add_argument(
        "-r", "--regulator", required=True, help="Locus tag for regulator"
    )
    parser.add_argument(
        "--u_target", type=int, default=150, help="Upstream bp for target"
    )
    parser.add_argument(
        "--u_regulator", type=int, default=300, help="Upstream bp for regulator"
    )
    parser.add_argument("-k", "--kmer", type=int, default=6, help="K-mer length")
    parser.add_argument("--top", type=int, default=20, help="Top N k-mers to report")
    return parser.parse_args()


def main() -> None:
    args = get_args()

    # Validation
    if args.kmer < 1:
        sys.exit("[!] --kmer must be at least 1.")
    if args.u_target < 1 or args.u_regulator < 1:
        sys.exit("[!] Upstream values must be positive integers.")

    output_path = Path(args.output) if args.output else None

    try:
        # Extract sequences
        t_seq, _, _, _ = extract_upstream_sequence(
            args.input, args.target, args.u_target
        )
        r_seq, _, _, _ = extract_upstream_sequence(
            args.input, args.regulator, args.u_regulator
        )

        # Count k-mers
        t_counts = get_kmer_counts(t_seq, args.kmer)
        r_counts = get_kmer_counts(r_seq, args.kmer)

        all_kmers = sorted(set(t_counts) | set(r_counts))

        # Output preparation
        if output_path:
            with open(output_path, "w", encoding="utf-8") as f:
                f.write(
                    "Kmer\tTarget_Count\tRegulator_Count\tTarget_CPK\tRegulator_CPK\n"
                )
                for kmer in all_kmers:
                    t_c = t_counts[kmer]
                    r_c = r_counts[kmer]
                    t_cpk = (t_c / len(t_seq)) * 1000
                    r_cpk = (r_c / len(r_seq)) * 1000
                    f.write(f"{kmer}\t{t_c}\t{r_c}\t{t_cpk:.2f}\t{r_cpk:.2f}\n")
            print(f"[*] Success! Analysis written to {output_path.name}")

        else:
            print(f"{'Kmer':<10} | {'Target':<10} | {'Reg':<10}")
            print("-" * 35)
            # Just print the top N sorted by Target frequency for terminal
            top_kmers = sorted(t_counts.items(), key=lambda x: x[1], reverse=True)[
                : args.top
            ]
            for kmer, count in top_kmers:
                print(f"{kmer:<10} | {count:<10} | {r_counts[kmer]:<10}")

    except (FileNotFoundError, ValueError) as e:
        # Cleanup partial file on error
        if output_path and output_path.exists():
            output_path.unlink()
        sys.exit(f"\n[!] Pipeline Error: {e}")
    except KeyboardInterrupt:
        if output_path and output_path.exists():
            output_path.unlink()
        sys.exit("\n[!] Pipeline interrupted by user.")


if __name__ == "__main__":
    main()
