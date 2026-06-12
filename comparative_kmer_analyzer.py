"""
Comparative K-mer Frequency Regulatory Analyzer

Computes and compares normalized k-mer frequencies in the upstream regulatory
regions of target and regulator genes to identify potential transcription
factor binding sites (TFBS).

This tool extracts upstream sequences for two genes, computes k-mer frequency
distributions, normalizes them by sequence length (CPK - Counts per Kilobase),
and identifies motifs enriched in one region versus the other.

License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).

Example Usage:
    # Basic run: Compare two genes with default k=6, show top 20 k-mers
    $ python3 comparative_kmer_analyzer.py -i genome.gbk -t ctg1_50 -r ctg1_74 -o analysis.tsv

    # Custom k-mer size: Use 8-mers instead of 6-mers
    $ python3 comparative_kmer_analyzer.py -i genome.gbk -t ctg1_50 -r ctg1_74 -k 8 -o analysis.tsv

    # Custom upstream windows: Different for target (100bp) vs regulator (200bp)
    $ python3 comparative_kmer_analyzer.py -i genome.gbk -t ctg1_50 -r ctg1_74 --u_target 100 --u_regulator 200 -o analysis.tsv

    # Terminal output: Show top 10 k-mers sorted by CPK difference (no file)
    $ python3 comparative_kmer_analyzer.py -i genome.gbk -t ctg1_50 -r ctg1_74 --top 10

    # All custom: Eukaryotic enhancer analysis with 7-mers, large upstream windows
    $ python3 comparative_kmer_analyzer.py -i yeast_genome.gbff -t YAL001C -r YAL003W -k 7 --u_target 2000 --u_regulator 3000 --top 15 -o enhancer_analysis.tsv
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.1.0"

import sys
import argparse

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

    # args.output is already a Path from base_parser (type=Path); no re-wrapping needed
    output_path = args.output

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

        # The correct denominator for CPK is L - k + 1 (the actual number of
        # sliding windows), NOT the raw sequence length L.
        # Using L inflates the denominator and deflates CPK density.
        # For example, a 150bp sequence with k=6 has 145 windows, not 150.
        total_t_windows = max(1, len(t_seq) - args.kmer + 1)
        total_r_windows = max(1, len(r_seq) - args.kmer + 1)

        def calc_cpk(count: int, total_windows: int) -> float:
            """Counts Per Kilobase, normalized by actual k-mer window count."""
            return (count / total_windows) * 1000

        if output_path:
            with open(output_path, "w", encoding="utf-8") as f:
                f.write(
                    "Kmer\tTarget_Count\tRegulator_Count\t"
                    "Target_CPK\tRegulator_CPK\tCPK_Diff\n"
                )
                for kmer in all_kmers:
                    t_c = t_counts.get(kmer, 0)
                    r_c = r_counts.get(kmer, 0)
                    t_cpk = calc_cpk(t_c, total_t_windows)
                    r_cpk = calc_cpk(r_c, total_r_windows)
                    cpk_diff = abs(t_cpk - r_cpk)
                    f.write(
                        f"{kmer}\t{t_c}\t{r_c}\t"
                        f"{t_cpk:.2f}\t{r_cpk:.2f}\t{cpk_diff:.2f}\n"
                    )
            print(
                f"[*] Success! Analysis written to {output_path.resolve()}",
                file=sys.stderr,
            )

        else:
            # Terminal output: sort by absolute CPK difference to surface the most
            # biologically distinct k-mers. Sorting by raw count is wrong here because
            # target and regulator windows may have different lengths — only CPK-
            # normalised values are fairly comparable across the two sequences.
            print(
                f"[*] Showing top {args.top} k-mers by absolute CPK difference\n",
                file=sys.stderr,
            )
            print(f"{'Kmer':<10} | {'Target CPK':<12} | {'Reg CPK':<12} | {'|Diff|'}")
            print("-" * 52)

            top_kmers = sorted(
                all_kmers,
                key=lambda k: abs(
                    calc_cpk(t_counts.get(k, 0), total_t_windows)
                    - calc_cpk(r_counts.get(k, 0), total_r_windows)
                ),
                reverse=True,
            )[: args.top]

            for kmer in top_kmers:
                t_cpk = calc_cpk(t_counts.get(kmer, 0), total_t_windows)
                r_cpk = calc_cpk(r_counts.get(kmer, 0), total_r_windows)
                diff = abs(t_cpk - r_cpk)
                print(f"{kmer:<10} | {t_cpk:<12.2f} | {r_cpk:<12.2f} | {diff:.2f}")

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
