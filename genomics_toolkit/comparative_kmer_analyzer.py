#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

"""Comparative K-mer Frequency Regulatory Analyzer

Computes and compares normalized k-mer frequencies in the upstream regulatory
regions of target and regulator genes to identify potential transcription
factor binding sites (TFBS).

This tool extracts upstream sequences for two genes, computes k-mer frequency
distributions, normalizes them by sequence length (CPK - Counts per Kilobase),
and identifies motifs enriched in one region versus the other.

CANONICAL K-MER ANALYSIS (strand-aware):
    Transcription factors bind double-stranded DNA and do not distinguish which
    strand the genome annotator labeled the coding strand. A TF binding GATA
    on the template strand appears as TATC in your extracted sequence. Without
    canonical k-mers, you count GATA and TATC as separate entities and see
    each at half the real frequency.

    This script uses canonical k-mers: for each k-mer extracted from the
    sequence, it computes the reverse complement and keeps whichever is
    lexicographically smaller. The counts of a k-mer and its reverse complement
    are merged into a single canonical count. This correctly represents TF
    binding affinity regardless of strand orientation.

    Example: GATA and TATC (its reverse complement) both become GATA (if GATA
    <= TATC lexicographically), so all occurrences on either strand are counted
    together.

    Note: upstream sequences produced by universal_promoter_extractor.py are
    already strand-corrected (5'→3' relative to the gene). Canonical k-mers
    additionally handle palindromic TF binding sites and sequences supplied
    from external sources.

ENRICHMENT METRIC — LOG2 FOLD CHANGE (L2FC):
    Raw CPK difference is a misleading enrichment metric. Consider:
        K-mer A: Target CPK = 1010, Regulator CPK = 1000  |diff| = 10
        K-mer B: Target CPK = 10,   Regulator CPK = 0     |diff| = 10
    Raw difference ranks these equally, but K-mer B is infinitely enriched
    while K-mer A is background noise.

    This script computes Log2 Fold Change (L2FC) using Haldane-Anscombe
    pseudo-count correction (+0.5 to counts, +1 to window total) to prevent
    log(0) and stabilize estimates for rare k-mers:

        t_freq = (t_count + 0.5) / (t_windows + 1)
        r_freq = (r_count + 0.5) / (r_windows + 1)
        L2FC   = log2(t_freq / r_freq)

    Positive L2FC = enriched in target; negative = enriched in regulator.
    Terminal output sorts by |L2FC| to surface the most biologically
    distinct k-mers. CPK values are retained in the TSV for reference.

Note:
    Associated with ongoing, unpublished research (manuscript in
    preparation). Correct attribution is requested when used in
    derivative works.

Examples:
    # Basic run: Compare two genes with default k=6, show top 20 k-mers
    $ python3 comparative_kmer_analyzer.py -i genome.gbk -t ctg1_50 -r ctg1_74 -o analysis.tsv

    # Custom k-mer size: Use 8-mers instead of 6-mers
    $ python3 comparative_kmer_analyzer.py -i genome.gbk -t ctg1_50 -r ctg1_74 -k 8 -o analysis.tsv

    # Custom upstream windows: Different for target (100bp) vs regulator (200bp)
    $ python3 comparative_kmer_analyzer.py -i genome.gbk -t ctg1_50 -r ctg1_74 --u_target 100 --u_regulator 200 -o analysis.tsv

    # Terminal output: Show top 10 k-mers sorted by |L2FC| (no file)
    $ python3 comparative_kmer_analyzer.py -i genome.gbk -t ctg1_50 -r ctg1_74 --top 10

    # All custom: Eukaryotic enhancer analysis with 7-mers, large upstream windows
    $ python3 comparative_kmer_analyzer.py -i yeast_genome.gbff -t YAL001C -r YAL003W -k 7 --u_target 2000 --u_regulator 3000 --top 15 -o enhancer_analysis.tsv
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.2.0"

import math
import sys
import argparse

from collections import Counter
from utils import base_parser, extract_upstream_sequence

# ── Canonical k-mer helpers ───────────────────────────────────────────────────

_RC_TABLE = str.maketrans("ACGT", "TGCA")


def _revcomp(seq: str) -> str:
    """Return the reverse complement of an uppercase DNA string."""
    return seq.translate(_RC_TABLE)[::-1]


def _canonical(kmer: str) -> str:
    """Return the canonical form of a k-mer (lexicographically smaller of kmer/revcomp).

    Canonical k-mers merge the counts of a k-mer and its reverse complement
    into a single entity, correctly representing TF binding affinity regardless
    of which DNA strand the site appears on.

    Args:
        kmer: Uppercase DNA k-mer string (ACGT only).

    Returns:
        The k-mer or its reverse complement, whichever sorts first.
    """
    rc = _revcomp(kmer)
    return kmer if kmer <= rc else rc


# ── Enrichment metric ─────────────────────────────────────────────────────────


def calc_l2fc(t_count: int, r_count: int, t_windows: int, r_windows: int) -> float:
    """Log2 fold change of a k-mer's frequency in target vs regulator.

    Uses Haldane-Anscombe pseudo-count correction (+0.5 to counts, +1 to
    window totals) to handle zero-count k-mers without log(0) and to
    stabilize fold-change estimates for rare k-mers.

    A positive value indicates enrichment in the target; negative indicates
    enrichment in the regulator.

    This is superior to raw CPK difference, which conflates effect size with
    absolute frequency: two k-mers with |CPK diff| = 10 could represent
    1010 vs 1000 (noise) or 10 vs 0 (infinite enrichment). L2FC correctly
    ranks the second case far above the first.

    Args:
        t_count:   Raw k-mer count in target sequence.
        r_count:   Raw k-mer count in regulator sequence.
        t_windows: Total k-mer windows in target (sequence_length - k + 1).
        r_windows: Total k-mer windows in regulator.

    Returns:
        Log2 fold change (positive = enriched in target).
    """
    t_freq = (t_count + 0.5) / (t_windows + 1)
    r_freq = (r_count + 0.5) / (r_windows + 1)
    return math.log2(t_freq / r_freq)


# ── K-mer counting ────────────────────────────────────────────────────────────


def get_kmer_counts(sequence: str, k: int) -> Counter:
    """
    Returns the canonical k-mer frequency count for a sequence.

    Canonical k-mers merge the count of each k-mer with its reverse complement
    into a single entity (the lexicographically smaller of the two). This
    correctly captures TF binding sites regardless of which strand they appear
    on, because a TF binding GATA on the template strand reads as TATC on the
    coding strand — without canonicalization, you see each at half frequency.

    Args:
        sequence: A nucleotide string.
        k: The length of each k-mer. Must be >= 1.

    Returns:
        A Counter mapping canonical k-mer strings to their integer counts.
        Each count reflects occurrences on BOTH strands combined.

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
    kmers = [_canonical(seq[i : i + k]) for i in range(len(seq) - k + 1)]
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
                    "Target_CPK\tRegulator_CPK\tCPK_Diff\tL2FC\n"
                )
                for kmer in all_kmers:
                    t_c = t_counts.get(kmer, 0)
                    r_c = r_counts.get(kmer, 0)
                    t_cpk = calc_cpk(t_c, total_t_windows)
                    r_cpk = calc_cpk(r_c, total_r_windows)
                    cpk_diff = abs(t_cpk - r_cpk)
                    l2fc = calc_l2fc(t_c, r_c, total_t_windows, total_r_windows)
                    f.write(
                        f"{kmer}\t{t_c}\t{r_c}\t"
                        f"{t_cpk:.2f}\t{r_cpk:.2f}\t{cpk_diff:.2f}\t{l2fc:.3f}\n"
                    )
            print(
                f"[*] Success! Analysis written to {output_path.resolve()}",
                file=sys.stderr,
            )

        else:
            # Terminal output: sort by |L2FC| to surface the most biologically
            # distinct k-mers. Raw CPK difference is misleading because it
            # conflates effect size with absolute frequency — two k-mers at
            # CPK 1010 vs 1000 and CPK 10 vs 0 both give diff=10, but only
            # the second is biologically enriched. L2FC correctly ranks the
            # infinitely enriched case above background noise.
            print(
                f"[*] Showing top {args.top} canonical k-mers by |L2FC|\n",
                file=sys.stderr,
            )
            print(f"{'Kmer':<10} | {'Target CPK':<12} | {'Reg CPK':<12} | {'L2FC':>8}")
            print("-" * 56)

            top_kmers = sorted(
                all_kmers,
                key=lambda k: abs(
                    calc_l2fc(
                        t_counts.get(k, 0),
                        r_counts.get(k, 0),
                        total_t_windows,
                        total_r_windows,
                    )
                ),
                reverse=True,
            )[: args.top]

            for kmer in top_kmers:
                t_cpk = calc_cpk(t_counts.get(kmer, 0), total_t_windows)
                r_cpk = calc_cpk(r_counts.get(kmer, 0), total_r_windows)
                l2fc = calc_l2fc(
                    t_counts.get(kmer, 0),
                    r_counts.get(kmer, 0),
                    total_t_windows,
                    total_r_windows,
                )
                print(f"{kmer:<10} | {t_cpk:<12.2f} | {r_cpk:<12.2f} | {l2fc:>8.3f}")

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
