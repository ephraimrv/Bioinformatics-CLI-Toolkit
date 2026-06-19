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

    v1.3.0 fixes: (1) ``get_kmer_counts()`` now excludes any window
    containing a non-ACGT character (e.g. 'N' from assembly gaps or
    contig-boundary truncation) from BOTH the k-mer counts and the
    window-count denominator used for CPK normalization. Previously every
    window was counted regardless of content — a sequence with an N-gap
    silently had a meaningless "k-mer" (e.g. 'ACNNNN') counted as real
    biological signal, and that same window still inflated the CPK
    denominator even though it contributed nothing real to any k-mer's
    count. motif_discovery.py already excluded N-windows from seed
    scoring; this brings get_kmer_counts() in line with that. The
    function now returns ``(counts, n_valid_windows)`` instead of just
    ``counts``, and raises ValueError if zero valid windows remain.
    (2) Target/regulator upstream extraction now uses
    ``extract_upstream_sequence_with_length()`` (utils.py v1.3.0) and
    warns to stderr if either side's extracted window is shorter than
    requested (contig-boundary truncation) — previously this was
    silently undetectable, meaning two windows of different real length
    could be compared under the assumption they were both the requested
    length.

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
__version__ = "1.3.0"

import math
import sys
import argparse

from collections import Counter
from utils import base_parser, extract_upstream_sequence_with_length, revcomp


# ── Canonical k-mer helpers ───────────────────────────────────────────────────
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
    rc = revcomp(kmer)
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


def get_kmer_counts(sequence: str, k: int) -> tuple[Counter, int]:
    """
    Returns the canonical k-mer frequency count for a sequence, excluding
    any window that touches a non-ACGT character.

    Canonical k-mers merge the count of each k-mer with its reverse complement
    into a single entity (the lexicographically smaller of the two). This
    correctly captures TF binding sites regardless of which strand they appear
    on, because a TF binding GATA on the template strand reads as TATC on the
    coding strand — without canonicalization, you see each at half frequency.

    Windows containing 'N' (or any other non-ACGT character — assembly gaps,
    contig-boundary truncation, ambiguity codes) are excluded entirely, from
    both the count and the valid-window total. An 'N' window represents an
    unknown base, not a real k-mer; counting it both dilutes true frequencies
    with meaningless entries (e.g. 'ACNNNN') and — if included in the window
    total used for CPK normalization while excluded from the numerator —
    would systematically deflate every CPK value. motif_discovery.py already
    excludes N-windows from seed scoring; this keeps the two scripts'
    treatment of ambiguous bases consistent.

    Args:
        sequence: A nucleotide string.
        k: The length of each k-mer. Must be >= 1.

    Returns:
        A tuple of (counts, n_valid_windows):
          - counts: Counter mapping canonical k-mer strings to their integer
            counts. Each count reflects occurrences on BOTH strands combined.
          - n_valid_windows: number of windows that were ACGT-only and
            therefore counted. Use this (not raw sequence length) as the
            CPK normalization denominator.

    Raises:
        ValueError: If k < 1, if the sequence is shorter than k, or if zero
            valid (ACGT-only) windows remain after exclusion.
    """
    if k < 1:
        raise ValueError(f"k must be at least 1, got {k}.")
    if len(sequence) < k:
        raise ValueError(
            f"Sequence length ({len(sequence)} bp) is shorter than k ({k}). "
            "Reduce k or increase the upstream window."
        )

    seq = sequence.upper()
    valid_bases = frozenset("ACGT")
    kmers = []
    n_valid = 0
    for i in range(len(seq) - k + 1):
        window = seq[i : i + k]
        if not set(window) <= valid_bases:
            continue
        n_valid += 1
        kmers.append(_canonical(window))

    if n_valid == 0:
        raise ValueError(
            f"No valid ACGT-only {k}-mer windows found in a {len(seq)}bp "
            "sequence (all windows touch 'N' or another non-ACGT character). "
            "Reduce k, increase the upstream window, or check the input "
            "for excessive assembly gaps."
        )

    return Counter(kmers), n_valid


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
        # Extract sequences. Using the *_with_length variant so we can
        # detect and warn about contig-boundary truncation — otherwise a
        # gene sitting near the edge of its contig would silently return
        # a shorter-than-requested window with no indication in the output,
        # and the two windows being compared could be of meaningfully
        # different real length without the user ever knowing.
        t_seq, _, _, _, t_actual_upstream = extract_upstream_sequence_with_length(
            args.input, args.target, args.u_target
        )
        r_seq, _, _, _, r_actual_upstream = extract_upstream_sequence_with_length(
            args.input, args.regulator, args.u_regulator
        )

        if t_actual_upstream < args.u_target:
            print(
                f"[!] Warning: target '{args.target}' is only "
                f"{t_actual_upstream}bp from its contig's edge — requested "
                f"--u_target {args.u_target}bp, got {t_actual_upstream}bp.",
                file=sys.stderr,
            )
        if r_actual_upstream < args.u_regulator:
            print(
                f"[!] Warning: regulator '{args.regulator}' is only "
                f"{r_actual_upstream}bp from its contig's edge — requested "
                f"--u_regulator {args.u_regulator}bp, got {r_actual_upstream}bp.",
                file=sys.stderr,
            )

        # Count k-mers. total_*_windows now comes directly from
        # get_kmer_counts()'s n_valid_windows (ACGT-only windows actually
        # counted), not from raw sequence length — keeping the CPK
        # normalization denominator consistent with what's in the numerator.
        t_counts, total_t_windows = get_kmer_counts(t_seq, args.kmer)
        r_counts, total_r_windows = get_kmer_counts(r_seq, args.kmer)

        all_kmers = sorted(set(t_counts) | set(r_counts))

        def calc_cpk(count: int, total_windows: int) -> float:
            """Counts Per Kilobase, normalized by actual valid k-mer window count."""
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
