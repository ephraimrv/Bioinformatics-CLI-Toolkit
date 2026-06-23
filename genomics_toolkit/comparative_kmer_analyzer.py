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

PROKARYOTE-ONLY ANCHOR — NOT A WINDOW-SIZE PROBLEM:
    Both the target and regulator upstream windows are extracted via
    utils.extract_upstream_sequence_with_length(), which anchors on CDS
    start (the translation start / ATG) — not on the Transcription Start
    Site (TSS). In prokaryotes these coincide, since there is no 5' UTR
    separating them. In eukaryotes they do not: the TSS sits upstream of
    the CDS start, often separated by a 5' UTR that itself contains
    introns.

    Increasing --u_target/--u_regulator on a eukaryotic genome does NOT
    fix this — it just extracts a longer stretch of 5' UTR/intron
    sequence anchored at the wrong coordinate for BOTH genes being
    compared, not the actual promoter/enhancer region. An earlier version
    of this docstring's Examples section showed a "eukaryotic enhancer
    analysis" run against real yeast locus tags with large windows, with
    no such caveat — that example has been removed as of v1.3.1 (see
    changelog below) because it implied this script already supports
    eukaryotic use, which it does not. This script's k-mer/L2FC math
    itself (canonicalization, CPK normalization, Haldane-Anscombe
    correction) is fully organism-agnostic and needs no eukaryotic
    rework — only the upstream-extraction step it depends on is
    currently prokaryote-only. There is no eukaryote mode here, the same
    way there is none in the sibling script regulon_scanner.py, which
    uses the identical CDS-anchored mechanism and already documents this
    limitation. A future version could accept pre-extracted FASTA
    sequences directly (e.g. TSS-anchored output from
    universal_promoter_extractor.py) instead of locus-tag lookup, letting
    the k-mer/L2FC logic run unchanged on eukaryotic input — not yet
    implemented.

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

    v1.3.1: Corrected misleading eukaryotic guidance, found while
    auditing this script against its sibling regulon_scanner.py. Both
    scripts ultimately anchor on a CDS-start coordinate via utils.py's
    CDS-first feature resolution — but regulon_scanner.py already
    documented itself as PROKARYOTE-ONLY for exactly that reason, while
    this script's own Examples section demonstrated a "eukaryotic
    enhancer analysis" against real yeast locus tags as if it were a
    supported use case, with no caveat at all. Removed that example,
    added the PROKARYOTE-ONLY ANCHOR docstring section above (mirroring
    regulon_scanner.py's wording, since it is the same root cause), and
    added a one-time runtime warning — via the new shared
    utils.looks_eukaryotic() heuristic (mRNA-feature detection), checked
    once before either gene's extraction — pointing to
    universal_promoter_extractor.py / target_promoter_pipeline.py for
    actual eukaryotic TSS-anchored extraction. No change to the k-mer
    counting, L2FC, or CPK logic, and no change to prokaryote behavior or
    output format.

    v1.4.0: Five additions found during review, none changing default
    behavior unless their new flag is explicitly passed.
    (1) STATISTICAL SIGNIFICANCE (``--permutations N``, ``permutation_test()``):
    L2FC alone is an effect size, not a significance test — a k-mer
    occurring 2x in target vs 0x in regulator,
    at this script's own default window sizes (143 vs 293 windows at
    k=8), produces L2FC ~3.35 (a reported ~10x enrichment) from a
    difference of two raw counts. ``--permutations`` runs N
    nucleotide-shuffles of both sequences and reports what fraction of
    shuffles produced an equally-or-more-extreme L2FC by chance, given
    each sequence's own base composition. This does NOT by itself fix
    the deeper n=1-vs-n=1 problem (see permutation_test()'s own CAVEAT) —
    it only distinguishes "small-sample noise" from "more extreme than
    this specific composition's chance baseline."
    (2) LOW-COMPLEXITY FILTERING (``--min-entropy``, ``_shannon_entropy()``):
    overlapping-window k-mer counting inflates repetitive stretches — a
    10bp run of 'AAAAAAAAAA' produces 5
    overlapping 'AAAAAA' windows at k=6, counted as 5 occurrences of one
    actual low-complexity region. ``--min-entropy`` excludes windows
    below a Shannon-entropy threshold the same way N-windows already are.
    Default 0.0 applies no filtering.
    (3) REVERSE-COMPLEMENT REPORTING: canonicalization merges a k-mer
    and its reverse complement into one count, correctly, but previously
    only ever displayed the canonical (arbitrarily-chosen,
    lexicographically-smaller) form — a user could easily not realize
    the literal sequence they should search for elsewhere might be the
    OTHER orientation. Both TSV and terminal output now include a
    Reverse_Complement column alongside Kmer.
    (4) MOTIF CLUSTERING (``--cluster-distance N``, ``cluster_kmers()``):
    a real TFBS is often k bp wide but reported as several different
    near-identical k-mers across slightly different window offsets or
    single-base variants (e.g. TTGACA/TTGACG/TTGACC/TTGACT) — without
    grouping, these look like 4 independent hits rather than 1 motif
    with variation. Greedy Hamming-distance clustering (checked in both
    strand orientations via ``_kmer_distance()``) groups them; this is
    explicitly NOT a consensus-motif builder, just a grouping aid.
    (5) MEMORY: ``get_kmer_counts()`` previously built an intermediate
    list of every valid k-mer before constructing its Counter — now
    increments the Counter directly per window, never materializing that
    list. Irrelevant at promoter scale; matters at megabase scale.
    None of the above change k-mer counting, L2FC, or CPK behavior for
    any caller that doesn't pass the new flags — every new flag defaults
    to "off, identical to v1.3.1."

Examples:
    # Basic run: Compare two genes with default k=6, show top 20 k-mers
    $ python3 comparative_kmer_analyzer.py -i genome.gbk -t ctg1_50 -r ctg1_74 -o analysis.tsv

    # Custom k-mer size: Use 8-mers instead of 6-mers
    $ python3 comparative_kmer_analyzer.py -i genome.gbk -t ctg1_50 -r ctg1_74 -k 8 -o analysis.tsv

    # Custom upstream windows: Different for target (100bp) vs regulator (200bp)
    $ python3 comparative_kmer_analyzer.py -i genome.gbk -t ctg1_50 -r ctg1_74 --u_target 100 --u_regulator 200 -o analysis.tsv

    # Terminal output: Show top 10 k-mers sorted by |L2FC| (no file)
    $ python3 comparative_kmer_analyzer.py -i genome.gbk -t ctg1_50 -r ctg1_74 --top 10
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.4.0"

import math
import random
import sys
import argparse

from collections import Counter
from utils import (
    base_parser,
    extract_upstream_sequence_with_length,
    revcomp,
    looks_eukaryotic,
)


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


def _shannon_entropy(window: str) -> float:
    """Computes the Shannon entropy (bits) of a k-mer's base composition.

    Low entropy flags low-complexity windows. A homopolymer run like
    'AAAAAA' has entropy 0.0 (single symbol, no information content); a
    window using all 4 bases in roughly equal proportion approaches 2.0,
    the maximum possible for a 4-symbol alphabet.

    Used by ``get_kmer_counts()``'s optional ``min_entropy`` filter — see
    that function's LOW-COMPLEXITY FILTERING note for why this matters:
    overlapping-window k-mer counting otherwise inflates repetitive
    stretches (a 10bp run of 'AAAAAAAAAA' produces
    5 overlapping 'AAAAAA' windows at k=6, counted as if 5 independent
    occurrences of a real motif existed, when only one repetitive region
    does).

    Args:
        window: A k-mer string (ACGT only, by the time this is called —
                this function does not itself validate the alphabet).

    Returns:
        Entropy in bits, range [0.0, 2.0].
    """
    length = len(window)
    counts = Counter(window)
    return -sum((c / length) * math.log2(c / length) for c in counts.values())


def get_kmer_counts(
    sequence: str, k: int, min_entropy: float = 0.0
) -> tuple[Counter, int]:
    """
    Returns the canonical k-mer frequency count for a sequence, excluding
    any window that touches a non-ACGT character or (optionally) falls
    below a minimum complexity threshold.

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

    LOW-COMPLEXITY FILTERING (v1.4.0):
        Overlapping-window k-mer counting inherently inflates repetitive
        stretches: a 10bp homopolymer run like 'AAAAAAAAAA' produces 5
        overlapping 'AAAAAA' windows at k=6 — so
        one biological low-complexity region is counted as if 5
        independent occurrences of a real motif existed. ``min_entropy``
        (Shannon entropy in bits, computed per-window via
        ``_shannon_entropy()``) lets windows below the threshold be
        excluded the same way N-windows already are. Default 0.0 applies
        NO filtering — behavior for any existing caller that doesn't pass
        this argument is completely unchanged from pre-v1.4.0.

    MEMORY (v1.4.0): previously built an intermediate Python list of
    every valid canonical k-mer before constructing the ``Counter`` from
    it — O(n_valid) extra memory just to hand the list to ``Counter()``.
    Now increments the ``Counter`` directly per window in a single pass,
    never materializing that list. Irrelevant at promoter scale (a few
    hundred bp); matters if this is ever fed megabase-scale sequences.

    Args:
        sequence:    A nucleotide string.
        k:           The length of each k-mer. Must be >= 1.
        min_entropy: Minimum Shannon entropy (bits, range 0.0-2.0) a
                     window must have to be counted. Default 0.0 (no
                     filtering). Try 1.0-1.5 to exclude low-complexity
                     runs.

    Returns:
        A tuple of (counts, n_valid_windows):
          - counts: Counter mapping canonical k-mer strings to their integer
            counts. Each count reflects occurrences on BOTH strands combined.
          - n_valid_windows: number of windows that passed both the ACGT-only
            and (if set) min_entropy checks. Use this (not raw sequence
            length) as the CPK normalization denominator.

    Raises:
        ValueError: If k < 1, if min_entropy is outside [0.0, 2.0], if the
            sequence is shorter than k, or if zero valid windows remain
            after exclusion.
    """
    if k < 1:
        raise ValueError(f"k must be at least 1, got {k}.")
    if not (0.0 <= min_entropy <= 2.0):
        raise ValueError(f"min_entropy must be between 0.0 and 2.0, got {min_entropy}.")
    if len(sequence) < k:
        raise ValueError(
            f"Sequence length ({len(sequence)} bp) is shorter than k ({k}). "
            "Reduce k or increase the upstream window."
        )

    seq = sequence.upper()
    valid_bases = frozenset("ACGT")
    counts: Counter = Counter()
    n_valid = 0

    for i in range(len(seq) - k + 1):
        window = seq[i : i + k]
        if not set(window) <= valid_bases:
            continue
        if min_entropy > 0.0 and _shannon_entropy(window) < min_entropy:
            continue
        counts[_canonical(window)] += 1
        n_valid += 1

    if n_valid == 0:
        reason = (
            "all windows touch 'N' or another non-ACGT character, or fall "
            "below --min-entropy"
            if min_entropy > 0.0
            else "all windows touch 'N' or another non-ACGT character"
        )
        raise ValueError(
            f"No valid {k}-mer windows found in a {len(seq)}bp sequence "
            f"({reason}). Reduce k, increase the upstream window, lower "
            f"--min-entropy, or check the input for excessive assembly gaps."
        )

    return counts, n_valid


# ── Motif clustering ──────────────────────────────────────────────────────────


def _kmer_distance(a: str, b: str) -> int:
    """Minimum Hamming distance between two canonical k-mers, in EITHER
    relative orientation.

    Canonicalization (``_canonical()``) arbitrarily picks one strand's
    representation for each k-mer — so two biologically related k-mers
    could end up canonicalized to opposite strands and look unrelated by
    naive same-orientation Hamming distance alone. Checking both the
    direct comparison and the comparison against ``b``'s reverse
    complement catches that case.

    Args:
        a: First canonical k-mer.
        b: Second canonical k-mer.

    Returns:
        The smaller of the two Hamming distances (direct, and against
        ``revcomp(b)``). K-mers of different length are never considered
        close (returns the longer length, guaranteed to exceed any
        realistic ``max_distance`` threshold).
    """
    if len(a) != len(b):
        return max(len(a), len(b))
    direct = sum(1 for x, y in zip(a, b) if x != y)
    flipped = sum(1 for x, y in zip(a, revcomp(b)) if x != y)
    return min(direct, flipped)


def cluster_kmers(kmers: list[str], max_distance: int = 1) -> list[list[str]]:
    """Greedily clusters k-mers that likely represent variants of one motif.

    Two k-mers within ``max_distance`` Hamming distance of each other
    (checked in both direct and reverse-complement orientation — see
    ``_kmer_distance()``) are grouped into the same cluster. This is
    deliberately simple greedy single-linkage clustering, NOT a true
    consensus-motif builder: it groups likely-related k-mers (e.g.
    TTGACA/TTGACG/TTGACC/TTGACT, each one mismatch from the others) so a
    reader can see they probably represent one motif rather than four
    independent ones — it does not attempt to reconstruct an IUPAC
    consensus or claim which member is the "true" motif.

    Args:
        kmers:        List of canonical k-mer strings to cluster (e.g.
                       the top N by |L2FC|). Order matters for greedy
                       assignment — earlier k-mers seed clusters that
                       later k-mers may join.
        max_distance: Maximum Hamming distance (in either orientation)
                      for two k-mers to be grouped together. Default: 1.

    Returns:
        A list of clusters, each a list of k-mers, in the order clusters
        were first created. Singletons (no close match found) form their
        own one-member cluster.
    """
    clusters: list[list[str]] = []
    for kmer in kmers:
        placed = False
        for cluster in clusters:
            if any(_kmer_distance(kmer, member) <= max_distance for member in cluster):
                cluster.append(kmer)
                placed = True
                break
        if not placed:
            clusters.append([kmer])
    return clusters


# ── Permutation significance testing ──────────────────────────────────────────


def permutation_test(
    t_seq: str,
    r_seq: str,
    k: int,
    min_entropy: float,
    observed_l2fc: dict[str, float],
    n_permutations: int = 1000,
    seed: int | None = None,
) -> dict[str, float]:
    """Computes a two-tailed permutation p-value for each k-mer's observed L2FC.

    RATIONALE: ``calc_l2fc()`` reports effect size, not statistical
    significance. A k-mer occurring just 2x in a
    143-window target vs 0x in a 293-window regulator (this script's own
    default --u_target/--u_regulator at k=8) produces an L2FC of ~3.35 —
    a reported ~10x enrichment from a difference of two raw counts. With
    only ONE target sequence and ONE regulator sequence being compared
    (n=1 per condition), there is no way to know whether an observed L2FC
    reflects real biology or pure chance without some null model.

    METHOD: each permutation independently shuffles the NUCLEOTIDES of
    both the target and regulator sequences (preserving each sequence's
    own single-base composition, but destroying any real positional
    signal), recomputes k-mer counts and L2FC for every k-mer under that
    random shuffle, and checks whether the shuffled |L2FC| meets or
    exceeds the REAL observed |L2FC|. The p-value is the fraction of
    permutations where that happens (two-tailed: enrichment in either
    direction counts as meeting the bar). This is a standard Monte-Carlo
    permutation approach for testing whether an enrichment statistic is
    more extreme than chance, given each sequence's own base composition
    as the null model. Implemented as a running exceedance count per
    k-mer rather than storing full null distributions, to keep memory
    flat regardless of ``n_permutations`` or k-mer vocabulary size.

    CAVEAT — this does NOT fully solve the n=1-vs-n=1 problem: this tests
    "is this L2FC more extreme than chance, given THESE two sequences'
    own base composition," not "would this hold across biological
    replicates of target/regulator-type promoters." It controls for the
    small-window-size effect-size-inflation problem, but true biological
    replication requires comparing SETS of promoters, not one sequence
    per side — see this module's discussion of a possible future
    multi-promoter mode.

    PERFORMANCE: each permutation re-runs ``get_kmer_counts()`` on both
    shuffled sequences (cheap at promoter scale — a few hundred bp) and
    then does one dict lookup per observed k-mer. The default of 1000
    permutations is fast at default window sizes; very large k or very
    large ``--u_target``/``--u_regulator`` values will scale this up
    proportionally.

    Args:
        t_seq:          Target upstream sequence (already extracted).
        r_seq:          Regulator upstream sequence (already extracted).
        k:              K-mer length (must match the real analysis).
        min_entropy:    Same low-complexity filter as the real analysis.
        observed_l2fc:  {kmer: real_l2fc} from the actual (unshuffled)
                        sequences — the values being tested.
        n_permutations: Number of random shuffles. Default 1000.
        seed:           Optional RNG seed for reproducibility.

    Returns:
        {kmer: p_value} for every kmer in ``observed_l2fc``.
    """
    rng = random.Random(seed)
    exceed_counts: Counter = Counter()
    t_chars = list(t_seq.upper())
    r_chars = list(r_seq.upper())
    valid_permutations = 0

    for _ in range(n_permutations):
        rng.shuffle(t_chars)
        rng.shuffle(r_chars)
        shuffled_t = "".join(t_chars)
        shuffled_r = "".join(r_chars)

        try:
            t_counts, t_windows = get_kmer_counts(shuffled_t, k, min_entropy)
        except ValueError:
            continue
        try:
            r_counts, r_windows = get_kmer_counts(shuffled_r, k, min_entropy)
        except ValueError:
            continue

        valid_permutations += 1
        for kmer, obs_l2fc in observed_l2fc.items():
            null_l2fc = calc_l2fc(
                t_counts.get(kmer, 0), r_counts.get(kmer, 0), t_windows, r_windows
            )
            if abs(null_l2fc) >= abs(obs_l2fc):
                exceed_counts[kmer] += 1

    if valid_permutations == 0:
        # Every shuffle happened to produce zero valid windows (e.g.
        # min_entropy set so high nothing ever passes) — report 1.0
        # (no evidence of significance) rather than dividing by zero.
        return {kmer: 1.0 for kmer in observed_l2fc}

    return {kmer: exceed_counts[kmer] / valid_permutations for kmer in observed_l2fc}


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
        "--u_target",
        type=int,
        default=150,
        help=(
            "Upstream bp for target. Default: 150. NOTE: extraction is "
            "CDS-anchored, not TSS-anchored — increasing this value does "
            "NOT adapt this script for eukaryotic use (see the "
            "PROKARYOTE-ONLY ANCHOR note in this script's docstring)."
        ),
    )
    parser.add_argument(
        "--u_regulator",
        type=int,
        default=300,
        help=(
            "Upstream bp for regulator. Default: 300. Same CDS-anchored "
            "caveat as --u_target applies."
        ),
    )
    parser.add_argument("-k", "--kmer", type=int, default=6, help="K-mer length")
    parser.add_argument("--top", type=int, default=20, help="Top N k-mers to report")
    parser.add_argument(
        "--min-entropy",
        type=float,
        default=0.0,
        metavar="BITS",
        help=(
            "Exclude k-mer windows below this Shannon entropy (bits, "
            "range 0.0-2.0). Default: 0.0 (no filtering — identical to "
            "pre-v1.4.0 behavior). A homopolymer run like 'AAAAAA' has "
            "entropy 0.0; try 1.0-1.5 to exclude low-complexity runs that "
            "would otherwise inflate apparent k-mer counts through pure "
            "window overlap, with no real biological signal."
        ),
    )
    parser.add_argument(
        "--permutations",
        type=int,
        default=0,
        metavar="N",
        help=(
            "Run N nucleotide-shuffle permutations to compute a p-value "
            "per k-mer, testing whether its L2FC is more extreme than "
            "chance given each sequence's own base composition. Default: "
            "0 (disabled — L2FC alone is an effect size, not a "
            "significance test; see this script's docstring). Does NOT "
            "by itself solve the deeper n=1-vs-n=1 problem (one target "
            "sequence vs one regulator sequence) — only that this "
            "specific L2FC value isn't simply small-sample noise."
        ),
    )
    parser.add_argument(
        "--perm-seed",
        type=int,
        default=None,
        metavar="N",
        help="Optional RNG seed for --permutations, for reproducible p-values.",
    )
    parser.add_argument(
        "--cluster-distance",
        type=int,
        default=0,
        metavar="N",
        help=(
            "Group displayed/output k-mers within N Hamming distance of "
            "each other (checked in both strand orientations) into "
            "labeled clusters — e.g. TTGACA/TTGACG/TTGACC/TTGACT likely "
            "represent variants of one motif, not four independent ones. "
            "Default: 0 (disabled, no clustering). This is simple greedy "
            "grouping, not a consensus-motif builder."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = get_args()

    # Validation
    if args.kmer < 1:
        sys.exit("[!] --kmer must be at least 1.")
    if args.u_target < 1 or args.u_regulator < 1:
        sys.exit("[!] Upstream values must be positive integers.")
    if not (0.0 <= args.min_entropy <= 2.0):
        sys.exit("[!] --min-entropy must be between 0.0 and 2.0.")
    if args.permutations < 0:
        sys.exit("[!] --permutations must be >= 0.")
    if args.cluster_distance < 0:
        sys.exit("[!] --cluster-distance must be >= 0.")

    # args.output is already a Path from base_parser (type=Path); no re-wrapping needed
    output_path = args.output

    try:
        # One-time heuristic warning: mRNA features are a strong signal of
        # eukaryotic annotation (Prokka/Bakta prokaryote output never emits
        # them). Both extraction calls below are CDS-anchored, not
        # TSS-anchored (see the PROKARYOTE-ONLY ANCHOR docstring section)
        # — checked once here rather than once per gene, since both genes
        # come from the same input file.
        if looks_eukaryotic(args.input):
            print(
                "[!] Warning: mRNA features detected — this looks like a "
                "eukaryotic genome. Both the target and regulator upstream "
                "windows are anchored on CDS start, not the transcription "
                "start site (TSS), so the true promoter/enhancer region "
                "for either gene will likely be missed. See the "
                "PROKARYOTE-ONLY ANCHOR note in this script's docstring. "
                "For eukaryotic regulatory analysis, extract upstream "
                "regions with universal_promoter_extractor.py or "
                "target_promoter_pipeline.py instead.",
                file=sys.stderr,
            )

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
        # get_kmer_counts()'s n_valid_windows (ACGT-only, and (if set)
        # min-entropy-passing windows actually counted), not from raw
        # sequence length — keeping the CPK normalization denominator
        # consistent with what's in the numerator.
        t_counts, total_t_windows = get_kmer_counts(t_seq, args.kmer, args.min_entropy)
        r_counts, total_r_windows = get_kmer_counts(r_seq, args.kmer, args.min_entropy)

        all_kmers = sorted(set(t_counts) | set(r_counts))

        def calc_cpk(count: int, total_windows: int) -> float:
            """Counts Per Kilobase, normalized by actual valid k-mer window count."""
            return (count / total_windows) * 1000

        # Precompute L2FC for every k-mer once, used by output, clustering,
        # and (if enabled) as the observed values for the permutation test.
        l2fc_map = {
            kmer: calc_l2fc(
                t_counts.get(kmer, 0),
                r_counts.get(kmer, 0),
                total_t_windows,
                total_r_windows,
            )
            for kmer in all_kmers
        }

        # Clustering is computed once, ranked by |L2FC| descending regardless
        # of output mode — the most-enriched k-mer should seed each cluster,
        # not whichever happens to come first alphabetically (the TSV's
        # natural sort order). cluster_id_map is then applied to whichever
        # rows actually get written.
        cluster_id_map: dict[str, int] = {}
        if args.cluster_distance > 0:
            ranked_for_clustering = sorted(
                all_kmers, key=lambda km: abs(l2fc_map[km]), reverse=True
            )
            clusters = cluster_kmers(
                ranked_for_clustering, max_distance=args.cluster_distance
            )
            for cluster_idx, cluster in enumerate(clusters, 1):
                for kmer in cluster:
                    cluster_id_map[kmer] = cluster_idx
            multi_member = sum(1 for c in clusters if len(c) > 1)
            print(
                f"[*] Clustering (max distance {args.cluster_distance}): "
                f"{len(clusters)} cluster(s), {multi_member} with more than "
                f"one member.",
                file=sys.stderr,
            )

        if output_path:
            # Permutation p-values, if requested, computed for every k-mer
            # written to the TSV (i.e. all_kmers, not just a displayed
            # subset — the TSV is the complete-record output mode).
            p_value_map: dict[str, float] = {}
            if args.permutations > 0:
                print(
                    f"[*] Running {args.permutations:,} permutations for "
                    f"{len(all_kmers):,} k-mer(s)...",
                    file=sys.stderr,
                )
                p_value_map = permutation_test(
                    t_seq,
                    r_seq,
                    args.kmer,
                    args.min_entropy,
                    l2fc_map,
                    n_permutations=args.permutations,
                    seed=args.perm_seed,
                )

            header_cols = [
                "Kmer",
                "Reverse_Complement",
                "Target_Count",
                "Regulator_Count",
                "Target_CPK",
                "Regulator_CPK",
                "CPK_Diff",
                "L2FC",
            ]
            if args.permutations > 0:
                header_cols.append("P_Value")
            if args.cluster_distance > 0:
                header_cols.append("Cluster_ID")

            with open(output_path, "w", encoding="utf-8") as f:
                f.write("\t".join(header_cols) + "\n")
                for kmer in all_kmers:
                    t_c = t_counts.get(kmer, 0)
                    r_c = r_counts.get(kmer, 0)
                    t_cpk = calc_cpk(t_c, total_t_windows)
                    r_cpk = calc_cpk(r_c, total_r_windows)
                    cpk_diff = abs(t_cpk - r_cpk)
                    l2fc = l2fc_map[kmer]
                    row = [
                        kmer,
                        revcomp(kmer),
                        str(t_c),
                        str(r_c),
                        f"{t_cpk:.2f}",
                        f"{r_cpk:.2f}",
                        f"{cpk_diff:.2f}",
                        f"{l2fc:.3f}",
                    ]
                    if args.permutations > 0:
                        row.append(f"{p_value_map[kmer]:.4f}")
                    if args.cluster_distance > 0:
                        row.append(str(cluster_id_map.get(kmer, "")))
                    f.write("\t".join(row) + "\n")
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

            top_kmers = sorted(
                all_kmers, key=lambda km: abs(l2fc_map[km]), reverse=True
            )[: args.top]

            # Permutation p-values, if requested, computed only for the
            # k-mers actually being displayed here — no point spending
            # 1000 shuffles per k-mer on entries that never make the cut.
            p_value_map: dict[str, float] = {}
            if args.permutations > 0:
                print(
                    f"[*] Running {args.permutations:,} permutations for "
                    f"{len(top_kmers):,} displayed k-mer(s)...",
                    file=sys.stderr,
                )
                observed_subset = {km: l2fc_map[km] for km in top_kmers}
                p_value_map = permutation_test(
                    t_seq,
                    r_seq,
                    args.kmer,
                    args.min_entropy,
                    observed_subset,
                    n_permutations=args.permutations,
                    seed=args.perm_seed,
                )

            header = f"{'Kmer':<10} | {'RevComp':<10} | {'Target CPK':<12} | {'Reg CPK':<12} | {'L2FC':>8}"
            if args.permutations > 0:
                header += f" | {'P_Value':>8}"
            if args.cluster_distance > 0:
                header += f" | {'Cluster':>7}"
            print(header)
            print("-" * len(header))

            for kmer in top_kmers:
                t_cpk = calc_cpk(t_counts.get(kmer, 0), total_t_windows)
                r_cpk = calc_cpk(r_counts.get(kmer, 0), total_r_windows)
                l2fc = l2fc_map[kmer]
                line = (
                    f"{kmer:<10} | {revcomp(kmer):<10} | {t_cpk:<12.2f} | "
                    f"{r_cpk:<12.2f} | {l2fc:>8.3f}"
                )
                if args.permutations > 0:
                    line += f" | {p_value_map[kmer]:>8.4f}"
                if args.cluster_distance > 0:
                    line += f" | {cluster_id_map.get(kmer, ''):>7}"
                print(line)

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
