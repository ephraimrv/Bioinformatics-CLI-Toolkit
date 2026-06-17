#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

r"""EM Motif Discovery

Discovers transcription factor binding sites and regulatory motifs in upstream
sequences using Expectation-Maximization (MEME-style OOPS model).

Background and algorithm:
    The classic sliding window + IC approach finds the most absolutely conserved
    windows. This is wrong for regulatory motif discovery because:

    1. Background matters. Your upstream sequences are AT-rich (~65% A+T).
       An AT-rich conserved block (e.g. AACATTTCAACATGA) scores 30 bits IC
       but its log-odds vs. background is LOWER than a GC-rich regulatory
       motif, because A and T are already expected by chance. MEME adjusts
       for this; sliding window IC does not.

    2. Motifs float. A transcription factor binding site does not occur at
       the same absolute position in every sequence. MEME's EM assigns each
       sequence independently to the position where the motif best fits.
       Sliding window uses the same absolute position in all sequences.

    This script implements a simplified MEME OOPS model (One Occurrence Per
    Sequence, Bailey & Elkan 1994):

      Seed     → For each distinct k-mer in the input sequences, build an
                 initial Position Weight Matrix (PWM) and score it.
      Cluster  → Before advancing seeds to EM, filter out seeds within
                 --min-hamming Hamming distance of a higher-scoring seed.
                 This prevents redundant EM runs on shifted variants of the
                 same underlying motif (the "seed collision" problem).
      Quick EM → Run 25 EM iterations for the top --seeds diverse candidates.
      Full EM  → Run --iter EM iterations for the top --refine candidates.
      Report   → The best converged motif is reported. Its instances are
                 masked and the process repeats for the next motif.

    E-step: P(motif at position j | sequence, PWM) ∝ exp(log-odds score)
    M-step: Update PWM using the fractional expected counts (soft EM),
            using whichever strand (forward or reverse complement) gave
            the higher log-odds score at each position.
    Score:  Sum of log-likelihood across all sequences.

BIDIRECTIONAL STRAND SCANNING:
    Transcription factors bind double-stranded DNA. They do not care which
    strand the genome annotators labeled the "coding" strand. A TF binding
    GATA on the template strand appears as TATC in your extracted sequence.

    This script scans BOTH strands at every position during the E-step and
    M-step. At each position j, the forward window score and the reverse
    complement window score are both computed; only the higher of the two
    contributes to the probability and PWM update. The strand that produced
    the best score for the best-position instance in each sequence is reported
    in the output as '+' or '-'.

    Note: upstream sequences produced by universal_promoter_extractor.py are
    already strand-corrected (5'→3' relative to the gene). Bidirectional
    scanning additionally handles palindromic TF binding sites and sequences
    supplied from external sources.

SEED CLUSTERING (diversity-aware seeding):
    All unique W-mers from the input sequences are scored and ranked. Without
    clustering, the top --seeds seeds may all be shifted variants of the same
    dominant motif (e.g., ATCGATCG, TCGATCGA, CGATCGAT), wasting all --seeds
    EM runs on convergence to the same PWM.

    After scoring, seeds are greedily clustered by Hamming distance. A seed
    is only advanced to EM if it differs by at least --min-hamming positions
    from every already-selected seed. This ensures the seed pool is physically
    diverse, so EM has a real chance of discovering multiple distinct motifs.
    Default min-hamming is 3 (configurable via --min-hamming).

MASKING THRESHOLD (prevents corruption of multi-motif discovery):
    After each motif is found, the OOPS model forces a "best" position in
    every sequence to be masked. If a sequence does not actually contain
    Motif 1, the forced match may overwrite a real binding site for Motif 2.

    The --mask-threshold flag (default 0.0) sets a log-odds score floor.
    A position is only masked if its score exceeds the threshold, meaning
    the subsequence looks more like the motif than background. Sequences
    with no plausible Motif 1 site are left unmasked, protecting their
    Motif 2 sites. A threshold of 0.0 corresponds to break-even with the
    background model; increase it (e.g., 2.0) for stricter masking.

PERFORMANCE — zip+Counter COLUMN PROFILING:
    The profile builder in _build_profile_from_instances() uses
    zip(*sequences) + Counter for column-wise counting rather than a
    nested Python loop. Counter runs in C (CPython), substantially reducing
    the number of Python-level operations for large sequence sets.

Note:
    Associated with ongoing, unpublished research (manuscript in
    preparation). Correct attribution is requested when used in
    derivative works.

Examples:
    # Find top 3 motifs of default 15bp width
    $ python3 motif_discovery.py -i upstream.fasta --top 3 -o motifs.tsv

    # Custom width, more seeding, verbose
    $ python3 motif_discovery.py -i upstream.fasta --top 3 --width 12 \
      --seeds 80 --refine 30 -o motifs.tsv

    # Diverse seeds: reject seeds within 4 positions of a higher-scoring seed
    $ python3 motif_discovery.py -i upstream.fasta --top 3 --min-hamming 4

    # Stricter masking: only mask positions scoring >2.0 log-odds
    $ python3 motif_discovery.py -i upstream.fasta --top 3 --mask-threshold 2.0

Caveats:
    - Input sequences do NOT need to be aligned (raw upstream FASTA is fine).
    - All sequences must be at least --width bp long.
    - Results are supporting calculations for MEME Suite, not a replacement.
      Always confirm discovered motifs with full MEME analysis.
    - For TSV output: open with File > Open in Excel (not Data > Get Data)
      to avoid the 'Transform Data' wizard.
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.1.0"

import math
import sys
from collections import Counter

try:
    from Bio import SeqIO
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )
from utils import base_parser

# ── Reverse complement ────────────────────────────────────────────────────────

_RC_TABLE = str.maketrans("ACGTN", "TGCAN")


def _revcomp(seq: str) -> str:
    """Return the reverse complement of a DNA string. N maps to N."""
    return seq.translate(_RC_TABLE)[::-1]


# ── Background model ──────────────────────────────────────────────────────────


def _compute_background(sequences: list[str]) -> dict[str, float]:
    """Compute empirical nucleotide background frequencies.

    Uses observed base frequencies across all input sequences combined.
    A minimum floor of 0.01 prevents log(0) in log-odds calculations.

    Args:
        sequences: List of uppercase DNA strings.

    Returns:
        Dict mapping "ACGT" → frequency, summing to ~1.0.
    """
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    for seq in sequences:
        for ch in seq:
            if ch in counts:
                counts[ch] += 1
    total = sum(counts.values()) or 1
    return {b: max(counts[b] / total, 0.01) for b in "ACGT"}


# ── PWM initialisation ────────────────────────────────────────────────────────


def _init_pwm(kmer: str, pseudocount: float) -> dict[str, list[float]]:
    """Build an initial PWM from a single k-mer seed.

    Adds pseudocounts to all bases at all positions to avoid zero probabilities,
    then places an additional count of 1.0 at the seed's observed base.

    Args:
        kmer:        The seed k-mer string (uppercase ACGT).
        pseudocount: Smoothing value added to every cell before normalisation.

    Returns:
        PWM dict mapping "ACGT" → list of probabilities (one per position).
    """
    W = len(kmer)
    counts = {b: [pseudocount] * W for b in "ACGT"}
    for pos, base in enumerate(kmer):
        if base in counts:
            counts[base][pos] += 1.0
    pwm: dict[str, list[float]] = {}
    for pos in range(W):
        total = sum(counts[b][pos] for b in "ACGT")
        for b in "ACGT":
            if b not in pwm:
                pwm[b] = [0.0] * W
            pwm[b][pos] = counts[b][pos] / total
    return pwm


# ── EM core ───────────────────────────────────────────────────────────────────


def _score_window(window: str, pwm: dict, bg: dict) -> float:
    """Log-odds score for a pre-extracted window string against a PWM.

    Core scoring primitive used by both the forward and reverse complement
    paths in the E-step, M-step, and instance extraction. Returns -inf if
    any character in the window is not in the PWM (e.g. masked 'N').

    Args:
        window: DNA string of exactly the motif width.
        pwm:    Current PWM dict (``"ACGT"`` → list of probabilities).
        bg:     Background nucleotide frequencies.

    Returns:
        Sum of position-wise log(PWM[base] / background[base]), or -inf.
    """
    score = 0.0
    for k, base in enumerate(window):
        if base not in pwm:
            return float("-inf")
        p = max(pwm[base][k], 1e-300)
        b = max(bg.get(base, 0.25), 1e-300)
        score += math.log(p / b)
    return score


def _pos_log_odds(seq: str, start: int, pwm: dict, bg: dict, width: int) -> float:
    """Log-odds score for the motif at position ``start`` in ``seq`` (forward strand).

    Thin wrapper around ``_score_window`` for positional indexing.
    For bidirectional scoring, call ``_score_window`` on both the forward
    window and its reverse complement and take the max.
    """
    return _score_window(seq[start : start + width], pwm, bg)


def _e_step(sequences: list[str], pwm: dict, bg: dict, width: int) -> list[list[float]]:
    """E-step: compute P(motif at position j | sequence, PWM) via softmax.

    Scans BOTH strands at every position. At each start position j, the
    forward window and its reverse complement are both scored; the higher
    log-odds score is used as the position's score. This correctly captures
    TF binding sites regardless of which strand they appear on.

    For each sequence, all valid start positions are scored, then normalised
    to a probability distribution via softmax (numerically stable: subtract
    max before exponentiating).

    Args:
        sequences: Input sequences (may contain 'N' at masked positions).
        pwm:       Current PWM (dict ``"ACGT"`` → list[float]).
        bg:        Background model.
        width:     Motif width.

    Returns:
        List of probability lists, one per sequence.
    """
    all_probs = []
    for seq in sequences:
        n_pos = len(seq) - width + 1
        if n_pos <= 0:
            all_probs.append([1.0])
            continue

        log_scores = []
        for j in range(n_pos):
            window = seq[j : j + width]
            fwd = _score_window(window, pwm, bg)
            rev = _score_window(_revcomp(window), pwm, bg)
            log_scores.append(max(fwd, rev))

        # Softmax with -inf handling
        valid = [s for s in log_scores if s > float("-inf")]
        if not valid:
            all_probs.append([1.0 / n_pos] * n_pos)
            continue

        max_s = max(valid)
        exp_s = [math.exp(s - max_s) if s > float("-inf") else 0.0 for s in log_scores]
        total = sum(exp_s) or 1.0
        all_probs.append([e / total for e in exp_s])

    return all_probs


def _m_step(
    sequences: list[str],
    probs: list[list[float]],
    width: int,
    pseudocount: float,
    pwm: dict,
    bg: dict,
) -> dict[str, list[float]]:
    """M-step: update PWM using fractional expected counts.

    For each position j with probability z, determines whether the forward
    window or its reverse complement scored higher under the current PWM,
    then uses the better-scoring strand's sequence to update the counts.
    This ensures the PWM is updated toward the actual strand orientation
    of the motif, rather than always accumulating forward-strand windows
    even when the true binding site is on the opposite strand.

    Args:
        sequences:   Input sequences.
        probs:       E-step output (position probabilities per sequence).
        width:       Motif width.
        pseudocount: Smoothing value.
        pwm:         Current PWM (used to determine which strand scores higher).
        bg:          Background model.

    Returns:
        Updated PWM.
    """
    counts = {b: [pseudocount] * width for b in "ACGT"}

    for seq, seq_probs in zip(sequences, probs):
        for j, z in enumerate(seq_probs):
            if z < 1e-15:
                continue
            window = seq[j : j + width]
            if len(window) < width:
                break

            # Use whichever strand's window scored higher under the current PWM
            rc_window = _revcomp(window)
            fwd_score = _score_window(window, pwm, bg)
            rev_score = _score_window(rc_window, pwm, bg)
            use_window = rc_window if rev_score > fwd_score else window

            for k, base in enumerate(use_window):
                if base in counts:
                    counts[base][k] += z

    pwm_new: dict[str, list[float]] = {}
    for pos in range(width):
        total = sum(counts[b][pos] for b in "ACGT") or 1.0
        for b in "ACGT":
            if b not in pwm_new:
                pwm_new[b] = [0.0] * width
            pwm_new[b][pos] = counts[b][pos] / total

    return pwm_new


def _log_likelihood(sequences: list[str], pwm: dict, bg: dict, width: int) -> float:
    """Total log-likelihood of all sequences under the current PWM.

    Scores both forward and reverse complement at each position (taking the
    max), consistent with the bidirectional E-step. Uses log-sum-exp for
    numerical stability.
    """
    total = 0.0
    for seq in sequences:
        n_pos = len(seq) - width + 1
        if n_pos <= 0:
            continue
        log_scores = []
        for j in range(n_pos):
            window = seq[j : j + width]
            fwd = _score_window(window, pwm, bg)
            rev = _score_window(_revcomp(window), pwm, bg)
            log_scores.append(max(fwd, rev))
        valid = [s for s in log_scores if s > float("-inf")]
        if not valid:
            continue
        max_s = max(valid)
        total += max_s + math.log(sum(math.exp(s - max_s) for s in valid))
    return total


def _run_em(
    sequences: list[str],
    initial_pwm: dict,
    bg: dict,
    width: int,
    max_iter: int,
    tol: float,
    pseudocount: float,
) -> tuple[dict, float]:
    """Run EM to convergence from an initial PWM.

    Args:
        sequences:   Input (possibly masked) sequences.
        initial_pwm: Starting PWM.
        bg:          Background model.
        width:       Motif width.
        max_iter:    Maximum number of E/M cycles.
        tol:         Convergence threshold (change in log-likelihood).
        pseudocount: M-step smoothing value.

    Returns:
        (converged_pwm, final_log_likelihood)
    """
    pwm = {b: initial_pwm[b][:] for b in "ACGT"}
    prev_ll = float("-inf")

    for _ in range(max_iter):
        probs = _e_step(sequences, pwm, bg, width)
        pwm = _m_step(sequences, probs, width, pseudocount, pwm, bg)
        ll = _log_likelihood(sequences, pwm, bg, width)
        if ll - prev_ll < tol:
            break
        prev_ll = ll

    return pwm, _log_likelihood(sequences, pwm, bg, width)


# ── Post-EM: extract instances and mask ───────────────────────────────────────


def _extract_instances(
    sequences: list[str],
    seq_ids: list[str],
    pwm: dict,
    bg: dict,
    width: int,
) -> list[tuple[str, str, int, str]]:
    """Find the best-scoring position of the motif in each sequence.

    Scores both forward and reverse complement at every position and picks
    the overall best. Reports which strand produced the winning score.

    Returns:
        List of (seq_id, extracted_sequence, start_position_1indexed, strand)
        where strand is '+' (forward) or '-' (reverse complement).
    """
    instances = []
    for seq, sid in zip(sequences, seq_ids):
        n_pos = len(seq) - width + 1
        if n_pos <= 0:
            instances.append((sid, seq[:width].ljust(width, "N"), 1, "+"))
            continue

        best_j = 0
        best_score = float("-inf")
        best_strand = "+"

        for j in range(n_pos):
            window = seq[j : j + width]
            fwd_score = _score_window(window, pwm, bg)
            rev_score = _score_window(_revcomp(window), pwm, bg)
            if rev_score > fwd_score:
                score, strand = rev_score, "-"
            else:
                score, strand = fwd_score, "+"
            if score > best_score:
                best_score = score
                best_j = j
                best_strand = strand

        instances.append((sid, seq[best_j : best_j + width], best_j + 1, best_strand))
    return instances


def _mask_sequences(
    sequences: list[str],
    instances: list[tuple],
    pwm: dict,
    bg: dict,
    width: int,
    threshold: float = 0.0,
) -> tuple[list[str], int]:
    """Replace confirmed motif instances with 'N' to mask for subsequent searches.

    Only masks a position if its forward-strand log-odds score exceeds
    ``threshold``. This prevents the OOPS masking danger: if a sequence does
    not actually contain Motif 1, the forced "best" position may coincide with
    a real Motif 2 binding site. Without a threshold, that Motif 2 site gets
    masked unconditionally, corrupting subsequent motif discovery. With a
    threshold, sequences where the best Motif 1 score is below background
    (score ≤ threshold) are left unmasked.

    A threshold of 0.0 means "only mask if the subsequence looks more like
    the motif than background." Increase (e.g., 2.0) for stricter masking.

    'N' characters score -inf in ``_score_window``, making them invisible
    to the EM in subsequent motif searches.

    Args:
        sequences: Current (possibly already masked) working sequences.
        instances: Output of ``_extract_instances`` — 4-tuples.
        pwm:       Converged PWM of the motif just found.
        bg:        Background model.
        width:     Motif width.
        threshold: Minimum log-odds score required to mask (default: 0.0).

    Returns:
        (masked_sequences, n_masked) — the updated sequence list and a count
        of how many sequences were actually masked.
    """
    masked = []
    n_masked = 0
    for seq, (_, _, start_1, _) in zip(sequences, instances):
        start_0 = start_1 - 1
        score = _pos_log_odds(seq, start_0, pwm, bg, width)
        if score > threshold:
            s = list(seq)
            for k in range(width):
                if start_0 + k < len(s):
                    s[start_0 + k] = "N"
            masked.append("".join(s))
            n_masked += 1
        else:
            masked.append(seq)  # score too low — likely not a real site; leave unmasked
    return masked, n_masked


# ── Profile from instances (for display) ──────────────────────────────────────


def _build_profile_from_instances(
    extracted: list[str],
) -> tuple[str, dict[str, list[float]], list[float]]:
    """Build consensus, PPM, and IC from a list of extracted motif sequences.

    Uses zip(*extracted) + Counter for column-wise counting (C-backed,
    faster than a Python nested loop for large instance sets).

    Identical in output to alignment_conservation_profiler._build_profile
    but inlined here so this script is standalone.
    """
    seq_len = len(extracted[0])
    n = len(extracted)
    valid = "ACGT-"
    counts = {c: [0] * seq_len for c in valid}
    max_bits = math.log2(4.0)

    for i, column in enumerate(zip(*extracted)):
        col_counts = Counter(column)
        for c in valid:
            counts[c][i] = col_counts[c]  # Counter returns 0 for missing keys

    consensus_list = []
    ppm = {c: [0.0] * seq_len for c in valid}
    ic = [0.0] * seq_len

    for i in range(seq_len):
        best, bch = -1, "-"
        ent = 0.0
        acgt = sum(counts[c][i] for c in "ACGT")
        for c in valid:
            cnt = counts[c][i]
            ppm[c][i] = cnt / n
            if cnt > best:
                best, bch = cnt, c
            if c in "ACGT" and acgt > 0:
                p = cnt / acgt
                if p > 0:
                    ent += p * math.log2(p)
        consensus_list.append(bch)
        fp = acgt / n
        ic[i] = max(0.0, (max_bits + ent) * fp)

    return "".join(consensus_list), ppm, ic


# ── Seed clustering ───────────────────────────────────────────────────────────


def _hamming(a: str, b: str) -> int:
    """Hamming distance between two equal-length strings."""
    return sum(x != y for x, y in zip(a, b))


def _cluster_seeds(
    scored_seeds: list[tuple[float, str]],
    n_seeds: int,
    min_hamming: int,
) -> list[str]:
    """Greedily select diverse seeds by Hamming distance.

    Iterates through seeds in descending score order. A seed is only
    selected if it differs from every already-selected seed by at least
    ``min_hamming`` positions. This prevents --seeds EM runs from all
    converging to shifted variants of the same dominant motif.

    Args:
        scored_seeds: List of (score, kmer) sorted descending by score.
        n_seeds:      Maximum number of diverse seeds to return.
        min_hamming:  Minimum Hamming distance between any two selected seeds.

    Returns:
        List of up to ``n_seeds`` diverse seed k-mers.
    """
    diverse: list[str] = []
    for _, kmer in scored_seeds:
        if all(_hamming(kmer, kept) >= min_hamming for kept in diverse):
            diverse.append(kmer)
            if len(diverse) >= n_seeds:
                break
    return diverse


# ── Main discovery pipeline ───────────────────────────────────────────────────


def discover_motifs(
    sequences: list[str],
    seq_ids: list[str],
    width: int,
    top_n: int,
    n_seeds: int,
    n_refine: int,
    pseudocount: float,
    min_hamming: int = 3,
    mask_threshold: float = 0.0,
) -> list[dict]:
    """Run EM motif discovery: find top N non-overlapping motifs.

    Strategy:
      1. Collect all unique W-mers from all input sequences as candidate seeds.
      2. Score each seed by initial log-likelihood (bidirectional).
      3. Cluster seeds by Hamming distance — only advance seeds that differ
         by at least ``min_hamming`` from all higher-scoring seeds already
         selected. This ensures the seed pool is physically diverse.
      4. Quick EM (25 iterations) on the top ``n_seeds`` diverse seeds.
      5. Full EM (``--iter`` iterations) on the top ``n_refine`` from step 4.
      6. Report the best converged motif.
      7. Mask confirmed instances (score > ``mask_threshold``) with 'N' and
         repeat for the next motif.

    Args:
        sequences:      Input sequences (>= width bp each).
        seq_ids:        Sequence identifiers (same order).
        width:          Motif width in bp.
        top_n:          Number of motifs to find.
        n_seeds:        Diverse seeds to advance to quick EM.
        n_refine:       Quick-EM results to advance to full EM.
        pseudocount:    PWM smoothing value.
        min_hamming:    Minimum Hamming distance between selected seeds.
        mask_threshold: Log-odds threshold for masking (default 0.0 = above background).

    Returns:
        List of motif dicts, each containing:
          consensus (str), ppm (dict), ic (list), total_ic (float),
          log_likelihood (float),
          instances (list of (seq_id, seq, start_1indexed, strand)).
    """
    bg = _compute_background(sequences)
    working_seqs = list(sequences)
    results = []

    print(
        f"  [*] Background: {', '.join(f'{b}={v:.3f}' for b, v in bg.items())}",
        file=sys.stderr,
    )

    for motif_num in range(1, top_n + 1):
        print(f"\n  [*] Searching for Motif {motif_num}/{top_n}...", file=sys.stderr)

        # ── Collect and score all unique k-mer seeds ──────────────────────
        all_kmers: list[str] = []
        for seq in working_seqs:
            for start in range(len(seq) - width + 1):
                kmer = seq[start : start + width]
                if "N" not in kmer:
                    all_kmers.append(kmer)

        if not all_kmers:
            print(
                f"  [!] No valid {width}bp windows remain after masking. "
                f"Stopping at motif {motif_num - 1}.",
                file=sys.stderr,
            )
            break

        seen_kmers: set[str] = set()
        seed_scores: list[tuple[float, str]] = []
        for kmer in all_kmers:
            if kmer in seen_kmers:
                continue
            seen_kmers.add(kmer)
            pwm0 = _init_pwm(kmer, pseudocount)
            ll0 = _log_likelihood(working_seqs, pwm0, bg, width)
            seed_scores.append((ll0, kmer))

        seed_scores.sort(reverse=True)

        # ── Cluster seeds by Hamming distance ─────────────────────────────
        top_seeds = _cluster_seeds(seed_scores, n_seeds, min_hamming)
        n_collapsed = max(0, len(seed_scores) - len(top_seeds))
        print(
            f"      Unique seeds: {len(seed_scores)}  |  "
            f"Diverse (>={min_hamming} Hamming): {len(top_seeds)}  |  "
            f"Collapsed: {n_collapsed}",
            file=sys.stderr,
        )

        # ── Quick EM (25 iterations) to narrow the field ──────────────────
        quick_results: list[tuple[float, dict]] = []
        for kmer in top_seeds:
            pwm0 = _init_pwm(kmer, pseudocount)
            pwm_q, ll_q = _run_em(
                working_seqs,
                pwm0,
                bg,
                width,
                max_iter=25,
                tol=1e-4,
                pseudocount=pseudocount,
            )
            quick_results.append((ll_q, pwm_q))

        quick_results.sort(key=lambda x: x[0], reverse=True)
        refine_pwms = [pwm for _, pwm in quick_results[:n_refine]]
        print(f"      Advancing top {len(refine_pwms)} to full EM", file=sys.stderr)

        # ── Full EM to convergence ────────────────────────────────────────
        best_ll: float = float("-inf")
        best_pwm: dict = {}
        for pwm0 in refine_pwms:
            pwm_f, ll_f = _run_em(
                working_seqs,
                pwm0,
                bg,
                width,
                max_iter=200,
                tol=1e-6,
                pseudocount=pseudocount,
            )
            if ll_f > best_ll:
                best_ll = ll_f
                best_pwm = pwm_f

        if not best_pwm:
            print(
                f"  [!] EM failed to converge for motif {motif_num}.", file=sys.stderr
            )
            break

        # ── Build profile from best-position instances ────────────────────
        instances = _extract_instances(working_seqs, seq_ids, best_pwm, bg, width)
        extracted_seqs = [seq for _, seq, _, _ in instances]
        consensus, ppm, ic = _build_profile_from_instances(extracted_seqs)

        results.append(
            {
                "consensus": consensus,
                "ppm": ppm,
                "ic": ic,
                "total_ic": sum(ic),
                "log_likelihood": best_ll,
                "instances": instances,
            }
        )

        print(
            f"      Converged: LL={best_ll:.3f} | "
            f"Consensus: {consensus} | Total IC: {sum(ic):.3f} bits",
            file=sys.stderr,
        )

        # ── Mask confirmed sites before next search ───────────────────────
        working_seqs, n_masked = _mask_sequences(
            working_seqs, instances, best_pwm, bg, width, threshold=mask_threshold
        )
        print(
            f"      Masked {n_masked}/{len(working_seqs)} sequence(s) "
            f"(threshold={mask_threshold:.1f}).",
            file=sys.stderr,
        )

    return results


# ── CLI ───────────────────────────────────────────────────────────────────────


def main() -> None:
    parser = base_parser("EM Motif Discovery (MEME-style OOPS)")
    parser.add_argument(
        "--top",
        type=int,
        required=True,
        metavar="N",
        help="Number of motifs to discover (required).",
    )
    parser.add_argument(
        "--width",
        type=int,
        default=15,
        metavar="BP",
        help="Motif width in base pairs (Default: 15).",
    )
    parser.add_argument(
        "--seeds",
        type=int,
        default=50,
        metavar="K",
        help="Top K seeds to advance to quick EM (Default: 50).",
    )
    parser.add_argument(
        "--refine",
        type=int,
        default=20,
        metavar="K",
        help="Top K quick-EM results to advance to full EM (Default: 20).",
    )
    parser.add_argument(
        "--iter",
        type=int,
        default=200,
        metavar="N",
        help="Max EM iterations per seed in full EM (Default: 200).",
    )
    parser.add_argument(
        "--pseudocount",
        type=float,
        default=0.1,
        metavar="F",
        help="PWM pseudocount for smoothing (Default: 0.1).",
    )
    parser.add_argument(
        "--min-hamming",
        type=int,
        default=3,
        metavar="D",
        help=(
            "Minimum Hamming distance between any two selected seeds. "
            "Seeds within this distance of a higher-scoring seed are discarded "
            "before EM to prevent redundant runs on shifted variants of the same "
            "motif (seed collision). Default: 3 (~20%% of a 15bp width). "
            "Increase for wider motifs; decrease for very short motifs."
        ),
    )
    parser.add_argument(
        "--mask-threshold",
        type=float,
        default=0.0,
        metavar="F",
        help=(
            "Log-odds score threshold for masking (Default: 0.0). "
            "A position is only masked if its score exceeds this value. "
            "Score > 0 means the subsequence looks more like the motif than "
            "background. Sequences with no plausible site are left unmasked, "
            "protecting their binding sites for subsequent motif discovery. "
            "Increase (e.g., 2.0) for stricter masking."
        ),
    )
    args = parser.parse_args()

    if args.top < 1:
        sys.exit("[!] --top must be a positive integer.")
    if args.width < 1:
        sys.exit("[!] --width must be a positive integer.")

    print(f"[*] EM Motif Discovery", file=sys.stderr)
    print(f"[*] Input          : {args.input.name}", file=sys.stderr)
    print(f"[*] Width          : {args.width}bp", file=sys.stderr)
    print(f"[*] Motifs         : {args.top}", file=sys.stderr)
    print(
        f"[*] Seeds          : {args.seeds} quick -> {args.refine} full",
        file=sys.stderr,
    )
    print(
        f"[*] Min Hamming    : {args.min_hamming} (seed diversity filter)",
        file=sys.stderr,
    )
    print(
        f"[*] Mask threshold : {args.mask_threshold} (log-odds floor for masking)",
        file=sys.stderr,
    )
    print(f"[*] Max iter       : {args.iter}", file=sys.stderr)

    try:
        records = list(SeqIO.parse(args.input, "fasta"))
        if not records:
            sys.exit("[!] Pipeline Halted: FASTA file is empty.")

        sequences = [str(r.seq).upper().replace("U", "T") for r in records]
        seq_ids = [r.id for r in records]

        for seq in sequences:
            if len(seq) < args.width:
                sys.exit(
                    f"[!] Sequence shorter than --width ({len(seq)} < {args.width}bp). "
                    f"Reduce --width."
                )

        print(f"[*] Sequences : {len(sequences)}", file=sys.stderr)

        motifs = discover_motifs(
            sequences,
            seq_ids,
            width=args.width,
            top_n=args.top,
            n_seeds=args.seeds,
            n_refine=args.refine,
            pseudocount=args.pseudocount,
            min_hamming=args.min_hamming,
            mask_threshold=args.mask_threshold,
        )

    except ValueError as e:
        sys.exit(f"\n[!] Pipeline Halted: {e}")
    except FileNotFoundError:
        sys.exit(f"\n[!] Pipeline Halted: Could not find {args.input}")
    except KeyboardInterrupt:
        sys.exit("\n[!] Pipeline interrupted by user.")

    if not motifs:
        sys.exit("[!] No motifs found.")

    n_found = len(motifs)
    col_w = 7
    border = "═" * 64

    # ── Terminal output ───────────────────────────────────────────────────────
    print(f"\n{'=' * 64}")
    print(f"  RESULTS: {n_found} motif(s) discovered")
    print(f"{'=' * 64}")

    for idx, m in enumerate(motifs, 1):
        W = len(m["consensus"])
        print(f"\n{border}")
        print(
            f"  MOTIF {idx}/{n_found}  |  "
            f"Total IC: {m['total_ic']:.3f} bits  |  "
            f"LL: {m['log_likelihood']:.3f}"
        )
        print(border)
        print(f"  Consensus  :  {m['consensus']}\n")

        print("  Position Probability Matrix (PPM):")
        pos_hdr = "".join(f"{i:<{col_w}}" for i in range(1, W + 1))
        print(f"  {'pos':<{col_w}}{pos_hdr}")
        for ch in "ACGT-":
            row = "".join(f"{m['ppm'][ch][i]:<{col_w}.3f}" for i in range(W))
            print(f"  {ch + ':':<{col_w}}{row}")
        ic_row = "".join(f"{v:<{col_w}.3f}" for v in m["ic"])
        print(f"  {'IC:':<{col_w}}{ic_row}")

        print(f"\n  Instances (best position per sequence):")
        max_id = max(len(sid) for sid, _, _, _ in m["instances"])
        for sid, inst_seq, start_1, strand in m["instances"]:
            print(f"    {sid:<{max_id}}  pos {start_1:>4}  strand {strand}  {inst_seq}")

    # ── TSV output ────────────────────────────────────────────────────────────
    if args.output:
        tsv_blocks: list[str] = []

        for idx, m in enumerate(motifs, 1):
            W = len(m["consensus"])
            tsv_blocks.append(
                f"## Motif {idx}/{n_found}"
                f"\tConsensus: {m['consensus']}"
                f"\tTotal IC: {m['total_ic']:.3f} bits"
                f"\tLog-likelihood: {m['log_likelihood']:.3f}"
            )
            pos_row = "\t".join(str(i) for i in range(1, W + 1))
            tsv_blocks.append(f"Position\t{pos_row}")
            tsv_blocks.append(f"Consensus\t{chr(9).join(m['consensus'])}")
            for ch in "ACGT-":
                row = "\t".join(f"{m['ppm'][ch][i]:.3f}" for i in range(W))
                tsv_blocks.append(f"{ch}\t{row}")
            ic_row = "\t".join(f"{v:.3f}" for v in m["ic"])
            tsv_blocks.append(f"IC\t{ic_row}")
            tsv_blocks.append("## Instances")
            max_id = max(len(sid) for sid, _, _, _ in m["instances"])
            for sid, inst_seq, start_1, strand in m["instances"]:
                tsv_blocks.append(f"{sid}\t{start_1}\t{strand}\t{inst_seq}")
            tsv_blocks.append("")

        try:
            args.output.write_text("\n".join(tsv_blocks) + "\n", encoding="utf-8-sig")
            print(
                f"\n[*] TSV written to: {args.output.resolve()}",
                file=sys.stderr,
            )
        except OSError as e:
            sys.exit(f"[!] Could not write to {args.output.name}: {e}")
    else:
        print("\n[*] Note: No -o specified. TSV not saved.", file=sys.stderr)


if __name__ == "__main__":
    main()
