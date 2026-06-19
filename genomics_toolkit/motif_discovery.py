#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

r"""EM Motif Discovery

Discovers transcription factor binding sites and regulatory motifs in upstream
sequences using Expectation-Maximization (MEME-style OOPS/ZOOPS-lite hybrid).

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

    3. Not every sequence necessarily contains the motif. This is the
       difference between OOPS (One Occurrence Per Sequence) and ZOOPS
       (Zero or One Occurrence Per Sequence, Bailey & Elkan's later MEME
       work). Pure OOPS forces every sequence to contribute a full unit of
       "occurrence" probability regardless of whether it has a real site —
       fine for tight prokaryotic operons, but on datasets where many
       sequences genuinely lack the motif (e.g. eukaryotic upstream
       regions with distal/combinatorial regulation), that forced
       contribution is pure noise that can wash out the PWM.

    This script implements a simplified MEME OOPS/ZOOPS-lite hybrid
    (Bailey & Elkan 1994 for OOPS; ZOOPS-lite is a threshold-based
    approximation of the full textbook ZOOPS, which instead learns a
    mixture weight — see ``--zoops-threshold`` and ``_e_step``'s docstring):

      Seed     → For each distinct k-mer in the input sequences, build an
                 initial Position Weight Matrix (PWM) and score it
                 (ZOOPS-lite-aware by default — see --zoops-threshold).
      Cluster  → Before advancing seeds to EM, filter out seeds within
                 --min-hamming Hamming distance of a higher-scoring seed.
                 This prevents redundant EM runs on shifted variants of the
                 same underlying motif (the "seed collision" problem).
      Quick EM → Run 25 EM iterations for the top --seeds diverse candidates
                 (ZOOPS-lite by default — see ZOOPS-LITE section below for
                 why this applies here too, not just in full EM).
      Full EM  → Run --iter EM iterations for the top --refine candidates
                 (ZOOPS-lite by default — see --zoops-threshold/--no-zoops).
      Report   → The best converged motif is reported. Instances whose
                 score didn't clear --zoops-threshold are flagged
                 low-confidence and excluded from the displayed profile.
                 Confirmed instances are masked and the process repeats
                 for the next motif.

    E-step: P(motif at position j | sequence, PWM) ∝ exp(log-odds score),
            or an all-zero vector under ZOOPS-lite if the sequence's best
            score doesn't clear --zoops-threshold.
    M-step: Update PWM using the fractional expected counts (soft EM),
            using whichever strand (forward or reverse complement) gave
            the higher log-odds score at each position.
    Score:  Sum of log-likelihood across all sequences (ZOOPS-lite-aware
            when applicable, for consistency with the E/M-step model).

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

ZOOPS-LITE (protects the PWM on datasets where the motif isn't universal):
    Pure OOPS (One Occurrence Per Sequence) forces every sequence to
    contribute a full unit of "occurrence" probability to the PWM, even
    one with no real site — its softmax over scores still sums to 1.0,
    concentrated on whatever window happens to score relatively highest,
    even if that's just noise. This is a reasonable assumption for tight
    prokaryotic operons, but not for eukaryotic upstream regions, where a
    real fraction of sequences may lack the motif entirely (distal
    enhancers, combinatorial regulation, etc.) — that forced contribution
    is then pure noise diluting the true signal.

    The --zoops-threshold flag (default 0.0) sets a log-odds floor applied
    from seed scoring onward — seed scoring, quick EM, AND full EM all use
    it. (An earlier design held it back to full EM only, reasoning that
    an unrefined PWM couldn't yet reliably tell a real-but-divergent site
    from noise; that turned out to under-protect the result, since pure
    OOPS quick EM still let decoys soften the PWM into a starting point
    full EM's local optimization couldn't escape — see "Strategy" above
    for the empirical detail.) A sequence whose best score across both
    strands and all positions does not clear the threshold contributes
    nothing to that round's PWM update, instead of a forced distribution.
    Reported instances below the threshold are tagged low-confidence and
    excluded from the displayed/exported profile too. Pass --no-zoops to
    disable this and reproduce pre-v1.3.0 pure-OOPS behaviour throughout.

    This is a threshold-based approximation of textbook ZOOPS (Bailey &
    Elkan), not the full version, which instead learns a mixture weight λ
    (the prior probability any given sequence contains the motif at all)
    via its own EM update rule. The threshold approach is simpler and
    reuses the same mental model as --mask-threshold, at the cost of not
    being a formally Bayesian treatment of "does this sequence have a
    site."

PERFORMANCE — zip+Counter COLUMN PROFILING:
    The profile builder in _build_profile_from_instances() uses
    zip(*sequences) + Counter for column-wise counting rather than a
    nested Python loop. Counter runs in C (CPython), substantially reducing
    the number of Python-level operations for large sequence sets.

Note:
    Associated with ongoing, unpublished research (manuscript in
    preparation). Correct attribution is requested when used in
    derivative works.

    v1.2.0 fixes: (1) ``_mask_sequences`` now scores both strands before
    deciding whether to mask, matching the bidirectional scoring used
    everywhere else in this script — previously a reverse-strand-only
    instance could stay permanently unmasked and resurface as a later
    "motif"; (2) sequences with zero valid (non-N) windows now contribute
    zero probability to the M-step instead of a uniform fallback, which
    was injecting artificial PWM weight from invalid, partially-masked
    windows; (3) ``_build_profile_from_instances`` now applies the same
    WebLogo small-sample correction as alignment_conservation_profiler.py,
    keeping the two scripts' IC output identical as documented.

    v1.3.0: Added ZOOPS-lite, a threshold-based approximation of full
    textbook ZOOPS for datasets where a real fraction of sequences may
    lack the motif entirely (e.g. eukaryotic upstream regions). This went
    through several rounds of empirical correction before landing on the
    final design below — earlier iterations seemed reasonable but were
    each disproven by direct testing; see the module-level ZOOPS-LITE
    section and ``discover_motifs``'s docstring for the specifics.

    (1) ``--zoops-threshold``/``--no-zoops``: a sequence whose best score
    doesn't clear the threshold is excluded from that round's PWM update
    entirely, instead of being forced to contribute a full unit of OOPS
    "occurrence" probability built from noise. Applied uniformly from
    seed scoring onward (not held back to full EM only — an earlier
    design that did so was empirically disproven: pure-OOPS quick EM
    still let decoys soften the PWM into a starting point full EM
    couldn't escape).
    (2) ``--zoops-no-calibrate``/``--zoops-seed``: by default, the flat
    threshold is auto-raised per round via permutation calibration
    (shuffle sequences, score against the current PWM, take the 95th
    percentile of the resulting noise distribution). A flat threshold
    alone does not scale with sequence length: pure noise was confirmed
    to clear a flat 0.0 floor 90% of the time on realistic-length
    sequences, since longer sequences offer more candidate windows for
    noise to spuriously exceed any fixed bar by chance.
    (3) Reported instances are tagged low-confidence when they don't
    clear the (calibrated) threshold, and are now also excluded from the
    displayed/exported PPM, consensus, and IC — not just the underlying
    learned PWM — so the profile the user actually reads doesn't get
    silently diluted by decoys even when the model itself is protected.
    (4) Fixed ``--iter`` being silently ignored: full EM was hardcoded to
    200 iterations regardless of this flag.

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

    # Eukaryotic dataset where many sequences may lack the motif: stricter
    # ZOOPS-lite floor so only confident sites inform the PWM
    $ python3 motif_discovery.py -i upstream.fasta --top 3 --zoops-threshold 1.0

    # Reproduce pre-v1.3.0 pure-OOPS behaviour
    $ python3 motif_discovery.py -i upstream.fasta --top 3 --no-zoops

Caveats:
    - Input sequences do NOT need to be aligned (raw upstream FASTA is fine).
    - All sequences must be at least --width bp long.
    - Results are supporting calculations for MEME Suite, not a replacement.
      Always confirm discovered motifs with full MEME analysis.
    - For TSV output: open with File > Open in Excel (not Data > Get Data)
      to avoid the 'Transform Data' wizard.
    - ZOOPS-lite is a threshold-based approximation, not textbook ZOOPS
      (which learns a mixture weight λ via EM rather than using a fixed
      score cutoff). It only applies during full EM, not seed scoring or
      quick EM.
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.3.0"

import math
import random
import sys
from collections import Counter

try:
    from Bio import SeqIO
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )
from utils import base_parser, revcomp

# ── Reverse complement ────────────────────────────────────────────────────────

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


def _e_step(
    sequences: list[str],
    pwm: dict,
    bg: dict,
    width: int,
    zoops_threshold: float = float("-inf"),
) -> list[list[float]]:
    """E-step: compute P(motif at position j | sequence, PWM) via softmax.

    Scans BOTH strands at every position. At each start position j, the
    forward window and its reverse complement are both scored; the higher
    log-odds score is used as the position's score. This correctly captures
    TF binding sites regardless of which strand they appear on.

    For each sequence, all valid start positions are scored, then normalised
    to a probability distribution via softmax (numerically stable: subtract
    max before exponentiating).

    ZOOPS-lite: this is OOPS (One Occurrence Per Sequence) by default
    (``zoops_threshold`` left at ``-inf``, its every-sequence probabilities
    always sum to 1.0). Passing a real ``zoops_threshold`` makes it
    ZOOPS-lite (Zero or One Occurrence Per Sequence): if a sequence's best
    score (across both strands, across every position) does not exceed
    the threshold, the entire sequence is given an all-zero probability
    vector instead of a forced softmax distribution — it contributes
    nothing to the following M-step, rather than injecting a full unit of
    "occurrence" mass built entirely from background-level noise. This
    protects the PWM on datasets (e.g. eukaryotic upstream regions) where
    a real fraction of sequences genuinely lack the motif. This is a
    threshold-based approximation of textbook ZOOPS, not the full
    EM-learned-mixture-weight (λ) version — see the module docstring.

    Args:
        sequences:       Input sequences (may contain 'N' at masked positions).
        pwm:             Current PWM (dict ``"ACGT"`` → list[float]).
        bg:              Background model.
        width:           Motif width.
        zoops_threshold: Log-odds floor for ZOOPS-lite participation. Default
                          ``-inf`` disables it (pure OOPS — every sequence
                          always contributes). Pass a real value (e.g. 0.0,
                          break-even with background) to enable ZOOPS-lite.

    Returns:
        List of probability lists, one per sequence. A sequence whose
        every window touches a masked 'N' (no valid window at all), or
        whose best score does not clear ``zoops_threshold``, gets an
        all-zero list rather than a forced distribution — so it
        contributes nothing to the following M-step.
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
            rev = _score_window(revcomp(window), pwm, bg)
            log_scores.append(max(fwd, rev))

        # Softmax with -inf handling. If literally every window in this
        # sequence touches a masked 'N' (possible when the sequence is
        # short relative to `width` and heavily masked from prior motif
        # rounds), there is no legitimate position to assign probability
        # to. Assigning a uniform fallback here (as opposed to zero) would
        # inject artificial M-step weight onto whichever non-N bases
        # happen to sit in these otherwise-invalid windows, contaminating
        # the PWM with noise from a sequence that has no real signal this
        # round. Zero ensures `_m_step`'s `if z < 1e-15: continue` guard
        # correctly skips this sequence entirely.
        valid = [s for s in log_scores if s > float("-inf")]
        if not valid:
            all_probs.append([0.0] * n_pos)
            continue

        max_s = max(valid)

        # ZOOPS-lite: a sequence whose single best window doesn't even
        # beat background-level odds is treated as not containing the
        # motif at all, rather than forced to contribute a full unit of
        # "occurrence" probability built from noise.
        if max_s <= zoops_threshold:
            all_probs.append([0.0] * n_pos)
            continue

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
            rc_window = revcomp(window)
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


def _log_likelihood(
    sequences: list[str],
    pwm: dict,
    bg: dict,
    width: int,
    zoops_threshold: float = float("-inf"),
) -> float:
    """Total log-likelihood of all sequences under the current PWM.

    Scores both forward and reverse complement at each position (taking the
    max), consistent with the bidirectional E-step. Uses log-sum-exp for
    numerical stability.

    ZOOPS-lite (see ``_e_step``): a sequence whose best score does not
    clear ``zoops_threshold`` contributes 0 rather than the OOPS
    `max_s + log(sum(exp))` term, keeping this evaluation metric
    consistent with the E/M-step model used to reach this PWM. Default
    ``-inf`` disables this (pure OOPS).
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
            rev = _score_window(revcomp(window), pwm, bg)
            log_scores.append(max(fwd, rev))
        valid = [s for s in log_scores if s > float("-inf")]
        if not valid:
            continue
        max_s = max(valid)
        if max_s <= zoops_threshold:
            continue
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
    zoops_threshold: float = float("-inf"),
) -> tuple[dict, float]:
    """Run EM to convergence from an initial PWM.

    Args:
        sequences:       Input (possibly masked) sequences.
        initial_pwm:     Starting PWM.
        bg:              Background model.
        width:           Motif width.
        max_iter:        Maximum number of E/M cycles.
        tol:             Convergence threshold (change in log-likelihood).
        pseudocount:     M-step smoothing value.
        zoops_threshold: ZOOPS-lite participation floor, forwarded to
                          ``_e_step``/``_log_likelihood`` every iteration.
                          Default ``-inf`` = pure OOPS (recommended for
                          quick/seed-ranking EM, where the PWM hasn't
                          refined enough yet to reliably distinguish real
                          sites from noise — see ``discover_motifs``).

    Returns:
        (converged_pwm, final_log_likelihood)
    """
    pwm = {b: initial_pwm[b][:] for b in "ACGT"}
    prev_ll = float("-inf")

    for _ in range(max_iter):
        probs = _e_step(sequences, pwm, bg, width, zoops_threshold=zoops_threshold)
        pwm = _m_step(sequences, probs, width, pseudocount, pwm, bg)
        ll = _log_likelihood(sequences, pwm, bg, width, zoops_threshold=zoops_threshold)
        if ll - prev_ll < tol:
            break
        prev_ll = ll

    return pwm, _log_likelihood(
        sequences, pwm, bg, width, zoops_threshold=zoops_threshold
    )


# ── Post-EM: extract instances and mask ───────────────────────────────────────


def _extract_instances(
    sequences: list[str],
    seq_ids: list[str],
    pwm: dict,
    bg: dict,
    width: int,
) -> list[tuple[str, str, int, str, float]]:
    """Find the best-scoring position of the motif in each sequence.

    Scores both forward and reverse complement at every position and picks
    the overall best. Reports which strand produced the winning score, and
    the score itself (needed by callers to flag low-confidence instances
    under ZOOPS-lite — a sequence with no real site still gets a "best"
    position reported here, since this function always reports *something*,
    but its score will typically sit at or below background).

    Returns:
        List of (seq_id, extracted_sequence, start_position_1indexed,
        strand, score) where strand is '+' (forward) or '-' (reverse
        complement). A sequence too short for even one window gets a
        score of ``-inf``.
    """
    instances = []
    for seq, sid in zip(sequences, seq_ids):
        n_pos = len(seq) - width + 1
        if n_pos <= 0:
            instances.append(
                (sid, seq[:width].ljust(width, "N"), 1, "+", float("-inf"))
            )
            continue

        best_j = 0
        best_score = float("-inf")
        best_strand = "+"

        for j in range(n_pos):
            window = seq[j : j + width]
            fwd_score = _score_window(window, pwm, bg)
            rev_score = _score_window(revcomp(window), pwm, bg)
            if rev_score > fwd_score:
                score, strand = rev_score, "-"
            else:
                score, strand = fwd_score, "+"
            if score > best_score:
                best_score = score
                best_j = j
                best_strand = strand

        instances.append(
            (sid, seq[best_j : best_j + width], best_j + 1, best_strand, best_score)
        )
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

    Scores BOTH strands at the instance position and only masks if the
    better of the two exceeds ``threshold`` — consistent with the
    bidirectional scoring used everywhere else in this script (E-step,
    M-step, instance extraction). Checking only the forward strand here
    would silently fail to mask any instance whose best score came from
    the reverse complement (a common case, since TF binding sites are not
    strand-specific): the forward-strand score at that same position is
    often far below threshold even when the true site, on the other
    strand, scored very high. An unmasked reverse-strand site can then
    resurface and get reported again as a later motif.

    This prevents the OOPS masking danger: if a sequence does not actually
    contain Motif 1, the forced "best" position may coincide with a real
    Motif 2 binding site. Without a threshold, that Motif 2 site gets
    masked unconditionally, corrupting subsequent motif discovery. With a
    threshold, sequences where the best Motif 1 score (either strand) is
    below background (score ≤ threshold) are left unmasked.

    A threshold of 0.0 means "only mask if the subsequence looks more like
    the motif than background." Increase (e.g., 2.0) for stricter masking.

    'N' characters score -inf in ``_score_window``, making them invisible
    to the EM in subsequent motif searches.

    Args:
        sequences: Current (possibly already masked) working sequences.
        instances: Tuples whose first 4 elements are (seq_id, seq,
                   start_1indexed, strand) — either raw ``_extract_instances``
                   output or ``discover_motifs``'s ZOOPS-lite-tagged
                   6-tuples; any trailing elements are ignored here.
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
    for seq, (_, _, start_1, _, _, _) in zip(sequences, instances):
        start_0 = start_1 - 1
        window = seq[start_0 : start_0 + width]
        fwd_score = _score_window(window, pwm, bg)
        rev_score = _score_window(revcomp(window), pwm, bg)
        score = max(fwd_score, rev_score)
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
    but inlined here so this script is standalone — including the WebLogo
    small-sample correction term e_n = (s-1)/(2 ln(2) n); the two must stay
    in sync, since this docstring promises identical output.
    """
    seq_len = len(extracted[0])
    n = len(extracted)
    valid = "ACGT-"
    counts = {c: [0] * seq_len for c in valid}
    max_bits = math.log2(4.0)
    small_sample_numerator = 3.0  # (s - 1) for s = 4 DNA symbols

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
        e_n = small_sample_numerator / (2 * math.log(2) * acgt) if acgt > 0 else 0.0
        ic[i] = max(0.0, (max_bits + ent - e_n) * fp)

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


# ── ZOOPS-lite threshold calibration ──────────────────────────────────────────

# Internal tuning constants for _calibrate_zoops_threshold. Not exposed as
# CLI flags to keep --zoops-threshold/--no-zoops the only knobs a user
# needs to reason about; change here if you need to retune the
# speed/precision tradeoff.
_ZOOPS_CALIBRATION_PERCENTILE = 95.0
_ZOOPS_CALIBRATION_SHUFFLES_PER_SEQ = 3
_ZOOPS_CALIBRATION_MAX_SEQS = 60


def _calibrate_zoops_threshold(
    sequences: list[str],
    pwm: dict,
    bg: dict,
    width: int,
    floor: float,
    rng: random.Random,
) -> float:
    """Estimates a length- and PWM-appropriate ZOOPS-lite floor by permutation.

    A flat log-odds threshold does not account for how many candidate
    windows a sequence actually offers: a 200bp sequence gives a width-15
    PWM roughly 372 chances (186 positions x 2 strands) to find SOME
    window that happens to score well purely by chance, even with no real
    site present. Confirmed empirically: pure-random 85bp sequences
    (81 positions x 2 strands = 162 chances) cleared a flat 0.0 threshold
    90% of the time against a realistic PWM. A fixed threshold tuned for
    short sequences is therefore not just imprecise but actively
    non-functional on realistic-length upstream sequences.

    This shuffles each sampled sequence (composition-preserving,
    structure-destroying — the standard null model: "what would this
    sequence's best score be if its bases were in random order"),
    scores the shuffled sequence's best window bidirectionally against
    ``pwm``, and returns the requested percentile of that null-score
    distribution — i.e., "the score level that pure noise of this exact
    length only exceeds 5% of the time" — floored at the user-supplied
    ``--zoops-threshold`` so an explicit user choice is never lowered,
    only potentially raised.

    For datasets with many sequences, only a random subsample (capped at
    ``_ZOOPS_CALIBRATION_MAX_SEQS``) is shuffled, to bound runtime; this
    sample size is still large enough for a stable percentile estimate at
    realistic dataset sizes.

    Args:
        sequences: The current working sequences (post-masking) to draw
                   length/composition statistics from.
        pwm:       Reference PWM to score shuffled sequences against (the
                   top seed's raw PWM before quick EM, or the top
                   quick-EM PWM before full EM — see ``discover_motifs``).
        bg:        Background model.
        width:     Motif width.
        floor:     The user-supplied ``--zoops-threshold``; the result
                   never goes below this value.
        rng:       A ``random.Random`` instance (caller controls the seed
                   for reproducibility — see ``--zoops-seed``).

    Returns:
        ``max(floor, calibrated_percentile)``. Falls back to ``floor``
        alone if no sequence offered even one scoreable window.
    """
    sample_seqs = (
        sequences
        if len(sequences) <= _ZOOPS_CALIBRATION_MAX_SEQS
        else rng.sample(sequences, _ZOOPS_CALIBRATION_MAX_SEQS)
    )

    null_scores: list[float] = []
    for seq in sample_seqs:
        chars = list(seq)
        for _ in range(_ZOOPS_CALIBRATION_SHUFFLES_PER_SEQ):
            rng.shuffle(chars)
            shuffled = "".join(chars)
            n_pos = len(shuffled) - width + 1
            if n_pos <= 0:
                continue
            best = float("-inf")
            for j in range(n_pos):
                window = shuffled[j : j + width]
                if "N" in window:
                    continue
                fwd = _score_window(window, pwm, bg)
                rev = _score_window(revcomp(window), pwm, bg)
                s = max(fwd, rev)
                if s > best:
                    best = s
            if best > float("-inf"):
                null_scores.append(best)

    if not null_scores:
        return floor

    null_scores.sort()
    idx = min(
        len(null_scores) - 1,
        int(len(null_scores) * _ZOOPS_CALIBRATION_PERCENTILE / 100.0),
    )
    calibrated = null_scores[idx]
    return max(floor, calibrated)


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
    zoops_threshold: float | None = 0.0,
    n_iter: int = 200,
    zoops_calibrate: bool = True,
    zoops_seed: int = 0,
) -> list[dict]:
    """Run EM motif discovery: find top N non-overlapping motifs.

    Strategy:
      1. Collect all unique W-mers from all input sequences as candidate seeds.
      2. Score each seed by initial log-likelihood (bidirectional;
         ZOOPS-lite-aware if ``zoops_threshold`` is set).
      3. Cluster seeds by Hamming distance — only advance seeds that differ
         by at least ``min_hamming`` from all higher-scoring seeds already
         selected. This ensures the seed pool is physically diverse.
      4. Quick EM (25 iterations) on the top ``n_seeds`` diverse seeds.
      5. Full EM (``--iter`` iterations) on the top ``n_refine`` from step 4.
      6. Report the best converged motif. Instances whose best score does
         not clear ``zoops_threshold`` are flagged as low-confidence
         (sequence likely doesn't contain this motif at all), and excluded
         from the displayed/exported PPM and consensus too.
      7. Mask confirmed instances (score > ``mask_threshold``) with 'N' and
         repeat for the next motif.

    ZOOPS-lite is applied uniformly from step 2 onward (seed scoring,
    quick EM, AND full EM), not just full EM. An earlier design held it
    back to full EM only, reasoning that an unrefined single-seed PWM
    would be too coarse to safely separate a real-but-divergent site from
    noise this early. That reasoning turned out to be only half right:
    delaying ZOOPS-lite does NOT protect the final result, because pure
    OOPS quick EM still force-includes decoy sequences (those genuinely
    lacking the motif), softening the PWM it hands to full EM into a
    starting point where decoys already clear the threshold at iteration
    1 — a stable bad equilibrium that full EM's local optimization cannot
    escape on its own (confirmed empirically: a decoy scoring -6.5
    against a from-scratch ZOOPS-protected PWM still scored +3.6 against
    the same dataset's pure-OOPS-quick-EM result, and stayed there through
    200 more ZOOPS-lite full-EM iterations). Applying the threshold from
    the start was then tested against the original concern directly — a
    real site 2/5 bp diverged from the seed — and did not exhibit
    premature exclusion: even fairly divergent real sites scored well
    above a 0.0 floor (3.4-6.6 vs the floor of 0.0) and were retained.

    Args:
        sequences:       Input sequences (>= width bp each).
        seq_ids:         Sequence identifiers (same order).
        width:           Motif width in bp.
        top_n:           Number of motifs to find.
        n_seeds:         Diverse seeds to advance to quick EM.
        n_refine:        Quick-EM results to advance to full EM.
        pseudocount:     PWM smoothing value.
        min_hamming:     Minimum Hamming distance between selected seeds.
        mask_threshold:  Log-odds threshold for masking (default 0.0 = above
                          background). Independent of ``zoops_threshold`` —
                          this one decides whether to destroy a confirmed
                          site before searching for the next motif; that one
                          decides whether a sequence participates in
                          learning the *current* motif's PWM at all.
        zoops_threshold: Log-odds floor for ZOOPS-lite participation,
                          applied from seed scoring onward (default 0.0 =
                          break-even with background). Pass ``None`` to
                          disable ZOOPS-lite entirely and run pure OOPS
                          throughout, matching pre-v1.3.0 behaviour.
        n_iter:          Max EM iterations per candidate in full EM
                          (default 200). Quick EM is always a fixed 25
                          iterations regardless of this value, since its
                          purpose is fast triage, not convergence.
        zoops_calibrate: If True (default) and ``zoops_threshold`` is not
                          ``None``, automatically raise the effective
                          threshold per round via permutation calibration
                          — see ``_calibrate_zoops_threshold``. A flat
                          threshold alone does not scale with sequence
                          length (confirmed empirically: pure noise
                          cleared a flat 0.0 floor 90% of the time on
                          realistic-length sequences), so this is on by
                          default. ``zoops_threshold`` still acts as a
                          floor — calibration only ever raises it.
        zoops_seed:      Seed for the calibration shuffling RNG (default
                          0), so results are exactly reproducible across
                          runs on the same data. Only used when
                          ``zoops_calibrate`` is True.

    Returns:
        List of motif dicts, each containing:
          consensus (str), ppm (dict), ic (list), total_ic (float),
          log_likelihood (float),
          instances (list of (seq_id, seq, start_1indexed, strand, score,
          low_confidence)), where ``low_confidence`` is True if the
          instance's score did not clear ``zoops_threshold``.
          consensus/ppm/ic are built ONLY from high-confidence instances
          when ZOOPS-lite is enabled (falling back to every instance only
          if literally all of them were flagged) — this keeps the
          reported profile consistent with the learned PWM; including
          low-confidence (likely-decoy) instances in the displayed
          profile would silently re-dilute it even though the underlying
          EM never let them influence ``best_pwm``.
    """
    full_em_zoops = float("-inf") if zoops_threshold is None else zoops_threshold
    zoops_rng = random.Random(zoops_seed)
    bg = _compute_background(sequences)
    working_seqs = list(sequences)
    results = []

    print(
        f"  [*] Background: {', '.join(f'{b}={v:.3f}' for b, v in bg.items())}",
        file=sys.stderr,
    )
    if zoops_threshold is not None:
        print(
            f"  [*] ZOOPS-lite enabled (seed scoring + quick EM + full EM): "
            f"floor={zoops_threshold:.2f}"
            + (
                ", auto-calibrated per round above the noise floor"
                if zoops_calibrate
                else " (calibration disabled, using this flat value throughout)"
            ),
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
            ll0 = _log_likelihood(
                working_seqs, pwm0, bg, width, zoops_threshold=full_em_zoops
            )
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

        # ── Quick EM (25 iterations) to narrow the field ───────────────────
        # ZOOPS-lite is applied here too, not just in full EM. Empirically,
        # delaying it to full EM only does NOT protect the final PWM: pure
        # OOPS quick EM still force-includes decoy sequences, which softens
        # the PWM it hands to full EM into a starting point where decoys
        # already clear the threshold at iteration 1 — a stable bad
        # equilibrium that full EM's local optimization cannot escape from
        # on its own (verified by direct test: a decoy that scores -6.5
        # against a from-scratch ZOOPS-protected PWM still scored +3.6
        # against the same dataset's pure-OOPS quick-EM result, and stayed
        # there through 200 more ZOOPS-lite full-EM iterations). Applying
        # the threshold from iteration 1 instead was tested against a
        # harder case too (a real site 2/5 bp diverged from the seed) and
        # did not exhibit the bootstrapping failure mode that motivated
        # holding it back in the first place: even fairly divergent real
        # sites scored well above a 0.0 floor and were retained.
        #
        # Quick EM uses the FLAT floor here, not a calibrated one.
        # Calibrating against the raw single-kmer seed PWM was tried and
        # confirmed counterproductive: an unrefined PWM (pseudocount
        # smoothing only, ~0.79 vs ~0.07 per base) is not yet sharp enough
        # to separate true sites from noise AT ALL, so a calibrated
        # "noise ceiling" against it lands right on top of a true site's
        # own achievable score, rejecting real sites during quick EM
        # before they ever get a chance to refine the PWM. The flat floor
        # is lenient enough to admit real sites while still rejecting
        # clearly-negative noise (as confirmed in the original direct
        # test above). Calibration is deferred to full EM (below), once
        # quick EM has had a chance to actually sharpen the PWM.
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
                zoops_threshold=full_em_zoops,
            )
            quick_results.append((ll_q, pwm_q))

        quick_results.sort(key=lambda x: x[0], reverse=True)
        refine_pwms = [pwm for _, pwm in quick_results[:n_refine]]
        print(f"      Advancing top {len(refine_pwms)} to full EM", file=sys.stderr)

        # ── Full EM to convergence (ZOOPS-lite if zoops_threshold is set) ──
        # Calibrated here (not before quick EM) against the best
        # quick-EM-refined PWM, which by now is sharp enough for the
        # noise-vs-signal separation to be statistically meaningful —
        # confirmed empirically: real sites scored 12.7 against a
        # quick-EM-refined PWM with the calibrated threshold landing at
        # 4.3, comfortably between the true signal and the decoy scores
        # (-7.1 to 4.4).
        full_em_zoops_calibrated = full_em_zoops
        if zoops_threshold is not None and zoops_calibrate and refine_pwms:
            full_em_zoops_calibrated = _calibrate_zoops_threshold(
                working_seqs, refine_pwms[0], bg, width, full_em_zoops, zoops_rng
            )
            print(
                f"      ZOOPS-lite calibrated for full EM: "
                f"floor={full_em_zoops:.2f} \u2192 {full_em_zoops_calibrated:.2f}",
                file=sys.stderr,
            )

        best_ll: float = float("-inf")
        best_pwm: dict = {}
        for pwm0 in refine_pwms:
            pwm_f, ll_f = _run_em(
                working_seqs,
                pwm0,
                bg,
                width,
                max_iter=n_iter,
                tol=1e-6,
                pseudocount=pseudocount,
                zoops_threshold=full_em_zoops_calibrated,
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
        raw_instances = _extract_instances(working_seqs, seq_ids, best_pwm, bg, width)
        # Tag instances whose best score didn't clear the ZOOPS-lite
        # threshold ACTUALLY USED for this round's full EM (the
        # calibrated value, not the flat user-supplied floor) — this
        # sequence almost certainly does not contain the motif, and the
        # reported position is just "the best of a bad lot" rather than
        # a real site. (When zoops_threshold is None, full_em_zoops is
        # -inf and full_em_zoops_calibrated stays -inf too, so nothing
        # gets flagged.)
        instances = [
            (sid, seq, start_1, strand, score, score <= full_em_zoops_calibrated)
            for sid, seq, start_1, strand, score in raw_instances
        ]
        n_low_confidence = sum(1 for *_, low_conf in instances if low_conf)
        # Build the DISPLAYED profile from only the high-confidence
        # instances. Using every extracted instance here (including
        # sequences already excluded from the EM-learned best_pwm itself)
        # would silently re-introduce exactly the noise contamination
        # ZOOPS-lite exists to prevent — the underlying model would be
        # clean, but the reported PPM/consensus/IC the user actually
        # reads would still be diluted by decoys' meaningless "best
        # guess" positions. Falls back to every instance only in the
        # degenerate case where ALL of them were flagged (so there's
        # nothing else to build a profile from).
        high_conf_seqs = [
            seq
            for sid, seq, start_1, strand, score, low_conf in instances
            if not low_conf
        ]
        extracted_seqs = (
            high_conf_seqs
            if high_conf_seqs
            else [seq for sid, seq, start_1, strand, score, low_conf in instances]
        )
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
        if zoops_threshold is not None and n_low_confidence:
            print(
                f"      {n_low_confidence}/{len(instances)} instance(s) flagged "
                f"low-confidence (score did not clear ZOOPS-lite threshold "
                f"{full_em_zoops_calibrated:.2f}) \u2014 likely absent from this "
                f"sequence; excluded from both the learned PWM and the "
                f"profile/IC above.",
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
    parser = base_parser("EM Motif Discovery (MEME-style OOPS/ZOOPS-lite)")
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
    parser.add_argument(
        "--zoops-threshold",
        type=float,
        default=0.0,
        metavar="F",
        help=(
            "ZOOPS-lite log-odds floor (Default: 0.0 = break-even with "
            "background), applied to seed scoring, quick EM, and full EM "
            "alike. A sequence whose best score (either strand, any "
            "position) does not clear this value is excluded from the PWM "
            "update entirely that round, instead of being forced to "
            "contribute a full unit of OOPS 'occurrence' probability built "
            "from noise. Use this on datasets where many sequences may "
            "genuinely lack the motif (e.g. eukaryotic upstream regions "
            "with distal/combinatorial regulation). Independent from "
            "--mask-threshold, which governs a separate decision (whether "
            "to destroy a confirmed site before searching for the next "
            "motif). This is a threshold-based approximation of textbook "
            "ZOOPS, not the full EM-learned-mixture-weight version."
        ),
    )
    parser.add_argument(
        "--no-zoops",
        action="store_true",
        help=(
            "Disable ZOOPS-lite and run pure OOPS throughout (every "
            "sequence is always forced to contribute, matching pre-v1.3.0 "
            "behaviour). Overrides --zoops-threshold."
        ),
    )
    parser.add_argument(
        "--zoops-no-calibrate",
        action="store_true",
        help=(
            "Use --zoops-threshold as a flat, literal value with no "
            "per-dataset auto-calibration. Calibration is on by default "
            "because a flat threshold does not scale with sequence length "
            "— pure noise was confirmed to clear a flat 0.0 floor 90%% of "
            "the time on realistic-length sequences, since longer "
            "sequences offer more candidate windows for noise to "
            "spuriously exceed any fixed bar by chance. Only disable this "
            "if you have already determined an appropriate threshold for "
            "your specific data (e.g. via your own permutation testing) "
            "and want exact control."
        ),
    )
    parser.add_argument(
        "--zoops-seed",
        type=int,
        default=0,
        metavar="N",
        help=(
            "Random seed for ZOOPS-lite calibration shuffling (Default: 0). "
            "Fixed by default so results are exactly reproducible across "
            "runs on the same data. Only used when calibration is active."
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
    if args.no_zoops:
        print("[*] ZOOPS-lite     : disabled (pure OOPS, --no-zoops)", file=sys.stderr)
    elif args.zoops_no_calibrate:
        print(
            f"[*] ZOOPS-lite     : enabled, flat threshold={args.zoops_threshold} "
            f"(calibration disabled, seeding + quick EM + full EM)",
            file=sys.stderr,
        )
    else:
        print(
            f"[*] ZOOPS-lite     : enabled, floor={args.zoops_threshold}, "
            f"auto-calibrated per round (seed={args.zoops_seed})",
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
            zoops_threshold=None if args.no_zoops else args.zoops_threshold,
            n_iter=args.iter,
            zoops_calibrate=not args.zoops_no_calibrate,
            zoops_seed=args.zoops_seed,
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
        max_id = max(len(sid) for sid, _, _, _, _, _ in m["instances"])
        for sid, inst_seq, start_1, strand, score, low_conf in m["instances"]:
            tag = "  [low-confidence: likely no real site]" if low_conf else ""
            print(
                f"    {sid:<{max_id}}  pos {start_1:>4}  strand {strand}  "
                f"{inst_seq}  (score {score:.2f}){tag}"
            )

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
            tsv_blocks.append(
                "Seq_ID\tPosition\tStrand\tSequence\tScore\tLow_Confidence"
            )
            for sid, inst_seq, start_1, strand, score, low_conf in m["instances"]:
                tsv_blocks.append(
                    f"{sid}\t{start_1}\t{strand}\t{inst_seq}\t{score:.3f}\t"
                    f"{'yes' if low_conf else 'no'}"
                )
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
