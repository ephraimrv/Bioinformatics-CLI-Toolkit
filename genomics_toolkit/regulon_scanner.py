#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

"""Genome-Wide Regulon Mapper

Maps transcriptional networks by identifying operator motifs in upstream regions.

The pipeline isolates the upstream sequence (default 150bp) of every CDS in
a GenBank assembly. It performs a regex-based motif search for a provided
IUPAC/Regex operator footprint and compiles matches into a genomic matrix
suitable for network analysis.

Both strands are scanned independently. Motif positions are reported as
negative integers relative to the Translation Start Site (TSS), following
standard molecular biology convention (e.g., -10 and -35 boxes). The motif
strand column (+/-) indicates which DNA strand the binding site was found on.

STATISTICAL MODEL — Regex-Based vs PWM-Based Scoring:
    This script uses combinatorial probability: P-values are calculated as the
    probability of a random sequence matching the defined regex pattern, given
    empirical genomic nucleotide frequencies. This contrasts with MEME/FIMO,
    which use Position-Weight Matrices (PWMs) and report continuous match scores
    for every position in the sequence.

    Regex model (this script):
    - Binary classification: match or no-match
    - P-value = product of character class probabilities (fixed base or [ACGT])
    - Reports only hits with statistically significant q-values (FDR control)
    - Implicit threshold: q-value < 0.05 (or user-specified α)

    PWM model (FIMO):
    - Continuous scoring: reports match scores for every position
    - P-value = tail probability of the score distribution (requires training)
    - Reports all matches, ranked by significance
    - User manually filters results by p-value threshold (no default)

WHY RESULTS DIFFER:
    FIMO returns MANY results (often thousands) because it reports every position
    with a match score, even weak/insignificant ones. The user must manually
    filter by p-value. This script returns only hits with q-values passing
    Benjamini-Hochberg FDR correction at α=0.05, automatically discarding noise
    from the massive genome-wide search space.

    If you compare this script's results to FIMO's unfiltered output:
    - FIMO will have ~100-1000× more hits
    - This script's hits should be a subset of FIMO's significant hits (p < 0.05)
    - Occasional mismatches (different motif, position off-by-one) are expected
      due to different scoring models

INTENDED USE:
    For prokaryotic regulatory network discovery where the motif footprint is
    known or suspected (lacticin operator boxes, SigmaA boxes, etc.) and binary
    classification (matches the pattern or doesn't) is appropriate.

PROKARYOTE-ONLY — NOT A WINDOW-SIZE PROBLEM:
    This script anchors every upstream window on the CDS start (the
    translation start / ATG), not on the Transcription Start Site (TSS).
    In prokaryotes these coincide, since there is no 5' UTR separating
    them. In eukaryotes they do not: the TSS sits upstream of the CDS
    start, often separated by a 5' UTR that itself contains introns.

    Increasing --upstream on a eukaryotic genome does NOT fix this — it
    just extracts a longer stretch of 5' UTR/intron sequence anchored at
    the wrong coordinate, not the actual promoter. There is currently no
    eukaryote mode here (unlike universal_promoter_extractor.py, which
    resolves the TSS from mRNA features across isoforms). For eukaryotic
    regulatory motif discovery: extract upstream regions with
    universal_promoter_extractor.py or target_promoter_pipeline.py (both
    TSS-anchored), then search that output with MEME/FIMO directly.

Note:
    Associated with ongoing, unpublished research (manuscript in
    preparation). Correct attribution is requested when used in
    derivative works.

    v1.4.0: ``_compute_background``, ``_count_total_windows``, and
    ``stream_regulon_hits`` each independently reimplemented the same
    strand-aware upstream-slicing arithmetic (forward:
    ``slice_start = max(0, start - upstream_bp)``; reverse:
    ``slice_end = min(len(seq), end + upstream_bp)`` + reverse-complement)
    — three copies in this file alone, on top of further copies in
    ``universal_promoter_extractor.py`` and
    ``utils.extract_upstream_sequence``. All three now delegate to
    ``utils.resolve_upstream_window()`` / ``utils.extract_upstream_window()``
    (see ``utils.py``'s v1.2.0 changelog note for the full rationale).
    ``_count_total_windows`` uses the arithmetic-only
    ``resolve_upstream_window()`` since it never needs the actual bases;
    the other two use ``extract_upstream_window()``. No behavior change —
    same slice boundaries, same truncation handling, same genomic
    coordinates for motif hits.

    v1.5.0: Fixed two issues found while validating criticism of the
    sibling script gbk_promoter_finder.py, which advertises the same
    "IUPAC/regex motif" support and had the identical bug.
    (1) Raw IUPAC ambiguity codes (W, R, Y, S, K, M, B, D, H, V, N) were
    passed directly to re.compile(), which has no concept of them —
    "TATAWAW" matched zero times against a
    sequence containing the valid instance "TATAAAA". Now translated via
    utils.translate_iupac_to_regex() before compiling.
    (2) A second, independent bug was found layered on top of the first:
    _motif_pvalue() already correctly handles explicit "[TC]"-style
    bracket notation, but did not recognize raw IUPAC letters at all —
    it fell through to its "complex token, no constraint" branch for
    them. The same motif spelled "TATAWAW" vs
    "TATA[AT]A[AT]" produced two different p-values (0.0024 vs 0.00088)
    purely from spelling, with the raw-IUPAC spelling understating
    significance. Fixed by translating the motif ONCE in main() and using
    that single translated string for the search, the motif-length
    estimate, and the p-value calculation — _motif_pvalue() itself didn't
    need to change, since it already does the right thing given explicit
    bracket notation.
    Also replaced the shared "UNKNOWN" default for CDS features lacking
    /locus_tag with a coordinate-based fallback (matching the format used
    in pairwise_homolog_finder.py and universal_promoter_extractor.py) —
    this script aggregates no data by locus_tag, so there was no data-loss
    risk, but multiple untagged genes previously showed up identically as
    "UNKNOWN" in the output with no way to distinguish them.
    These fixes are unrelated to and fully compatible with v1.4.0's
    slicing-arithmetic consolidation above — v1.5.0 builds on top of it
    rather than replacing it.

Example:
    $ python3 regulon_scanner.py -i C5_genome.gbk -u 200 -m "GCGCAG[CT]G[GT]T[TA]AAAT" -o regulon.tsv
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.5.0"

import re
import sys
import argparse
import traceback
from pathlib import Path
from typing import Iterator
from collections import Counter

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )

from utils import (
    base_parser,
    resolve_upstream_window,
    extract_upstream_window,
    translate_iupac_to_regex,
)


def _compute_background(gbk_path: Path, upstream_bp: int) -> dict[str, float]:
    """Compute empirical nucleotide background frequencies from all upstream sequences.

    Scans every CDS upstream region in the genome and tallies base frequencies.
    Used to compute positional p-values for motif matches. A minimum floor of
    0.01 is applied to prevent division-by-zero on low-frequency bases.

    Uses collections.Counter for the per-sequence tally (C-backed, faster than
    a manual character loop at genome scale). Counter.update() will also tally
    any non-ACGT characters present (e.g. 'N' runs in draft assemblies), but
    the final total and per-base frequencies below restrict to A/C/G/T only —
    matching the original behaviour of silently excluding ambiguous bases from
    both the numerator and the denominator.

    The strand-aware upstream slice itself is delegated to
    ``utils.extract_upstream_window()`` rather than reimplemented here —
    see the module changelog (v1.4.0).

    Args:
        gbk_path:    Path to the GenBank file.
        upstream_bp: Upstream window size (must match the scan window).

    Returns:
        Dict mapping "ACGT" → frequency, summing to ~1.0.
    """
    counts: Counter = Counter()
    for record in SeqIO.parse(gbk_path, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue
            start = int(feature.location.start)
            end = int(feature.location.end)
            strand = feature.location.strand
            upstream, _actual_upstream, _slice_start, _slice_end = (
                extract_upstream_window(record, start, end, strand, upstream_bp)
            )
            counts.update(upstream.upper())
    total = sum(counts[b] for b in "ACGT") or 1
    return {b: max(counts[b] / total, 0.01) for b in "ACGT"}


def _count_total_windows(gbk_path: Path, upstream_bp: int, motif_len: int) -> int:
    """Count the total number of sliding window positions evaluated genome-wide.

    This is the true N for Benjamini-Hochberg correction — every position
    on both strands of every CDS upstream region, regardless of whether a
    motif was found there. Using the number of HITS as N instead would
    collapse the correction down to only the successful matches, completely
    erasing the statistical penalty for searching a massive genome and
    producing q-values that are falsely tiny.

    Accounts for contig-boundary truncation: genes near the start of a
    contig have shorter actual upstream regions, and therefore fewer windows.

    Only the window LENGTH is needed here, never the actual bases, so this
    uses ``utils.resolve_upstream_window()`` (arithmetic only) rather than
    ``extract_upstream_window()`` — no sequence slicing or reverse-complement
    work is wasted just to count positions. See the module changelog (v1.4.0).

    Args:
        gbk_path:    Path to the GenBank file.
        upstream_bp: Upstream window size used in the scan.
        motif_len:   Motif length in bp (windows shorter than this contribute 0).

    Returns:
        Total integer count of sliding window positions tested (both strands).
    """
    total = 0
    for record in SeqIO.parse(gbk_path, "genbank"):
        seq_len = len(record.seq)
        for feature in record.features:
            if feature.type != "CDS":
                continue
            start = int(feature.location.start)
            end = int(feature.location.end)
            strand = feature.location.strand
            _slice_start, _slice_end, actual_up = resolve_upstream_window(
                seq_len, start, end, strand, upstream_bp
            )
            windows = max(0, actual_up - motif_len + 1)
            total += windows * 2  # both strands scanned
    return total


def _motif_pvalue(regex_pattern: str, bg: dict[str, float]) -> float:
    """Compute the p-value of a motif match at a single random position.

    DEFINITION (regex-based model):
        P-value = the probability that a random sequence of the same length as the
        motif matches the motif pattern, given background nucleotide frequencies.
        This is the standard definition (FIMO Bailey et al.): "the probability of a
        random sequence matching this position with as good or better a score."

        For a binary regex scanner (match / no-match), "as good or better" means any
        sequence satisfying the regex. The p-value is therefore the product of
        per-position match probabilities:
          - Fixed base (e.g. A):     bg[A]
          - Character class [TC]:    bg[T] + bg[C]  (sum of disjoint options)
          - Unknown / complex token: 1.0            (unconstrained, conservative)

    IMPORTANT CAVEAT — Instance vs Motif P-Value:
        This gives the probability of the MOTIF (the regex pattern itself), not the
        specific observed instance. For a degenerate position [TC], both the T-instance
        and the C-instance have the same p-value — the probability of ANY match to the
        motif pattern. This is the correct definition for a binary regex scanner, but
        it differs from PWM-based p-values (FIMO), which assign per-instance scores
        based on a trained weight matrix.

    CONTRAST WITH FIMO:
        FIMO computes p-values per position from a PWM, producing a different (often
        higher-precision) score for every location. This script computes one p-value
        per motif pattern, applied uniformly to all matches. FIMO is more granular;
        this script is simpler but adequate for operator footprints with known
        consensus structure.

    Args:
        regex_pattern: The IUPAC/regex motif string used in the scan
                       (e.g., "ATCG[TC]TGCGCAGCGG").
        bg:            Background frequencies from _compute_background().

    Returns:
        Float p-value in the range (0, 1].
    """
    p = 1.0
    i = 0
    pattern = regex_pattern.upper()
    while i < len(pattern):
        ch = pattern[i]
        if ch == "[":
            # Character class: sum background probs of all matching bases
            j = pattern.index("]", i)
            bases = set(pattern[i + 1 : j]) & set("ACGT")
            if bases:
                p *= sum(bg.get(b, 0.25) for b in bases)
            i = j + 1
        elif ch in "ACGT":
            p *= bg.get(ch, 0.25)
            i += 1
        else:
            # Complex token (regex quantifier, wildcard, etc.) — no constraint
            i += 1
    return p


def _bh_qvalues(pvalues: list[float], total_tests: int) -> list[float]:
    """Compute Benjamini-Hochberg FDR-corrected q-values.

    Applies the standard BH procedure (Benjamini & Hochberg 1995) to control the
    False Discovery Rate across all motif hits from the entire genome scan.

    INTERPRETATION:
        q-value = the false discovery rate (FDR) if you accept this match as real.
        At a threshold of q < 0.05, you expect ~5% of accepted matches to be false
        positives (random background). Matches with q > 0.05 are noise and should
        be discarded.

        Only matches passing this threshold are included in the final TSV output
        (implicitly filtered, q_value < 0.05 or user-specified α).

    WHY TOTAL_TESTS MATTERS:
        The BH formula requires N = the total number of hypotheses tested, not just
        the successful ones. For a genome-wide scan, this is the number of sliding
        windows evaluated across all CDS on both strands (~544,000 for a typical
        bacterial genome). Using only the hit count as N collapses the correction
        to near-zero and produces q-values that are statistically meaningless.

    CONTRAST WITH FIMO:
        FIMO reports q-values for every position without automatic filtering. The
        user manually selects a threshold (often q < 0.05) to discard noise. This
        script applies the threshold automatically, so the TSV contains only
        significant hits. Both approaches arrive at similar conclusions when FIMO's
        results are filtered to the same significance level.

    Args:
        pvalues:     List of p-values (one per motif match, any order).
        total_tests: Total sliding window positions evaluated genome-wide
                    (from _count_total_windows()). This is the true N for BH —
                    the denominator that scales the correction to the search space.

    Returns:
        List of q-values in the same order as the input p-values.
    """
    n = len(pvalues)
    if n == 0:
        return []
    indexed = sorted(enumerate(pvalues), key=lambda x: x[1])
    qvalues = [1.0] * n
    running_min = 1.0
    for rank_offset, (orig_idx, p) in enumerate(reversed(indexed)):
        rank = n - rank_offset  # rank among hits: 1 (best) to n (worst)
        q = p * total_tests / rank  # BH: p × (total genome-wide tests) / rank
        running_min = min(running_min, q)
        qvalues[orig_idx] = min(1.0, running_min)
    return qvalues


def stream_regulon_hits(
    gbk_path: Path, regex_pattern: str, upstream_bp: int
) -> Iterator[dict]:
    """Scans every CDS upstream region for a motif on both DNA strands.

    Motif positions are returned as negative integers relative to the
    Translation Start Site (TSS), following standard molecular biology
    convention (e.g., the -10 and -35 boxes in prokaryotic promoters).

    Both the coding strand (+) and the template strand (-) are scanned.
    This ensures Transcription Factor binding sites in either orientation
    are detected, including palindromic and non-palindromic motifs.

    The strand-aware upstream slice (and the ``slice_start``/``slice_end``
    needed to map a motif hit back to genomic coordinates) is delegated to
    ``utils.extract_upstream_window()`` rather than reimplemented here —
    see the module changelog (v1.4.0).

    Args:
        gbk_path:       Path to the GenBank file.
        regex_pattern:  IUPAC/regex motif string. IUPAC ambiguity codes
                        (W, R, Y, S, K, M, B, D, H, V, N) are translated to
                        regex character classes via
                        ``utils.translate_iupac_to_regex()`` before
                        compiling — raw codes passed directly to
                        ``re.compile()`` are matched as literal characters
                        and never match real DNA
                        (e.g. "TATAWAW" found zero hits against a sequence
                        containing the valid instance "TATAAAA").
                        Matching is case-insensitive
                        by uppercasing the pattern (once, at compile time) and
                        every extracted upstream sequence (below) rather than
                        passing re.IGNORECASE through the regex engine for
                        every character compared across genome-wide windows.
        upstream_bp:    Number of bases upstream of each CDS start to extract.

    Yields:
        A dict per CDS with at least one motif hit, containing locus_tag,
        product, contig, gene strand, and a sorted list of
        (rel_pos, matched_seq, motif_strand) tuples.
    """
    try:
        safe_pattern = re.compile(f"(?=({translate_iupac_to_regex(regex_pattern)}))")
    except re.error as e:
        raise ValueError(f"Invalid regex pattern: '{regex_pattern}'") from e

    warned_eukaryote = False

    try:
        for record in SeqIO.parse(gbk_path, "genbank"):
            for feature in record.features:
                if feature.type == "mRNA" and not warned_eukaryote:
                    print(
                        "[!] Warning: mRNA features detected — this looks like "
                        "a eukaryotic genome. This script anchors upstream "
                        "windows on CDS start, not the transcription start "
                        "site (TSS), so the true promoter will likely be "
                        "missed. See the module docstring (PROKARYOTE-ONLY) "
                        "for details.",
                        file=sys.stderr,
                    )
                    warned_eukaryote = True
                if feature.type == "CDS":

                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    # Coordinate-based fallback instead of a shared "UNKNOWN"
                    # string: this script is prokaryote-only and prokaryote
                    # files always carry /locus_tag in practice, so this is a
                    # low-probability edge case (e.g. raw Prodigal output
                    # without Prokka/Bakta annotation) — but since results
                    # here are yielded per-feature rather than aggregated
                    # into a dict (no data-loss risk the way there was in
                    # pairwise_homolog_finder.py), the only real cost of the
                    # old default was traceability: every unlabeled gene in
                    # a report would show the identical string "UNKNOWN" with
                    # no way to tell them apart. Same coordinate-fallback
                    # format used throughout the toolkit for consistency.
                    locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                    if not locus_tag:
                        locus_tag = f"UNANNOTATED_{record.id}_{start}_{end}"
                    product = feature.qualifiers.get(
                        "product", ["hypothetical protein"]
                    )[0]
                    strand = feature.location.strand

                    # Extract upstream with boundary tracking. slice_start/
                    # slice_end are needed below to convert motif hit positions
                    # back to genomic coordinates (gstart/gend).
                    upstream_seq, actual_upstream, slice_start, slice_end = (
                        extract_upstream_window(record, start, end, strand, upstream_bp)
                    )
                    upstream_seq = upstream_seq.upper()

                    if actual_upstream < upstream_bp:
                        print(
                            f"    [!] Warning: {locus_tag} upstream truncated to "
                            f"{actual_upstream}bp (contig boundary).",
                            file=sys.stderr,
                        )

                    matches = []

                    # Forward (coding) strand scan
                    # Position reported as negative distance from TSS:
                    #   match.start() 0  → -(actual_upstream)  (farthest from ATG)
                    #   match.start() L-1 → -1                 (one base before ATG)
                    for match in safe_pattern.finditer(upstream_seq):
                        rel_pos = -(actual_upstream - match.start())
                        W = len(match.group(1))
                        if strand == 1:
                            gstart = slice_start + match.start() + 1
                            gend = slice_start + match.start() + W
                        else:
                            j = match.start()
                            gstart = slice_end - j - W + 1
                            gend = slice_end - j
                        matches.append((rel_pos, match.group(1), "+", gstart, gend))

                    # Reverse complement (template) strand scan
                    # The RC of the upstream sequence is scanned with the same pattern.
                    # The coordinate is mapped back to the forward-strand TSS origin:
                    #   forward position of motif 5' end = len(upstream_seq) - true_match_end
                    #   biological coord = -(actual_upstream - forward_position)
                    #
                    # IMPORTANT: Because we use a zero-width lookahead assertion (?=(...)),
                    # match.end() always equals match.start() — the outer match consumes
                    # zero characters. Using match.end() directly would place every RC hit
                    # at the wrong position (off by the motif length). We must calculate
                    # the true end from the captured group's length instead.
                    rc_seq = str(Seq(upstream_seq).reverse_complement())
                    for match in safe_pattern.finditer(rc_seq):
                        true_match_end = match.start() + len(match.group(1))
                        fwd_pos = len(upstream_seq) - true_match_end
                        rel_pos = -(actual_upstream - fwd_pos)
                        W = len(match.group(1))
                        if strand == 1:
                            gstart = slice_start + fwd_pos + 1
                            gend = slice_start + fwd_pos + W
                        else:
                            gstart = slice_end - fwd_pos - W + 1
                            gend = slice_end - fwd_pos
                        matches.append((rel_pos, match.group(1), "-", gstart, gend))

                    if matches:
                        # Sort biologically: 5' → 3' relative to TSS (most negative first)
                        matches.sort(key=lambda x: x[0])
                        yield {
                            "locus_tag": locus_tag,
                            "product": product,
                            "contig": record.id,
                            "strand": strand,
                            "matches": matches,
                        }

    except (FileNotFoundError, OSError, UnicodeDecodeError, ValueError) as e:
        raise ValueError(f"GenBank Parsing Error in {gbk_path.name}: {e}") from e


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Genome-Wide Regulon Scanner\n"
        "Regex-Based Statistical Model\n"
        "Prokaryote-only (CDS-anchored, not TSS-anchored)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-u",
        "--upstream",
        type=int,
        default=150,
        help=(
            "Bases upstream of each CDS to extract and scan. "
            "Default: 150 (appropriate for prokaryotes). This script is "
            "prokaryote-only — it anchors on CDS start, not the TSS, so "
            "increasing this value does not make it usable on eukaryotic "
            "genomes. See universal_promoter_extractor.py for TSS-anchored "
            "eukaryotic upstream extraction."
        ),
    )
    parser.add_argument(
        "-m",
        "--motif",
        required=True,
        help="Regex/IUPAC motif to search for on both strands",
    )
    args = parser.parse_args()

    # Translate IUPAC ambiguity codes to regex character classes ONCE here,
    # and use this single translated string everywhere downstream (search,
    # motif-length estimate, and p-value calculation). This matters beyond
    # just the search: _motif_pvalue() already correctly handles explicit
    # "[TC]"-style character classes, but does NOT recognize raw IUPAC
    # letters (W, R, Y, ...) — it falls through to its "complex token, no
    # constraint" branch for them, silently treating an ambiguous position
    # as if it could be ANY base for free. The same
    # biological motif spelled as "TATAWAW" vs "TATA[AT]A[AT]" produced two
    # different p-values (0.0024 vs 0.00088) purely from spelling, with the
    # raw-IUPAC spelling UNDERSTATING significance. Translating once here,
    # before either consumer sees the pattern, makes both correct and
    # consistent without needing to change _motif_pvalue() itself — it
    # already does the right thing once given explicit bracket notation.
    translated_motif = translate_iupac_to_regex(args.motif)

    print(f"[*] Scanning genome    : {args.input.name}", file=sys.stderr)
    print(f"[*] Upstream region    : {args.upstream}bp", file=sys.stderr)
    print(f"[*] Motif              : {args.motif}", file=sys.stderr)
    print(f"[*] Strands scanned    : Both (+) coding and (-) template", file=sys.stderr)
    print(f"[*] Position reference : TSS (negative = upstream of ATG)", file=sys.stderr)

    try:
        # ── Collect all hits ──────────────────────────────────────────────────
        # Results are collected into memory first so that BH q-values can be
        # computed across ALL matches from ALL genes before writing the TSV.
        print(f"[*] Collecting hits...", file=sys.stderr)
        all_hits = list(
            stream_regulon_hits(args.input, translated_motif, args.upstream)
        )

        total_genes_hit = len(all_hits)
        total_motifs_found = sum(len(h["matches"]) for h in all_hits)

        for hit in all_hits:
            print(
                f"    -> Regulon member found: {hit['locus_tag']} "
                f"({hit['product'][:40]}...)",
                file=sys.stderr,
            )

        # ── Compute background + total tests + p-values + q-values ──────────────
        print(f"[*] Computing background frequencies...", file=sys.stderr)
        bg = _compute_background(args.input, args.upstream)
        print(
            f"    Background: {', '.join(f'{b}={v:.3f}' for b,v in bg.items())}",
            file=sys.stderr,
        )

        # Approximate motif length from regex: replace each [...] group with a
        # single placeholder character so each degenerate position counts as 1bp.
        # (Using '' instead of 'X' would undercount by one per degenerate position.)
        # Uses translated_motif so a raw IUPAC code (which becomes a bracket
        # group after translation) is correctly counted as 1bp too, exactly
        # like a hand-written bracket group already was.
        approx_motif_len = len(re.sub(r"\[.*?\]", "X", translated_motif))

        print(f"[*] Counting total windows tested...", file=sys.stderr)
        total_tests = _count_total_windows(args.input, args.upstream, approx_motif_len)
        print(
            f"    Approx motif length : {approx_motif_len}bp",
            file=sys.stderr,
        )
        print(
            f"    Total windows tested: {total_tests:,}  "
            f"(all CDS × both strands × upstream positions)",
            file=sys.stderr,
        )

        # P-value = probability of ANY sequence matching the motif at a random
        # position (motif footprint probability, not the specific instance).
        # All hits from the same motif scan share this single p-value — correct
        # for a binary regex scanner where every match satisfies the same pattern.
        # Uses translated_motif so raw IUPAC codes are correctly constrained
        # (see the note above translated_motif's assignment).
        motif_p = _motif_pvalue(translated_motif, bg)
        print(f"    Motif p-value       : {motif_p:.3e}", file=sys.stderr)

        all_pvalues: list[float] = [motif_p] * sum(len(h["matches"]) for h in all_hits)

        # Q-values use total_tests (genome-wide windows) as N — not len(hits).
        all_qvalues = _bh_qvalues(all_pvalues, total_tests)

        # Re-attach p-values and q-values to each match
        pval_idx = 0
        for hit in all_hits:
            enriched = []
            for match in hit["matches"]:
                # match = (rel_pos, matched_seq, motif_strand, gstart, gend)
                p = all_pvalues[pval_idx]
                q = all_qvalues[pval_idx]
                enriched.append((*match, p, q))
                pval_idx += 1
            hit["matches"] = enriched
        # enriched match tuple:
        # (rel_pos, matched_seq, motif_strand, gstart, gend, p_value, q_value)

        # ── Write output ──────────────────────────────────────────────────────
        if args.output:
            with open(args.output, "w", encoding="utf-8") as tsv:
                tsv.write(
                    "Locus_Tag\tContig\tGene_Strand\tMotif_Count\t"
                    "Positions_Relative_to_TSS\tGenomic_Start\tGenomic_End\t"
                    "P_Value\tQ_Value\tMatched_Sequences\tProduct\n"
                )
                for hit in all_hits:
                    m = hit["matches"]
                    positions = ",".join(f"{x[0]}({x[2]})" for x in m)
                    gstarts = ",".join(str(x[3]) for x in m)
                    gends = ",".join(str(x[4]) for x in m)
                    pvals = ",".join(f"{x[5]:.3e}" for x in m)
                    qvals = ",".join(f"{x[6]:.3e}" for x in m)
                    sequences = ",".join(x[1] for x in m)
                    tsv.write(
                        f"{hit['locus_tag']}\t{hit['contig']}\t{hit['strand']}\t"
                        f"{len(m)}\t{positions}\t{gstarts}\t{gends}\t"
                        f"{pvals}\t{qvals}\t{sequences}\t{hit['product']}\n"
                    )
        else:
            print(
                "\n[*] Note: No output file specified (-o). Results printed to terminal only.",
                file=sys.stderr,
            )

        print("\n" + "=" * 40, file=sys.stderr)
        print("          PIPELINE COMPLETE", file=sys.stderr)
        print("=" * 40, file=sys.stderr)
        print(f"Total Genes in Regulon : {total_genes_hit}", file=sys.stderr)
        print(f"Total Motifs Bound     : {total_motifs_found}", file=sys.stderr)
        if args.output:
            print(f"Results written to     : {args.output.resolve()}", file=sys.stderr)
        print("=" * 40, file=sys.stderr)

    except (ValueError, FileNotFoundError, PermissionError) as e:
        sys.exit(f"\n[!] Pipeline Halted: {e}")
    except KeyboardInterrupt:
        sys.exit("\n[!] Pipeline gracefully interrupted by user.")
    except Exception:
        print("\n[!] UNEXPECTED BUG ENCOUNTERED:", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
