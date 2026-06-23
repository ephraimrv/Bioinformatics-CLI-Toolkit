#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

r"""Universal Promoter Extractor

Extracts upstream regulatory regions (promoter sequences) from GenBank files
for MEME motif discovery. Works on both prokaryotic and eukaryotic genomes
without requiring separate scripts.

ORGANISM MODES:
    auto        Detects organism type from the file (default). Prokaryotic
                GenBank files (Prokka, Bakta, NCBI prokaryote) never contain
                mRNA features. Eukaryotic files (Ensembl, RefSeq eukaryote,
                Augustus) always contain mRNA features. Detection is reliable
                and fires per file, so mixed directories work correctly.

    prokaryote  Forces CDS-based extraction. Upstream sequence is extracted
                relative to the CDS start coordinate (= translation start = ATG),
                which is correct for prokaryotes because there is no 5' UTR
                separating the promoter from the coding sequence.

    eukaryote   Forces mRNA-based extraction. Upstream sequence is extracted
                relative to the mRNA start coordinate (= Transcription Start
                Site = TSS), which is the biologically correct anchor for
                promoter analysis in eukaryotes. Extracting upstream of the
                CDS start in eukaryotes is wrong because the CDS begins at
                the ATG, which is separated from the TSS by a 5' UTR (and
                possibly introns), so the "upstream" region would actually be
                inside a UTR or intron rather than the true promoter.

                REQUIREMENT: The GenBank file must contain mRNA features with
                /locus_tag qualifiers. Files that only have CDS features will
                yield no results in eukaryote mode.

KEYWORD SEARCH IN EUKARYOTE MODE:
    mRNA features often lack /product qualifiers. This script performs a
    two-pass scan per record in eukaryote mode:
        Pass 1 — CDS features: collect {locus_tag: product} for all CDS
                 whose /product annotation matches a keyword.
        Pass 2 — mRNA features: for locus tags found in pass 1, resolve
                 the 5'-most TSS across all isoforms, then extract upstream.
    This guarantees that keyword matching always uses the /product annotation
    regardless of which feature type carries it, while coordinate extraction
    always uses the biologically correct mRNA anchor.

ALTERNATIVE SPLICING (ISOFORM HANDLING):
    In eukaryotes, a single locus_tag can have multiple mRNA features
    representing alternative splice variants. Each isoform may have a
    different Transcription Start Site (TSS). By default, the script
    resolves this by finding the 5'-most TSS across all isoforms of a
    locus:
        Forward strand (+): 5'-most TSS = smallest mRNA start coordinate.
        Reverse strand (-): 5'-most TSS = largest mRNA end coordinate.
    When multiple isoforms are detected, a [Multi-isoform] message is printed
    to stderr showing how many isoforms were found and which TSS was used.

    CAVEAT: the 5'-most TSS is a conservative, deterministic choice — not
    necessarily the biologically dominant one. A rare, minimally-expressed
    transcript variant with an unusually far-upstream TSS would be chosen
    over a well-supported, highly-expressed canonical isoform with a
    closer-in TSS, since this script has no expression data and most
    GenBank files carry no "canonical transcript" annotation to prefer
    instead. Pass ``--all-isoforms`` (v1.7.0) to extract a separate
    upstream region per isoform instead of merging them — see that flag's
    help text below.

CAVEAT — applies to every mode of this script, not fixable in code:
    (1) DETECTION HEURISTIC: ``--mode auto`` calls a file eukaryotic if and
        only if an ``mRNA`` feature key appears before ``ORIGIN``. This is
        reliable for standard NCBI/Ensembl/Augustus output, but is not a
        universal biological test — a GFF-to-GenBank conversion, a
        partial/ncRNA-focused annotation, or a custom pipeline could
        produce a genuinely eukaryotic file with no ``mRNA`` feature at
        all, which this heuristic would misclassify as prokaryotic. Pass
        ``--mode eukaryote`` explicitly if you know the organism and
        auto-detection seems wrong.
    (2) TSS ACCURACY: an mRNA feature's start coordinate is only as
        reliable as the annotation that produced it. For NCBI RefSeq
        entries built from cap-trapped/CAGE-seq transcript evidence, this
        is a genuine experimentally-supported TSS. For computational gene
        predictions (Augustus, MAKER, and many draft eukaryotic
        annotations), "mRNA start" usually means "transcript MODEL start"
        — the predictor's best guess, not a direct experimental
        measurement. This script cannot tell which kind of annotation it
        has; treat extracted "promoters" from predicted annotations as a
        best-available approximation, not a confirmed TSS, especially
        before reporting them in a manuscript.
    (3) PROMOTER LENGTH IS NOT ONE-SIZE-FITS-ALL: ``--upstream N`` extracts
        a fixed window for every gene. In bacteria and fungi this is
        usually a reasonable approximation (operons/promoters are
        physically compact). In plants and animals it frequently is not:
        a real regulatory element can sit hundreds of bp to tens of kb
        from the TSS, and no single ``--upstream`` value is correct for
        every gene in such a genome. There is no universal promoter
        length to default to — this is a biological reality, not a
        parameter this script can tune around. See the equivalent CAVEAT
        in find_gbk_features.py's ``--context`` documentation for the
        same underlying point applied to coordinate-proximity search.

ALL-ISOFORMS MODE (v1.7.0):
    Pass ``--all-isoforms`` to extract one upstream region PER ISOFORM
    instead of merging all isoforms of a locus into a single 5'-most-TSS
    region (the default — see the CAVEAT above for why that default isn't
    always biologically representative). Each isoform gets a disambiguated
    identifier: ``{locus_tag}#{transcript_id}`` when the mRNA feature
    carries a ``/transcript_id`` qualifier, else ``{locus_tag}#isoform{N}``.
    This mirrors ``pairwise_homolog_finder.py``'s ``--keep-all-isoforms``
    flag and identifier format exactly, for consistency across the
    toolkit. Off by default — existing output is completely unchanged
    unless you opt in. Increases output size proportionally to isoform
    count; fine for a handful of isoforms per locus, more verbose at
    whole-genome scale with many alternatively-spliced genes.

Note:
    Associated with ongoing, unpublished research (manuscript in
    preparation). Correct attribution is requested when used in
    derivative works.

    v1.5.0: The strand-aware upstream-slicing arithmetic in
    ``_extract_upstream_seq`` (below) now delegates to
    ``utils.extract_upstream_window()`` instead of reimplementing the
    forward/reverse slice math locally. This was one of several
    independent copies of the same coordinate logic spread across this
    script, ``regulon_scanner.py``, and ``utils.extract_upstream_sequence``
    — see ``utils.py``'s v1.2.0 changelog note for the full consolidation.
    This function's own signature, return contract, and truncation-warning
    behavior are unchanged; only the internal slicing math moved.

    v1.6.0: Fixed a critical data-loss bug, mirroring the one just fixed
    in pairwise_homolog_finder.py. Every locus_tag extraction site in this
    file previously defaulted features lacking ``/locus_tag`` to the
    literal string ``"UNKNOWN"``. On a file where ``/locus_tag`` is
    absent entirely (common for eukaryotic assemblies, GFF3-to-GenBank
    conversions, draft genomes), every such feature collided on that one
    key. I tested this directly: 50 keyword-matching CDS collapsed to 1
    entry in Pass 1, and 5,000 unannotated mRNA features across an entire
    simulated genome collapsed into a single fabricated "super-gene"
    spanning the most extreme min/max coordinates in the file. All six
    extraction sites (both branches of ``extract_regulatory_regions()``
    and ``extract_by_loci()``) now use ``_resolve_identifier()`` — the
    same locus_tag -> protein_id -> gene -> coordinate fallback hierarchy
    as ``pairwise_homolog_finder.py``, using an identical coordinate-fallback
    string format so the two scripts always agree on an ID for the same
    feature (critical for ``target_promoter_pipeline.py``, which passes
    IDs from one script's output into the other's ``extract_by_loci()``
    lookup).

    A second, subtler issue turned up in the same pass: simply
    swapping in ``_resolve_identifier()`` is not sufficient for the
    eukaryotic two-pass design when a gene has NO locus_tag/protein_id/
    gene at all. The CDS pass and the mRNA pass call the resolver on
    DIFFERENT features, and CDS coordinates routinely differ from mRNA
    coordinates for the same gene (UTRs) — so the two passes' fallback
    IDs would never string-match even after the basic fix, silently
    dropping every fully-unannotated eukaryotic gene. Both eukaryote
    branches now also track each fallback-keyed CDS's coordinate span and,
    when an mRNA's own fallback ID doesn't directly match, correlate it
    via coordinate containment instead (``_find_overlapping_fallback_id()``)
    — does this mRNA's span fully contain a known fallback CDS's span?
    Genes with real locus_tag/protein_id/gene annotations are unaffected;
    they still match directly as before.

    v1.7.0: Two additions found during review, neither changing default
    behavior unless explicitly requested.
    (1) Added CAVEAT documentation (module docstring) covering four
    pre-existing, non-code limitations I'd left real
    but undocumented: the mRNA-presence eukaryote-detection heuristic
    isn't universal (GFF-to-GenBank conversions, ncRNA-only annotations
    could be eukaryotic with no mRNA feature at all); mRNA start is only
    as reliable as the annotation that produced it (an experimentally-
    supported TSS for NCBI RefSeq, but often just a transcript MODEL
    start for Augustus/MAKER predictions); ``--upstream N`` is necessarily
    a fixed-window approximation, biologically reasonable for bacteria/
    fungi but not for plants/animals, where regulatory elements can sit
    much farther from the TSS than any single window value captures; and
    the default 5'-most-TSS isoform choice is conservative/deterministic,
    not necessarily the biologically dominant transcript.
    (2) Added ``--all-isoforms`` (both ``extract_regulatory_regions()``
    and ``extract_by_loci()``): extracts one upstream region per mRNA
    isoform instead of merging all isoforms of a locus to a single
    5'-most-TSS region. Directly addresses the limitation in (1) above —
    mirrors ``pairwise_homolog_finder.py``'s ``--keep-all-isoforms`` flag
    and identifier format (``{locus_tag}#{transcript_id}``, or
    ``{locus_tag}#isoform{N}`` without one) exactly, for consistency
    across the toolkit. Off by default; existing callers and existing
    output are completely unaffected unless this flag is explicitly
    passed.

    v1.7.1: The isoform-disambiguation
    logic added in v1.7.0 (``{locus_tag}#{transcript_id}`` /
    ``{locus_tag}#isoform{N}``) was duplicated inline in two places in
    this file, AND duplicated again, separately, in
    pairwise_homolog_finder.py's ``--keep-all-isoforms`` — four copies of
    essentially the same pattern across two scripts, differing only in
    which qualifier (``transcript_id`` vs ``protein_id``) each looked up.
    Both call sites here now use the new shared
    ``utils.disambiguate_isoform_id()`` instead. No behavior change —
    I re-ran this file's existing ``--all-isoforms`` tests
    before and after the swap to be sure.

Examples:
    # Prokaryote (auto-detected)
    $ python3 universal_promoter_extractor.py \
        -i C5_genome.gbk -o C5_promoters.fasta -u 150 \
        -k bacteriocin lactobin cerein

    # Eukaryote (auto-detected)
    $ python3 universal_promoter_extractor.py \
        -i Arabidopsis.gbff -o arab_promoters.fasta -u 500 \
        -k "WRKY transcription factor" "disease resistance"

    # Force mode explicitly
    $ python3 universal_promoter_extractor.py \
        -i genome.gbk -o out.fasta -u 150 -k bacteriocin \
        --mode prokaryote

    # Scan a directory of mixed genomes (auto-detects each file)
    $ python3 universal_promoter_extractor.py \
        -i references/ -o all_promoters.fasta -u 150 \
        -k bacteriocin lactobin

    # Eukaryotic file, keeping every splice isoform separately
    $ python3 universal_promoter_extractor.py \
        -i Arabidopsis.gbff -o arab_all_isoforms.fasta -u 500 \
        -k "WRKY transcription factor" --all-isoforms
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.7.1"

import sys
import argparse
import re
from pathlib import Path
from typing import Iterator

try:
    from Bio import SeqIO
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )
from utils import (
    stream_reference_files,
    extract_upstream_window,
    disambiguate_isoform_id,
)

# ── Mode constants ─────────────────────────────────────────────────────────────

MODE_AUTO = "auto"
MODE_PROKARYOTE = "prokaryote"
MODE_EUKARYOTE = "eukaryote"
_VALID_MODES = (MODE_AUTO, MODE_PROKARYOTE, MODE_EUKARYOTE)


# ── Auto-detection ─────────────────────────────────────────────────────────────


def _detect_organism_mode(gbk_path: Path) -> str:
    """Detect whether a GenBank file is prokaryotic or eukaryotic.

    Uses raw text streaming rather than BioPython's SeqIO.parse, which
    is orders of magnitude faster — especially for large eukaryotic
    GenBank records where SeqIO would instantiate thousands of Python
    objects just to check one feature type.

    Strategy: stream lines until the first ``mRNA`` feature key is found
    (eukaryote) or until the ``ORIGIN`` keyword is reached (end of the
    feature table — prokaryote). Only the first record is scanned since
    annotation style is uniform across all records in a standard assembly.

    GenBank format note: feature keys are indented with exactly 5 spaces.
    Qualifier lines are indented with 21 spaces. Therefore
    ``line.startswith("     mRNA")`` safely matches feature declarations
    without false-positives from qualifier values like
    ``/product="mRNA processing factor"``.

    Args:
        gbk_path: Path to the GenBank file to inspect.

    Returns:
        ``"eukaryote"`` if an mRNA feature is found before ORIGIN,
        ``"prokaryote"`` otherwise.
    """
    try:
        with open(gbk_path, "r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("     mRNA"):
                    return MODE_EUKARYOTE
                # ORIGIN marks the end of the feature table and start of
                # the raw sequence — no point reading further.
                if line.startswith("ORIGIN"):
                    break
    except Exception:
        pass
    return MODE_PROKARYOTE


# ── Shared coordinate helper ───────────────────────────────────────────────────


def _extract_upstream_seq(
    record,
    start: int,
    end: int,
    strand: int,
    upstream_bp: int,
    locus_tag: str,
) -> tuple[str, int]:
    """Extract the upstream sequence and return (sequence, actual_length).

    Handles both strands and contig boundary truncation. Prints a warning
    to stderr when the requested window is truncated.

    The strand-aware slicing arithmetic itself lives in
    ``utils.extract_upstream_window()`` (shared with ``regulon_scanner.py``
    and ``utils.extract_upstream_sequence``); this function only adds the
    truncation-warning behavior specific to this script's output.

    Args:
        record:      BioPython SeqRecord containing the genome sequence.
        start:       0-based feature start coordinate (BioPython convention).
        end:         0-based feature end coordinate (BioPython convention).
        strand:      1 for forward strand, -1 for reverse strand.
        upstream_bp: Requested upstream window in base pairs.
        locus_tag:   Locus tag string used only for the truncation warning.

    Returns:
        A 2-tuple of (upstream_dna_sequence, actual_extracted_length).
        The actual length may be shorter than upstream_bp when the feature
        is within upstream_bp bases of a contig boundary.
    """
    upstream_seq, actual_upstream, _slice_start, _slice_end = extract_upstream_window(
        record, start, end, strand, upstream_bp
    )

    if actual_upstream < upstream_bp:
        print(
            f"      [!] Warning: {locus_tag} upstream truncated to "
            f"{actual_upstream}bp (contig boundary — requested {upstream_bp}bp).",
            file=sys.stderr,
        )

    return upstream_seq, actual_upstream


# ── Identifier resolution (fixes the "UNKNOWN" collision data-loss bug) ────────


def _resolve_identifier(feature, record_id: str) -> str:
    """Resolves a unique grouping key for a feature, with safe fallbacks.

    Mirrors ``pairwise_homolog_finder._resolve_identifier()`` exactly — same
    fallback order, same coordinate-fallback string format — so that the
    two scripts always agree on an identifier for the same feature. This
    matters specifically for ``extract_by_loci()``, the bridge function
    ``target_promoter_pipeline.py`` uses to look up loci that
    ``pairwise_homolog_finder.py`` already found: if the two scripts computed
    fallback IDs differently, a homolog found under one script's ID would
    never be found by the other's lookup.

    Previously this script (like pairwise_homolog_finder.py before its own
    fix) defaulted every feature lacking ``/locus_tag`` to the literal
    string ``"UNKNOWN"``. On a file where ``/locus_tag`` is absent
    entirely, every such feature collided on that one key — collapsing
    an entire eukaryotic genome's worth of mRNA
    features into a single fabricated "super-gene" spanning the most
    extreme min/max coordinates in the file.

    Fallback order (each step only used if the previous one is empty):
      1. ``/locus_tag``
      2. ``/protein_id``
      3. ``/gene``
      4. ``record_id`` + genomic coordinates (guaranteed unique; includes
         the contig/record ID specifically because coordinates alone are
         only unique within one contig).

    Args:
        feature:   A Biopython SeqFeature to resolve an identifier for.
        record_id: The parent record's ``.id`` (contig/chromosome name),
                   used only in the final coordinate-based fallback.

    Returns:
        A non-empty string suitable for use as a deduplication/lookup key.
        Never ``"UNKNOWN"`` or any other shared constant. Coordinate-
        fallback identifiers are prefixed ``"UNANNOTATED_"`` so callers can
        detect them (see ``_find_overlapping_fallback_id()`` below).
    """
    identifier = feature.qualifiers.get("locus_tag", [""])[0]
    if identifier:
        return identifier

    identifier = feature.qualifiers.get("protein_id", [""])[0]
    if identifier:
        return identifier

    identifier = feature.qualifiers.get("gene", [""])[0]
    if identifier:
        return identifier

    start = int(feature.location.start)
    end = int(feature.location.end)
    return f"UNANNOTATED_{record_id}_{start}_{end}"


def _find_overlapping_fallback_id(
    mrna_start: int,
    mrna_end: int,
    fallback_cds_spans: dict[str, tuple[int, int]],
) -> str | None:
    """Correlates an unannotated mRNA to its CDS via coordinate containment.

    Why this exists: when a gene has NO ``/locus_tag``, ``/protein_id``, or
    ``/gene`` qualifier anywhere, ``_resolve_identifier()`` falls back to a
    coordinate-based string built from whichever feature it was called on.
    The CDS pass and the mRNA pass call it on DIFFERENT features — and CDS
    coordinates routinely differ from mRNA coordinates for the same gene,
    because the mRNA spans the full transcript (5' UTR + exons + 3' UTR)
    while the CDS spans only the coding portion (ATG to stop). For
    example: a CDS at 1500-2800 and its own mRNA at 1200-3100 (same
    gene, no shared qualifiers) resolve to two DIFFERENT fallback strings
    even though they're the same gene — so a literal string match between
    the CDS pass's identifier and the mRNA pass's identifier would never
    succeed for this case, silently dropping every fully-unannotated
    eukaryotic gene even after the basic UNKNOWN-collision fix.

    This function instead checks which (if any) already-known fallback-ID
    CDS span is fully CONTAINED within the mRNA's span — i.e. the mRNA's
    transcript structurally includes this CDS's coding region, which is
    exactly the UTR relationship described above. This correctly
    correlates the two without relying on matching qualifier text that
    doesn't exist for either feature.

    Args:
        mrna_start:         0-based mRNA feature start coordinate.
        mrna_end:           0-based mRNA feature end coordinate.
        fallback_cds_spans: Dict mapping fallback CDS identifier ->
                             (cds_start, cds_end), built during the CDS
                             pass for every CDS whose identifier was
                             itself a coordinate fallback (i.e. starts
                             with ``"UNANNOTATED_"``).

    Returns:
        The matching CDS's fallback identifier if exactly one fallback CDS
        span is contained within the mRNA's span, else ``None``. In the
        rare case of multiple contained spans (e.g. an unusual nested
        annotation), the first one found is used — this is a known,
        accepted simplification rather than an attempt to disambiguate
        genuinely ambiguous nested annotations.
    """
    for cds_id, (cds_start, cds_end) in fallback_cds_spans.items():
        if mrna_start <= cds_start and cds_end <= mrna_end:
            return cds_id
    return None


# ── Core extraction functions ──────────────────────────────────────────────────


def extract_regulatory_regions(
    gbk_path: Path,
    keywords: list[str],
    upstream_bp: int,
    mode: str = MODE_AUTO,
    all_isoforms: bool = False,
) -> Iterator[tuple[str, str, str, str, int]]:
    """Scan a GenBank file for keyword-matching genes and extract their upstream regions.

    Supports both prokaryotic and eukaryotic GenBank files via the ``mode``
    parameter. When ``mode="auto"`` (default), organism type is detected
    automatically from the file.

    Prokaryote mode:
        Scans ``CDS`` features for keyword matches in /product annotations.
        Extracts upstream sequence relative to the CDS start coordinate
        (= ATG / translation start).

    Eukaryote mode:
        Uses a two-pass strategy per record:
        Pass 1 — ``CDS`` features: builds a {locus_tag: product} map for all
                 CDS whose /product matches a keyword. This is necessary because
                 mRNA features often lack /product qualifiers.
        Pass 2 — ``mRNA`` features: for locus_tags found in pass 1, extracts
                 upstream sequence relative to the mRNA start coordinate (= TSS).
                 By default, isoforms of one locus are merged to a single
                 5'-most-TSS region; pass ``all_isoforms=True`` to instead
                 yield one region per individual isoform (see
                 ``ALL-ISOFORMS MODE`` in the module docstring).

    Tracks actual extracted length separately from the requested window.
    These differ when a gene is within upstream_bp bases of a contig boundary.

    Args:
        gbk_path:     Path to the target .gbk or .gbff file.
        keywords:     Keywords to match against /product annotations (case-insensitive).
        upstream_bp:  Number of base pairs to extract upstream of the anchor coordinate.
        mode:         Organism mode: ``"auto"``, ``"prokaryote"``, or ``"eukaryote"``.
        all_isoforms: If ``True``, yield one upstream region per mRNA isoform
                      instead of merging to the 5'-most TSS per locus. No
                      effect in prokaryote mode (no isoforms to merge).

    Yields:
        A 5-item tuple:
        ``(seq_id, locus_tag, product, upstream_seq, actual_upstream)``

    Raises:
        ValueError: If the GenBank file is malformed or unreadable.
    """
    resolved_mode = _detect_organism_mode(gbk_path) if mode == MODE_AUTO else mode

    try:
        with open(gbk_path, "r", encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, "genbank"):

                if resolved_mode == MODE_PROKARYOTE:
                    # ── Prokaryote: single pass over CDS features ──────────────
                    for feature in record.features:
                        if feature.type != "CDS":
                            continue

                        product = feature.qualifiers.get("product", [""])[0]
                        if not any(k.lower() in product.lower() for k in keywords):
                            continue

                        locus_tag = _resolve_identifier(feature, record.id)
                        start = int(feature.location.start)
                        end = int(feature.location.end)
                        strand = feature.location.strand

                        upstream_seq, actual_upstream = _extract_upstream_seq(
                            record, start, end, strand, upstream_bp, locus_tag
                        )
                        yield record.id, locus_tag, product, upstream_seq, actual_upstream

                else:
                    # ── Eukaryote: two-pass per record ────────────────────────
                    # Pass 1 — collect keyword-matching locus_tags from CDS features.
                    # mRNA features often lack /product, so keyword matching is always
                    # done against CDS annotations regardless of extraction mode.
                    #
                    # fallback_cds_spans tracks the CDS coordinate span for every
                    # entry whose identifier came from the coordinate fallback
                    # (no locus_tag/protein_id/gene at all) — needed in Pass 2 to
                    # correlate that gene's mRNA feature, since the mRNA's OWN
                    # fallback identifier would be computed from different
                    # coordinates (UTRs) and would never string-match this one.
                    # See _find_overlapping_fallback_id()'s docstring.
                    keyword_loci: dict[str, str] = {}  # {locus_tag: product}
                    fallback_cds_spans: dict[str, tuple[int, int]] = {}
                    for feature in record.features:
                        if feature.type != "CDS":
                            continue
                        product = feature.qualifiers.get("product", [""])[0]
                        if any(k.lower() in product.lower() for k in keywords):
                            lt = _resolve_identifier(feature, record.id)
                            keyword_loci[lt] = product
                            if lt.startswith("UNANNOTATED_"):
                                fallback_cds_spans[lt] = (
                                    int(feature.location.start),
                                    int(feature.location.end),
                                )

                    if not keyword_loci:
                        continue

                    if all_isoforms:
                        # ── ALL-ISOFORMS MODE (v1.7.0): yield one region per
                        # individual mRNA feature, not merged to the 5'-most
                        # TSS — see the module docstring's ALL-ISOFORMS MODE
                        # section. Each isoform gets a disambiguated
                        # identifier, mirroring pairwise_homolog_finder.py's
                        # --keep-all-isoforms format exactly.
                        isoform_counters: dict[str, int] = {}
                        for feature in record.features:
                            if feature.type != "mRNA":
                                continue
                            start = int(feature.location.start)
                            end = int(feature.location.end)
                            strand = feature.location.strand

                            locus_tag = _resolve_identifier(feature, record.id)
                            if locus_tag not in keyword_loci:
                                if locus_tag.startswith("UNANNOTATED_"):
                                    matched = _find_overlapping_fallback_id(
                                        start, end, fallback_cds_spans
                                    )
                                    if matched is None:
                                        continue
                                    locus_tag = matched
                                else:
                                    continue

                            isoform_key = disambiguate_isoform_id(
                                locus_tag,
                                feature,
                                isoform_counters,
                                id_qualifier="transcript_id",
                            )

                            product = keyword_loci[locus_tag]
                            upstream_seq, actual_upstream = _extract_upstream_seq(
                                record, start, end, strand, upstream_bp, isoform_key
                            )
                            yield (
                                record.id,
                                isoform_key,
                                product,
                                upstream_seq,
                                actual_upstream,
                            )
                        continue

                    # Pass 2 — find the 5'-most TSS per locus across all isoforms.
                    # In eukaryotes, a single locus_tag can have multiple mRNA
                    # features (alternative splice variants / isoforms). Each isoform
                    # may have a different Transcription Start Site. Keeping the first
                    # one seen is wrong — it picks an arbitrary isoform. The default
                    # anchor is the 5'-most TSS, which represents the full extent of
                    # the regulatory region (see the module docstring's CAVEAT on why
                    # this is a conservative, not necessarily biologically dominant,
                    # choice — pass all_isoforms=True above to avoid this merge
                    # entirely instead).
                    #   Forward strand (+): 5'-most TSS = smallest start coordinate.
                    #   Reverse strand (-): 5'-most TSS = largest end coordinate.
                    locus_tss: dict[str, dict] = (
                        {}
                    )  # {locus_tag: {start, end, strand, n_isoforms}}

                    for feature in record.features:
                        if feature.type != "mRNA":
                            continue
                        start = int(feature.location.start)
                        end = int(feature.location.end)
                        strand = feature.location.strand

                        locus_tag = _resolve_identifier(feature, record.id)
                        if locus_tag not in keyword_loci:
                            # Direct match failed. If this mRNA itself has no
                            # locus_tag/protein_id/gene either, its fallback ID
                            # was computed from ITS OWN coordinates — which
                            # differ from the CDS's coordinates whenever the
                            # gene has UTRs, so it will never string-match the
                            # CDS pass's fallback ID for the same gene. Try
                            # coordinate containment instead: does this mRNA's
                            # span fully contain a known fallback CDS's span?
                            if locus_tag.startswith("UNANNOTATED_"):
                                matched = _find_overlapping_fallback_id(
                                    start, end, fallback_cds_spans
                                )
                                if matched is None:
                                    continue
                                locus_tag = matched
                            else:
                                continue

                        if locus_tag not in locus_tss:
                            locus_tss[locus_tag] = {
                                "start": start,
                                "end": end,
                                "strand": strand,
                                "n_isoforms": 1,
                            }
                        else:
                            entry = locus_tss[locus_tag]
                            entry["n_isoforms"] += 1
                            if strand == 1:
                                entry["start"] = min(entry["start"], start)
                            else:
                                entry["end"] = max(entry["end"], end)

                    # Yield one upstream region per locus using the resolved TSS
                    for locus_tag, bounds in locus_tss.items():
                        if bounds["n_isoforms"] > 1:
                            print(
                                f"      [Multi-isoform] {locus_tag}: "
                                f"{bounds['n_isoforms']} mRNA isoforms found — "
                                f"using 5'-most TSS.",
                                file=sys.stderr,
                            )
                        product = keyword_loci[locus_tag]
                        upstream_seq, actual_upstream = _extract_upstream_seq(
                            record,
                            bounds["start"],
                            bounds["end"],
                            bounds["strand"],
                            upstream_bp,
                            locus_tag,
                        )
                        yield record.id, locus_tag, product, upstream_seq, actual_upstream

    except Exception as e:
        raise ValueError(f"Failed to parse {gbk_path.name}: {e}") from e


def extract_by_loci(
    gbk_path: Path,
    locus_tags: list[str],
    upstream_bp: int,
    mode: str = MODE_AUTO,
    warn_missing: bool = True,
    all_isoforms: bool = False,
) -> Iterator[tuple[str, str, str, str, int, int, str]]:
    """Extract upstream regions for a specific list of locus tags.

    The programmatic counterpart to ``extract_regulatory_regions``. While
    ``extract_regulatory_regions`` discovers targets by keyword, this function
    extracts targets whose locus tags are already known — making it the correct
    tool for bridge scripts that receive locus tags from ``pairwise_homolog_finder``.

    Supports both prokaryotic and eukaryotic GenBank files via the ``mode``
    parameter. When ``mode="auto"`` (default), organism type is detected
    automatically from the file.

    Prokaryote mode:
        Searches ``CDS`` features for matching locus tags. Extracts upstream
        sequence relative to the CDS start coordinate (= ATG).

    Eukaryote mode:
        Uses a two-pass strategy per record:
        Pass 1 — ``CDS`` features: builds a {locus_tag: product} map for
                 all target loci. This populates product names since mRNA
                 features often lack /product qualifiers.
        Pass 2 — ``mRNA`` features: resolves the 5'-most TSS across all
                 isoforms, then extracts upstream of that coordinate. Pass
                 ``all_isoforms=True`` to instead yield one region per
                 individual isoform — see ``ALL-ISOFORMS MODE`` in the
                 module docstring.

    Args:
        gbk_path:     Path to the GenBank file to scan.
        locus_tags:   Locus tags to extract. Duplicates are silently removed.
        upstream_bp:  Number of bases to extract upstream of the anchor coordinate.
        mode:         Organism mode: ``"auto"``, ``"prokaryote"``, or ``"eukaryote"``.
        warn_missing: If ``True`` (default), print a warning to stderr listing
                      any locus tags not found in this file after scanning
                      completes. Set to ``False`` when calling from a pipeline
                      that scans multiple files — locus tags absent from one
                      file are expected if they belong to a different genome.
                      The pipeline should track missing tags globally across
                      all files and report once at the end.
        all_isoforms: If ``True``, yield one upstream region per mRNA isoform
                      instead of merging to the 5'-most TSS per locus. No
                      effect in prokaryote mode (no isoforms to merge).

    Yields:
        A 7-item tuple:
        ``(seq_id, locus_tag, product, upstream_seq, actual_upstream, strand, genome_label)``

        - seq_id:          Contig/record ID from the GenBank file.
        - locus_tag:       The matched locus tag (or, with ``all_isoforms=True``,
                           the disambiguated ``{locus_tag}#{transcript_id}``
                           / ``{locus_tag}#isoformN`` identifier).
        - product:         The /product annotation (from CDS in eukaryote mode).
        - upstream_seq:    Strand-corrected upstream DNA sequence.
        - actual_upstream: Actual extracted length; may be < upstream_bp near
                           contig boundaries.
        - strand:          1 for forward strand, -1 for reverse strand.
        - genome_label:    Stem of the GenBank filename for FASTA header use.

    Raises:
        ValueError: If the GenBank file cannot be parsed.

    Notes:
        - If ``warn_missing=True`` and any requested locus tag is not found,
          a warning is printed to stderr after scanning, listing all missing tags.
        - Duplicate locus tags in the input list are silently deduplicated.
        - If the same locus tag appears on multiple CDS features in prokaryote
          mode (malformed annotation), only the first occurrence is yielded
          and a warning is printed.
        - In eukaryote mode, multiple mRNA features per locus (alternative
          isoforms) are handled by selecting the 5'-most TSS by default, or
          yielded individually when ``all_isoforms=True`` — never just the
          first feature seen.
        - Scanning stops early once all requested locus tags have been found.
    """
    resolved_mode = _detect_organism_mode(gbk_path) if mode == MODE_AUTO else mode

    target_set: set[str] = set(locus_tags)
    remaining: set[str] = set(locus_tags)
    already_yielded: set[str] = set()
    genome_label = gbk_path.stem

    try:
        with open(gbk_path, "r", encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, "genbank"):

                if not remaining:
                    break

                if resolved_mode == MODE_PROKARYOTE:
                    # ── Prokaryote: scan CDS features ─────────────────────────
                    for feature in record.features:
                        if not remaining:
                            break
                        if feature.type != "CDS":
                            continue

                        locus_tag = _resolve_identifier(feature, record.id)
                        if locus_tag not in target_set:
                            continue

                        if locus_tag in already_yielded:
                            print(
                                f"      [!] Warning: {locus_tag} appears on multiple "
                                f"CDS features in {gbk_path.name} — skipping duplicate.",
                                file=sys.stderr,
                            )
                            continue

                        product = feature.qualifiers.get(
                            "product", ["Unknown product"]
                        )[0]
                        start = int(feature.location.start)
                        end = int(feature.location.end)
                        strand = feature.location.strand

                        upstream_seq, actual_upstream = _extract_upstream_seq(
                            record, start, end, strand, upstream_bp, locus_tag
                        )

                        already_yielded.add(locus_tag)
                        remaining.discard(locus_tag)

                        yield (
                            record.id,
                            locus_tag,
                            product,
                            upstream_seq,
                            actual_upstream,
                            strand,
                            genome_label,
                        )

                else:
                    # ── Eukaryote: two-pass per record ────────────────────────
                    # Pass 1 — collect product names from CDS features.
                    # mRNA features often lack /product, so we populate product
                    # names here and use them when yielding from mRNA features.
                    #
                    # fallback_cds_spans tracks the CDS coordinate span for every
                    # requested ID that came from the coordinate fallback (no
                    # locus_tag/protein_id/gene at all) — needed in Pass 2 to
                    # correlate that gene's mRNA feature, whose own fallback ID
                    # would be computed from different coordinates (UTRs) and
                    # would never string-match this one. This is the exact
                    # mechanism that lets this function find a locus that
                    # pairwise_homolog_finder.py resolved via its own coordinate
                    # fallback — see _find_overlapping_fallback_id()'s docstring.
                    cds_products: dict[str, str] = {}
                    fallback_cds_spans: dict[str, tuple[int, int]] = {}
                    for feature in record.features:
                        if feature.type != "CDS":
                            continue
                        lt = _resolve_identifier(feature, record.id)
                        if lt in remaining:
                            cds_products[lt] = feature.qualifiers.get(
                                "product", ["Unknown product"]
                            )[0]
                            if lt.startswith("UNANNOTATED_"):
                                fallback_cds_spans[lt] = (
                                    int(feature.location.start),
                                    int(feature.location.end),
                                )

                    # Pass 2 — find the 5'-most TSS per locus across all isoforms.
                    # In eukaryotes, a single locus_tag can have multiple mRNA
                    # features (alternative splice variants). Each may have a
                    # different TSS. We must find the 5'-most TSS — the one
                    # farthest upstream — which represents the full regulatory region.
                    #   Forward strand (+): 5'-most TSS = smallest start coordinate.
                    #   Reverse strand (-): 5'-most TSS = largest end coordinate.
                    # The already_yielded guard below prevents re-processing a locus
                    # that was found in a previous record (should not happen in
                    # well-formed files, but guards against malformed ones).
                    if all_isoforms:
                        # ── ALL-ISOFORMS MODE (v1.7.0): yield one region per
                        # individual mRNA feature instead of merging to the
                        # 5'-most TSS — see the module docstring's
                        # ALL-ISOFORMS MODE section. Each isoform gets a
                        # disambiguated identifier, mirroring
                        # pairwise_homolog_finder.py's --keep-all-isoforms
                        # format exactly. The ORIGINAL locus_tag is only
                        # marked found/removed from `remaining` after every
                        # isoform in THIS record has been yielded, so a
                        # multi-isoform locus split across feature order
                        # within one record is never partially missed.
                        isoform_counters: dict[str, int] = {}
                        found_this_locus: set[str] = set()
                        for feature in record.features:
                            if feature.type != "mRNA":
                                continue
                            start = int(feature.location.start)
                            end = int(feature.location.end)
                            strand = feature.location.strand

                            locus_tag = _resolve_identifier(feature, record.id)
                            if locus_tag not in target_set:
                                if locus_tag.startswith("UNANNOTATED_"):
                                    matched = _find_overlapping_fallback_id(
                                        start, end, fallback_cds_spans
                                    )
                                    if matched is None:
                                        continue
                                    locus_tag = matched
                                else:
                                    continue
                            if locus_tag in already_yielded:
                                continue  # processed in a prior record

                            isoform_key = disambiguate_isoform_id(
                                locus_tag,
                                feature,
                                isoform_counters,
                                id_qualifier="transcript_id",
                            )

                            product = cds_products.get(locus_tag, "Unknown product")
                            upstream_seq, actual_upstream = _extract_upstream_seq(
                                record, start, end, strand, upstream_bp, isoform_key
                            )
                            found_this_locus.add(locus_tag)
                            yield (
                                record.id,
                                isoform_key,
                                product,
                                upstream_seq,
                                actual_upstream,
                                strand,
                                genome_label,
                            )

                        for lt in found_this_locus:
                            already_yielded.add(lt)
                            remaining.discard(lt)
                        continue

                    locus_tss: dict[str, dict] = {}

                    for feature in record.features:
                        if feature.type != "mRNA":
                            continue
                        start = int(feature.location.start)
                        end = int(feature.location.end)
                        strand = feature.location.strand

                        locus_tag = _resolve_identifier(feature, record.id)
                        if locus_tag not in target_set:
                            # Direct match failed. If this mRNA itself has no
                            # locus_tag/protein_id/gene either, try coordinate
                            # containment against a known fallback CDS span
                            # for one of the requested target IDs.
                            if locus_tag.startswith("UNANNOTATED_"):
                                matched = _find_overlapping_fallback_id(
                                    start, end, fallback_cds_spans
                                )
                                if matched is None:
                                    continue
                                locus_tag = matched
                            else:
                                continue
                        if locus_tag in already_yielded:
                            continue  # processed in a prior record

                        if locus_tag not in locus_tss:
                            locus_tss[locus_tag] = {
                                "start": start,
                                "end": end,
                                "strand": strand,
                                "n_isoforms": 1,
                            }
                        else:
                            entry = locus_tss[locus_tag]
                            entry["n_isoforms"] += 1
                            if strand == 1:
                                entry["start"] = min(entry["start"], start)
                            else:
                                entry["end"] = max(entry["end"], end)

                    # Yield one upstream region per locus using the resolved TSS
                    for locus_tag, bounds in locus_tss.items():
                        if bounds["n_isoforms"] > 1:
                            print(
                                f"      [Multi-isoform] {locus_tag}: "
                                f"{bounds['n_isoforms']} mRNA isoforms found — "
                                f"using 5'-most TSS.",
                                file=sys.stderr,
                            )
                        product = cds_products.get(locus_tag, "Unknown product")
                        upstream_seq, actual_upstream = _extract_upstream_seq(
                            record,
                            bounds["start"],
                            bounds["end"],
                            bounds["strand"],
                            upstream_bp,
                            locus_tag,
                        )
                        already_yielded.add(locus_tag)
                        remaining.discard(locus_tag)
                        yield (
                            record.id,
                            locus_tag,
                            product,
                            upstream_seq,
                            actual_upstream,
                            bounds["strand"],
                            genome_label,
                        )

    except Exception as e:
        raise ValueError(f"Failed to parse {gbk_path.name}: {e}") from e

    if remaining and warn_missing:
        print(
            f"\n      [!] Warning: {len(remaining)} locus tag(s) not found in "
            f"{gbk_path.name}:\n"
            + "\n".join(f"              - {tag}" for tag in sorted(remaining)),
            file=sys.stderr,
        )


# ── CLI ────────────────────────────────────────────────────────────────────────


def get_args() -> argparse.Namespace:
    """Configure the CLI and return parsed arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Extract upstream regulatory regions from prokaryotic or eukaryotic "
            "GenBank files for MEME motif discovery."
        ),
    )

    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        default=Path("."),
        help="Input GenBank file OR a directory to scan (Default: current directory)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("upstream_MEME.fasta"),
        help="Output FASTA file name (Default: upstream_MEME.fasta)",
    )
    parser.add_argument(
        "-u",
        "--upstream",
        type=int,
        default=150,
        help=(
            "Number of upstream base pairs to extract. "
            "Default: 150 (prokaryotes). For eukaryotes, 500-2000bp is typical."
        ),
    )
    parser.add_argument(
        "-k",
        "--keywords",
        type=str,
        nargs="+",
        required=True,
        help=(
            "Keywords to search for in /product annotations (case-insensitive). "
            "Example: -k bacteriocin lactobin cerein"
        ),
    )
    parser.add_argument(
        "--mode",
        choices=_VALID_MODES,
        default=MODE_AUTO,
        help=(
            "Organism mode. "
            "'auto' (default): detects organism type from each file automatically. "
            "'prokaryote': extracts upstream of CDS start (ATG). "
            "'eukaryote': extracts upstream of mRNA start (TSS). "
            "Use explicit mode when auto-detection is unreliable."
        ),
    )
    parser.add_argument(
        "--all-isoforms",
        action="store_true",
        default=False,
        help=(
            "Extract one upstream region PER ISOFORM instead of merging all "
            "isoforms of a locus into a single 5'-most-TSS region (the "
            "default). The default is a conservative, deterministic choice "
            "but not necessarily the biologically dominant isoform — a "
            "rare transcript variant with an unusually far-upstream TSS "
            "would be chosen over a well-supported canonical isoform with "
            "a closer TSS. Each isoform gets a disambiguated identifier "
            "({locus_tag}#{transcript_id}, or {locus_tag}#isoformN without "
            "a transcript_id). No effect on prokaryotic files. Off by "
            "default — existing output is unchanged unless you opt in."
        ),
    )

    return parser.parse_args()


# ── Main ───────────────────────────────────────────────────────────────────────


def main() -> None:
    """Coordinate file routing, deduplication, header sanitization, and output."""
    args = get_args()

    print(f"[*] Scanning target    : {args.input}", file=sys.stderr)
    print(f"[*] Upstream window    : {args.upstream}bp", file=sys.stderr)
    print(f"[*] Keywords           : {args.keywords}", file=sys.stderr)
    print(f"[*] Mode               : {args.mode}", file=sys.stderr)
    if args.all_isoforms:
        print(
            f"[*] Isoforms           : ALL (one region per isoform, not merged)",
            file=sys.stderr,
        )
    print(file=sys.stderr)

    hits_found = 0
    duplicates_skipped = 0
    seen_loci = set()

    try:
        with open(args.output, "w", encoding="utf-8") as out_file:

            for file_path in stream_reference_files(args.input):

                if file_path.suffix.lower() in (".fasta", ".fa", ".faa"):
                    print(
                        f"  [!] Skipping {file_path.name}: "
                        f"Cannot extract upstream DNA from FASTA format.",
                        file=sys.stderr,
                    )
                    continue

                # Resolve and report detected mode when auto
                if args.mode == MODE_AUTO:
                    detected = _detect_organism_mode(file_path)
                    print(
                        f"  -> Parsing {file_path.name}... "
                        f"[auto-detected: {detected}]",
                        file=sys.stderr,
                    )
                    file_mode = detected
                else:
                    print(f"  -> Parsing {file_path.name}...", file=sys.stderr)
                    file_mode = args.mode

                for seq_id, locus, prod, seq, actual_up in extract_regulatory_regions(
                    file_path,
                    args.keywords,
                    args.upstream,
                    mode=file_mode,
                    all_isoforms=args.all_isoforms,
                ):
                    dedup_key = (file_path.stem, locus)
                    if dedup_key in seen_loci:
                        duplicates_skipped += 1
                        continue

                    seen_loci.add(dedup_key)
                    hits_found += 1

                    clean_prod = re.sub(r"[^\w\-]", "_", prod)
                    fasta_header = f">{seq_id}_{locus}_{clean_prod}_up{actual_up}"

                    out_file.write(f"{fasta_header}\n{seq}\n")
                    print(f"      [Hit] {locus} | {prod[:50]}", file=sys.stderr)

        print("\n" + "=" * 50, file=sys.stderr)
        print(
            f"[*] SUCCESS: {hits_found} unique regulatory regions extracted.",
            file=sys.stderr,
        )
        if duplicates_skipped > 0:
            print(
                f"[*] WARNING: {duplicates_skipped} duplicate(s) skipped.",
                file=sys.stderr,
            )
        print(f"[*] Output saved to: {args.output.resolve()}", file=sys.stderr)
        print("=" * 50, file=sys.stderr)

    except ValueError as e:
        sys.exit(f"\n[!] Pipeline Error: {e}")
    except KeyboardInterrupt:
        sys.exit("\n[!] Pipeline gracefully interrupted by user.")


if __name__ == "__main__":
    main()
