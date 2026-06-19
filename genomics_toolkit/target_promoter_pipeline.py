#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

r"""Targeted Promoter Pipeline

A bridge script that connects pairwise_homolog_finder.py and
universal_promoter_extractor.py into a single automated workflow.

Replaces homology_extractor.py, which used exact substring matching
(if peptide in translation) — a method that silently misses distant
orthologs with even a single amino acid substitution.

This pipeline uses real sequence alignment (Smith-Waterman / BLOSUM62)
to identify true evolutionary orthologs, then extracts their upstream
regulatory regions for MEME motif discovery. Works on both prokaryotic
and eukaryotic genomes without requiring separate scripts.

WORKFLOW:
    Step 1  Load query proteins from a GenBank file (always treated as prokaryotic
            for protein extraction; query itself can be eukaryotic, only the protein
            extraction logic matters, not the coordinate system).
    Step 2  Align them against a target genome using Smith-Waterman.
    Step 3  Collect the locus tags of all hits above identity/coverage thresholds.
    Step 4  Extract upstream regions for those loci using the target genome's mode
            (auto-detected from each reference file).
    Step 5  (Optional, --dedup-identity) Collapse near-identical upstream DNA
            sequences down to one representative each, so paralogs/recent gene
            duplications don't masquerade as independent observations.

WHY SEQUENCE-LEVEL REDUNDANCY MATTERS (--dedup-identity):
    Step 3's locus-tag dedup only removes hits that point to the exact same
    locus_tag (e.g. two query proteins both matching the same reference gene).
    It does NOT catch the more common case: two or three DIFFERENT loci whose
    upstream DNA is still 90%+ identical because they are paralogs from a
    recent duplication, or near-identical strain variants pulled in from
    multiple reference files.

    If those near-duplicate sequences are fed into motif_discovery.py (or
    MEME) as if they were independent observations, the same regulatory
    signal gets counted multiple times. This inflates apparent motif
    confidence/significance without adding real evidence — a problem if
    you intend to report that significance in a manuscript.

    --dedup-identity THRESHOLD performs greedy single-linkage clustering on
    the extracted upstream DNA sequences themselves (local pairwise
    alignment, identity = matches / alignment_length, BLAST/EMBOSS
    convention — same definition used by pairwise_homolog_finder.py).
    Sequences arrive already ordered by ortholog-hit identity (best matches
    to the query first), so the first sequence in each redundant cluster —
    the best-supported one — is kept as the representative; the rest are
    logged to <stem>.redundancy.tsv (locus, matched representative, percent
    identity) for full traceability, not silently discarded.

    Off by default — existing output is unchanged unless you opt in.

ORGANISM SUPPORT:
    Prokaryotic target genomes:
        Extracts upstream of CDS start (= ATG / translation start).
        Auto-detected: CDS-only GenBank files = prokaryote.

    Eukaryotic target genomes:
        Extracts upstream of mRNA start (= Transcription Start Site = TSS).
        Auto-detected: GenBank files with mRNA features = eukaryote.
        Keyword matching still uses /product from CDS features (which always
        have annotations). Coordinate extraction uses mRNA features (which
        have correct TSS coordinates but often lack /product).

OUTPUTS:
    <stem>.fasta       MEME-ready upstream sequences (wrapped at 60 chars).
                       Contains only kept (non-redundant) sequences when
                       --dedup-identity is used.
    <stem>.tsv         Extraction results table: one row per upstream region,
                       with columns for locus_tag, genome, contig, product,
                       strand, and upstream length.
    <stem>.hits.tsv    Ortholog alignment hits table: one row per hit from
                       Step 2, with query/ref locus, identity, alignment
                       length, and source genome. Full version of the
                       truncated terminal display.
    <stem>.redundancy.tsv  (Only with --dedup-identity) One row per DROPPED
                       locus: which representative locus it matched and at
                       what percent identity. Not written if no sequences
                       were dropped.

WHY THIS REPLACES homology_extractor.py:
    The old script used:
        if core_peptide in translation:
    This is exact substring matching — one amino acid substitution and
    the ortholog is silently missed. Real homology requires alignment,
    which is what pairwise_homolog_finder.py provides.

Note:
    Associated with ongoing, unpublished research (manuscript in
    preparation). Correct attribution is requested when used in
    derivative works.

Examples:
    # Prokaryotic target (auto-detected)
    $ python3 target_promoter_pipeline.py \
        -q C5_genome.gbk \
        -r GCF_000014445_1_genomic.gbff \
        -o atcc8293_targeted_promoters.fasta \
        -u 150 \
        --identity 0.35

    # Eukaryotic target (auto-detected)
    $ python3 target_promoter_pipeline.py \
        -q Arabidopsis_proteins.gbk \
        -r Zea_mays.gbff \
        -o maize_promoters.fasta \
        -u 1000 \
        --identity 0.40 \
        --coverage 0.70

    # Explicit mode (when auto-detection is unreliable)
    $ python3 target_promoter_pipeline.py \
        -q query.gbk -r target.gbff --mode eukaryote -u 1500

    # Mixed directory of genomes (each auto-detected independently)
    $ python3 target_promoter_pipeline.py \
        -q C5_genome.gbk \
        -r references/ \
        -o all_promoters.fasta \
        --mature --max-length 150 --identity 0.40

    # Collapse near-identical paralog promoters before motif discovery
    # (recommended whenever output will feed motif_discovery.py / MEME)
    $ python3 target_promoter_pipeline.py \
        -q C5_genome.gbk \
        -r references/ \
        -o all_promoters.fasta \
        --dedup-identity 0.90
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.4.0"

import sys
import argparse
from pathlib import Path
from dataclasses import dataclass

try:
    from Bio.Align import PairwiseAligner
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )

# ── Import from existing toolkit scripts ──────────────────────────────────────
# These are your two production scripts. This bridge adds NO new alignment or
# extraction logic — it only connects them. Any improvements made to either
# script are automatically inherited here.

try:
    from pairwise_homolog_finder import extract_proteins_from_gbk, find_homologs
except ImportError as e:
    sys.exit(
        f"[!] Cannot import from pairwise_homolog_finder.py.\n"
        f"    Ensure it is in the same directory as this script.\n"
        f"    Details: {e}"
    )

try:
    from universal_promoter_extractor import extract_by_loci
except ImportError as e:
    sys.exit(
        f"[!] Cannot import from universal_promoter_extractor.py.\n"
        f"    Ensure it is in the same directory as this script.\n"
        f"    Details: {e}"
    )

try:
    from utils import wrap_fasta
except ImportError as e:
    sys.exit(
        f"[!] Cannot import from utils.py.\n"
        f"    Ensure it is in the same directory as this script.\n"
        f"    Details: {e}"
    )

# ── Module-level DNA aligner (created ONCE, reused for all dedup comparisons) ──
# Local alignment (Smith-Waterman) with simple match/mismatch scoring — there is
# no standard substitution matrix for raw nucleotides the way BLOSUM62 exists
# for proteins, so a flat +1/-1 scheme is the conventional choice (same
# approach EMBOSS/water uses for DNA-DNA identity). Only used when
# --dedup-identity is supplied; otherwise never instantiated/called.
_DNA_ALIGNER = PairwiseAligner()
_DNA_ALIGNER.mode = "local"
_DNA_ALIGNER.match_score = 1
_DNA_ALIGNER.mismatch_score = -1
_DNA_ALIGNER.open_gap_score = -2
_DNA_ALIGNER.extend_gap_score = -1


# ─────────────────────────────────────────────────────────────────────────────
# REDUNDANCY FILTERING (--dedup-identity)
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class PromoterRecord:
    """One extracted upstream region, buffered before output so the optional
    redundancy filter can compare all of them against each other first."""

    locus: str
    seq_id: str
    genome_label: str
    product: str
    seq: str
    actual_up: int
    strand_symbol: str


def _dna_identity(seq_a: str, seq_b: str) -> float:
    """Computes percent identity between two DNA sequences via local alignment.

    Mirrors the identity definition in pairwise_homolog_finder.calculate_identity()
    (identical positions / alignment_length, the BLAST/EMBOSS convention) so
    --dedup-identity thresholds are directly comparable to --identity thresholds
    elsewhere in this toolkit.

    Args:
        seq_a: First DNA sequence.
        seq_b: Second DNA sequence.

    Returns:
        Percent identity as a 0.0-1.0 fraction. Returns 0.0 for empty input
        or if no alignment is found.
    """
    if not seq_a or not seq_b:
        return 0.0

    alignments = _DNA_ALIGNER.align(seq_a, seq_b)
    best = next(iter(alignments), None)
    if best is None:
        return 0.0

    aligned_a, aligned_b = best[0], best[1]
    alignment_length = len(aligned_a)
    if alignment_length == 0:
        return 0.0

    identical = sum(1 for a, b in zip(aligned_a, aligned_b) if a == b and a != "-")
    return identical / alignment_length


def deduplicate_promoters(
    records: list[PromoterRecord], threshold: float
) -> tuple[list[PromoterRecord], list[tuple[PromoterRecord, str, float]]]:
    """Greedy single-linkage clustering on extracted upstream DNA sequences.

    Processes records in the order given (which, as called from main(), is
    ortholog-hit-identity order — best matches to the query first). For each
    record, compares it against every representative kept so far; if its
    identity to ANY existing representative meets or exceeds the threshold,
    it is dropped as redundant and the first (best-supported) representative
    is kept. This means cluster representatives are always the most
    query-relevant member of their redundancy group, not an arbitrary pick.

    This is O(n^2) pairwise alignments. Fine for the dozens-to-low-hundreds
    of loci typical of a targeted promoter search; would need a smarter
    approach (e.g. k-mer pre-filtering like pairwise_homolog_finder.py uses)
    if ever applied to thousands of sequences.

    Args:
        records: Extracted promoter records, in priority order.
        threshold: Minimum identity (0.0-1.0) to consider two sequences
            redundant.

    Returns:
        A tuple of (kept_records, dropped) where dropped is a list of
        (record, matched_representative_locus, identity) for every record
        that was removed as redundant.
    """
    kept: list[PromoterRecord] = []
    dropped: list[tuple[PromoterRecord, str, float]] = []

    for rec in records:
        best_identity = 0.0
        best_rep: PromoterRecord | None = None
        for rep in kept:
            identity = _dna_identity(rec.seq, rep.seq)
            if identity > best_identity:
                best_identity = identity
                best_rep = rep

        if best_rep is not None and best_identity >= threshold:
            dropped.append((rec, best_rep.locus, best_identity))
        else:
            kept.append(rec)

    return kept, dropped


# ── Argument parsing ───────────────────────────────────────────────────────────


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="target_promoter_pipeline.py",
        description=(
            "Find true orthologs via Smith-Waterman alignment and extract "
            "their upstream regulatory regions in one automated pass."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  python3 target_promoter_pipeline.py -q C5_genome.gbk -r ATCC8293.gbff
  python3 target_promoter_pipeline.py -q C5_genome.gbk -r references/ --mature --max-length 150
  python3 target_promoter_pipeline.py -q C5_genome.gbk -r NBRC.gbff --identity 0.70 --coverage 0.80
        """,
    )

    parser.add_argument(
        "-q",
        "--query",
        type=Path,
        required=True,
        metavar="GBK",
        help="Query GenBank file containing the target proteins to search for.",
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=Path,
        required=True,
        metavar="GBK_OR_DIR",
        help=(
            "Target genome GenBank file (or directory of .gbk/.gbff files) "
            "to search inside."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        metavar="FASTA",
        help=(
            "Output FASTA file path. "
            "If omitted, saves to 'targeted_promoters.fasta' in the current directory. "
            "Two TSV files are automatically generated alongside it: "
            "<stem>.tsv (extraction results) and <stem>.hits.tsv (ortholog hits)."
        ),
    )
    parser.add_argument(
        "-u",
        "--upstream",
        type=int,
        default=150,
        metavar="BP",
        help="Upstream base pairs to extract per locus. Default: 150",
    )
    parser.add_argument(
        "--identity",
        type=float,
        default=0.35,
        metavar="FLOAT",
        help=(
            "Minimum Smith-Waterman alignment identity to report a hit (0.0-1.0). "
            "Default: 0.35 (35%%). Increase for stricter ortholog definition."
        ),
    )
    parser.add_argument(
        "--coverage",
        type=float,
        default=0.65,
        metavar="FLOAT",
        help=(
            "Minimum alignment coverage fraction (0.0-1.0). Default: 0.65 (65%%). "
            "Uses 'min' mode: shorter sequence as denominator (correct for "
            "bacteriocin domain searches)."
        ),
    )
    parser.add_argument(
        "--mature",
        action="store_true",
        default=False,
        help=(
            "Apply mature core trimming (removes signal/leader peptides) before "
            "alignment. Recommended for bacteriocins. Uses calculate_mature_core() "
            "from utils.py."
        ),
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=None,
        metavar="AA",
        help=(
            "Skip query proteins longer than this many amino acids. "
            "Useful to filter out large unrelated proteins when searching "
            "for short peptides like bacteriocins (e.g. --max-length 150)."
        ),
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=10,
        metavar="AA",
        help="Skip query proteins shorter than this many amino acids. Default: 10",
    )
    parser.add_argument(
        "--mode",
        choices=["auto", "prokaryote", "eukaryote"],
        default="auto",
        help=(
            "Organism mode for upstream extraction. "
            "'auto' (default): detects from each reference file automatically. "
            "'prokaryote': extracts upstream of CDS start (ATG). "
            "'eukaryote': extracts upstream of mRNA start (TSS). "
            "The query file is always treated as prokaryotic for protein extraction."
        ),
    )
    parser.add_argument(
        "--dedup-identity",
        type=float,
        default=None,
        metavar="FLOAT",
        help=(
            "Collapse near-identical upstream DNA sequences before writing "
            "output (0.0-1.0). If a promoter is >= this fraction identical "
            "to an already-kept promoter — e.g. a recent paralog or "
            "duplicated gene — it is dropped and logged to "
            "<stem>.redundancy.tsv; only the best-supported representative "
            "is kept. Recommended before motif discovery (MEME / "
            "motif_discovery.py), since near-duplicate sequences are not "
            "statistically independent observations and can inflate "
            "apparent motif significance. Off by default — suggested "
            "value: 0.90."
        ),
    )

    return parser


# ── Main ──────────────────────────────────────────────────────────────────────


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    # Maximum number of hits/loci printed to the terminal.
    # Beyond this, a truncation message directs the user to the TSV.
    _MAX_DISPLAY = 20

    output = args.output or Path("targeted_promoters.fasta")
    tsv_output = output.with_suffix(".tsv")
    hits_tsv = output.with_stem(output.stem + ".hits").with_suffix(".tsv")
    redundancy_tsv = output.with_stem(output.stem + ".redundancy").with_suffix(".tsv")

    # Validate inputs
    if not args.query.exists():
        sys.exit(f"[!] Query file not found: {args.query}")
    if not args.reference.exists():
        sys.exit(f"[!] Reference path not found: {args.reference}")
    if args.dedup_identity is not None and not (0.0 <= args.dedup_identity <= 1.0):
        sys.exit(
            f"[!] --dedup-identity must be between 0.0 and 1.0 "
            f"(got {args.dedup_identity})."
        )

    print("=" * 60, file=sys.stderr)
    print("  TARGETED PROMOTER PIPELINE", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    print(f"  Query           : {args.query.name}", file=sys.stderr)
    print(f"  Reference       : {args.reference}", file=sys.stderr)
    print(f"  FASTA output    : {output}", file=sys.stderr)
    print(f"  TSV output      : {tsv_output}", file=sys.stderr)
    print(f"  Hits TSV        : {hits_tsv}", file=sys.stderr)
    print(f"  Upstream        : {args.upstream}bp", file=sys.stderr)
    print(f"  Min identity    : {args.identity * 100:.0f}%", file=sys.stderr)
    print(f"  Min coverage    : {args.coverage * 100:.0f}% (min mode)", file=sys.stderr)
    print(
        f"  Mature core     : {'YES — leader peptides trimmed' if args.mature else 'NO'}",
        file=sys.stderr,
    )
    print(f"  Organism mode   : {args.mode}", file=sys.stderr)
    print(
        f"  Dedup identity  : "
        f"{f'{args.dedup_identity * 100:.0f}% (DNA, local alignment)' if args.dedup_identity is not None else 'OFF'}",
        file=sys.stderr,
    )
    print("=" * 60, file=sys.stderr)

    # ── Step 1: Load query proteins ───────────────────────────────────────────
    print(
        f"\n[Step 1] Loading query proteins from {args.query.name}...", file=sys.stderr
    )

    query_proteins = extract_proteins_from_gbk(
        gbk_path=args.query,
        apply_mature=args.mature,
        min_length=args.min_length,
        max_length=args.max_length,
        verbose=False,
    )

    if not query_proteins:
        sys.exit(
            "[!] No proteins extracted from query.\n"
            "    Check that the file is a valid GenBank with /translation= fields.\n"
            "    If using --mature, check that mature cores are not zero-length."
        )

    print(f"    {len(query_proteins)} protein(s) loaded.", file=sys.stderr)

    # ── Step 2: Find orthologs via Smith-Waterman ─────────────────────────────
    print(
        f"\n[Step 2] Searching for orthologs via Smith-Waterman alignment...",
        file=sys.stderr,
    )

    # find_homologs() in pairwise_homolog_finder.py calls SeqIO.parse() directly
    # and expects a single file path — it does not handle directories.
    # We resolve the directory here and call find_homologs() per file,
    # then aggregate all hits. This mirrors how Step 3 handles directories.
    ref_files: list[Path] = []
    if args.reference.is_dir():
        for ext in ("*.gbk", "*.gbff"):
            ref_files.extend(sorted(args.reference.rglob(ext)))
        if not ref_files:
            sys.exit(
                f"[!] No GenBank files (.gbk / .gbff) found in: {args.reference}\n"
                f"    Check the directory path and file extensions."
            )
        print(
            f"    Reference is a directory — found {len(ref_files)} file(s).",
            file=sys.stderr,
        )
    else:
        ref_files = [args.reference]

    all_hits = []
    for ref_file in ref_files:
        print(f"    Searching {ref_file.name}...", file=sys.stderr)
        file_hits = find_homologs(
            query_proteins=query_proteins,
            ref_path=ref_file,
            min_identity=args.identity,
            use_mature=args.mature,
            min_coverage=args.coverage,
            coverage_mode="min",  # correct for bacteriocin/domain searches
        )
        if file_hits:
            print(f"      -> {len(file_hits)} hit(s) found.", file=sys.stderr)
            all_hits.extend(file_hits)
        else:
            print(f"      -> No hits above threshold.", file=sys.stderr)

    hits = all_hits

    if not hits:
        sys.exit(
            f"[!] No orthologs found above {args.identity*100:.0f}% identity / "
            f"{args.coverage*100:.0f}% coverage in any reference file.\n"
            f"    Try lowering --identity or --coverage, or check your query file."
        )

    print(
        f"\n    {len(hits)} total ortholog hit(s) across all reference files.",
        file=sys.stderr,
    )

    # Write all hits to .hits.tsv — full record for downstream review.
    # Terminal display is truncated to _MAX_DISPLAY; TSV always has everything.
    HITS_TSV_HEADER = "\t".join(
        [
            "query_locus",
            "query_product",
            "ref_locus",
            "ref_product",
            "identity_pct",
            "alignment_length",
            "mismatches",
            "query_length",
            "ref_length",
            "ref_file",
        ]
    )
    sorted_hits = sorted(hits, key=lambda h: (-h.identity, h.ref_locus))
    with open(hits_tsv, "w", encoding="utf-8") as hf:
        hf.write(HITS_TSV_HEADER + "\n")
        for hit in sorted_hits:
            hf.write(
                "\t".join(
                    [
                        hit.query_locus,
                        hit.query_product,
                        hit.ref_locus,
                        hit.ref_product,
                        f"{hit.identity * 100:.1f}",
                        str(hit.alignment_length),
                        str(hit.mismatches),
                        str(hit.query_length),
                        str(hit.ref_length),
                        hit.ref_file,
                    ]
                )
                + "\n"
            )

    # Terminal: show first _MAX_DISPLAY hits, then a truncation notice
    for hit in sorted_hits[:_MAX_DISPLAY]:
        print(
            f"      {hit.ref_locus:<20} {hit.identity*100:.1f}% identity  "
            f"{hit.ref_product[:45]}",
            file=sys.stderr,
        )
    if len(sorted_hits) > _MAX_DISPLAY:
        hidden = len(sorted_hits) - _MAX_DISPLAY
        print(
            f"      ... and {hidden} more hit(s) not shown. "
            f"See {hits_tsv.name} for the full list.",
            file=sys.stderr,
        )

    # Deduplicate reference locus tags.
    # The same ref_locus can appear in multiple hits when more than one query
    # protein matches it (paralogs in the query genome). We only need to extract
    # upstream of each reference locus once, regardless of how many query
    # proteins pointed to it.
    target_loci = list(dict.fromkeys(hit.ref_locus for hit in sorted_hits))
    print(
        f"\n    {len(target_loci)} unique reference locus tag(s) to extract "
        f"({len(hits) - len(target_loci)} redundant hit(s) collapsed):",
        file=sys.stderr,
    )

    # Terminal: show first _MAX_DISPLAY loci, then a truncation notice
    for tag in target_loci[:_MAX_DISPLAY]:
        print(f"      {tag}", file=sys.stderr)
    if len(target_loci) > _MAX_DISPLAY:
        hidden = len(target_loci) - _MAX_DISPLAY
        print(
            f"      ... and {hidden} more. " f"See {hits_tsv.name} for the full list.",
            file=sys.stderr,
        )

    # ── Step 3: Extract upstream regions ─────────────────────────────────────
    print(
        f"\n[Step 3] Extracting {args.upstream}bp upstream of identified loci...",
        file=sys.stderr,
    )

    # Extraction is buffered into memory (rather than streamed straight to
    # disk, as in prior versions) so that Step 4's optional redundancy filter
    # can compare every extracted sequence against every other one before
    # anything is written. When --dedup-identity is not used, the buffer is
    # written out unfiltered and in the same order as before — output is
    # byte-identical to pre-v1.4.0 behavior.
    extracted: list[PromoterRecord] = []
    found_loci: set[str] = set()  # tracks which loci were found across ALL files

    try:
        for ref_file in ref_files:
            file_count = 0

            for (
                seq_id,
                locus,
                product,
                seq,
                actual_up,
                strand,
                genome_label,
            ) in extract_by_loci(
                ref_file,
                target_loci,
                args.upstream,
                mode=args.mode,
                warn_missing=False,  # suppressed: loci absent from one
            ):  # genome are expected in multi-file runs

                strand_symbol = "+" if strand == 1 else "-"
                extracted.append(
                    PromoterRecord(
                        locus=locus,
                        seq_id=seq_id,
                        genome_label=genome_label,
                        product=product,
                        seq=seq,
                        actual_up=actual_up,
                        strand_symbol=strand_symbol,
                    )
                )
                found_loci.add(locus)
                file_count += 1

            # Per-file summary — one clean line instead of one line per locus
            status = f"{file_count} region(s) extracted" if file_count else "no matches"
            print(f"      {ref_file.name:<45} {status}", file=sys.stderr)

    except ValueError as e:
        sys.exit(f"[!] Extraction error: {e}")

    extracted_count = len(extracted)

    # Consolidated missing-loci report — only loci never found in ANY file
    never_found = sorted(set(target_loci) - found_loci)
    if never_found:
        print(
            f"\n  [!] {len(never_found)} locus tag(s) not found in any reference file:",
            file=sys.stderr,
        )
        for tag in never_found:
            print(f"        - {tag}", file=sys.stderr)

    # ── Step 4: Redundancy filter (optional) ─────────────────────────────────
    dropped: list[tuple[PromoterRecord, str, float]] = []
    if args.dedup_identity is not None:
        print(
            f"\n[Step 4] Checking {extracted_count} extracted sequence(s) for "
            f">= {args.dedup_identity * 100:.0f}% DNA identity redundancy...",
            file=sys.stderr,
        )
        kept, dropped = deduplicate_promoters(extracted, args.dedup_identity)
        print(
            f"    {len(kept)} kept, {len(dropped)} dropped as redundant.",
            file=sys.stderr,
        )
        for rec, rep_locus, identity in dropped[:_MAX_DISPLAY]:
            print(
                f"      {rec.locus:<20} -> redundant with {rep_locus} "
                f"({identity*100:.1f}% identity)",
                file=sys.stderr,
            )
        if len(dropped) > _MAX_DISPLAY:
            print(
                f"      ... and {len(dropped) - _MAX_DISPLAY} more. "
                f"See {redundancy_tsv.name} for the full list.",
                file=sys.stderr,
            )
        extracted = kept

    # ── Write outputs ─────────────────────────────────────────────────────────
    TSV_HEADER = "\t".join(
        [
            "locus_tag",
            "genome_label",
            "contig_id",
            "product",
            "strand",
            "upstream_requested",
            "upstream_extracted",
        ]
    )

    with (
        open(output, "w", encoding="utf-8") as fasta_out,
        open(tsv_output, "w", encoding="utf-8") as tsv_out,
    ):
        tsv_out.write(TSV_HEADER + "\n")
        for rec in extracted:
            header = (
                f">{rec.locus} | {rec.seq_id} | {rec.genome_label} | "
                f"{rec.actual_up}bp upstream | strand {rec.strand_symbol} | "
                f"{rec.product[:45]}"
            )
            fasta_out.write(f"{header}\n{wrap_fasta(rec.seq)}\n")
            tsv_out.write(
                "\t".join(
                    [
                        rec.locus,
                        rec.genome_label,
                        rec.seq_id,
                        rec.product,
                        rec.strand_symbol,
                        str(args.upstream),
                        str(rec.actual_up),
                    ]
                )
                + "\n"
            )

    if dropped:
        REDUNDANCY_HEADER = "\t".join(
            ["dropped_locus", "matched_representative_locus", "identity_pct"]
        )
        with open(redundancy_tsv, "w", encoding="utf-8") as rf:
            rf.write(REDUNDANCY_HEADER + "\n")
            for rec, rep_locus, identity in dropped:
                rf.write(f"{rec.locus}\t{rep_locus}\t{identity*100:.1f}\n")

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"\n{'=' * 60}", file=sys.stderr)
    print(f"  PIPELINE COMPLETE", file=sys.stderr)
    print(f"{'=' * 60}", file=sys.stderr)
    print(f"  Orthologs found    : {len(hits)}", file=sys.stderr)
    print(f"  Unique loci        : {len(target_loci)}", file=sys.stderr)
    print(f"  Regions extracted  : {extracted_count}", file=sys.stderr)
    if args.dedup_identity is not None:
        print(f"  Dropped (redundant): {len(dropped)}", file=sys.stderr)
        print(f"  Written to output  : {len(extracted)}", file=sys.stderr)
    print(f"  FASTA written to   : {output.resolve()}", file=sys.stderr)
    print(f"  TSV written to     : {tsv_output.resolve()}", file=sys.stderr)
    print(f"  Hits TSV written   : {hits_tsv.resolve()}", file=sys.stderr)
    if dropped:
        print(f"  Redundancy TSV     : {redundancy_tsv.resolve()}", file=sys.stderr)
    print(f"{'=' * 60}", file=sys.stderr)

    if extracted_count == 0:
        print(
            "\n[!] WARNING: Orthologs were found but no upstream sequences were extracted.\n"
            "    This can happen if the locus tags in the ortholog hits do not\n"
            "    match the locus tags in the reference GenBank file.\n"
            "    Check that you are using the same GenBank file for both steps.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit("\n[!] Pipeline interrupted by user.")
