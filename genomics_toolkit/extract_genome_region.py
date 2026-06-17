#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

"""Genomic region extractor from GenBank files.

Extracts a user-defined genomic region from a GBFF or GBK file and writes
it in one or more output formats. Works with NCBI GBFF files, antiSMASH
region GBK files, and Prokka-annotated GBK files.

Three extraction modes are supported:

    Coordinate mode (``--c1`` + ``--c2``):
        Provide explicit 1-based genomic start and end positions. These can
        be obtained from ``find_gbk_features.py --context`` or from
        antiSMASH region output.

    Locus tag mode (``--locus FIRST_TAG LAST_TAG``):
        The script locates both locus tags in the file and extracts the
        genomic span between them. Order does not matter. To extract a
        single gene, provide the same tag twice: ``--locus TAG TAG``.

    Whole-sequence mode (``--seq`` with no ``--c1``/``--c2`` or ``--locus``):
        Extracts the entire contig.

CIRCULAR GENOME SUPPORT (``--circular``):
    Bacterial chromosomes and plasmids are circular. If a feature of
    interest spans the sequencer-assigned replication origin, a standard
    linear extraction will either fail (coordinate mode) or extract the
    entire genome minus the gap (locus mode).

    Use ``--circular`` to enable origin-spanning extraction:

    - Coordinate mode: c1 may be greater than c2. The region is assembled
      by concatenating the two flanking segments:
        target[c1-1 : contig_len] + target[0 : c2]

    - Locus mode: both paths (direct and cross-origin) are computed; the
      shorter path is used automatically.

    Note: GBK output in circular mode will have locally adjusted
    coordinates starting at 0 for the first segment. FASTA headers
    always show the original global genome coordinates.

TRUNCATED BOUNDARY FEATURE HANDLING:
    When a region boundary slices through a CDS, Biopython clips the
    feature's DNA coordinates to the extraction window but does NOT
    update the cached ``/translation`` qualifier. This causes a critical
    desync: ``--gene-fna`` contains truncated DNA while ``--faa`` would
    write the original full-length protein.

    This script detects truncated features by comparing the ``/translation``
    string length against the extracted coding sequence length. Truncated
    CDS features are skipped in both ``--faa`` and ``--gene-fna`` output
    with a clear warning. ``--gbk`` retains them with their clipped
    coordinates so the region boundary is visible in genome browsers.

GLOBAL COORDINATE TRACEABILITY:
    FASTA headers always show the feature's original genomic coordinates
    from the source file, NOT the local coordinates of the extracted
    region. This ensures that when a researcher finds an interesting
    protein in a ``.faa`` file, they can immediately locate it in the
    full genome.

    Example: a gene at global position 53,508..54,911 extracted as part
    of a region starting at 53,317 will appear in ``.faa`` and
    ``--gene-fna`` as ``[global_location=53508..54911]``, not
    ``[location=192..1594]``.

PIPELINE INTEGRATION (``--genes-file``):
    The ``--genes`` flag accepts locus tags on the command line. For large
    gene lists (from scanner TSV output), use ``--genes-file`` instead.
    It accepts either a plain text file (one locus tag per line) or the
    TSV files produced by ``conserved_annotation_scanner.py`` and
    ``cross_genome_keyword_scanner.py`` — the ``Locus_Tag`` column is
    automatically detected and extracted.

    This avoids shell ARG_MAX limits and enables direct pipeline integration:
        python3 cross_genome_keyword_scanner.py ... -o hits.tsv
        python3 extract_genome_region.py -i genome.gbk --genes-file hits.tsv --faa out.faa

PERFORMANCE:
    A ``{locus_tag: feature}`` index is built once per genome record in
    O(N). All subsequent locus tag lookups (``--locus``, ``--genes``,
    ``--genes-file``) are O(1), making the script suitable for extracting
    thousands of genes from a single record without redundant iteration.

OUTPUT FORMAT REFERENCE:
    ``--faa FILE.faa``       Protein FASTA (one entry per non-truncated CDS)
    ``--fna FILE.fna``       Region DNA FASTA (one entry, includes intergenic)
    ``--gene-fna FILE.fna``  Gene DNA FASTA (one entry per non-truncated CDS)
    ``--gbk FILE.gbk``       Annotated GenBank with local coordinates

Note:
    Associated with ongoing, unpublished research (manuscript in
    preparation). Correct attribution is requested when used in
    derivative works.

Examples:
    Check what sequences are in the file::

        python3 extract_genome_region.py -i genome.gbff --list-sequences

    Extract by coordinates::

        python3 extract_genome_region.py -i genome.gbff \\
            --seq NZ_CP134351.1 --c1 53317 --c2 78823 \\
            --gbk region.gbk --faa region.faa --fna region.fna

    Extract by locus tag range::

        python3 extract_genome_region.py -i genome.gbff \\
            --seq NZ_CP134351.1 \\
            --locus RHP56_RS00340 RHP56_RS00455 \\
            --faa region.faa --fna region.fna --gene-fna region_genes.fna

    Extract cross-origin region (circular genome)::

        python3 extract_genome_region.py -i genome.gbff \\
            --seq NZ_CP134351.1 --c1 4999000 --c2 100 --circular \\
            --faa wrap_region.faa --gene-fna wrap_genes.fna

    Extract a whole contig::

        python3 extract_genome_region.py -i genome.gbff \\
            --seq NZ_CP134351.1 \\
            --faa all_proteins.faa --fna full_contig.fna

    Extract genes from scanner TSV output (pipeline integration)::

        python3 extract_genome_region.py -i genome.gbk \\
            --genes-file conserved_hits.tsv --faa core_proteins.faa
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "2.2.0"

import sys
import argparse
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )
from utils import wrap_fasta

# ── Feature index ─────────────────────────────────────────────────────────────


def _build_feature_index(record: SeqRecord) -> dict:
    """Build a ``{locus_tag: feature}`` index for O(1) locus tag lookup.

    Constructed once per genome record. All subsequent locus tag lookups
    (``--locus``, ``--genes``, ``--genes-file``) use this index rather than
    iterating over the feature list, reducing O(N×M) to O(N + M).

    Args:
        record: A BioPython SeqRecord.

    Returns:
        Dict mapping locus tag strings to their CDS SeqFeature objects.
    """
    index = {}
    for feature in record.features:
        if feature.type == "CDS":
            lt = feature.qualifiers.get("locus_tag", [""])[0]
            if lt:
                index[lt] = feature
    return index


def _build_coord_index(record: SeqRecord) -> dict:
    """Build a ``{locus_tag: (global_start_1based, global_end_1based, strand)}`` index.

    Records the original pre-slice genomic coordinates for every CDS.
    Passed to FASTA header builders so headers always show global genome
    coordinates, not local coordinates of the extracted region.

    Args:
        record: A BioPython SeqRecord (must be the ORIGINAL, unsliced record).

    Returns:
        Dict mapping locus tag strings to (start, end, strand) tuples.
        ``start`` and ``end`` are 1-based (matching NCBI/GenBank convention).
        ``strand`` is 1 (forward) or -1 (reverse).
    """
    index = {}
    for feature in record.features:
        if feature.type == "CDS":
            lt = feature.qualifiers.get("locus_tag", [""])[0]
            if lt:
                index[lt] = (
                    int(feature.location.start) + 1,  # 1-based
                    int(feature.location.end),  # 1-based end
                    feature.location.strand,
                )
    return index


# ── Truncation detection ──────────────────────────────────────────────────────


def _find_truncated_loci(region: SeqRecord) -> set[str]:
    """Return locus tags of CDS features truncated by slice boundaries.

    After Biopython slices a region, features that crossed the boundary are
    retained but have their DNA coordinates clipped to the slice window.
    Critically, the ``/translation`` qualifier is NOT updated — it still
    holds the original full-length protein sequence. This causes a dangerous
    desync between ``--gene-fna`` (truncated DNA) and ``--faa`` (full protein).

    Detection: compare the ``/translation`` string length against the number
    of codons in the extracted coding sequence. A mismatch > 1 codon (to
    allow for stop codon conventions) indicates truncation.

    Args:
        region: A sliced BioPython SeqRecord.

    Returns:
        Set of locus tag strings for truncated CDS features. These are
        skipped in ``write_faa`` and ``write_gene_fna`` with a warning.
    """
    truncated: set[str] = set()
    for feature in region.features:
        if feature.type != "CDS":
            continue
        translation = feature.qualifiers.get("translation", [""])[0]
        if not translation:
            continue
        coding_seq = str(feature.location.extract(region.seq))
        expected_aa = len(coding_seq) // 3
        if abs(len(translation) - expected_aa) > 1:
            lt = feature.qualifiers.get("locus_tag", ["?"])[0]
            truncated.add(lt)
    return truncated


# ── CLI ───────────────────────────────────────────────────────────────────────


def get_args() -> argparse.Namespace:
    """Configures the CLI parser and returns parsed arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Extract a genomic region from a GBFF or GBK file in one or more "
            "output formats. Three extraction modes: coordinate (--c1/--c2), "
            "locus tag (--locus FIRST LAST), or whole-sequence (--seq alone). "
            "Use --circular for origin-spanning regions on circular genomes."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # ── Input ─────────────────────────────────────────────────────────────────
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Input GenBank file (.gbff or .gbk).",
    )
    parser.add_argument(
        "--seq",
        type=str,
        default=None,
        help=(
            "Sequence ID to extract from (e.g. 'NZ_CP134351.1'). "
            "Required for multi-record GBFF files; optional for single-record files. "
            "Run --list-sequences to find available IDs."
        ),
    )
    parser.add_argument(
        "--circular",
        action="store_true",
        default=False,
        help=(
            "Enable circular genome mode. Allows c1 > c2 in coordinate mode "
            "(cross-origin extraction). In locus mode, automatically selects "
            "the shorter of the two possible paths around the chromosome. "
            "The region is assembled as target[c1-1:contig_len] + target[0:c2]."
        ),
    )

    # ── Extraction modes ──────────────────────────────────────────────────────
    mode = parser.add_argument_group(
        "Extraction Modes",
        "Choose one mode. If none is given, the entire contig is extracted.",
    )
    mode.add_argument(
        "--c1",
        type=int,
        default=None,
        metavar="START_BP",
        help="Coordinate-mode start position (1-based, inclusive). Requires --c2.",
    )
    mode.add_argument(
        "--c2",
        type=int,
        default=None,
        metavar="END_BP",
        help=(
            "Coordinate-mode end position (1-based, inclusive). Requires --c1. "
            "May be less than --c1 when --circular is set (cross-origin extraction)."
        ),
    )
    mode.add_argument(
        "--locus",
        type=str,
        nargs=2,
        metavar=("FIRST_TAG", "LAST_TAG"),
        default=None,
        help=(
            "Locus tag mode: extract the genomic span from FIRST_TAG to LAST_TAG. "
            "Tag order does not matter. Same tag twice extracts a single gene. "
            "Cannot be combined with --c1/--c2 or --genes/--genes-file."
        ),
    )
    mode.add_argument(
        "--genes",
        type=str,
        nargs="+",
        metavar="TAG",
        default=None,
        help=(
            "Extract multiple individual genes by locus tag. "
            "Provide one or more tags (e.g. --genes ctg1_47 ctg1_58 ctg1_74). "
            "Cannot be combined with --c1/--c2 or --locus. "
            "For large gene lists, prefer --genes-file."
        ),
    )
    mode.add_argument(
        "--genes-file",
        type=Path,
        default=None,
        metavar="FILE",
        help=(
            "Read locus tags from a file instead of the command line. "
            "Accepts either a plain text file (one locus tag per line) or a TSV "
            "from conserved_annotation_scanner.py / cross_genome_keyword_scanner.py "
            "(the 'Locus_Tag' column is detected and extracted automatically). "
            "Duplicates are silently removed while preserving order. "
            "Cannot be combined with --c1/--c2 or --locus."
        ),
    )

    # ── Output formats ────────────────────────────────────────────────────────
    out = parser.add_argument_group(
        "Output Formats",
        "Specify at least one. Multiple flags can be used in one run.",
    )
    out.add_argument(
        "--gbk",
        type=Path,
        default=None,
        metavar="FILE.gbk",
        help="Annotated GenBank with coordinates adjusted to local origin.",
    )
    out.add_argument(
        "--faa",
        type=Path,
        default=None,
        metavar="FILE.faa",
        help=(
            "Protein FASTA. One entry per non-truncated CDS. "
            "Headers show global genome coordinates for traceability."
        ),
    )
    out.add_argument(
        "--fna",
        type=Path,
        default=None,
        metavar="FILE.fna",
        help=(
            "Region DNA FASTA. One entry for the entire region including "
            "intergenic spacers."
        ),
    )
    out.add_argument(
        "--gene-fna",
        type=Path,
        default=None,
        metavar="FILE.fna",
        help=(
            "Gene DNA FASTA. One entry per non-truncated CDS. "
            "Minus-strand genes are auto reverse-complemented."
        ),
    )

    # ── Utility ───────────────────────────────────────────────────────────────
    parser.add_argument(
        "--list-sequences",
        action="store_true",
        help="List all sequence IDs, lengths, and organisms in the file, then exit.",
    )

    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    """Validates argument combinations and exits with a clear error message."""
    if args.list_sequences:
        return

    errors: list[str] = []

    if not any([args.gbk, args.faa, args.fna, args.gene_fna]):
        errors.append(
            "  No output format specified. Use at least one of:\n"
            "    --gbk FILE.gbk       annotated GenBank\n"
            "    --faa FILE.faa       protein FASTA\n"
            "    --fna FILE.fna       region DNA FASTA\n"
            "    --gene-fna FILE.fna  individual gene DNA FASTA\n"
        )

    if (args.c1 is None) != (args.c2 is None):
        errors.append("  --c1 and --c2 must always be used together.")

    if args.c1 is not None and args.c2 is not None:
        if not args.circular and args.c1 >= args.c2:
            errors.append(
                f"  --c1 ({args.c1:,}) must be less than --c2 ({args.c2:,}).\n"
                f"  For cross-origin extraction on circular genomes, add --circular."
            )
        if args.circular and args.c1 == args.c2:
            errors.append(f"  --c1 and --c2 must differ even in --circular mode.")

    # Count mutually exclusive extraction modes
    has_genes = (args.genes is not None) or (args.genes_file is not None)
    extraction_modes = sum(
        [
            args.locus is not None,
            args.c1 is not None,
            has_genes,
        ]
    )
    if extraction_modes > 1:
        errors.append(
            "  --locus, --c1/--c2, and --genes/--genes-file cannot be combined.\n"
            "  Choose one extraction mode."
        )
    if args.genes and args.genes_file:
        errors.append("  --genes and --genes-file cannot be used together.")

    if has_genes and args.fna:
        errors.append(
            "  --fna (region DNA) is not supported in --genes/--genes-file mode.\n"
            "  Use --gene-fna (individual gene DNA) instead."
        )

    if errors:
        sys.exit("\n[!] Argument errors:\n" + "\n".join(errors) + "\n")


# ── Gene tag loading ──────────────────────────────────────────────────────────


def _load_gene_tags_from_file(path: Path) -> list[str]:
    """Load locus tags from a plain text or TSV file.

    Auto-detects format:
    - If the first line contains tabs → TSV mode: extracts the ``Locus_Tag``
      column. Compatible with output from ``conserved_annotation_scanner.py``
      and ``cross_genome_keyword_scanner.py``.
    - Otherwise → plain text mode: one locus tag per line (# comments ignored).

    Duplicates are removed while preserving the order of first occurrence.

    Args:
        path: Path to the input file.

    Returns:
        Deduplicated list of locus tag strings.
    """
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    if not lines:
        sys.exit(f"[!] --genes-file '{path.name}' is empty.")

    tags: list[str] = []

    if "\t" in lines[0]:
        # TSV mode: find the Locus_Tag column
        header = lines[0].split("\t")
        try:
            col_idx = header.index("Locus_Tag")
        except ValueError:
            sys.exit(
                f"[!] TSV file '{path.name}' has no 'Locus_Tag' column.\n"
                f"    Available columns: {', '.join(header)}"
            )
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) > col_idx:
                tag = parts[col_idx].strip()
                if tag:
                    tags.append(tag)
    else:
        # Plain text mode
        for line in lines:
            tag = line.strip()
            if tag and not tag.startswith("#"):
                tags.append(tag)

    # Deduplicate preserving order
    return list(dict.fromkeys(tags))


# ── Target record resolution ──────────────────────────────────────────────────


def _find_target_record(args: argparse.Namespace) -> SeqRecord:
    """Finds and returns the target SeqRecord from the input file."""
    target = None
    found_second = False

    for record in SeqIO.parse(args.input, "genbank"):
        if args.seq is None or record.id == args.seq:
            if target is None:
                target = record
                if args.seq is not None:
                    break
            else:
                found_second = True
                break
        else:
            found_second = True

    if target is None:
        if args.seq:
            sys.exit(
                f"\n[!] Sequence '{args.seq}' not found in '{args.input.name}'.\n"
                f"    Run --list-sequences to see available IDs.\n"
            )
        else:
            sys.exit(f"\n[!] No records found in '{args.input.name}'.\n")

    if args.seq is None and found_second:
        sys.exit(
            f"\n[!] '{args.input.name}' contains multiple sequences.\n"
            f"    Use --seq SEQ_ID to specify which contig to extract.\n"
            f"    Run --list-sequences to see available IDs.\n"
        )

    return target


# ── Coordinate resolution ─────────────────────────────────────────────────────


def _resolve_locus_range(
    feature_index: dict,
    first_tag: str,
    last_tag: str,
    contig_len: int,
    is_circular: bool = False,
) -> tuple[int, int]:
    """Look up two locus tags and return their genomic coordinate span.

    For linear genomes (``is_circular=False``): returns the minimum start
    and maximum end of the two locus tags.

    For circular genomes (``is_circular=True``): computes both the direct
    path and the cross-origin path, and returns the shorter one. If the
    cross-origin path is shorter, returns ``(max_end, min_start)`` where
    ``c1 > c2`` — the caller must handle this by concatenating two slices.

    Args:
        feature_index: ``{locus_tag: feature}`` index from ``_build_feature_index``.
        first_tag:     Locus tag of the first boundary gene.
        last_tag:      Locus tag of the second boundary gene.
        contig_len:    Total length of the contig in bp.
        is_circular:   Whether to consider the cross-origin path.

    Returns:
        A tuple ``(c1, c2)`` of 1-based genomic coordinates.
        If ``c1 > c2``, the region wraps across the replication origin.

    Raises:
        SystemExit: If either locus tag is not found in the index.
    """
    missing = [t for t in {first_tag, last_tag} if t not in feature_index]
    if missing:
        sys.exit(
            f"\n[!] Locus tag(s) not found: {', '.join(missing)}\n"
            f"    Use find_gbk_features.py to browse available tags.\n"
        )

    def _coords(tag: str) -> tuple[int, int]:
        f = feature_index[tag]
        return int(f.location.start) + 1, int(f.location.end)

    start_a, end_a = _coords(first_tag)
    start_b, end_b = _coords(last_tag)

    min_start = min(start_a, start_b)
    max_end = max(end_a, end_b)

    if not is_circular:
        return min_start, max_end

    # Circular: compare direct path vs cross-origin path
    direct_span = max_end - min_start
    cross_origin_span = contig_len - max_end + min_start

    if cross_origin_span < direct_span:
        # Cross-origin is shorter: return (max_end, min_start) → c1 > c2
        return max_end, min_start
    else:
        return min_start, max_end


def _resolve_coordinates(
    record: SeqRecord,
    feature_index: dict,
    args: argparse.Namespace,
) -> tuple[int, int]:
    """Returns the (c1, c2) extraction coordinates for the active mode."""
    if args.locus:
        return _resolve_locus_range(
            feature_index,
            args.locus[0],
            args.locus[1],
            len(record.seq),
            args.circular,
        )
    if args.c1 is not None:
        return args.c1, args.c2
    # Whole-sequence mode
    return 1, len(record.seq)


def _resolve_gene_tags(
    feature_index: dict,
    gene_tags: list[str],
) -> list[tuple[int, int]]:
    """Look up multiple locus tags and return their coordinates in order.

    Uses the pre-built feature index for O(1) lookup per tag.

    Args:
        feature_index: ``{locus_tag: feature}`` from ``_build_feature_index``.
        gene_tags:     Locus tags to look up (in the order to process them).

    Returns:
        List of ``(c1, c2)`` 1-based coordinate tuples, one per tag.

    Raises:
        SystemExit: If any locus tag is absent from the index.
    """
    missing = [t for t in gene_tags if t not in feature_index]
    if missing:
        sys.exit(f"\n[!] Locus tag(s) not found: {', '.join(missing)}\n")
    return [
        (
            int(feature_index[t].location.start) + 1,
            int(feature_index[t].location.end),
        )
        for t in gene_tags
    ]


# ── FASTA header builders ─────────────────────────────────────────────────────


def _location_string(start_1: int, end_1: int, strand: int) -> str:
    """Returns a 1-based location string matching GenBank convention.

    Args:
        start_1: 1-based start coordinate.
        end_1:   1-based end coordinate.
        strand:  1 for forward strand, -1 for reverse.

    Returns:
        A string such as ``'53508..54911'`` or ``'complement(77319..78044)'``.
    """
    coords = f"{start_1}..{end_1}"
    return f"complement({coords})" if strand == -1 else coords


def _protein_header(
    feature,
    source_seq_id: str,
    organism: str,
    coord_index: dict | None = None,
) -> str:
    """Builds an NCBI-style protein FASTA header with global coordinates.

    When ``coord_index`` is provided, FASTA headers show the original
    genomic coordinates from the source file rather than the local
    coordinates of the extracted region, ensuring traceability back to
    the full genome.

    Args:
        feature:       A BioPython SeqFeature.
        source_seq_id: ID of the parent sequence record.
        organism:      Organism name from record annotations.
        coord_index:   ``{locus_tag: (start_1, end_1, strand)}`` from
                       ``_build_coord_index()``. If provided, global
                       coordinates are used; otherwise local coordinates
                       are used as a fallback.

    Returns:
        A complete FASTA header string starting with ``>``.
    """
    protein_id = feature.qualifiers.get("protein_id", [""])[0]
    locus_tag = feature.qualifiers.get("locus_tag", ["?"])[0]
    product = feature.qualifiers.get("product", ["unknown product"])[0]
    seq_id = protein_id if protein_id else locus_tag

    if coord_index and locus_tag in coord_index:
        g_start, g_end, g_strand = coord_index[locus_tag]
        loc_str = _location_string(g_start, g_end, g_strand)
        loc_label = "global_location"
    else:
        loc_str = _location_string(
            int(feature.location.start) + 1,
            int(feature.location.end),
            feature.location.strand,
        )
        loc_label = "location"

    header = f">{seq_id} {product}"
    if organism:
        header += f" [{organism}]"
    header += f" [locus_tag={locus_tag}] [{loc_label}={loc_str}]"
    if protein_id:
        header += f" [source={source_seq_id}]"
    return header


def _gene_dna_header(
    feature,
    source_seq_id: str,
    organism: str,
    coord_index: dict | None = None,
) -> str:
    """Builds a gene DNA FASTA header with global coordinates.

    Args:
        feature:       A BioPython SeqFeature.
        source_seq_id: ID of the parent sequence record.
        organism:      Organism name.
        coord_index:   ``{locus_tag: (start_1, end_1, strand)}`` for global coords.

    Returns:
        A complete FASTA header string starting with ``>``.
    """
    locus_tag = feature.qualifiers.get("locus_tag", ["?"])[0]
    product = feature.qualifiers.get("product", ["unknown product"])[0]
    protein_id = feature.qualifiers.get("protein_id", [""])[0]

    if coord_index and locus_tag in coord_index:
        g_start, g_end, g_strand = coord_index[locus_tag]
        loc_str = _location_string(g_start, g_end, g_strand)
        loc_label = "global_location"
    else:
        loc_str = _location_string(
            int(feature.location.start) + 1,
            int(feature.location.end),
            feature.location.strand,
        )
        loc_label = "location"

    header = f">{locus_tag} {product}"
    if organism:
        header += f" [{organism}]"
    header += f" [{loc_label}={loc_str}]"
    if protein_id:
        header += f" [protein_id={protein_id}]"
    return header


# ── Output writers ────────────────────────────────────────────────────────────


def write_gbk(region: SeqRecord, output_path: Path) -> int:
    """Writes the BioPython region record as a GenBank file."""
    with open(output_path, "w", encoding="utf-8") as f:
        SeqIO.write(region, f, "genbank")
    return sum(1 for feat in region.features if feat.type == "CDS")


def write_fna(
    region: SeqRecord,
    c1: int,
    c2: int,
    source_seq_id: str,
    organism: str,
    output_path: Path,
    is_circular: bool = False,
) -> None:
    """Writes the entire region as a single wrapped DNA FASTA entry."""
    seq_str = str(region.seq)
    if is_circular and c1 > c2:
        coord_label = f"wrap:{c1}-end+start-{c2}"
    else:
        coord_label = f"{c1}-{c2}"
    header = f">{source_seq_id}:{coord_label} {len(seq_str):,} bp"
    if organism:
        header += f" {organism} genomic region"
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(f"{header}\n{wrap_fasta(seq_str)}\n")


def write_faa(
    region: SeqRecord,
    source_seq_id: str,
    organism: str,
    output_path: Path,
    coord_index: dict | None = None,
    truncated_loci: set | None = None,
) -> int:
    """Writes all non-truncated CDS protein sequences as wrapped FASTA.

    Truncated features (detected by translation-vs-coding-length mismatch)
    are skipped with a warning. Global coordinates from ``coord_index`` are
    embedded in headers for full genomic traceability.
    """
    truncated_loci = truncated_loci or set()
    written = 0
    with open(output_path, "w", encoding="utf-8") as f:
        for feature in region.features:
            if feature.type != "CDS":
                continue
            lt = feature.qualifiers.get("locus_tag", ["?"])[0]
            if lt in truncated_loci:
                print(
                    f"  [!] Skipping truncated CDS at boundary: {lt} "
                    f"(translation/sequence length mismatch — boundary feature).",
                    file=sys.stderr,
                )
                continue
            translation = feature.qualifiers.get("translation", [""])[0]
            if not translation:
                print(
                    f"  [!] Skipping {lt}: no /translation= qualifier.",
                    file=sys.stderr,
                )
                continue
            f.write(
                f"{_protein_header(feature, source_seq_id, organism, coord_index)}\n"
            )
            f.write(f"{wrap_fasta(translation)}\n")
            written += 1
    return written


def write_gene_fna(
    region: SeqRecord,
    source_seq_id: str,
    organism: str,
    output_path: Path,
    coord_index: dict | None = None,
    truncated_loci: set | None = None,
) -> int:
    """Writes each non-truncated CDS coding sequence as a wrapped FASTA entry.

    Truncated features are skipped with a warning (consistent with write_faa).
    Minus-strand genes are auto reverse-complemented so every entry reads
    5' to 3' from the start codon.
    """
    truncated_loci = truncated_loci or set()
    written = 0
    with open(output_path, "w", encoding="utf-8") as f:
        for feature in region.features:
            if feature.type != "CDS":
                continue
            lt = feature.qualifiers.get("locus_tag", ["?"])[0]
            if lt in truncated_loci:
                print(
                    f"  [!] Skipping truncated CDS at boundary: {lt}.",
                    file=sys.stderr,
                )
                continue
            coding_seq: Seq = feature.location.extract(region.seq)
            seq_str = str(coding_seq)
            if not seq_str:
                print(
                    f"  [!] Skipping {lt}: could not extract sequence.",
                    file=sys.stderr,
                )
                continue
            f.write(
                f"{_gene_dna_header(feature, source_seq_id, organism, coord_index)}\n"
            )
            f.write(f"{wrap_fasta(seq_str)}\n")
            written += 1
    return written


# ── Utility ───────────────────────────────────────────────────────────────────


def list_sequences(input_path: Path) -> None:
    """Lists all sequence IDs, lengths, and organisms in a GenBank file."""
    print(f"\n[*] Sequences in '{input_path.name}':\n")
    count = 0
    for record in SeqIO.parse(input_path, "genbank"):
        org = record.annotations.get("organism", "")
        print(f"  {record.id:<35}  {len(record.seq):>12,} bp  {org[:50]}")
        count += 1
    print(f"\n[*] {count} sequence(s). Use one of these IDs with --seq.\n")


# ── Core extraction ───────────────────────────────────────────────────────────


def extract_and_write(args: argparse.Namespace) -> None:
    """Resolves the extraction mode and dispatches to the appropriate handler."""
    target = _find_target_record(args)
    organism = target.annotations.get("organism", "")

    # Build indices ONCE from the original unsliced record — O(N) total.
    # All subsequent locus tag lookups are O(1) from these dicts.
    feature_index = _build_feature_index(target)
    coord_index = _build_coord_index(target)

    # Resolve gene tags from CLI or file
    gene_tags: list[str] | None = None
    if args.genes:
        gene_tags = args.genes
    elif args.genes_file:
        gene_tags = _load_gene_tags_from_file(args.genes_file)
        print(
            f"[*] Loaded {len(gene_tags)} locus tag(s) from '{args.genes_file.name}'.",
            file=sys.stderr,
        )

    if gene_tags is not None:
        extract_genes_mode(
            target, feature_index, coord_index, gene_tags, organism, args
        )
    else:
        extract_region_mode(target, feature_index, coord_index, organism, args)


def extract_region_mode(
    target: SeqRecord,
    feature_index: dict,
    coord_index: dict,
    organism: str,
    args: argparse.Namespace,
) -> None:
    """Extracts a single contiguous region (coordinate, locus, or whole-sequence mode)."""
    c1, c2 = _resolve_coordinates(target, feature_index, args)
    contig_len = len(target.seq)
    is_circular = args.circular and c1 > c2

    if args.locus:
        mode_str = f"locus tag range  ({args.locus[0]} \u2192 {args.locus[1]})"
    elif args.c1 is not None:
        mode_str = "coordinate range"
    else:
        mode_str = "whole sequence"

    if is_circular:
        mode_str += "  [CIRCULAR — cross-origin]"

    # Bounds check (only for non-circular linear extractions)
    if not is_circular:
        if c1 < 1:
            sys.exit(f"\n[!] Resolved start ({c1:,}) must be >= 1.\n")
        if c2 > contig_len:
            sys.exit(
                f"\n[!] Resolved end ({c2:,}) exceeds '{target.id}' "
                f"length ({contig_len:,} bp).\n"
            )

    print(f"[*] Source  : {target.id}  ({contig_len:,} bp)")
    print(f"[*] Organism: {organism or '(not annotated)'}")
    print(f"[*] Mode    : {mode_str}")
    if is_circular:
        print(
            f"[*] Range   : {c1:,} \u2013 end + start \u2013 {c2:,}  (cross-origin wrap)"
        )
    else:
        print(f"[*] Range   : {c1:,} \u2013 {c2:,} bp  (1-based, inclusive)")

    # ── Slice ─────────────────────────────────────────────────────────────────
    if is_circular:
        # Concatenate the two flanking segments around the origin
        region = target[c1 - 1 : contig_len] + target[0:c2]
    else:
        region = target[c1 - 1 : c2]

    cds_count = sum(1 for f in region.features if f.type == "CDS")
    region_len = len(region.seq)

    # Detect truncated boundary features BEFORE writing any output
    truncated_loci = _find_truncated_loci(region)
    if truncated_loci:
        print(
            f"  [!] {len(truncated_loci)} boundary CDS feature(s) detected as truncated "
            f"and will be excluded from FAA and gene-FNA output.",
            file=sys.stderr,
        )

    print(
        f"[*] Extracted {region_len:,} bp  |  "
        f"{cds_count} CDS features  |  "
        f"{len(truncated_loci)} truncated (excluded from FAA/gene-FNA)"
    )
    print()

    # ── Write outputs ─────────────────────────────────────────────────────────
    generated_outputs = []

    if args.gbk:
        region.id = f"{target.id}:{c1}-{c2}"
        region.name = region.id[:16]
        region.description = (
            f"Region {c1:,}-{c2:,} from {target.id}"
            + (f" {organism}" if organism else "")
            + (" [circular wrap]" if is_circular else "")
        )
        n = write_gbk(region, args.gbk)
        generated_outputs.append(args.gbk)
        print(f"  [+] GBK      \u2192 {args.gbk}  ({n} CDS features)")

    if args.fna:
        write_fna(region, c1, c2, target.id, organism, args.fna, is_circular)
        generated_outputs.append(args.fna)
        print(f"  [+] FNA      \u2192 {args.fna}  ({region_len:,} bp, 1 entry)")

    if args.faa:
        n = write_faa(
            region,
            target.id,
            organism,
            args.faa,
            coord_index=coord_index,
            truncated_loci=truncated_loci,
        )
        generated_outputs.append(args.faa)
        print(f"  [+] FAA      \u2192 {args.faa}  ({n} proteins)")

    if args.gene_fna:
        n = write_gene_fna(
            region,
            target.id,
            organism,
            args.gene_fna,
            coord_index=coord_index,
            truncated_loci=truncated_loci,
        )
        generated_outputs.append(args.gene_fna)
        print(f"  [+] gene-FNA \u2192 {args.gene_fna}  ({n} gene sequences)")

    if generated_outputs:
        label = (
            "OUTPUT FILE GENERATED:"
            if len(generated_outputs) == 1
            else "ALL OUTPUT FILES GENERATED:"
        )
        print(f"\n{'=' * 40}\n{label}")
        for path in generated_outputs:
            print(f"  - {Path(path).resolve()}")
        print("=" * 40)


def extract_genes_mode(
    target: SeqRecord,
    feature_index: dict,
    coord_index: dict,
    gene_tags: list[str],
    organism: str,
    args: argparse.Namespace,
) -> None:
    """Extracts multiple individual genes by locus tag using the pre-built index.

    Each gene is sliced individually. FAA and gene-FNA accumulate all genes
    in one file. GBK files are written separately per gene.
    """
    contig_len = len(target.seq)

    print(f"[*] Source  : {target.id}  ({contig_len:,} bp)")
    print(f"[*] Organism: {organism or '(not annotated)'}")
    print(f"[*] Mode    : multiple genes  ({len(gene_tags)} tag(s))")
    print(
        f"[*] Tags    : {', '.join(gene_tags[:10])}"
        + (f"  ...and {len(gene_tags)-10} more" if len(gene_tags) > 10 else "")
    )
    print()

    coords_list = _resolve_gene_tags(feature_index, gene_tags)

    total_bp = 0
    total_genes = 0
    regions = []

    for tag, (c1, c2) in zip(gene_tags, coords_list):
        if c1 < 1 or c2 > contig_len:
            sys.exit(
                f"\n[!] Gene '{tag}': coordinates {c1:,}-{c2:,} are out of bounds.\n"
            )
        region = target[c1 - 1 : c2]
        n_cds = sum(1 for f in region.features if f.type == "CDS")
        total_bp += len(region.seq)
        total_genes += n_cds

        # Detect truncated features for this individual gene slice
        truncated = _find_truncated_loci(region)
        regions.append((tag, c1, c2, region, n_cds, truncated))

    print(
        f"[*] Extracted {len(regions)} gene(s)  |  "
        f"{total_bp:,} total bp  |  {total_genes} total CDS"
    )
    print()

    generated_outputs = []

    if args.gbk:
        for tag, c1, c2, region, n_cds, _ in regions:
            gbk_path = Path(str(args.gbk).replace(".gbk", f"_{tag}.gbk"))
            region.id = f"{target.id}:{c1}-{c2}"
            region.name = region.id[:16]
            region.description = f"Gene {tag} from {target.id}" + (
                f" {organism}" if organism else ""
            )
            write_gbk(region, gbk_path)
            generated_outputs.append(str(gbk_path))
            print(f"  [+] GBK      \u2192 {gbk_path}  ({n_cds} CDS)")

    if args.faa:
        total_written = 0
        with open(args.faa, "w", encoding="utf-8") as f_out:
            for tag, c1, c2, region, _, truncated in regions:
                written = _append_faa_to_file(
                    region, target.id, organism, f_out, coord_index, truncated
                )
                total_written += written
        generated_outputs.append(args.faa)
        print(f"  [+] FAA      \u2192 {args.faa}  ({total_written} proteins)")

    if args.gene_fna:
        total_written = 0
        with open(args.gene_fna, "w", encoding="utf-8") as f_out:
            for tag, c1, c2, region, _, truncated in regions:
                written = _append_gene_fna_to_file(
                    region, target.id, organism, f_out, coord_index, truncated
                )
                total_written += written
        generated_outputs.append(args.gene_fna)
        print(
            f"  [+] gene-FNA \u2192 {args.gene_fna}  ({total_written} gene sequences)"
        )

    if generated_outputs:
        label = (
            "OUTPUT FILE GENERATED:"
            if len(generated_outputs) == 1
            else "ALL OUTPUT FILES GENERATED:"
        )
        print(f"\n{'=' * 40}\n{label}")
        for path in generated_outputs:
            print(f"  - {Path(path).resolve()}")
        print("=" * 40)


def _append_faa_to_file(
    region: SeqRecord,
    source_seq_id: str,
    organism: str,
    f_out,
    coord_index: dict | None = None,
    truncated_loci: set | None = None,
) -> int:
    """Appends non-truncated CDS protein sequences to an open file handle."""
    truncated_loci = truncated_loci or set()
    written = 0
    for feature in region.features:
        if feature.type != "CDS":
            continue
        lt = feature.qualifiers.get("locus_tag", ["?"])[0]
        if lt in truncated_loci:
            print(
                f"  [!] Skipping truncated CDS at boundary: {lt}.",
                file=sys.stderr,
            )
            continue
        translation = feature.qualifiers.get("translation", [""])[0]
        if not translation:
            print(f"  [!] Skipping {lt}: no /translation= qualifier.", file=sys.stderr)
            continue
        f_out.write(
            f"{_protein_header(feature, source_seq_id, organism, coord_index)}\n"
        )
        f_out.write(f"{wrap_fasta(translation)}\n")
        written += 1
    return written


def _append_gene_fna_to_file(
    region: SeqRecord,
    source_seq_id: str,
    organism: str,
    f_out,
    coord_index: dict | None = None,
    truncated_loci: set | None = None,
) -> int:
    """Appends non-truncated CDS coding sequences to an open file handle."""
    truncated_loci = truncated_loci or set()
    written = 0
    for feature in region.features:
        if feature.type != "CDS":
            continue
        lt = feature.qualifiers.get("locus_tag", ["?"])[0]
        if lt in truncated_loci:
            print(
                f"  [!] Skipping truncated CDS at boundary: {lt}.",
                file=sys.stderr,
            )
            continue
        coding_seq: Seq = feature.location.extract(region.seq)
        seq_str = str(coding_seq)
        if not seq_str:
            print(f"  [!] Skipping {lt}: could not extract sequence.", file=sys.stderr)
            continue
        f_out.write(
            f"{_gene_dna_header(feature, source_seq_id, organism, coord_index)}\n"
        )
        f_out.write(f"{wrap_fasta(seq_str)}\n")
        written += 1
    return written


# ── Entry point ───────────────────────────────────────────────────────────────


def main() -> None:
    """Parses arguments and runs the extraction pipeline."""
    args = get_args()
    validate_args(args)

    if not args.input.exists():
        sys.exit(f"\n[!] File not found: {args.input}\n")

    if args.list_sequences:
        list_sequences(args.input)
        return

    extract_and_write(args)


if __name__ == "__main__":
    main()
