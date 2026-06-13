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
        genomic span between them. Order does not matter — the leftmost
        start and rightmost end are always used. To extract a single gene,
        provide the same tag twice: ``--locus TAG TAG``.

    Whole-sequence mode (``--seq`` with no ``--c1``/``--c2`` or ``--locus``):
        Extracts the entire contig. No coordinate lookup is needed. Useful
        for generating a full protein or gene inventory of one assembly
        sequence.

Output Format Reference:
    ``--faa FILE.faa``  — Protein FASTA
        One entry per CDS containing the translated amino acid sequence
        from the GenBank ``/translation=`` qualifier. Intergenic regions
        are excluded. Use for protein-function tools: InterPro (domain and
        family classification), BLASTp (protein homology), HMMER.

    ``--fna FILE.fna``  — Region DNA FASTA
        A single FASTA entry containing the full extracted region as one
        continuous nucleotide sequence. Unlike ``--gene-fna``, intergenic
        spacers between genes are included. Use for whole-region analysis:
        MEME Suite motif discovery (including cross-gene motifs), BLASTn
        against other genome regions.

    ``--gene-fna FILE.fna``  — Gene DNA FASTA
        One entry per CDS containing only the coding nucleotide sequence
        of each individual gene; intergenic regions are excluded. Minus-
        strand genes are automatically reverse-complemented so every entry
        reads 5' to 3' from the start codon. Use for per-gene DNA analysis:
        BLASTn gene homology, codon usage, per-gene GC content.

    ``--gbk FILE.gbk``  — Annotated GenBank
        The extracted region as a GenBank file with all feature coordinates
        adjusted to the new local origin. Opens in Geneious, SnapGene,
        Artemis, and UGENE.

Note:
    This script is part of ongoing research and is associated with an upcoming
    publication. Correct attribution is requested when used in derivative works.
    Released under the MIT License. See LICENSE in the repository root.

Example:
    Check what sequences are in the file::

        python3 extract_genome_region.py -i genome.gbff --list-sequences

    Extract by coordinates::

        python3 extract_genome_region.py -i genome.gbff \\
            --seq NZ_CP134351.1 --c1 53317 --c2 78823 \\
            --gbk region.gbk --faa region.faa --fna region.fna

    Extract by locus tag range (avoids manual coordinate lookup)::

        python3 extract_genome_region.py -i genome.gbff \\
            --seq NZ_CP134351.1 \\
            --locus RHP56_RS00340 RHP56_RS00455 \\
            --faa region.faa --fna region.fna --gene-fna region_genes.fna

    Extract a single gene by locus tag::

        python3 extract_genome_region.py -i genome.gbff \\
            --seq NZ_CP134351.1 \\
            --locus RHP56_RS00345 RHP56_RS00345 \\
            --faa gene.faa --gene-fna gene.fna

    Extract a whole contig (no coordinates needed)::

        python3 extract_genome_region.py -i genome.gbff \\
            --seq NZ_CP134351.1 \\
            --faa all_proteins.faa --fna full_contig.fna
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "2.1.0"

import sys
import argparse
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )
from utils import wrap_fasta

# ── CLI ───────────────────────────────────────────────────────────────────────


def get_args() -> argparse.Namespace:
    """Configures the CLI parser and returns parsed arguments.

    Returns:
        An ``argparse.Namespace`` object containing all parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Extract a genomic region from a GBFF or GBK file in one or more "
            "output formats. Three extraction modes: coordinate (--c1/--c2), "
            "locus tag (--locus FIRST LAST), or whole-sequence (--seq alone). "
            "Works with NCBI GBFF, antiSMASH GBK, and Prokka GBK files."
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

    # ── Extraction modes ──────────────────────────────────────────────────────
    mode = parser.add_argument_group(
        "Extraction Modes",
        "Choose one mode. If none is given, the entire contig is extracted "
        "(requires --seq for multi-record files).",
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
        help="Coordinate-mode end position (1-based, inclusive). Requires --c1.",
    )
    mode.add_argument(
        "--locus",
        type=str,
        nargs=2,
        metavar=("FIRST_TAG", "LAST_TAG"),
        default=None,
        help=(
            "Locus tag mode: extract the genomic span from FIRST_TAG to LAST_TAG. "
            "The script looks up both locus tags and computes the coordinates "
            "automatically, so no manual coordinate lookup is needed. "
            "Tag order does not matter — the leftmost start and rightmost end "
            "are always used. "
            "To extract a single gene: --locus TAG TAG (same tag twice). "
            "Cannot be combined with --c1/--c2 or --genes."
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
            "Provide one or more locus tags (e.g. --genes ctg1_47 ctg1_58 ctg1_74). "
            "Each gene is extracted separately and written as individual FASTA entries "
            "to the output file(s). Order is preserved as given on the command line. "
            "Cannot be combined with --c1/--c2 or --locus. Note: --fna is not "
            "supported in --genes mode (use --gene-fna instead)."
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
        help=(
            "Annotated GenBank with all feature coordinates adjusted to the "
            "new local origin. Opens in Geneious, SnapGene, Artemis, UGENE."
        ),
    )
    out.add_argument(
        "--faa",
        type=Path,
        default=None,
        metavar="FILE.faa",
        help=(
            "Protein FASTA (amino acid sequences). One entry per CDS, "
            "using the /translation= qualifier. Intergenic regions excluded. "
            "Use for: InterPro, BLASTp, HMMER."
        ),
    )
    out.add_argument(
        "--fna",
        type=Path,
        default=None,
        metavar="FILE.fna",
        help=(
            "Region DNA FASTA (nucleotide). One entry for the entire extracted "
            "region including intergenic spacers. "
            "Use for: MEME Suite whole-region motif discovery, BLASTn."
        ),
    )
    out.add_argument(
        "--gene-fna",
        type=Path,
        default=None,
        metavar="FILE.fna",
        help=(
            "Gene DNA FASTA (nucleotide). One entry per CDS containing only "
            "the coding sequence; intergenic regions excluded. Minus-strand "
            "genes are auto reverse-complemented (reads 5' to 3'). "
            "Use for: BLASTn gene homology, codon usage, per-gene GC content."
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
    """Validates argument combinations and exits with a clear error message.

    Args:
        args: Parsed argument namespace from ``get_args()``.

    Raises:
        SystemExit: If a required argument combination is missing or invalid.
    """
    if args.list_sequences:
        return

    errors: list[str] = []

    # At least one output format is always required
    if not any([args.gbk, args.faa, args.fna, args.gene_fna]):
        errors.append(
            "  No output format specified. Use at least one of:\n"
            "    --gbk FILE.gbk       annotated GenBank\n"
            "    --faa FILE.faa       protein FASTA   (InterPro / BLASTp)\n"
            "    --fna FILE.fna       region DNA FASTA (MEME / BLASTn)\n"
            "    --gene-fna FILE.fna  individual gene DNA FASTA\n"
        )

    # --c1 and --c2 must appear together
    if (args.c1 is None) != (args.c2 is None):
        errors.append("  --c1 and --c2 must always be used together.")

    if args.c1 is not None and args.c2 is not None and args.c1 >= args.c2:
        errors.append(f"  --c1 ({args.c1:,}) must be less than --c2 ({args.c2:,}).")

    # --locus, --c1/--c2, and --genes are mutually exclusive
    extraction_modes = sum(
        [
            args.locus is not None,
            args.c1 is not None,
            args.genes is not None,
        ]
    )
    if extraction_modes > 1:
        errors.append(
            "  --locus, --c1/--c2, and --genes cannot be combined.\n"
            "  Choose one extraction mode."
        )

    # --fna is not supported in --genes mode
    if args.genes and args.fna:
        errors.append(
            "  --fna (region DNA) is not supported in --genes mode.\n"
            "  Use --gene-fna (individual gene DNA) instead."
        )

    if errors:
        sys.exit("\n[!] Argument errors:\n" + "\n".join(errors) + "\n")


# ── Target record resolution ──────────────────────────────────────────────────


def _find_target_record(args):
    """Finds and returns the target SeqRecord from the input file.

    Handles single-record files (no ``--seq`` needed) and multi-record files
    (``--seq`` required). For whole-sequence mode without ``--seq``, exits
    with a helpful error if the file has more than one sequence.

    Args:
        args: Parsed argument namespace.

    Returns:
        The matching BioPython ``SeqRecord``.

    Raises:
        SystemExit: If the sequence is not found or ``--seq`` is required
            but not provided for a multi-record file.
    """
    target = None
    found_second = False

    for record in SeqIO.parse(args.input, "genbank"):
        if args.seq is None or record.id == args.seq:
            if target is None:
                target = record
                if args.seq is not None:
                    break  # Found the named sequence; no need to read further
            else:
                found_second = True
                break
        else:
            found_second = True  # There is at least one record that didn't match

    if target is None:
        if args.seq:
            sys.exit(
                f"\n[!] Sequence '{args.seq}' not found in '{args.input.name}'.\n"
                f"    Run --list-sequences to see available IDs.\n"
            )
        else:
            sys.exit(f"\n[!] No records found in '{args.input.name}'.\n")

    # If no --seq was given and there are multiple records, we cannot safely
    # guess which contig the user wants.
    if args.seq is None and found_second:
        sys.exit(
            f"\n[!] '{args.input.name}' contains multiple sequences.\n"
            f"    Use --seq SEQ_ID to specify which contig to extract.\n"
            f"    Run --list-sequences to see available IDs.\n"
        )

    return target


# ── Coordinate resolution ─────────────────────────────────────────────────────


def _resolve_locus_range(record, first_tag: str, last_tag: str) -> tuple[int, int]:
    """Looks up two locus tags and returns their genomic coordinate span.

    Searches all CDS features in ``record`` for ``first_tag`` and
    ``last_tag``. Returns the 1-based coordinates that span both genes,
    taking the minimum start and maximum end so that the order of the
    supplied tags does not matter.

    If the same tag is supplied twice, returns the start and end of that
    single gene (single-gene extraction mode).

    Args:
        record:    A BioPython ``SeqRecord`` to search.
        first_tag: Locus tag of the first boundary gene.
        last_tag:  Locus tag of the second boundary gene.

    Returns:
        A tuple ``(c1, c2)`` of 1-based genomic coordinates.

    Raises:
        SystemExit: If either locus tag is not found in the record.
    """
    coords: dict[str, tuple[int, int]] = {}

    for feature in record.features:
        if feature.type != "CDS":
            continue
        lt = feature.qualifiers.get("locus_tag", [""])[0]
        if lt in (first_tag, last_tag):
            coords[lt] = (
                int(feature.location.start) + 1,  # 1-based start
                int(feature.location.end),  # 1-based end
            )

    missing = [t for t in {first_tag, last_tag} if t not in coords]
    if missing:
        sys.exit(
            f"\n[!] Locus tag(s) not found in '{record.id}': "
            f"{', '.join(missing)}\n"
            f"    Use find_gbk_features.py -i {record.id} --seq {record.id} "
            f"to browse available tags.\n"
        )

    all_starts = [coords[first_tag][0], coords[last_tag][0]]
    all_ends = [coords[first_tag][1], coords[last_tag][1]]

    return min(all_starts), max(all_ends)


def _resolve_coordinates(record, args) -> tuple[int, int]:
    """Returns the (c1, c2) extraction coordinates for the active mode.

    Dispatches to the appropriate coordinate source:
    - Locus mode:       looks up ``--locus`` tags in the record.
    - Coordinate mode:  returns ``--c1`` and ``--c2`` directly.
    - Whole-sequence:   returns ``(1, len(record.seq))``.

    Args:
        record: The target BioPython ``SeqRecord``.
        args:   Parsed argument namespace.

    Returns:
        A tuple ``(c1, c2)`` of 1-based genomic coordinates suitable for
        passing to BioPython record slicing as ``record[c1-1 : c2]``.
    """
    if args.locus:
        return _resolve_locus_range(record, args.locus[0], args.locus[1])

    if args.c1 is not None:
        return args.c1, args.c2

    # Whole-sequence mode
    return 1, len(record.seq)


def _resolve_gene_tags(record, gene_tags: list[str]) -> list[tuple[int, int]]:
    """Looks up multiple locus tags and returns their coordinates in order.

    Searches all CDS features in ``record`` for each locus tag in
    ``gene_tags``. Returns a list of (c1, c2) tuples preserving the
    order given by the user (not sorted by genomic position).

    Args:
        record:    A BioPython ``SeqRecord`` to search.
        gene_tags: List of locus tag strings to find.

    Returns:
        A list of tuples ``(c1, c2)`` of 1-based genomic coordinates,
        one per gene, in the order requested.

    Raises:
        SystemExit: If any locus tag is not found in the record.
    """
    coords: dict[str, tuple[int, int]] = {}

    for feature in record.features:
        if feature.type != "CDS":
            continue
        lt = feature.qualifiers.get("locus_tag", [""])[0]
        if lt in gene_tags:
            coords[lt] = (
                int(feature.location.start) + 1,  # 1-based start
                int(feature.location.end),  # 1-based end
            )

    missing = [t for t in gene_tags if t not in coords]
    if missing:
        sys.exit(
            f"\n[!] Locus tag(s) not found in '{record.id}': " f"{', '.join(missing)}\n"
        )

    # Return in the order requested by the user
    return [coords[t] for t in gene_tags]


# ── FASTA header builders ─────────────────────────────────────────────────────


def _location_string(start_0: int, end_0: int, strand: int) -> str:
    """Returns a 1-based location string matching GenBank convention.

    BioPython stores positions as 0-indexed; this function adds 1 to the
    start so displayed coordinates match what users see in GBFF files.

    Args:
        start_0: 0-indexed feature start from BioPython.
        end_0:   0-indexed feature end from BioPython (already 1-based
                 equivalent for the end position).
        strand:  1 for forward strand, -1 for reverse.

    Returns:
        A string such as ``'53508..54911'`` or ``'complement(77319..78044)'``.
    """
    coords = f"{start_0 + 1}..{end_0}"
    return f"complement({coords})" if strand == -1 else coords


def _protein_header(feature, source_seq_id: str, organism: str) -> str:
    """Builds an NCBI-style protein FASTA header.

    Format::

        >PROTEIN_ID PRODUCT [ORGANISM] [locus_tag=X] [location=LOCATION]

    Falls back to the locus tag as the sequence ID when no ``/protein_id=``
    is present (e.g. in Prokka-annotated files).

    Args:
        feature:       A BioPython ``SeqFeature`` object.
        source_seq_id: ID of the parent sequence record.
        organism:      Organism name from record annotations.

    Returns:
        A complete FASTA header string starting with ``>``.
    """
    protein_id = feature.qualifiers.get("protein_id", [""])[0]
    locus_tag = feature.qualifiers.get("locus_tag", ["?"])[0]
    product = feature.qualifiers.get("product", ["unknown product"])[0]
    seq_id = protein_id if protein_id else locus_tag
    loc_str = _location_string(
        int(feature.location.start),
        int(feature.location.end),
        feature.location.strand,
    )

    header = f">{seq_id} {product}"
    if organism:
        header += f" [{organism}]"
    header += f" [locus_tag={locus_tag}] [location={loc_str}]"
    if protein_id:
        header += f" [source={source_seq_id}]"
    return header


def _gene_dna_header(feature, source_seq_id: str, organism: str) -> str:
    """Builds a gene DNA FASTA header.

    Format::

        >LOCUS_TAG PRODUCT [ORGANISM] [location=LOCATION] [protein_id=X]

    Args:
        feature:       A BioPython ``SeqFeature`` object.
        source_seq_id: ID of the parent sequence record.
        organism:      Organism name from record annotations.

    Returns:
        A complete FASTA header string starting with ``>``.
    """
    locus_tag = feature.qualifiers.get("locus_tag", ["?"])[0]
    product = feature.qualifiers.get("product", ["unknown product"])[0]
    protein_id = feature.qualifiers.get("protein_id", [""])[0]
    loc_str = _location_string(
        int(feature.location.start),
        int(feature.location.end),
        feature.location.strand,
    )

    header = f">{locus_tag} {product}"
    if organism:
        header += f" [{organism}]"
    header += f" [location={loc_str}]"
    if protein_id:
        header += f" [protein_id={protein_id}]"
    return header


# ── Output writers ────────────────────────────────────────────────────────────


def write_gbk(region, output_path: Path) -> int:
    """Writes the BioPython region record as a GenBank file.

    Args:
        region:      A sliced BioPython ``SeqRecord`` with adjusted coordinates.
        output_path: Destination file path.

    Returns:
        The number of CDS features written.
    """
    with open(output_path, "w", encoding="utf-8") as f:
        SeqIO.write(region, f, "genbank")
    return sum(1 for feat in region.features if feat.type == "CDS")


def write_fna(
    region,
    c1: int,
    c2: int,
    source_seq_id: str,
    organism: str,
    output_path: Path,
) -> None:
    """Writes the entire region as a single wrapped DNA FASTA entry.

    The header encodes the original genomic coordinates so the file
    remains traceable without the accompanying GenBank. Intergenic regions
    between genes are included in the output.

    Args:
        region:        Sliced BioPython ``SeqRecord``.
        c1:            Original 1-based start coordinate.
        c2:            Original 1-based end coordinate.
        source_seq_id: ID of the source sequence.
        organism:      Organism name.
        output_path:   Destination file path.
    """
    seq_str = str(region.seq)
    header = f">{source_seq_id}:{c1}-{c2} {len(seq_str):,} bp"
    if organism:
        header += f" {organism} genomic region"
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(f"{header}\n{wrap_fasta(seq_str)}\n")


def write_faa(
    region,
    source_seq_id: str,
    organism: str,
    output_path: Path,
) -> int:
    """Writes all CDS protein sequences from the region as wrapped FASTA.

    Uses the ``/translation=`` qualifier stored in each CDS feature.
    Features without a translation qualifier are skipped with a warning.
    Intergenic regions are not included.

    Args:
        region:        Sliced BioPython ``SeqRecord``.
        source_seq_id: ID of the source sequence.
        organism:      Organism name.
        output_path:   Destination file path.

    Returns:
        The number of protein entries written.
    """
    written = 0
    with open(output_path, "w", encoding="utf-8") as f:
        for feature in region.features:
            if feature.type != "CDS":
                continue
            translation = feature.qualifiers.get("translation", [""])[0]
            if not translation:
                lt = feature.qualifiers.get("locus_tag", ["?"])[0]
                print(
                    f"  [!] Skipping {lt}: no /translation= qualifier.", file=sys.stderr
                )
                continue
            f.write(f"{_protein_header(feature, source_seq_id, organism)}\n")
            f.write(f"{wrap_fasta(translation)}\n")
            written += 1
    return written


def write_gene_fna(
    region,
    source_seq_id: str,
    organism: str,
    output_path: Path,
) -> int:
    """Writes each CDS coding sequence as a separate wrapped FASTA entry.

    Extracts sequences from the region DNA using
    ``feature.location.extract()``, which correctly handles complement-
    strand genes (auto reverse-complemented) and multi-exon join locations.
    Each sequence reads 5' to 3' from the start codon. Intergenic regions
    are excluded.

    Args:
        region:        Sliced BioPython ``SeqRecord``.
        source_seq_id: ID of the source sequence.
        organism:      Organism name.
        output_path:   Destination file path.

    Returns:
        The number of gene sequences written.
    """
    written = 0
    with open(output_path, "w", encoding="utf-8") as f:
        for feature in region.features:
            if feature.type != "CDS":
                continue
            coding_seq: Seq = feature.location.extract(region.seq)
            seq_str = str(coding_seq)
            if not seq_str:
                lt = feature.qualifiers.get("locus_tag", ["?"])[0]
                print(
                    f"  [!] Skipping {lt}: could not extract sequence.", file=sys.stderr
                )
                continue
            f.write(f"{_gene_dna_header(feature, source_seq_id, organism)}\n")
            f.write(f"{wrap_fasta(seq_str)}\n")
            written += 1
    return written


# ── Utility ───────────────────────────────────────────────────────────────────


def list_sequences(input_path: Path) -> None:
    """Lists all sequence IDs, lengths, and organisms in a GenBank file.

    Args:
        input_path: Path to the GenBank file.
    """
    print(f"\n[*] Sequences in '{input_path.name}':\n")
    count = 0
    for record in SeqIO.parse(input_path, "genbank"):
        org = record.annotations.get("organism", "")
        print(f"  {record.id:<35}  {len(record.seq):>12,} bp  {org[:50]}")
        count += 1
    print(f"\n[*] {count} sequence(s). Use one of these IDs with --seq.\n")


# ── Core extraction ───────────────────────────────────────────────────────────


def extract_and_write(args: argparse.Namespace) -> None:
    """Resolves the extraction mode and dispatches to the appropriate handler.

    Args:
        args: Parsed argument namespace from ``get_args()``.
    """
    target = _find_target_record(args)
    organism = target.annotations.get("organism", "")

    if args.genes:
        extract_genes_mode(target, args.genes, organism, args)
    else:
        extract_region_mode(target, organism, args)


def extract_region_mode(target, organism: str, args: argparse.Namespace) -> None:
    """Extracts a single contiguous region (coordinate, locus, or whole-sequence mode).

    Args:
        target:   The target BioPython ``SeqRecord``.
        organism: Organism name from record annotations.
        args:     Parsed argument namespace.
    """
    c1, c2 = _resolve_coordinates(target, args)
    contig_len = len(target.seq)

    # Determine which mode is active for the status display
    if args.locus:
        mode_str = f"locus tag range  ({args.locus[0]} \u2192 {args.locus[1]})"
    elif args.c1 is not None:
        mode_str = "coordinate range"
    else:
        mode_str = "whole sequence"

    # Coordinate bounds check
    if c1 < 1:
        sys.exit(f"\n[!] Resolved start ({c1:,}) must be \u2265 1.\n")
    if c2 > contig_len:
        sys.exit(
            f"\n[!] Resolved end ({c2:,}) exceeds '{target.id}' "
            f"length ({contig_len:,} bp).\n"
        )

    print(f"[*] Source  : {target.id}  ({contig_len:,} bp)")
    print(f"[*] Organism: {organism or '(not annotated)'}")
    print(f"[*] Mode    : {mode_str}")
    print(f"[*] Range   : {c1:,} \u2013 {c2:,} bp  (1-based, inclusive)")

    # ── Slice ─────────────────────────────────────────────────────────────────
    region = target[c1 - 1 : c2]
    cds_count = sum(1 for f in region.features if f.type == "CDS")
    region_len = len(region.seq)

    print(f"[*] Extracted {region_len:,} bp  |  {cds_count} CDS features")
    print()

    # ── Write outputs ─────────────────────────────────────────────────────────
    generated_outputs = []
    if args.gbk:
        region.id = f"{target.id}:{c1}-{c2}"
        region.name = region.id[:16]
        region.description = f"Region {c1:,}-{c2:,} from {target.id}" + (
            f" {organism}" if organism else ""
        )
        n = write_gbk(region, args.gbk)
        generated_outputs.append(args.gbk)
        print(f"  [+] GBK      \u2192 {args.gbk}  ({n} CDS features)")

    if args.fna:
        write_fna(region, c1, c2, target.id, organism, args.fna)
        generated_outputs.append(args.fna)
        print(f"  [+] FNA      \u2192 {args.fna}  ({region_len:,} bp, 1 entry)")

    if args.faa:
        n = write_faa(region, target.id, organism, args.faa)
        generated_outputs.append(args.faa)
        print(f"  [+] FAA      \u2192 {args.faa}  ({n} proteins)")

    if args.gene_fna:
        n = write_gene_fna(region, target.id, organism, args.gene_fna)
        generated_outputs.append(args.gene_fna)
        print(f"  [+] gene-FNA \u2192 {args.gene_fna}  ({n} gene sequences)")

    if generated_outputs:
        header = f"{'OUTPUT FILE' if len(generated_outputs) == 1 else 'ALL OUTPUT FILES'} GENERATED:"
        print(f"\n{'=' * 40}")
        print(header)
        for path in generated_outputs:
            print(f"  - {Path(path).resolve()}")
        print("=" * 40)


def extract_genes_mode(
    target, gene_tags: list[str], organism: str, args: argparse.Namespace
) -> None:
    """Extracts multiple individual genes specified by locus tags.

    Each gene is written separately. FAA and gene-FNA files accumulate all
    genes in one file. GBK files are written separately per gene.

    Args:
        target:     The target BioPython ``SeqRecord``.
        gene_tags:  List of locus tag strings to extract.
        organism:   Organism name from record annotations.
        args:       Parsed argument namespace.
    """
    contig_len = len(target.seq)

    print(f"[*] Source  : {target.id}  ({contig_len:,} bp)")
    print(f"[*] Organism: {organism or '(not annotated)'}")
    print(f"[*] Mode    : multiple genes  ({len(gene_tags)} tag(s))")
    print(f"[*] Tags    : {', '.join(gene_tags)}")
    print()

    # Resolve all locus tags to coordinates
    coords_list = _resolve_gene_tags(target, gene_tags)

    total_bp = 0
    total_genes = 0

    # Collect all regions for output
    regions = []
    for i, (c1, c2) in enumerate(coords_list):
        if c1 < 1 or c2 > contig_len:
            sys.exit(
                f"\n[!] Gene {i+1} ({gene_tags[i]}): "
                f"coordinates {c1:,}-{c2:,} are out of bounds.\n"
            )
        region = target[c1 - 1 : c2]
        cds_in_region = sum(1 for f in region.features if f.type == "CDS")
        total_bp += len(region.seq)
        total_genes += cds_in_region
        regions.append((gene_tags[i], c1, c2, region, cds_in_region))

    print(
        f"[*] Extracted {len(regions)} gene(s)  |  {total_bp:,} total bp  |  {total_genes} total CDS"
    )
    print()

    # ── Write outputs ─────────────────────────────────────────────────────────
    generated_outputs = []
    if args.gbk:
        for tag, c1, c2, region, n_cds in regions:
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
            for tag, c1, c2, region, _ in regions:
                written = _append_faa_to_file(region, target.id, organism, f_out)
                total_written += written
        generated_outputs.append(args.faa)
        print(f"  [+] FAA      \u2192 {args.faa}  ({total_written} proteins)")

    if args.gene_fna:
        total_written = 0
        with open(args.gene_fna, "w", encoding="utf-8") as f_out:
            for tag, c1, c2, region, _ in regions:
                written = _append_gene_fna_to_file(region, target.id, organism, f_out)
                total_written += written
        generated_outputs.append(args.gene_fna)
        print(
            f"  [+] gene-FNA \u2192 {args.gene_fna}  ({total_written} gene sequences)"
        )

    if generated_outputs:
        header = f"{'OUTPUT FILE' if len(generated_outputs) == 1 else 'ALL OUTPUT FILES'} GENERATED:"
        print("\n" + "=" * 40)
        print(f"{header}")
        for path in generated_outputs:
            print(f"  - {Path(path).resolve()}")
        print("=" * 40)


def _append_faa_to_file(region, source_seq_id: str, organism: str, f_out) -> int:
    """Appends all CDS protein sequences from region to an open file.

    Args:
        region:        Sliced BioPython ``SeqRecord``.
        source_seq_id: ID of the source sequence.
        organism:      Organism name.
        f_out:         Open file object in write mode.

    Returns:
        The number of protein entries written.
    """
    written = 0
    for feature in region.features:
        if feature.type != "CDS":
            continue
        translation = feature.qualifiers.get("translation", [""])[0]
        if not translation:
            lt = feature.qualifiers.get("locus_tag", ["?"])[0]
            print(f"  [!] Skipping {lt}: no /translation= qualifier.", file=sys.stderr)
            continue
        f_out.write(f"{_protein_header(feature, source_seq_id, organism)}\n")
        f_out.write(f"{wrap_fasta(translation)}\n")
        written += 1
    return written


def _append_gene_fna_to_file(region, source_seq_id: str, organism: str, f_out) -> int:
    """Appends all CDS coding sequences from region to an open file.

    Args:
        region:        Sliced BioPython ``SeqRecord``.
        source_seq_id: ID of the source sequence.
        organism:      Organism name.
        f_out:         Open file object in write mode.

    Returns:
        The number of gene sequences written.
    """
    written = 0
    for feature in region.features:
        if feature.type != "CDS":
            continue
        coding_seq: Seq = feature.location.extract(region.seq)
        seq_str = str(coding_seq)
        if not seq_str:
            lt = feature.qualifiers.get("locus_tag", ["?"])[0]
            print(f"  [!] Skipping {lt}: could not extract sequence.", file=sys.stderr)
            continue
        f_out.write(f"{_gene_dna_header(feature, source_seq_id, organism)}\n")
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
