"""
GBK Ortholog Finder (with custom pretty-print formatting)

Extracts protein-coding sequences directly from a query GBK/GBFF file,
optionally applies mature peptide trimming via calculate_mature_core(),
then compares them against reference genome(s) using pairwise sequence
alignment to identify orthologs by percent identity.

OUTPUT:
    - Terminal: Beautiful aligned columns (no external dependencies)
              [Unless using > output.txt, which includes protein sequences]
    - File (-o): Clean TSV format with protein sequences for detailed analysis

This is different from ortholog_extractor.py:
    ortholog_extractor.py  → takes a pre-made FASTA, uses EXACT substring match
    gbk_ortholog_finder.py → takes a GBK directly, uses SIMILARITY (% identity)

Example Usage:

    # Compare a region GBK against a single reference genome
    python3 gbk_ortholog_finder.py \\
        -q region001.gbk \\
        -r ATCC8293.gbff \\
        -o results.tsv

    # Compare a full genome GBK against a folder of reference genomes
    python3 gbk_ortholog_finder.py \\
        -q prokka_result.gbk \\
        -r references/ \\
        -o results.tsv \\
        --identity 0.40

    # Compare with mature core trimming enabled (for bacteriocins)
    python3 gbk_ortholog_finder.py \\
        -q region001.gbk \\
        -r references/ \\
        -o results.tsv \\
        --mature

    # Only compare small proteins (likely bacteriocins/peptides, <= 150 aa)
    python3 gbk_ortholog_finder.py \\
        -q region001.gbk \\
        -r references/ \\
        --mature \\
        --max-length 150 \\
        -o results.tsv

    # Output to text file via shell redirection (includes protein sequences)
    python3 gbk_ortholog_finder.py \\
        -q region001.gbk \\
        -r references/ \\
        --mature \\
        --max-length 150 > results.txt

License: MIT
Note: This code is part of ongoing research. Associated with upcoming publication.
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.2"

import sys
import argparse
import csv
from pathlib import Path
from dataclasses import dataclass

try:
    from Bio import SeqIO
    from Bio.Align import substitution_matrices, PairwiseAligner
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )

from utils import stream_reference_files, calculate_mature_core, smart_open


# ─────────────────────────────────────────────────────────────────────────────
# DATA STRUCTURES
# ─────────────────────────────────────────────────────────────────────────────
@dataclass
class Protein:
    """Holds a single protein extracted from a GBK file."""

    locus_tag: str
    product: str
    sequence: str  # Full sequence from GBK annotation
    mature_sequence: str  # After calculate_mature_core() (may equal sequence)
    source_file: str  # Which GBK file it came from
    length: int  # Length of the full sequence


@dataclass
class OrthoHit:
    """Holds a single ortholog comparison result."""

    query_locus: str
    query_product: str
    query_seq: str  # Sequence used for alignment (full or mature)
    ref_locus: str
    ref_product: str
    ref_seq: str  # Full reference protein sequence
    ref_file: str  # Name of the reference genome file
    identity: float  # Percent identity (0.0 to 1.0)
    alignment_length: int  # Length of the aligned region
    query_length: int  # Length of query sequence
    ref_length: int  # Length of reference sequence


# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: ARGUMENT PARSER
# ─────────────────────────────────────────────────────────────────────────────


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Find orthologs by extracting proteins from a GBK query and comparing to references.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Region vs single reference, save to TSV (includes protein sequences)
  python3 gbk_ortholog_finder.py -q region001.gbk -r ATCC8293.gbff -o hits.tsv

  # Full genome vs folder, print pretty table to terminal
  python3 gbk_ortholog_finder.py -q genome.gbff -r references/ --identity 0.40

  # Only bacteriocin-sized proteins, with mature core trimming
  python3 gbk_ortholog_finder.py -q region001.gbk -r refs/ --mature --max-length 150 -o hits.tsv

  # Redirect to text file (includes protein sequences)
  python3 gbk_ortholog_finder.py -q region001.gbk -r references/ > results.txt
        """,
    )

    parser.add_argument(
        "-q",
        "--query",
        type=Path,
        required=True,
        help="Query GBK/GBFF file to extract proteins from.",
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=Path,
        required=True,
        help="Reference GBK/GBFF file OR directory of reference files to compare against.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output TSV file with protein sequences. If not specified, prints table to terminal (no sequences).",
    )
    parser.add_argument(
        "--identity",
        type=float,
        default=0.30,
        help="Minimum percent identity threshold (0.0–1.0). Default: 0.30 (30%%)",
    )
    parser.add_argument(
        "--mature",
        action="store_true",
        default=False,
        help="Apply calculate_mature_core() to extract mature peptide before comparing.",
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=None,
        help="Only include proteins <= this length (aa). Useful for filtering bacteriocins.",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=10,
        help="Minimum protein length (aa) to include. Default: 10",
    )

    return parser.parse_args()


# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: EXTRACT PROTEINS FROM QUERY GBK
# ─────────────────────────────────────────────────────────────────────────────


def extract_proteins_from_gbk(
    gbk_path: Path,
    apply_mature: bool = False,
    min_length: int = 10,
    max_length: int | None = None,
) -> list[Protein]:
    """
    Reads a GBK/GBFF file and extracts all CDS protein sequences.

    For each CDS, it captures:
        - Locus tag (unique gene ID)
        - Product name (annotation)
        - Full protein sequence (from translation qualifier)
        - Mature sequence (if --mature is set, else same as full)

    Args:
        gbk_path:     Path to the query GBK/GBFF file.
        apply_mature: If True, runs calculate_mature_core() on each protein.
        min_length:   Skip proteins shorter than this.
        max_length:   Skip proteins longer than this (None = no limit).

    Returns:
        List of Protein objects extracted from the file.
    """
    proteins = []

    print(f"\n[*] Extracting proteins from query: {gbk_path.name}", file=sys.stderr)

    for record in SeqIO.parse(gbk_path, "genbank"):
        for feature in record.features:

            # Only process protein-coding sequences
            if feature.type != "CDS":
                continue

            # Get the translation (protein sequence)
            translation = feature.qualifiers.get("translation", [""])[0]
            if not translation:
                continue  # Skip CDS with no translation

            # Get annotation fields
            locus_tag = feature.qualifiers.get("locus_tag", ["UNKNOWN"])[0]
            product = feature.qualifiers.get("product", ["Unknown product"])[0]

            full_length = len(translation)

            # Apply length filters
            if full_length < min_length:
                continue
            if max_length is not None and full_length > max_length:
                continue

            # Apply mature core trimming if requested
            if apply_mature:
                mature_seq = calculate_mature_core(translation)
            else:
                mature_seq = translation  # Use full sequence

            protein = Protein(
                locus_tag=locus_tag,
                product=product,
                sequence=translation,
                mature_sequence=mature_seq,
                source_file=gbk_path.name,
                length=full_length,
            )
            proteins.append(protein)

            # Show what we extracted
            if apply_mature and mature_seq != translation:
                print(
                    f"   {locus_tag} ({full_length} aa) → mature: {len(mature_seq)} aa | {product[:47] + "..."}",
                    file=sys.stderr,
                )
            else:
                print(
                    f"   {locus_tag} ({full_length} aa) | {product[:47] + "..."}",
                    file=sys.stderr,
                )

    print(f"\n[*] Extracted {len(proteins)} protein(s) from query.\n", file=sys.stderr)
    return proteins


# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: PAIRWISE ALIGNMENT & IDENTITY CALCULATION
# ─────────────────────────────────────────────────────────────────────────────


def calculate_identity(seq_a: str, seq_b: str, blosum62) -> tuple[float, int]:
    """
    Performs a global pairwise alignment between two protein sequences
    and returns the percent identity and alignment length.

    Uses BioPython's PairwiseAligner with BLOSUM62 substitution matrix,
    which is the standard for protein sequence comparison.

    Args:
        seq_a:    First protein sequence (query).
        seq_b:    Second protein sequence (reference).
        blosum62: BLOSUM62 substitution matrix.

    Returns:
        A tuple of (percent_identity, alignment_length).
        percent_identity is between 0.0 and 1.0.
        alignment_length is the number of aligned positions.
    """
    # Skip if either sequence is empty
    if not seq_a or not seq_b:
        return 0.0, 0

    # Set up the aligner
    aligner = PairwiseAligner()
    aligner.mode = "global"  # Global = Needleman-Wunsch
    aligner.substitution_matrix = blosum62
    aligner.open_gap_score = -10  # Standard gap open penalty
    aligner.extend_gap_score = -0.5  # Standard gap extension penalty

    # Run alignment and take the best hit
    alignments = aligner.align(seq_a, seq_b)
    best = next(iter(alignments), None)

    if best is None:
        return 0.0, 0

    # Count identical positions in the alignment
    aligned_query = best[0]  # Query sequence with gaps
    aligned_ref = best[1]  # Reference sequence with gaps
    alignment_length = len(aligned_query)

    identical = sum(
        1 for q, r in zip(aligned_query, aligned_ref) if q == r and q != "-"
    )

    identity = identical / alignment_length if alignment_length > 0 else 0.0
    return identity, alignment_length


# ─────────────────────────────────────────────────────────────────────────────
# STEP 4: COMPARE QUERY PROTEINS AGAINST REFERENCE FILES
# ─────────────────────────────────────────────────────────────────────────────


def find_orthologs(
    query_proteins: list[Protein],
    ref_path: Path,
    min_identity: float,
    use_mature: bool,
    blosum62,
) -> list[OrthoHit]:
    """
    Compares each query protein against all CDS in a reference file.
    Returns all comparisons that meet the minimum identity threshold.

    Args:
        query_proteins: List of Protein objects from the query GBK.
        ref_path:       Path to the reference GBK/GBFF file.
        min_identity:   Minimum identity (0.0–1.0) to report a hit.
        use_mature:     If True, uses mature_sequence for alignment, otherwise full.
        blosum62:       BLOSUM62 substitution matrix.

    Returns:
        List of OrthoHit objects for all hits above threshold.
    """
    hits = []

    for record in SeqIO.parse(ref_path, "genbank"):
        for feature in record.features:

            # Only compare against CDS features
            if feature.type != "CDS":
                continue

            ref_translation = feature.qualifiers.get("translation", [""])[0]
            if not ref_translation:
                continue

            ref_locus = feature.qualifiers.get("locus_tag", ["UNKNOWN"])[0]
            ref_product = feature.qualifiers.get("product", ["Unknown product"])[0]

            # Compare each query protein against this reference protein
            for qprotein in query_proteins:

                # Choose sequence for comparison
                query_seq = (
                    qprotein.mature_sequence if use_mature else qprotein.sequence
                )

                # Quick pre-filter: skip if length ratio is too extreme
                ratio = len(query_seq) / max(len(ref_translation), 1)
                if ratio < 0.3 or ratio > 3.0:
                    continue

                # Run pairwise alignment
                ref_seq = (
                    calculate_mature_core(ref_translation)
                    if use_mature
                    else ref_translation
                )
                identity, aln_length = calculate_identity(
                    query_seq, ref_translation, blosum62
                )

                # Apply threshold
                if identity >= min_identity:
                    hit = OrthoHit(
                        query_locus=qprotein.locus_tag,
                        query_product=qprotein.product,
                        query_seq=query_seq,
                        ref_locus=ref_locus,
                        ref_product=ref_product,
                        ref_seq=ref_seq,
                        ref_file=ref_path.stem,
                        identity=identity,
                        alignment_length=aln_length,
                        query_length=len(query_seq),
                        ref_length=len(ref_translation),
                    )
                    hits.append(hit)

    return hits


# ─────────────────────────────────────────────────────────────────────────────
# STEP 5: PRETTY PRINT TABLE (NO EXTERNAL DEPENDENCIES)
# ─────────────────────────────────────────────────────────────────────────────


def print_hits_table(hits: list[OrthoHit]) -> None:
    """
    Prints ortholog hits in a beautifully formatted table with aligned columns.
    Uses the same approach as extract_contig_genes_v2.py - no external dependencies.

    Sorted by query locus tag, then by identity (highest first).

    Args:
        hits: List of OrthoHit objects to display.
    """
    if not hits:
        print("No hits found.")
        return

    # Define column widths (in characters)
    col_query_locus = 12
    col_query_product = 20
    col_ref_locus = 15
    col_ref_product = 25
    col_ref_file = 23
    col_identity = 10
    col_aln_length = 10
    col_query_len = 10
    col_ref_len = 10

    # Print header
    header = (
        f"{'Query Locus':<{col_query_locus}} | "
        f"{'Query Product':<{col_query_product}} | "
        f"{'Ref Locus':<{col_ref_locus}} | "
        f"{'Ref Product':<{col_ref_product}} | "
        f"{'Ref File':<{col_ref_file}} | "
        f"{'Identity %':<{col_identity}} | "
        f"{'Aln Len':<{col_aln_length}} | "
        f"{'Query Len':<{col_query_len}} | "
        f"{'Ref Len':<{col_ref_len}}"
    )
    print(header)

    # Print separator
    total_width = len(header)
    print("-" * total_width)

    # Sort by query locus, then identity descending
    sorted_hits = sorted(hits, key=lambda h: (h.query_locus, -h.identity))

    # Print each hit
    for hit in sorted_hits:
        # Truncate long product names to fit columns
        query_prod = (
            hit.query_product[:20] + "..."
            if len(hit.query_product) > 23
            else hit.query_product
        )
        ref_prod = (
            hit.ref_product[:20] + "..."
            if len(hit.ref_product) > 23
            else hit.ref_product
        )

        row = (
            f"{hit.query_locus:<{col_query_locus}} | "
            f"{query_prod:<{col_query_product}} | "
            f"{hit.ref_locus:<{col_ref_locus}} | "
            f"{ref_prod:<{col_ref_product}} | "
            f"{hit.ref_file:<{col_ref_file}} | "
            f"{hit.identity*100:>{col_identity-1}.1f}% | "
            f"{hit.alignment_length:>{col_aln_length}} | "
            f"{hit.query_length:>{col_query_len}} | "
            f"{hit.ref_length:>{col_ref_len}}"
        )
        print(row)

    # Print footer
    print("-" * total_width)
    print(f"Total hits: {len(hits)}\n")


# ─────────────────────────────────────────────────────────────────────────────
# STEP 6: WRITE TSV OUTPUT (WITH PROTEIN SEQUENCES)
# ─────────────────────────────────────────────────────────────────────────────

TSV_HEADERS = [
    "query_locus",
    "query_product",
    "ref_locus",
    "ref_product",
    "ref_file",
    "query_sequence",
    "ref_sequence",
    "identity_pct",
    "alignment_length",
    "query_length",
    "ref_length",
]


def write_tsv(hits: list[OrthoHit], out_handle) -> None:
    """
    Writes all ortholog hits to a TSV file with protein sequences.
    Sorted by query locus tag, then by identity (highest first).

    Args:
        hits:       List of OrthoHit results to write.
        out_handle: Open file handle or sys.stdout.
    """
    writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
    writer.writerow(TSV_HEADERS)

    # Sort: by query locus first, then identity descending
    sorted_hits = sorted(hits, key=lambda h: (h.query_locus, -h.identity))

    for hit in sorted_hits:
        writer.writerow(
            [
                hit.query_locus,
                hit.query_product,
                hit.ref_locus,
                hit.ref_product,
                hit.ref_file,
                hit.query_seq,
                hit.ref_seq,
                f"{hit.identity * 100:.1f}",
                hit.alignment_length,
                hit.query_length,
                hit.ref_length,
            ]
        )


# ─────────────────────────────────────────────────────────────────────────────
# STEP 7: LOAD BLOSUM62 MATRIX (once at module level)
# ─────────────────────────────────────────────────────────────────────────────

_BLOSUM62 = substitution_matrices.load("BLOSUM62")


# ─────────────────────────────────────────────────────────────────────────────
# STEP 8: MAIN PIPELINE
# ─────────────────────────────────────────────────────────────────────────────


def main() -> None:
    args = get_args()

    print("=" * 100, file=sys.stderr)
    print("GBK ORTHOLOG FINDER v1.2", file=sys.stderr)
    print("=" * 100, file=sys.stderr)
    print(f"  Query:       {args.query}", file=sys.stderr)
    print(f"  Reference:   {args.reference}", file=sys.stderr)
    print(f"  Min identity: {args.identity * 100:.0f}%", file=sys.stderr)
    print(f"  Mature core:  {'YES' if args.mature else 'NO'}", file=sys.stderr)
    if args.max_length:
        print(f"  Max length:  {args.max_length} aa", file=sys.stderr)
    print("=" * 100, file=sys.stderr)

    try:
        # STEP A: Extract proteins from query GBK
        query_proteins = extract_proteins_from_gbk(
            gbk_path=args.query,
            apply_mature=args.mature,
            min_length=args.min_length,
            max_length=args.max_length,
        )

        if not query_proteins:
            sys.exit("[!] No proteins extracted from query file. Check the file.")

        # STEP B: Scan all reference files
        all_hits: list[OrthoHit] = []
        ref_files = list(stream_reference_files(args.reference))

        if not ref_files:
            sys.exit("[!] No valid reference files found. Check path/extensions.")

        for i, ref_file in enumerate(ref_files, start=1):
            print(
                f"[*] Scanning reference {i}/{len(ref_files)}: {ref_file.name}...",
                file=sys.stderr,
            )

            file_hits = find_orthologs(
                query_proteins=query_proteins,
                ref_path=ref_file,
                min_identity=args.identity,
                use_mature=args.mature,
                blosum62=_BLOSUM62,
            )

            print(
                f"    → {len(file_hits)} hit(s) above {args.identity*100:.0f}% identity",
                file=sys.stderr,
            )
            all_hits.extend(file_hits)

        # STEP C: Output results
        print(f"\n[*] Total hits found: {len(all_hits)}\n", file=sys.stderr)

        if args.output:
            # Save to TSV file with protein sequences
            print(
                f"[*] Saving results with protein sequences to: {args.output.resolve()}",
                file=sys.stderr,
            )
            with smart_open(args.output) as out_handle:
                write_tsv(all_hits, out_handle)
        else:
            # Print pretty table to terminal (no sequences)
            print(
                "[*] NOTE: Protein sequences will be included if you save to a file:",
                file=sys.stderr,
            )
            print(
                "    python3 gbk_ortholog_finder.py -q ... -r ... -o results.tsv",
                file=sys.stderr,
            )
            print(
                "    or redirect output: python3 gbk_ortholog_finder.py -q ... -r ... > results.txt",
                file=sys.stderr,
            )
            print()  # Blank line before table
            print_hits_table(all_hits)

        print("=" * 100, file=sys.stderr)
        print("[*] Done.", file=sys.stderr)
        print("=" * 100, file=sys.stderr)

    except ValueError as e:
        sys.exit(f"\n[!] Error: {e}")
    except KeyboardInterrupt:
        sys.exit("\n[!] Interrupted by user.")


if __name__ == "__main__":
    main()
