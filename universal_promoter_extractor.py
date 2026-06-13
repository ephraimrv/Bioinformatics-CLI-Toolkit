"""
Universal Promoter Extractor

A flexible pipeline asset for extracting upstream promoter regions from
GenBank files. Supports three complementary extraction modes that can
be used independently or combined in a single run:

  1. Keyword mode (-k/--keywords):
     Scans CDS /product annotations for matching keywords.
     Use for broad discovery when you don't know exact locus tags.

  2. Locus mode (-l/--loci):
     Extracts upstream of specific, named locus tags (exact match).
     Accepts one or more non-contiguous locus tags anywhere in the genome.
     Use when you already know which genes you want.

  3. Range mode (--range LOCUS_START LOCUS_END):
     Extracts upstream of every CDS between two locus tags (inclusive),
     ordered by genomic coordinate. Operates within a single contig.
     Use for operons, genomic islands, or any contiguous gene cluster.

All results are automatically deduplicated. Outputs a single MEME-compatible
FASTA file. Truncation warnings are printed when a gene is too close to a
contig boundary to provide the full upstream window.

License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).

Example usage:
    # Keyword mode (broad search)
    $ python3 universal_promoter_extractor.py -i references/ \\
      -o promoters.fasta -u 150 -k bacteriocin lactobin cerein

    # Locus mode (specific, non-contiguous genes)
    $ python3 universal_promoter_extractor.py -i C5_genome.gbk \\
      -o promoters.fasta -u 150 -l ctg1_68 ctg1_50 ctg1_100

    # Range mode (all genes between two locus tags)
    $ python3 universal_promoter_extractor.py -i C5_genome.gbk \\
      -o promoters.fasta -u 150 --range ctg1_50 ctg1_75

    # Combined: keyword discovery + specific loci of interest
    $ python3 universal_promoter_extractor.py -i references/ \\
      -o promoters.fasta -u 150 -k bacteriocin -l ctg1_68 ctg1_50
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.2.0"

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
from utils import stream_reference_files


def get_args() -> argparse.Namespace:
    """Configures the CLI and returns parsed arguments.

    At least one extraction mode (-k, -l, or --range) is required.
    Multiple modes can be combined in a single run.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Extract upstream promoter regions from GenBank files. "
            "Three modes available: keyword search (-k), specific locus tags (-l), "
            "or coordinate range (--range). At least one mode is required."
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
        help="Number of upstream base pairs to extract (Default: 150)",
    )

    modes = parser.add_argument_group(
        "Extraction Modes",
        "At least one mode is required. Multiple modes can be combined.",
    )
    modes.add_argument(
        "-k",
        "--keywords",
        type=str,
        nargs="+",
        default=None,
        help=(
            "Keyword mode: case-insensitive substring search in CDS /product annotations. "
            "Example: -k bacteriocin lactobin cerein"
        ),
    )
    modes.add_argument(
        "-l",
        "--loci",
        type=str,
        nargs="+",
        default=None,
        help=(
            "Locus mode: extract upstream of specific locus tags (exact match). "
            "Accepts one or more non-contiguous locus tags. "
            "Example: -l ctg1_68 ctg1_50 ctg1_100"
        ),
    )
    modes.add_argument(
        "--range",
        type=str,
        nargs=2,
        metavar=("LOCUS_START", "LOCUS_END"),
        default=None,
        help=(
            "Range mode: extract upstream of ALL CDS features between two locus tags "
            "(inclusive), sorted by genomic coordinate. Operates within a single contig. "
            "Example: --range ctg1_50 ctg1_75"
        ),
    )

    args = parser.parse_args()

    if not args.keywords and not args.loci and not args.range:
        parser.error(
            "At least one extraction mode is required: "
            "-k/--keywords, -l/--loci, or --range"
        )

    return args


def _safe_locus_tag(feature) -> str:
    """Return locus_tag, or a coordinate-based fallback if absent.

    Using a plain "UNKNOWN" fallback collapses ALL unannotated CDS features
    in a file into a single deduplication key, causing every one after the
    first to be silently skipped. The coordinate-based fallback guarantees
    each unannotated feature gets its own unique key.
    """
    start = int(feature.location.start)
    fallback = f"UNKNOWN_CDS_{start}"
    return feature.qualifiers.get("locus_tag", [fallback])[0]


def _extract_upstream(record, feature, upstream_bp: int) -> tuple[str, int, int]:
    """Extract upstream DNA sequence from a CDS feature.

    Args:
        record:      BioPython SeqRecord containing the feature.
        feature:     BioPython SeqFeature (CDS type).
        upstream_bp: Requested upstream window in bp.

    Returns:
        (upstream_sequence, actual_length_extracted, strand)
        actual_length_extracted may be less than upstream_bp at contig edges.
    """
    start = int(feature.location.start)
    end = int(feature.location.end)
    strand = feature.location.strand

    if strand == 1:
        slice_start = max(0, start - upstream_bp)
        actual_upstream = start - slice_start
        upstream_seq = str(record.seq[slice_start:start])
    else:
        slice_end = min(len(record.seq), end + upstream_bp)
        actual_upstream = slice_end - end
        raw_seq = record.seq[end:slice_end]
        upstream_seq = str(raw_seq.reverse_complement())

    return upstream_seq, actual_upstream, strand


def extract_by_keywords(
    gbk_path: Path, keywords: list[str], upstream_bp: int
) -> Iterator[tuple[str, str, str, str, int, int]]:
    """Extract upstream sequences for CDS features matching product keywords.

    Args:
        gbk_path:    Path to GenBank file.
        keywords:    Keywords matched against /product qualifier (case-insensitive).
        upstream_bp: Upstream window in bp.

    Yields:
        (record_id, locus_tag, product, upstream_seq, actual_upstream_length, strand)
    """
    try:
        with open(gbk_path, "r", encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                for feature in record.features:
                    if feature.type != "CDS":
                        continue
                    product = feature.qualifiers.get("product", [""])[0]
                    if not any(k.lower() in product.lower() for k in keywords):
                        continue

                    locus_tag = _safe_locus_tag(feature)
                    upstream_seq, actual_upstream, strand = _extract_upstream(
                        record, feature, upstream_bp
                    )

                    if actual_upstream < upstream_bp:
                        print(
                            f"      [!] Warning: {locus_tag} upstream truncated to "
                            f"{actual_upstream}bp (contig boundary — requested {upstream_bp}bp).",
                            file=sys.stderr,
                        )

                    yield record.id, locus_tag, product, upstream_seq, actual_upstream, strand

    except Exception as e:
        raise ValueError(f"Failed to parse {gbk_path.name}: {e}") from e


def extract_by_loci(
    gbk_path: Path, locus_tags: list[str], upstream_bp: int
) -> Iterator[tuple[str, str, str, str, int, int]]:
    """Extract upstream sequences for specific locus tags (exact match).

    Accepts one or more non-contiguous locus tags anywhere in the genome.
    Reports any requested tags that were not found in the file.

    Args:
        gbk_path:    Path to GenBank file.
        locus_tags:  Locus tags to find and extract (exact match).
        upstream_bp: Upstream window in bp.

    Yields:
        (record_id, locus_tag, product, upstream_seq, actual_upstream_length, strand)
    """
    target_set = set(locus_tags)
    found = set()

    try:
        with open(gbk_path, "r", encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                for feature in record.features:
                    if feature.type != "CDS":
                        continue
                    locus_tag = _safe_locus_tag(feature)
                    if locus_tag not in target_set:
                        continue

                    product = feature.qualifiers.get(
                        "product", ["hypothetical protein"]
                    )[0]
                    upstream_seq, actual_upstream, strand = _extract_upstream(
                        record, feature, upstream_bp
                    )
                    found.add(locus_tag)

                    if actual_upstream < upstream_bp:
                        print(
                            f"      [!] Warning: {locus_tag} upstream truncated to "
                            f"{actual_upstream}bp (contig boundary — requested {upstream_bp}bp).",
                            file=sys.stderr,
                        )

                    yield record.id, locus_tag, product, upstream_seq, actual_upstream, strand

        # Report any requested loci that were absent in this file
        for missing in sorted(target_set - found):
            print(
                f"      [!] {missing} not found in {gbk_path.name}",
                file=sys.stderr,
            )

    except Exception as e:
        raise ValueError(f"Failed to parse {gbk_path.name}: {e}") from e


def extract_by_range(
    gbk_path: Path, locus_start: str, locus_end: str, upstream_bp: int
) -> Iterator[tuple[str, str, str, str, int, int]]:
    """Extract upstream sequences for all CDS between two locus tags.

    Collects all CDS features in the same contig, sorts by genomic coordinate,
    then yields upstream sequences for every CDS in the inclusive range
    [LOCUS_START, LOCUS_END]. Order of the two boundary tags does not matter.

    Note:
        Operates within a single contig. If the two boundary tags reside on
        different contigs (fragmented assembly), the range cannot be resolved
        and a warning is printed. This is intentional: cross-contig ranges have
        no defined genomic order.

    Args:
        gbk_path:    Path to GenBank file.
        locus_start: First boundary locus tag (inclusive).
        locus_end:   Second boundary locus tag (inclusive).
        upstream_bp: Upstream window in bp.

    Yields:
        (record_id, locus_tag, product, upstream_seq, actual_upstream_length, strand)
    """
    found_start = False
    found_end = False

    try:
        with open(gbk_path, "r", encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                # Collect all CDS in this record with coordinates and tags
                all_cds = []
                for feature in record.features:
                    if feature.type != "CDS":
                        continue
                    tag = _safe_locus_tag(feature)
                    all_cds.append((int(feature.location.start), tag, feature))

                # Sort by genomic coordinate
                all_cds.sort(key=lambda x: x[0])

                # Find boundary indices in this record
                start_idx = next(
                    (i for i, (_, tag, _) in enumerate(all_cds) if tag == locus_start),
                    None,
                )
                end_idx = next(
                    (i for i, (_, tag, _) in enumerate(all_cds) if tag == locus_end),
                    None,
                )

                if start_idx is None or end_idx is None:
                    # Both boundaries not in this record — skip, try next record
                    continue

                found_start = True
                found_end = True

                # Allow user to specify boundaries in any order
                lo, hi = sorted([start_idx, end_idx])
                range_count = hi - lo + 1

                print(
                    f"      [Range] {locus_start} → {locus_end}: "
                    f"{range_count} CDS features found in {record.id}",
                    file=sys.stderr,
                )

                for _, tag, feature in all_cds[lo : hi + 1]:
                    product = feature.qualifiers.get(
                        "product", ["hypothetical protein"]
                    )[0]
                    upstream_seq, actual_upstream, strand = _extract_upstream(
                        record, feature, upstream_bp
                    )

                    if actual_upstream < upstream_bp:
                        print(
                            f"      [!] Warning: {tag} upstream truncated to "
                            f"{actual_upstream}bp (contig boundary — requested {upstream_bp}bp).",
                            file=sys.stderr,
                        )

                    yield record.id, tag, product, upstream_seq, actual_upstream, strand

        if not found_start or not found_end:
            print(
                f"      [!] Warning: Range boundaries '{locus_start}' and '{locus_end}' "
                f"were not both found in the same contig of {gbk_path.name}. "
                f"Cross-contig ranges cannot be resolved.",
                file=sys.stderr,
            )

    except Exception as e:
        raise ValueError(f"Failed to parse {gbk_path.name}: {e}") from e


def main() -> None:
    """Routes to extraction modes, deduplicates, and writes MEME-ready FASTA."""
    args = get_args()

    print(f"[*] Target          : {args.input}", file=sys.stderr)
    print(f"[*] Upstream window : {args.upstream}bp", file=sys.stderr)
    if args.keywords:
        print(f"[*] Mode: Keywords  → {args.keywords}", file=sys.stderr)
    if args.loci:
        print(f"[*] Mode: Loci      → {args.loci}", file=sys.stderr)
    if args.range:
        print(
            f"[*] Mode: Range     → {args.range[0]} to {args.range[1]}", file=sys.stderr
        )
    print(file=sys.stderr)

    hits_found = 0
    duplicates_skipped = 0
    seen_loci: set[tuple[str, str]] = set()

    try:
        with open(args.output, "w", encoding="utf-8") as out_file:

            for file_path in stream_reference_files(args.input):

                # Skip nucleotide FASTA — no upstream DNA map available
                if file_path.suffix.lower() in (".fasta", ".fa", ".fna"):
                    print(
                        f"  [!] Skipping {file_path.name}: nucleotide FASTA cannot "
                        f"provide upstream DNA context.",
                        file=sys.stderr,
                    )
                    continue

                # Skip protein FASTA — no genomic coordinates
                if file_path.suffix.lower() in (".faa",):
                    print(
                        f"  [!] Skipping {file_path.name}: protein FASTA has no "
                        f"genomic coordinates.",
                        file=sys.stderr,
                    )
                    continue

                print(f"  -> Parsing {file_path.name}...", file=sys.stderr)

                # Build list of active iterators for this file
                iterators = []
                if args.keywords:
                    iterators.append(
                        extract_by_keywords(file_path, args.keywords, args.upstream)
                    )
                if args.loci:
                    iterators.append(
                        extract_by_loci(file_path, args.loci, args.upstream)
                    )
                if args.range:
                    iterators.append(
                        extract_by_range(
                            file_path, args.range[0], args.range[1], args.upstream
                        )
                    )

                for iterator in iterators:
                    for seq_id, locus, prod, seq, actual_up, strand in iterator:

                        # File-aware deduplication — unique per (file, locus_tag)
                        dedup_key = (file_path.stem, locus)
                        if dedup_key in seen_loci:
                            duplicates_skipped += 1
                            continue

                        seen_loci.add(dedup_key)
                        hits_found += 1

                        # Sanitize product name for FASTA header
                        clean_prod = re.sub(r"[^\w\-]", "_", prod)
                        strand_label = strand if strand in (1, -1) else "unknown"

                        fasta_header = (
                            f">{seq_id}_{locus}_{clean_prod}"
                            f"_upstream_{actual_up}bp_gene_strand_{strand_label}"
                        )
                        out_file.write(f"{fasta_header}\n{seq}\n")
                        print(
                            f"      [Hit] {locus} | strand {strand_label} | {prod[:50]}",
                            file=sys.stderr,
                        )

        print("\n" + "=" * 60, file=sys.stderr)
        print(
            f"[*] SUCCESS: {hits_found} unique upstream regions extracted.",
            file=sys.stderr,
        )
        if duplicates_skipped > 0:
            print(
                f"[*] Note: {duplicates_skipped} duplicate(s) skipped.",
                file=sys.stderr,
            )
        print(f"[*] Output: {args.output.resolve()}", file=sys.stderr)
        print("=" * 60, file=sys.stderr)

    except ValueError as e:
        sys.exit(f"\n[!] Pipeline Error: {e}")
    except KeyboardInterrupt:
        sys.exit("\n[!] Pipeline interrupted by user.")


if __name__ == "__main__":
    main()
