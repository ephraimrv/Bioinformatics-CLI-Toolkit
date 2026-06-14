"""
Universal Promoter Extractor

A flexible pipeline asset that dynamically scans either a single GenBank file
or an entire directory for target genes using dynamic keyword matching.
Extracts customizable upstream promoter regions and formats them into a
single MEME-compatible FASTA file for motif discovery. Includes automatic
deduplication to prevent skewed MEME statistical calculations.

License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).

Example usage:
    $ python3 universal_promoter_extractor.py -i C5_prokka_result.gbk \
      -o C5_promoters.fasta -u 150 -k bacteriocin lactobin cerein
      
    $ python3 universal_promoter_extractor.py -i references/ \
      -o C5_promoters.fasta -u 150 -k bacteriocin lactobin cerein
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.2.0"

import sys
import argparse
import re
from pathlib import Path
from typing import Iterator
from Bio import SeqIO
from utils import stream_reference_files


def get_args() -> argparse.Namespace:
    """
    Configures the CLI, handles input validation, and returns parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Extracts upstream regulatory regions based on keyword matching.",
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
    parser.add_argument(
        "-k",
        "--keywords",
        type=str,
        nargs="+",
        required=True,
        help="List of keywords to search for in product names (e.g., bacteriocin lactobin cerein)",
    )

    return parser.parse_args()


def extract_regulatory_regions(
    gbk_path: Path, keywords: list[str], upstream_bp: int
) -> Iterator[tuple[str, str, str, str, int]]:
    """
    Scans a GenBank file for specific keywords in CDS annotations and extracts
    the upstream DNA sequences (promoter regions) for MEME motif discovery.

    Tracks actual extracted length separately from the requested window —
    these differ when a gene is near a contig boundary. A warning is printed
    to stderr when truncation occurs so the user can exclude or flag those
    sequences before running MEME.

    Args:
        gbk_path:    Path object pointing to the target .gbk or .gbff file.
        keywords:    A list of string keywords to match against /product annotation.
        upstream_bp: The number of base pairs to extract upstream of the start codon.

    Yields:
        A 5-item tuple:
        (Sequence ID, Locus Tag, Product Annotation, Upstream DNA Sequence,
         Actual Extracted Length)

    Raises:
        ValueError: If the GenBank file is malformed, structurally invalid, or unreadable.
    """
    try:
        with open(gbk_path, "r", encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":

                        product = feature.qualifiers.get("product", [""])[0].lower()

                        if any(k.lower() in product for k in keywords):
                            locus_tag = feature.qualifiers.get(
                                "locus_tag", ["UNKNOWN"]
                            )[0]
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

                            # Warn if the upstream window was truncated by a contig boundary
                            if actual_upstream < upstream_bp:
                                print(
                                    f"      [!] Warning: {locus_tag} upstream truncated to "
                                    f"{actual_upstream}bp (contig boundary — requested {upstream_bp}bp).",
                                    file=sys.stderr,
                                )

                            yield record.id, locus_tag, product, upstream_seq, actual_upstream

    except Exception as e:
        raise ValueError(f"Failed to parse {gbk_path.name}: {e}") from e


def extract_by_loci(
    gbk_path: Path,
    locus_tags: list[str],
    upstream_bp: int,
) -> Iterator[tuple[str, str, str, str, int, int, str]]:
    """Extract upstream regions for a specific list of locus tags.

    This is the programmatic counterpart to ``extract_regulatory_regions``.
    While ``extract_regulatory_regions`` searches by keyword, this function
    searches by exact locus tag — making it suitable for programmatic pipelines
    where the target loci have already been identified (e.g., by
    ``gbk_ortholog_finder.find_orthologs()``).

    Designed to be imported and called by bridge scripts such as
    ``target_promoter_pipeline.py``. Does NOT require keyword input.

    PROKARYOTIC NOTE:
        In prokaryotes, the CDS start coordinate directly precedes the
        promoter, so extracting upstream of the CDS start is correct.
        For eukaryotic genomes, this logic would need updating to use
        the ``mRNA`` feature start (= TSS) instead of the ``CDS`` start,
        because eukaryotic CDS features begin at the ATG, which can be
        separated from the TSS by a 5' UTR or introns.

    Args:
        gbk_path:    Path to the GenBank file to scan.
        locus_tags:  List of locus tags to extract upstream sequences for.
        upstream_bp: Number of bases to extract upstream of the CDS start.

    Yields:
        A 7-item tuple:
        (seq_id, locus_tag, product, upstream_seq, actual_upstream, strand, genome_label)

        - seq_id:          Contig/record ID from the GenBank file
        - locus_tag:       The matched locus tag
        - product:         The /product annotation of the CDS
        - upstream_seq:    The extracted upstream DNA sequence (strand-corrected)
        - actual_upstream: Actual extracted length (may be < upstream_bp near
                           contig boundaries)
        - strand:          1 for forward strand, -1 for reverse strand
        - genome_label:    Stem of the GenBank filename (useful for FASTA headers
                           in multi-genome pipelines)

    Raises:
        ValueError: If the GenBank file cannot be parsed.
    """
    target_set = set(locus_tags)
    genome_label = gbk_path.stem

    try:
        with open(gbk_path, "r", encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                for feature in record.features:
                    if feature.type != "CDS":
                        continue

                    locus_tag = feature.qualifiers.get("locus_tag", ["UNKNOWN"])[0]
                    if locus_tag not in target_set:
                        continue

                    product = feature.qualifiers.get("product", ["Unknown product"])[0]
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

                    if actual_upstream < upstream_bp:
                        print(
                            f"      [!] Warning: {locus_tag} upstream truncated to "
                            f"{actual_upstream}bp (contig boundary — requested {upstream_bp}bp).",
                            file=sys.stderr,
                        )

                    yield (
                        record.id,
                        locus_tag,
                        product,
                        upstream_seq,
                        actual_upstream,
                        strand,
                        genome_label,
                    )

    except Exception as e:
        raise ValueError(f"Failed to parse {gbk_path.name}: {e}") from e


def main() -> None:
    """
    The Wrapper: Coordinates file routing, handles cross-file deduplication,
    sanitizes FASTA headers via Regex, and outputs MEME-ready sequences.
    """
    args = get_args()

    print(f"[*] Scanning target: {args.input}", file=sys.stderr)
    print(f"[*] Upstream extraction window: {args.upstream} bp", file=sys.stderr)
    print(f"[*] Active keywords: {args.keywords}\n", file=sys.stderr)

    hits_found = 0
    duplicates_skipped = 0
    seen_loci = set()

    try:
        with open(args.output, "w", encoding="utf-8") as out_file:

            for file_path in stream_reference_files(args.input):

                # Skip FASTA files because they don't contain upstream DNA maps
                if file_path.suffix.lower() in (".fasta", ".fa", ".faa"):
                    print(
                        f"  [!] Skipping {file_path.name}: Cannot extract upstream DNA from FASTA format.",
                        file=sys.stderr,
                    )
                    continue

                print(f"  -> Parsing {file_path.name}...", file=sys.stderr)

                for seq_id, locus, prod, seq, actual_up in extract_regulatory_regions(
                    file_path, args.keywords, args.upstream
                ):
                    # File-aware Deduplication Engine
                    dedup_key = (file_path.stem, locus)
                    if dedup_key in seen_loci:
                        duplicates_skipped += 1
                        continue

                    seen_loci.add(dedup_key)
                    hits_found += 1

                    # replace all non-word/non-hyphen characters with underscores
                    clean_prod = re.sub(r"[^\w\-]", "_", prod)

                    # Header uses actual extracted length — may differ from requested
                    # upstream if gene is near a contig boundary
                    fasta_header = f">{seq_id}_{locus}_{clean_prod}_up{actual_up}"

                    out_file.write(f"{fasta_header}\n{seq}\n")
                    print(f"      [Hit] {locus} | {prod[:40]}...", file=sys.stderr)

        print("\n" + "=" * 50, file=sys.stderr)
        print(
            f"[*] SUCCESS: {hits_found} unique regulatory regions extracted.",
            file=sys.stderr,
        )

        if duplicates_skipped > 0:
            print(
                f"[*] WARNING: {duplicates_skipped} duplicate sequences were safely skipped.",
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
