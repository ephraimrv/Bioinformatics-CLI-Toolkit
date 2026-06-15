r"""
Universal Promoter Extractor

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
        Pass 2 — mRNA features: for locus tags found in pass 1, extract
                 upstream of the mRNA start (TSS).
    This guarantees that keyword matching always uses the /product annotation
    regardless of which feature type carries it, while coordinate extraction
    always uses the biologically correct mRNA anchor.

License: MIT

Note:
    This module is part of ongoing research and is associated with an upcoming
    publication. Please cite appropriately when used in derivative works.
    See LICENSE file in the repository root for full license terms.

Example usage:
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
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.3.0"

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

# ── Mode constants ─────────────────────────────────────────────────────────────

MODE_AUTO = "auto"
MODE_PROKARYOTE = "prokaryote"
MODE_EUKARYOTE = "eukaryote"
_VALID_MODES = (MODE_AUTO, MODE_PROKARYOTE, MODE_EUKARYOTE)


# ── Auto-detection ─────────────────────────────────────────────────────────────


def _detect_organism_mode(gbk_path: Path) -> str:
    """Detect whether a GenBank file is prokaryotic or eukaryotic.

    Detection is based on the presence of ``mRNA`` features, which are
    absent from all prokaryotic GenBank files produced by standard
    annotation pipelines (Prokka, Bakta, NCBI prokaryote annotation)
    and always present in eukaryotic GenBank files (Ensembl, RefSeq
    eukaryote, Augustus, MAKER, etc.).

    Only the first record in the file is scanned for speed. This is
    reliable for all standard assemblies (single-chromosome files,
    multi-contig assemblies) because annotation style is uniform
    across records. For heterogeneous files, use ``--mode`` explicitly.

    Args:
        gbk_path: Path to the GenBank file to inspect.

    Returns:
        ``"eukaryote"`` if mRNA features are detected, ``"prokaryote"`` otherwise.
    """
    try:
        for record in SeqIO.parse(gbk_path, "genbank"):
            feature_types = {f.type for f in record.features}
            return MODE_EUKARYOTE if "mRNA" in feature_types else MODE_PROKARYOTE
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
    if strand == 1:
        slice_start = max(0, start - upstream_bp)
        actual_upstream = start - slice_start
        upstream_seq = str(record.seq[slice_start:start])
    else:
        slice_end = min(len(record.seq), end + upstream_bp)
        actual_upstream = slice_end - end
        upstream_seq = str(record.seq[end:slice_end].reverse_complement())

    if actual_upstream < upstream_bp:
        print(
            f"      [!] Warning: {locus_tag} upstream truncated to "
            f"{actual_upstream}bp (contig boundary — requested {upstream_bp}bp).",
            file=sys.stderr,
        )

    return upstream_seq, actual_upstream


# ── Core extraction functions ──────────────────────────────────────────────────


def extract_regulatory_regions(
    gbk_path: Path,
    keywords: list[str],
    upstream_bp: int,
    mode: str = MODE_AUTO,
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

    Tracks actual extracted length separately from the requested window.
    These differ when a gene is within upstream_bp bases of a contig boundary.

    Args:
        gbk_path:    Path to the target .gbk or .gbff file.
        keywords:    Keywords to match against /product annotations (case-insensitive).
        upstream_bp: Number of base pairs to extract upstream of the anchor coordinate.
        mode:        Organism mode: ``"auto"``, ``"prokaryote"``, or ``"eukaryote"``.

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

                        locus_tag = feature.qualifiers.get("locus_tag", ["UNKNOWN"])[0]
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
                    keyword_loci: dict[str, str] = {}  # {locus_tag: product}
                    for feature in record.features:
                        if feature.type != "CDS":
                            continue
                        product = feature.qualifiers.get("product", [""])[0]
                        if any(k.lower() in product.lower() for k in keywords):
                            lt = feature.qualifiers.get("locus_tag", ["UNKNOWN"])[0]
                            keyword_loci[lt] = product

                    if not keyword_loci:
                        continue

                    # Pass 2 — extract upstream of mRNA start (TSS) for matched loci.
                    for feature in record.features:
                        if feature.type != "mRNA":
                            continue

                        locus_tag = feature.qualifiers.get("locus_tag", ["UNKNOWN"])[0]
                        if locus_tag not in keyword_loci:
                            continue

                        product = keyword_loci[locus_tag]
                        start = int(feature.location.start)
                        end = int(feature.location.end)
                        strand = feature.location.strand

                        upstream_seq, actual_upstream = _extract_upstream_seq(
                            record, start, end, strand, upstream_bp, locus_tag
                        )
                        yield record.id, locus_tag, product, upstream_seq, actual_upstream

    except Exception as e:
        raise ValueError(f"Failed to parse {gbk_path.name}: {e}") from e


def extract_by_loci(
    gbk_path: Path,
    locus_tags: list[str],
    upstream_bp: int,
    mode: str = MODE_AUTO,
) -> Iterator[tuple[str, str, str, str, int, int, str]]:
    """Extract upstream regions for a specific list of locus tags.

    The programmatic counterpart to ``extract_regulatory_regions``. While
    ``extract_regulatory_regions`` discovers targets by keyword, this function
    extracts targets whose locus tags are already known — making it the correct
    tool for bridge scripts that receive locus tags from ``gbk_ortholog_finder``.

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
        Pass 2 — ``mRNA`` features: for matching locus_tags, extracts
                 upstream sequence relative to the mRNA start (= TSS).

    Args:
        gbk_path:    Path to the GenBank file to scan.
        locus_tags:  Locus tags to extract. Duplicates are silently removed.
        upstream_bp: Number of bases to extract upstream of the anchor coordinate.
        mode:        Organism mode: ``"auto"``, ``"prokaryote"``, or ``"eukaryote"``.

    Yields:
        A 7-item tuple:
        ``(seq_id, locus_tag, product, upstream_seq, actual_upstream, strand, genome_label)``

        - seq_id:          Contig/record ID from the GenBank file.
        - locus_tag:       The matched locus tag.
        - product:         The /product annotation (from CDS in eukaryote mode).
        - upstream_seq:    Strand-corrected upstream DNA sequence.
        - actual_upstream: Actual extracted length; may be < upstream_bp near
                           contig boundaries.
        - strand:          1 for forward strand, -1 for reverse strand.
        - genome_label:    Stem of the GenBank filename for FASTA header use.

    Raises:
        ValueError: If the GenBank file cannot be parsed.

    Notes:
        - If any requested locus tag is not found in the file, a warning is
          printed to stderr after scanning completes, listing all missing tags.
        - Duplicate locus tags in the input list are silently deduplicated.
        - If the same locus tag appears on multiple features of the same type
          in the GenBank file (malformed annotation), only the first occurrence
          is yielded and a warning is printed for subsequent duplicates.
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

                        locus_tag = feature.qualifiers.get("locus_tag", ["UNKNOWN"])[0]
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
                    cds_products: dict[str, str] = {}
                    for feature in record.features:
                        if feature.type != "CDS":
                            continue
                        lt = feature.qualifiers.get("locus_tag", ["UNKNOWN"])[0]
                        if lt in remaining:
                            cds_products[lt] = feature.qualifiers.get(
                                "product", ["Unknown product"]
                            )[0]

                    # Pass 2 — extract upstream of mRNA start (TSS).
                    for feature in record.features:
                        if not remaining:
                            break
                        if feature.type != "mRNA":
                            continue

                        locus_tag = feature.qualifiers.get("locus_tag", ["UNKNOWN"])[0]
                        if locus_tag not in target_set:
                            continue

                        if locus_tag in already_yielded:
                            print(
                                f"      [!] Warning: {locus_tag} appears on multiple "
                                f"mRNA features in {gbk_path.name} — skipping duplicate.",
                                file=sys.stderr,
                            )
                            continue

                        product = cds_products.get(locus_tag, "Unknown product")
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

    except Exception as e:
        raise ValueError(f"Failed to parse {gbk_path.name}: {e}") from e

    if remaining:
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

    return parser.parse_args()


# ── Main ───────────────────────────────────────────────────────────────────────


def main() -> None:
    """Coordinate file routing, deduplication, header sanitization, and output."""
    args = get_args()

    print(f"[*] Scanning target    : {args.input}", file=sys.stderr)
    print(f"[*] Upstream window    : {args.upstream}bp", file=sys.stderr)
    print(f"[*] Keywords           : {args.keywords}", file=sys.stderr)
    print(f"[*] Mode               : {args.mode}", file=sys.stderr)
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
                    file_path, args.keywords, args.upstream, mode=file_mode
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
