r"""
Targeted Promoter Pipeline

A bridge script that connects gbk_ortholog_finder.py and
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

WHY THIS REPLACES homology_extractor.py:
    The old script used:
        if core_peptide in translation:
    This is exact substring matching — one amino acid substitution and
    the ortholog is silently missed. Real homology requires alignment,
    which is what gbk_ortholog_finder.py provides.

License: MIT

Note:
    This module is part of ongoing research and is associated with an upcoming
    publication. Please cite appropriately when used in derivative works.
    See LICENSE file in the repository root for full license terms.

Example usage:
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
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.1.0"

import sys
import argparse
from pathlib import Path

# ── Import from existing toolkit scripts ──────────────────────────────────────
# These are your two production scripts. This bridge adds NO new alignment or
# extraction logic — it only connects them. Any improvements made to either
# script are automatically inherited here.

try:
    from gbk_ortholog_finder import extract_proteins_from_gbk, find_orthologs
except ImportError as e:
    sys.exit(
        f"[!] Cannot import from gbk_ortholog_finder.py.\n"
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
            "Output FASTA file path for MEME motif discovery. "
            "If omitted, saves to 'targeted_promoters.fasta' in the current directory."
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

    return parser


# ── Main ──────────────────────────────────────────────────────────────────────


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    output = args.output or Path("targeted_promoters.fasta")

    # Validate inputs
    if not args.query.exists():
        sys.exit(f"[!] Query file not found: {args.query}")
    if not args.reference.exists():
        sys.exit(f"[!] Reference path not found: {args.reference}")

    print("=" * 60, file=sys.stderr)
    print("  TARGETED PROMOTER PIPELINE", file=sys.stderr)
    print("=" * 60, file=sys.stderr)
    print(f"  Query           : {args.query.name}", file=sys.stderr)
    print(f"  Reference       : {args.reference}", file=sys.stderr)
    print(f"  Output          : {output}", file=sys.stderr)
    print(f"  Upstream        : {args.upstream}bp", file=sys.stderr)
    print(f"  Min identity    : {args.identity * 100:.0f}%", file=sys.stderr)
    print(f"  Min coverage    : {args.coverage * 100:.0f}% (min mode)", file=sys.stderr)
    print(
        f"  Mature core     : {'YES — leader peptides trimmed' if args.mature else 'NO'}",
        file=sys.stderr,
    )
    print(f"  Organism mode   : {args.mode}", file=sys.stderr)
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

    # find_orthologs() in gbk_ortholog_finder.py calls SeqIO.parse() directly
    # and expects a single file path — it does not handle directories.
    # We resolve the directory here and call find_orthologs() per file,
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
        file_hits = find_orthologs(
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

    # Report what was found
    for hit in sorted(hits, key=lambda h: (-h.identity, h.ref_locus)):
        print(
            f"      {hit.ref_locus:<20} {hit.identity*100:.1f}% identity  "
            f"{hit.ref_product[:45]}",
            file=sys.stderr,
        )

    # Deduplicate locus tags (one query can hit the same ref locus multiple times)
    target_loci = list(dict.fromkeys(hit.ref_locus for hit in hits))
    print(
        f"\n    {len(target_loci)} unique locus tag(s) to extract: {target_loci}",
        file=sys.stderr,
    )

    # ── Step 3: Extract upstream regions ─────────────────────────────────────
    print(
        f"\n[Step 3] Extracting {args.upstream}bp upstream of identified loci...",
        file=sys.stderr,
    )

    # ref_files already resolved in Step 2 — reused here directly.
    extracted_count = 0

    try:
        with open(output, "w", encoding="utf-8") as out_file:
            for ref_file in ref_files:
                for (
                    seq_id,
                    locus,
                    product,
                    seq,
                    actual_up,
                    strand,
                    genome_label,
                ) in extract_by_loci(
                    ref_file, target_loci, args.upstream, mode=args.mode
                ):

                    strand_symbol = "+" if strand == 1 else "-"
                    header = (
                        f">{locus} | {genome_label} | "
                        f"{actual_up}bp upstream | strand {strand_symbol} | "
                        f"{product[:45]}"
                    )
                    out_file.write(f"{header}\n{seq}\n")
                    extracted_count += 1
                    print(
                        f"      [Extracted] {locus} ({actual_up}bp, strand "
                        f"{strand_symbol}) from {genome_label}",
                        file=sys.stderr,
                    )

    except ValueError as e:
        sys.exit(f"[!] Extraction error: {e}")

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"\n{'=' * 60}", file=sys.stderr)
    print(f"  PIPELINE COMPLETE", file=sys.stderr)
    print(f"{'=' * 60}", file=sys.stderr)
    print(f"  Orthologs found    : {len(hits)}", file=sys.stderr)
    print(f"  Unique loci        : {len(target_loci)}", file=sys.stderr)
    print(f"  Regions extracted  : {extracted_count}", file=sys.stderr)
    print(f"  Output written to  : {output.resolve()}", file=sys.stderr)
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
