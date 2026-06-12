"""
Genomic CDS Region Extractor

Parses GenBank assemblies to extract and export coding sequences (CDS).

Supports two extraction modes:

  MODE 1 — PREFIX MODE (original behaviour):
    Extract by locus tag prefix and gene number range.
    Use this when you know the locus tag prefix and gene numbers.

  MODE 2 — COORDINATE MODE (new):
    Extract all CDS features that fall within a chromosomal coordinate range.
    Use this when you know the genomic coordinates but not the locus tag numbers.
    Works on ANY gbff/gbk file without needing to know the prefix first.

Example Usage:

  PREFIX MODE (original):
    $ python3 extract_gbk_region.py -i input.gbk -p "ctg1_" -s 46 -e 74 -o C5_locus
    $ python3 extract_gbk_region.py -i GCF_056530425.1.gbff -p "RHP56_RS" -s 340 -e 450 -o NCBI_locus

  COORDINATE MODE (new):
    $ python3 extract_gbk_region.py -i GCF_056530425.1.gbff --c1 53317 --c2 78823 -o Region1
    $ python3 extract_gbk_region.py -i any_genome.gbff --c1 100000 --c2 150000 -o my_region

  COORDINATE MODE with sequence filter (multi-record gbff):
    $ python3 extract_gbk_region.py -i GCF_056530425.1.gbff --c1 53317 --c2 78823 --seq NZ_CP134351.1 -o Region1

License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "2.0.0"

import sys
import argparse
from pathlib import Path
from Bio import SeqIO
from utils import base_parser, wrap_fasta


def get_args() -> argparse.Namespace:
    """Configures the CLI and returns parsed arguments."""
    parser = base_parser(
        description_text=(
            "Extracts bifunctional (DNA/Protein) CDS features from a GenBank file. "
            "Supports two modes: PREFIX MODE (by locus tag) and COORDINATE MODE (by genomic position)."
        ),
    )

    # ── PREFIX MODE arguments ─────────────────────────────────────────────────
    prefix_group = parser.add_argument_group(
        "PREFIX MODE",
        "Extract by locus tag prefix and gene number range. "
        "Use when you know the locus tag prefix (e.g. 'ctg1_' or 'RHP56_RS').",
    )
    prefix_group.add_argument(
        "-p",
        "--prefix",
        type=str,
        required=False,
        default=None,
        help="The locus tag prefix (e.g., 'ctg1_' or 'RHP56_RS')",
    )
    prefix_group.add_argument(
        "-s",
        "--start",
        type=int,
        required=False,
        default=None,
        help="Start gene number (Default: First available in prefix)",
    )
    prefix_group.add_argument(
        "-e",
        "--end",
        type=int,
        required=False,
        default=None,
        help="End gene number (Default: Last available in prefix)",
    )

    # ── COORDINATE MODE arguments ─────────────────────────────────────────────
    coord_group = parser.add_argument_group(
        "COORDINATE MODE",
        "Extract all CDS features within a chromosomal coordinate window. "
        "Use when you know the genomic positions but not the locus tag numbers. "
        "Works on any gbff/gbk file without needing a prefix.",
    )
    coord_group.add_argument(
        "--c1",
        type=int,
        required=False,
        default=None,
        metavar="START_BP",
        help="Chromosomal start coordinate in base pairs (e.g., 53317)",
    )
    coord_group.add_argument(
        "--c2",
        type=int,
        required=False,
        default=None,
        metavar="END_BP",
        help="Chromosomal end coordinate in base pairs (e.g., 78823)",
    )
    coord_group.add_argument(
        "--seq",
        type=str,
        required=False,
        default=None,
        metavar="SEQUENCE_ID",
        help=(
            "Sequence ID to target within a multi-record gbff file "
            "(e.g., 'NZ_CP134351.1'). If omitted, searches all records. "
            "Use this to target the chromosome vs the plasmid."
        ),
    )

    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> str:
    """
    Determines extraction mode and validates argument combinations.

    Returns:
        'prefix' or 'coord' indicating the extraction mode to use.

    Exits with an error message if arguments are invalid or contradictory.
    """
    using_prefix = args.prefix is not None
    using_coord = args.c1 is not None or args.c2 is not None

    # Must use at least one mode
    if not using_prefix and not using_coord:
        sys.exit(
            "\n[!] Error: You must specify either PREFIX MODE (-p) or COORDINATE MODE (--c1/--c2)."
            "\n    PREFIX:     python3 extract_gbk_region.py -i genome.gbff -p 'RHP56_RS' -s 340 -e 450 -o out"
            "\n    COORDINATE: python3 extract_gbk_region.py -i genome.gbff --c1 53317 --c2 78823 -o out"
        )

    # Cannot use both modes simultaneously
    if using_prefix and using_coord:
        sys.exit(
            "\n[!] Error: Cannot use PREFIX MODE and COORDINATE MODE together. "
            "Choose one or the other."
        )

    # PREFIX MODE validation
    if using_prefix:
        if args.start is not None and args.start < 1:
            sys.exit("\n[!] Error: --start must be a positive integer.")
        if args.end is not None and args.end < 1:
            sys.exit("\n[!] Error: --end must be a positive integer.")
        if args.start is not None and args.end is not None and args.start > args.end:
            sys.exit(
                f"\n[!] Error: Start gene ({args.start}) cannot be greater than "
                f"End gene ({args.end})."
            )
        return "prefix"

    # COORDINATE MODE validation
    if using_coord:
        if args.c1 is None or args.c2 is None:
            sys.exit(
                "\n[!] Error: COORDINATE MODE requires both --c1 (start) and --c2 (end). "
                "\n    Example: --c1 53317 --c2 78823"
            )
        if args.c1 < 0:
            sys.exit("\n[!] Error: --c1 must be a non-negative integer.")
        if args.c2 < 1:
            sys.exit("\n[!] Error: --c2 must be a positive integer.")
        if args.c1 >= args.c2:
            sys.exit(
                f"\n[!] Error: --c1 ({args.c1}) must be less than --c2 ({args.c2})."
            )
        return "coord"

    # Should never reach here
    sys.exit("\n[!] Error: Unknown argument state.")


def run_prefix_mode(args: argparse.Namespace, ffn_path: Path, faa_path: Path) -> None:
    """
    Extracts CDS features by locus tag prefix and gene number range.
    This is the original extraction behaviour.
    """
    start_label = args.start if args.start is not None else "[First Available]"
    end_label = args.end if args.end is not None else "[Last Available]"

    print(f"[*] MODE: PREFIX", file=sys.stderr)
    print(f"[*] Reading {args.input}...", file=sys.stderr)
    print(
        f"[*] Extracting {args.prefix}{start_label} to {args.prefix}{end_label}...\n",
        file=sys.stderr,
    )

    extracted_count = 0
    extracted_nums = []
    global_cds_counter = 1

    with open(ffn_path, "w", encoding="utf-8") as f_dna, open(
        faa_path, "w", encoding="utf-8"
    ) as f_prot:

        for record in SeqIO.parse(args.input, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":

                    locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                    product = feature.qualifiers.get(
                        "product", ["hypothetical protein"]
                    )[0]
                    translation = feature.qualifiers.get("translation", [""])[0]
                    protein_id = feature.qualifiers.get(
                        "protein_id", [f"SEQ_{global_cds_counter:05d}"]
                    )[0]
                    global_cds_counter += 1

                    if locus_tag.startswith(args.prefix) and translation:
                        raw_suffix = locus_tag[len(args.prefix) :]
                        digit_str = "".join(
                            char for char in raw_suffix if char.isdigit()
                        )

                        if not digit_str:
                            continue

                        gene_num = int(digit_str)

                        start_ok = args.start is None or gene_num >= args.start
                        end_ok = args.end is None or gene_num <= args.end

                        if start_ok and end_ok:
                            header = f">{locus_tag} | {protein_id} | {product}\n"
                            dna_seq = str(feature.extract(record.seq))

                            f_dna.write(f"{header}{wrap_fasta(dna_seq)}\n")
                            f_prot.write(f"{header}{wrap_fasta(translation)}\n")

                            extracted_count += 1
                            extracted_nums.append(gene_num)

                            print(
                                f"  -> Extracted {locus_tag} ({protein_id})",
                                file=sys.stderr,
                            )

    print(
        f"\n[*] Success! {extracted_count} features extracted.",
        file=sys.stderr,
    )
    print(f"[*] DNA file written to:     {ffn_path.resolve()}", file=sys.stderr)
    print(f"[*] Protein file written to: {faa_path.resolve()}", file=sys.stderr)

    if extracted_nums and (args.start is None or args.end is None):
        actual_start = min(extracted_nums)
        actual_end = max(extracted_nums)
        print(
            f"[*] Auto-detected bounds: {args.prefix}{actual_start} → "
            f"{args.prefix}{actual_end}",
            file=sys.stderr,
        )


def run_coord_mode(args: argparse.Namespace, ffn_path: Path, faa_path: Path) -> None:
    """
    Extracts all CDS features whose genomic coordinates fall within
    the user-specified chromosomal window (--c1 to --c2).

    A CDS is included if its start position >= c1 AND its end position <= c2.
    Partial overlaps are excluded to avoid extracting truncated sequences.

    If --seq is specified, only that sequence record is searched.
    """
    print(f"[*] MODE: COORDINATE", file=sys.stderr)
    print(f"[*] Reading {args.input}...", file=sys.stderr)
    print(
        f"[*] Extracting all CDS in range {args.c1:,} – {args.c2:,} bp",
        file=sys.stderr,
    )
    if args.seq:
        print(f"[*] Targeting sequence: {args.seq}", file=sys.stderr)
    print("", file=sys.stderr)

    extracted_count = 0
    skipped_records = 0
    global_cds_counter = 1

    with open(ffn_path, "w", encoding="utf-8") as f_dna, open(
        faa_path, "w", encoding="utf-8"
    ) as f_prot:

        for record in SeqIO.parse(args.input, "genbank"):

            # If --seq is specified, skip non-matching records
            if args.seq and record.id != args.seq:
                skipped_records += 1
                continue

            for feature in record.features:
                if feature.type == "CDS":

                    locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                    product = feature.qualifiers.get(
                        "product", ["hypothetical protein"]
                    )[0]
                    translation = feature.qualifiers.get("translation", [""])[0]
                    protein_id = feature.qualifiers.get(
                        "protein_id", [f"SEQ_{global_cds_counter:05d}"]
                    )[0]
                    global_cds_counter += 1

                    if not translation:
                        continue

                    # Get chromosomal coordinates (BioPython uses 0-based start)
                    feat_start = int(feature.location.start)
                    feat_end = int(feature.location.end)

                    # Include only features fully within the coordinate window
                    if feat_start >= args.c1 and feat_end <= args.c2:
                        header = (
                            f">{locus_tag} | {protein_id} | {product} "
                            f"| {feat_start + 1}-{feat_end}\n"
                        )
                        dna_seq = str(feature.extract(record.seq))

                        f_dna.write(f"{header}{wrap_fasta(dna_seq)}\n")
                        f_prot.write(f"{header}{wrap_fasta(translation)}\n")

                        extracted_count += 1

                        print(
                            f"  -> Extracted {locus_tag} ({protein_id}) "
                            f"[{feat_start + 1}-{feat_end}]",
                            file=sys.stderr,
                        )

    print(
        f"\n[*] Success! {extracted_count} features extracted.",
        file=sys.stderr,
    )
    print(f"[*] DNA file written to:     {ffn_path.resolve()}", file=sys.stderr)
    print(f"[*] Protein file written to: {faa_path.resolve()}", file=sys.stderr)

    if skipped_records > 0:
        print(
            f"[*] Skipped {skipped_records} record(s) not matching --seq '{args.seq}'.",
            file=sys.stderr,
        )

    if extracted_count == 0:
        print(
            f"\n[!] Warning: No CDS features found in range {args.c1:,}–{args.c2:,}.",
            file=sys.stderr,
        )
        if not args.seq:
            print(
                "[!] Tip: If your gbff has multiple records (chromosome + plasmid), "
                "use --seq to target the correct one.",
                file=sys.stderr,
            )
            print(
                "[!] Example: --seq NZ_CP134351.1",
                file=sys.stderr,
            )


def main() -> None:
    """
    Entry point. Validates arguments, determines extraction mode,
    routes to the appropriate extraction function.
    """
    args = get_args()

    # Mandatory output enforcement
    if not args.output:
        sys.exit(
            "\n[!] Error: This tool requires -o to name the output files "
            "(.ffn and .faa will be created)."
        )

    # Determine which mode to run
    mode = validate_args(args)

    # Safe Path Construction
    output_path = Path(args.output)
    out_base = output_path.parent / output_path.stem
    ffn_path = out_base.with_suffix(".ffn")
    faa_path = out_base.with_suffix(".faa")

    try:
        if mode == "prefix":
            run_prefix_mode(args, ffn_path, faa_path)
        elif mode == "coord":
            run_coord_mode(args, ffn_path, faa_path)

    except FileNotFoundError:
        sys.exit(f"\n[!] File not found: {args.input}")
    except (OSError, UnicodeDecodeError, ValueError) as e:
        sys.exit(f"\n[!] File parsing error: {e}")
    except KeyboardInterrupt:
        sys.exit("\n[!] Pipeline gracefully interrupted by user.")


if __name__ == "__main__":
    main()
