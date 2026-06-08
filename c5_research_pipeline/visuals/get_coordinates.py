"""
Genomic Coordinate Extractor

Extracts absolute start, end, and strand coordinates for a specific range
of locus tags from a GenBank file. Outputs a structured TSV file designed
for seamless ingestion by downstream visualization scripts.

License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).

Example usage:
    $ python3 get_coordinates.py -i input.gbk -o coords.tsv -p "ctg1_" -s 46 -e 74
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.0.1"

import sys
import argparse
import contextlib
from Bio import SeqIO
from utils import base_parser


def setup_cli() -> argparse.Namespace:
    parser = base_parser(
        description_text="Extracts genomic coordinates of CDS features to a TSV matrix.",
    )

    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        required=True,
        help="The locus tag prefix (e.g., 'ctg1_')",
    )

    parser.add_argument(
        "-s",
        "--start",
        type=int,
        required=False,
        default=None,
        help="Start gene number (Default: First available)",
    )
    parser.add_argument(
        "-e",
        "--end",
        type=int,
        required=False,
        default=None,
        help="End gene number (Default: Last available)",
    )
    return parser.parse_args()


def main() -> None:
    args = setup_cli()

    effective_start = args.start if args.start is not None else 0
    effective_end = args.end if args.end is not None else float("inf")

    extracted_count = 0
    extracted_nums = []

    if args.output:
        start_label = args.start if args.start is not None else "[First Available]"
        end_label = args.end if args.end is not None else "[Last Available]"
        print(f"[*] Reading {args.input.name}...", file=sys.stderr)
        print(
            f"[*] Extracting coordinates for {args.prefix}{start_label} to {args.prefix}{end_label}...\n",
            file=sys.stderr,
        )

    # Safely handle file contexts; fallback to stdout if no output specified
    ctx = (
        open(args.output, "w", encoding="utf-8")
        if args.output
        else contextlib.nullcontext(sys.stdout)
    )

    with ctx as out_handle:
        # Write TSV Header
        out_handle.write("locus_tag\tstart\tend\tstrand\tproduct\n")

        try:
            for record in SeqIO.parse(args.input, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":
                        locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                        product = feature.qualifiers.get("product", ["Unknown"])[0]

                        if locus_tag.startswith(args.prefix):
                            # Robust handling for NCBI RS-injected prefixes (e.g., ctg1_RS00046)
                            digit_str = "".join(
                                c for c in locus_tag[len(args.prefix) :] if c.isdigit()
                            )
                            if not digit_str:
                                continue

                            gene_num = int(digit_str)

                            if effective_start <= gene_num <= effective_end:
                                start_pos = int(feature.location.start)
                                end_pos = int(feature.location.end)
                                strand = int(feature.location.strand)

                                out_handle.write(
                                    f"{locus_tag}\t{start_pos}\t{end_pos}\t{strand}\t{product}\n"
                                )

                                extracted_count += 1
                                extracted_nums.append(gene_num)

        except FileNotFoundError:
            sys.exit(f"\n[!] File not found: {args.input}")
        except ValueError as e:
            sys.exit(f"\n[!] GenBank parsing error: {e}")

    if args.output:
        print(
            f"[✓] Success! {extracted_count} coordinate features saved to {args.output.name}",
            file=sys.stderr,
        )
        if extracted_nums and (args.start is None or args.end is None):
            actual_start = min(extracted_nums)
            actual_end = max(extracted_nums)
            print(
                f"[*] Automatically detected bounds: {args.prefix}{actual_start} to {args.prefix}{actual_end}",
                file=sys.stderr,
            )


if __name__ == "__main__":
    main()
