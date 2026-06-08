"""
Genomic CDS Region Extractor

Parses GenBank assemblies to extract and export coding sequences (CDS).

This tool iterates through GenBank (.gbk/.gbff) files to extract specific
genomic regions defined by locus tag prefixes. It handles alphanumeric
tag normalization (e.g., stripping NCBI 'RS' prefixes), performs coordinate
filtering, and writes both nucleotide and translated protein sequences
into synchronized FASTA files.

Author: Jan Ephraim R. Vallente (ephrvallente@gmail.com)
Date: 2026-06-07
License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).

Example Usage:
    # For a local Prokka assembly:
    $ python3 extract_gbk_region.py -i input.gbk -p "ctg1_" -s 46 -e 74 -o C5_locus

    # For an NCBI RefSeq assembly:
    $ python3 extract_gbk_region.py -i GCF_056530425.1.gbff -p "RHP56_" -s 5 -e 20 -o NCBI_locus
"""

__version__ = "1.0.2"

import sys
import argparse
from pathlib import Path
from Bio import SeqIO
from utils import base_parser, wrap_fasta


def get_args() -> argparse.Namespace:
    """Configures the CLI and returns parsed arguments."""
    parser = base_parser(
        description_text="Extracts bifunctional (DNA/Protein) CDS features from a GenBank file.",
    )

    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        required=True,
        help="The locus tag prefix (e.g., 'ctg1_' or 'RHP56_')",
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
    """
    Validates user bounds, ensures output destination exists, routes IO streams
    via context managers, and parses the GenBank file to extract, wrap, and
    output targeted bifunctional CDS features.
    """
    args = get_args()

    # Mandatory output enforcement
    if not args.output:
        sys.exit(
            "\n[!] Error: This tool requires -o to name the output files (.ffn and .faa will be created)."
        )

    # Strict Input Validation
    if args.start is not None and args.start < 1:
        sys.exit("\n[!] Error: --start must be a positive integer.")
    if args.end is not None and args.end < 1:
        sys.exit("\n[!] Error: --end must be a positive integer.")
    if args.start is not None and args.end is not None and args.start > args.end:
        sys.exit(
            f"\n[!] Error: Start gene ({args.start}) cannot be greater than End gene ({args.end})."
        )

    start_label = args.start if args.start is not None else "[First Available]"
    end_label = args.end if args.end is not None else "[Last Available]"

    print(f"[*] Reading {args.input}...", file=sys.stderr)
    print(
        f"[*] Extracting {args.prefix}{start_label} to {args.prefix}{end_label}...\n",
        file=sys.stderr,
    )

    extracted_count = 0
    extracted_nums = []
    global_cds_counter = 1

    # Safe Path Construction
    output_path = Path(args.output)
    out_base = output_path.parent / output_path.stem
    ffn_path = out_base.with_suffix(".ffn")
    faa_path = out_base.with_suffix(".faa")

    try:
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
                f"\n[*] Success! {extracted_count} dual-format features extracted.",
                file=sys.stderr,
            )
            print(f"[*] DNA file written to:     {ffn_path.resolve()}", file=sys.stderr)
            print(f"[*] Protein file written to: {faa_path.resolve()}", file=sys.stderr)

            if extracted_nums and (args.start is None or args.end is None):
                actual_start = min(extracted_nums)
                actual_end = max(extracted_nums)
                print(
                    f"[*] Automatically detected bounds: {args.prefix}{actual_start} to {args.prefix}{actual_end}",
                    file=sys.stderr,
                )

    except FileNotFoundError:
        sys.exit(f"\n[!] File not found: {args.input}")
    except (OSError, UnicodeDecodeError, ValueError) as e:
        sys.exit(f"\n[!] File parsing error: {e}")
    except KeyboardInterrupt:
        sys.exit("\n[!] Pipeline gracefully interrupted by user.")


if __name__ == "__main__":
    main()
