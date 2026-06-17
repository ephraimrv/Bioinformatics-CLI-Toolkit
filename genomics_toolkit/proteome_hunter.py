"""
GenBank Promoter Region Extractor

Scans GenBank assemblies for genes matching keywords and extracts upstream promoter sequences.

This tool performs keyword-based gene identification within GenBank assemblies. It
extracts a specified upstream base pair range for identified genes, performs
deduplication to ensure statistical validity for downstream motif analysis,
and exports the result to a FASTA file.

License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).

Example Usage:
    $ python3 universal_promoter_extractor.py -i C5_prokka.gbk \
      -o C5_promoters.fasta -u 150 -k bacteriocin lactobin cerein

    $ python3 universal_promoter_extractor.py -i references/ \
      -o C5_promoters.fasta -u 150 -k bacteriocin lactobin cerein
"""


__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.0.1"

import argparse
from pathlib import Path
from Bio import SeqIO
from utils import stream_reference_files, calculate_mature_core


def scan_for_peptide(gbk_path: Path, target_peptide: str) -> list[tuple[str, str]]:
    hits = []
    try:
        is_fasta = gbk_path.suffix.lower() in (".fasta", ".fa", ".faa")
        fmt = "fasta" if is_fasta else "genbank"

        with open(gbk_path, "r", encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, fmt):
                # Branch A: FASTA Reference
                if is_fasta:
                    translation = str(record.seq).upper()
                    if target_peptide in translation:
                        full_header_desc = record.description.replace(
                            record.id, ""
                        ).strip()
                        product = (
                            full_header_desc
                            if full_header_desc
                            else "Unannotated FASTA sequence"
                        )
                        hits.append((record.id, product))

                # Branch B: GenBank Reference
                else:
                    for feature in record.features:
                        if feature.type == "CDS":
                            translation = feature.qualifiers.get("translation", [""])[
                                0
                            ].upper()
                            if target_peptide in translation:
                                locus_tag = feature.qualifiers.get(
                                    "locus_tag", ["UNKNOWN"]
                                )[0]
                                product = feature.qualifiers.get(
                                    "product", ["Unknown product"]
                                )[0]
                                hits.append((locus_tag, product))
        return hits

    except (FileNotFoundError, OSError, UnicodeDecodeError, ValueError) as e:
        raise ValueError(f"Failed to parse or extract from {gbk_path.name}: {e}") from e


def main() -> None:
    parser = argparse.ArgumentParser(description="Bacteriocin Presence/Absence Hunter")
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
        required=False,
        help="Output TSV file for the presence/absence matrix",
    )
    args = parser.parse_args()

    print("=========================================")
    print(f"  Target Scope: {args.input.name}")
    if args.output:
        print(f"  Output File:  {args.output.resolve()}")
    print("=========================================")
    print("Paste a full protein sequence and press Enter.")
    print("Type 'quit' or press Ctrl+C to exit.\n")

    while True:
        try:
            user_input = input("Paste protein sequence: ").strip()

            if user_input.lower() in ["quit", "exit", "q"]:
                break
            if not user_input:
                continue

            print("\n  [*] Calculating structural core...")
            core_target = calculate_mature_core(user_input.upper())

            print(f"  [+] Core Extracted : {core_target}")
            print(f"  [i] Core Length    : {len(core_target)} amino acids\n")

            total_input_hits = 0

            if args.output:
                with open(args.output, "w", encoding="utf-8-sig") as tsv:
                    tsv.write("Genome_File\tLocus_Tag\tProduct\tStatus\n")
                    for file_path in stream_reference_files(args.input):
                        print(f"  [*] Scanning {file_path.name}...")
                        hits = scan_for_peptide(file_path, core_target)

                        if hits:
                            print(
                                f"      [!] ALERT: Found {len(hits)} match(es) in {file_path.name}"
                            )
                            total_input_hits += len(hits)
                            for locus, product in hits:
                                short_prod = (
                                    product[:45] + "..."
                                    if len(product) > 45
                                    else product
                                )
                                print(
                                    f"          -> Locus: {locus:<15} | Product: {short_prod}"
                                )
                                tsv.write(
                                    f"{file_path.name}\t{locus}\t{product}\tPRESENT\n"
                                )
                        else:
                            print(f"      [✓] ABSENT in {file_path.name}")
                            tsv.write(f"{file_path.name}\t-\t-\tABSENT\n")
            else:
                for file_path in stream_reference_files(args.input):
                    print(f"  [*] Scanning {file_path.name}...")
                    hits = scan_for_peptide(file_path, core_target)

                    if hits:
                        print(
                            f"      [!] ALERT: Found {len(hits)} match(es) in {file_path.name}"
                        )
                        total_input_hits += len(hits)
                        for locus, product in hits:
                            short_prod = (
                                product[:45] + "..." if len(product) > 45 else product
                            )
                            print(
                                f"          -> Locus: {locus:<15} | Product: {short_prod}"
                            )
                    else:
                        print(f"      [✓] ABSENT in {file_path.name}")
            if args.output:
                print(
                    f"\n  [=] BATCH SUMMARY: {total_input_hits} matches found. Matrix saved to {args.output.name}.\n"
                )
            print("-" * 60)

        except KeyboardInterrupt:
            print("\nForce quitting tool. Goodbye!")
            break
        except ValueError as e:
            print(f"\n  [X] Error: {e}\n")


if __name__ == "__main__":
    main()
