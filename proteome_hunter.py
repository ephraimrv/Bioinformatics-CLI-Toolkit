"""
Proteome Hunter — Interactive Protein Presence/Absence Scanner

An interactive REPL tool for scanning reference genomes for the presence
or absence of a protein of interest.

Paste a full protein sequence and press Enter. The tool automatically
calculates the mature peptide core (stripping signal peptides via
biochemical cleavage logic), then performs an exact substring search
across all reference genomes. Results are reported as a PRESENT/ABSENT
matrix per genome, with locus tag and product for each hit.

Useful for rapid cross-genome exploration during early-stage research,
before committing to a full alignment-based search with gbk_ortholog_finder.py.

Note:
    This tool uses exact substring matching. For homolog detection at
    lower identity, use gbk_ortholog_finder.py instead.

License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).

Example Usage:
    # Interactive scan against a directory of reference genomes
    $ python3 proteome_hunter.py -i references/ -o presence_matrix.tsv

    # Interactive scan against a single genome
    $ python3 proteome_hunter.py -i genome.gbff
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.1.0"

import argparse
from pathlib import Path

try:
    from Bio import SeqIO
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )
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

    print("=========================================", file=sys.stderr)
    print(f"  Target Scope: {args.input.name}", file=sys.stderr)
    if args.output:
        print(f"  Output File:  {args.output.resolve()}", file=sys.stderr)
    print("=========================================", file=sys.stderr)
    print("Paste a full protein sequence and press Enter.", file=sys.stderr)
    print("Type 'quit' or press Ctrl+C to exit.\n", file=sys.stderr)

    while True:
        try:
            user_input = input("Paste protein sequence: ").strip()

            if user_input.lower() in ["quit", "exit", "q"]:
                break
            if not user_input:
                continue

            print("\n  [*] Calculating structural core...", file=sys.stderr)
            core_target = calculate_mature_core(user_input.upper())

            # Guard: if calculate_mature_core returns an empty string
            # (signal peptide longer than the whole protein), an empty string
            # matches EVERY sequence — resulting in false positives across all genomes.
            if not core_target:
                print(
                    "  [!] Error: Mature core calculation returned an empty sequence.\n"
                    "      This can happen if the signal peptide spans the entire protein.\n"
                    "      Try a longer sequence or check the input.",
                    file=sys.stderr,
                )
                continue

            print(f"  [+] Core Extracted : {core_target}", file=sys.stderr)
            print(
                f"  [i] Core Length    : {len(core_target)} amino acids\n",
                file=sys.stderr,
            )

            total_input_hits = 0

            if args.output:
                with open(args.output, "w", encoding="utf-8-sig") as tsv:
                    tsv.write("Genome_File\tLocus_Tag\tProduct\tStatus\n")
                    for file_path in stream_reference_files(args.input):
                        print(f"  [*] Scanning {file_path.name}...")
                        hits = scan_for_peptide(file_path, core_target)

                        if hits:
                            print(
                                f"      [!] ALERT: Found {len(hits)} match(es) in {file_path.name}",
                                file=sys.stderr,
                            )
                            total_input_hits += len(hits)
                            for locus, product in hits:
                                short_prod = (
                                    product[:45] + "..."
                                    if len(product) > 45
                                    else product
                                )
                                print(
                                    f"          -> Locus: {locus:<15} | Product: {short_prod}",
                                    file=sys.stderr,
                                )
                                tsv.write(
                                    f"{file_path.name}\t{locus}\t{product}\tPRESENT\n"
                                )
                        else:
                            print(
                                f"      [✓] ABSENT in {file_path.name}", file=sys.stderr
                            )
                            tsv.write(f"{file_path.name}\t-\t-\tABSENT\n")
            else:
                for file_path in stream_reference_files(args.input):
                    print(f"  [*] Scanning {file_path.name}...")
                    hits = scan_for_peptide(file_path, core_target)

                    if hits:
                        print(
                            f"      [!] ALERT: Found {len(hits)} match(es) in {file_path.name}",
                            file=sys.stderr,
                        )
                        total_input_hits += len(hits)
                        for locus, product in hits:
                            short_prod = (
                                product[:45] + "..." if len(product) > 45 else product
                            )
                            print(
                                f"          -> Locus: {locus:<15} | Product: {short_prod}",
                                file=sys.stderr,
                            )
                    else:
                        print(f"      [✓] ABSENT in {file_path.name}", file=sys.stderr)
            if args.output:
                print(
                    f"\n  [=] BATCH SUMMARY: {total_input_hits} matches found. "
                    f"Matrix saved to {args.output.name}.\n",
                    file=sys.stderr,
                )
            print("-" * 60, file=sys.stderr)

        except KeyboardInterrupt:
            print("\nForce quitting tool. Goodbye!")
            break
        except ValueError as e:
            print(f"\n  [X] Error: {e}\n")


if __name__ == "__main__":
    main()
