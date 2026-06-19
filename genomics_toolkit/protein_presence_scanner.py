#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

"""Interactive Protein Presence/Absence Scanner

An interactive REPL tool for scanning reference genomes for the presence
or absence of a protein of interest.

Paste a full protein sequence and press Enter. The tool automatically
calculates the mature peptide core (stripping signal peptides via
biochemical cleavage logic), then performs an exact substring search
across all reference genomes. Results are reported as a PRESENT/ABSENT
matrix per genome, with locus tag and product for each hit.

Important:
    This tool uses EXACT substring matching. A single amino acid substitution
    will cause a miss. For divergent homolog detection, use pairwise_homolog_finder.py.

    Use --raw to skip bacteriocin core trimming and search with the full
    pasted sequence. Required for non-bacteriocin targets.

    SCOPE — bacterial cleavage logic only: calculate_mature_core() models
    bacterial Sec/Tat/RiPP-leader cleavage rules. It has no knowledge of
    eukaryotic secretory-pathway signal peptides (ER/Golgi-targeted),
    which use different cleavage motifs entirely. There is no reliable
    way to auto-detect this from a pasted sequence alone (it carries no
    organism metadata), so always pass --raw when the query peptide is
    not a bacterial bacteriocin/RiPP — otherwise the trimmed "core" may
    be wrong and the exact-match search will silently miss real hits.

    Only GenBank files (.gbk, .gbff) and protein FASTA files (.faa) are
    supported as references. Raw nucleotide FASTA files (.fa, .fasta) contain
    DNA strings — protein probes will never match them. These files are
    automatically detected and skipped with a warning.

    For eukaryotic genomes: use annotated GenBank files (NCBI .gbff) rather
    than raw genomic FASTA. Unannotated DNA files cannot be protein-searched
    without six-frame translation.

Note:
    Associated with ongoing, unpublished research (manuscript in
    preparation). Correct attribution is requested when used in
    derivative works.

Examples:
    # Interactive scan against a directory of reference genomes
    $ python3 protein_presence_scanner.py -i references/ -o presence_matrix.tsv

    # For non-bacteriocin proteins (skip core trimming)
    $ python3 protein_presence_scanner.py -i references/ --raw -o presence_matrix.tsv
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.2.1"

import sys
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
        suffix = gbk_path.suffix.lower()
        is_protein_fasta = suffix in (".faa", ".mpfa")
        is_nucleotide_fasta = suffix in (".fasta", ".fa", ".fna")

        # Nucleotide FASTA files contain DNA strings (ATCG...).
        # A protein probe (MKKTLV...) will never match them — silently
        # returning ABSENT for every genome. Catch this and warn explicitly.
        if is_nucleotide_fasta:
            print(
                f"  [!] Skipping {gbk_path.name}: nucleotide FASTA files cannot be "
                f"protein-searched. Provide a GenBank (.gbff/.gbk) or protein FASTA (.faa) instead.",
                file=sys.stderr,
            )
            return hits

        if is_protein_fasta:
            with open(gbk_path, "r", encoding="utf-8") as handle:
                for record in SeqIO.parse(handle, "fasta"):
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
        else:
            # GenBank format
            with open(gbk_path, "r", encoding="utf-8") as handle:
                for record in SeqIO.parse(handle, "genbank"):
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
    parser.add_argument(
        "--raw",
        action="store_true",
        help=(
            "Skip bacteriocin core trimming and use the full pasted sequence as the probe. "
            "Required for non-bacteriocin targets (housekeeping genes, TFs, kinases, etc.)."
        ),
    )
    args = parser.parse_args()

    print("=========================================", file=sys.stderr)
    print(f"  Target Scope: {args.input.name}", file=sys.stderr)
    if args.output:
        print(f"  Output File:  {args.output.resolve()}", file=sys.stderr)
    print("=========================================", file=sys.stderr)
    print("Paste a full protein sequence and press Enter.", file=sys.stderr)
    print("Type 'quit' or press Ctrl+C to exit.\n", file=sys.stderr)

    if args.raw:
        print(
            "[*] Mode: --raw (full sequences used, no core trimming)", file=sys.stderr
        )
    else:
        print(
            "[*] Mode: bacteriocin core trimming (bacterial Sec/Tat/RiPP "
            "cleavage rules; NOT valid for eukaryotic secretome proteins "
            "— use --raw for those)",
            file=sys.stderr,
        )

    # Track whether this is the first successful query in this session.
    # We open the output file in "w" mode for the first query (creating/overwriting)
    # and "a" mode for all subsequent queries (appending rows to same file).
    # This prevents the session eraser bug where "w" inside the loop would
    # wipe the previous query's results each time a new sequence is pasted.
    first_query = True

    while True:
        try:
            user_input = input("Paste protein sequence: ").strip()

            if user_input.lower() in ["quit", "exit", "q"]:
                break
            if not user_input:
                continue

            if args.raw:
                core_target = user_input.upper()
                print(
                    f"\n  [+] Using full sequence: {len(core_target)}aa\n",
                    file=sys.stderr,
                )
            else:
                print("\n  [*] Calculating structural core...", file=sys.stderr)
                core_target = calculate_mature_core(user_input.upper())

                if not core_target:
                    print(
                        "  [!] Error: Mature core calculation returned an empty sequence.\n"
                        "      This can happen if the signal peptide spans the entire protein,\n"
                        "      or if 'GG' appears at the very end of the sequence.\n"
                        "      Try a longer sequence, check the input, or use --raw.",
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
                file_mode = "w" if first_query else "a"
                with open(args.output, file_mode, encoding="utf-8-sig") as tsv:
                    if first_query:
                        # Write header only on the first query of this session
                        tsv.write(
                            "Query_Core\tGenome_File\tLocus_Tag\tProduct\tStatus\n"
                        )
                        first_query = False
                    for file_path in stream_reference_files(args.input):
                        print(f"  [*] Scanning {file_path.name}...", file=sys.stderr)
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
                                    f"{core_target}\t{file_path.name}\t{locus}\t{product}\tPRESENT\n"
                                )
                        else:
                            print(
                                f"      [✓] ABSENT in {file_path.name}", file=sys.stderr
                            )
                            tsv.write(
                                f"{core_target}\t{file_path.name}\t-\t-\tABSENT\n"
                            )
            else:
                for file_path in stream_reference_files(args.input):
                    print(f"  [*] Scanning {file_path.name}...", file=sys.stderr)
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
            print("\n[!] Force quitting. Goodbye!", file=sys.stderr)
            break
        except ValueError as e:
            print(f"\n  [X] Error: {e}\n", file=sys.stderr)


if __name__ == "__main__":
    main()
