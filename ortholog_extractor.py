"""
Bacteriocin Ortholog Extractor

Calculates mature peptide cores and identifies orthologs across reference genomes.

This script accepts pre-peptide sequences and calculates their mature,
membrane-inserting cores based on biochemical cleavage sites. It then
uses these core sequences to scan target reference genomes, identifying
and extracting orthologous protein sequences. If multiple targets map to
the same physical locus, it aggregates them into a single FASTA header
to prevent duplicate sequence outputs.

Example Usage:
    $ python3 ortholog_extractor.py -t target.faa -r references/ -o extracted_orthologs.faa

License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.0.1"

import sys
import argparse
from pathlib import Path
from typing import Iterator
from Bio import SeqIO
from utils import stream_reference_files, calculate_mature_core, smart_open


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Hunts for orthologs using auto-calculated core peptide homology."
    )
    parser.add_argument(
        "-t",
        "--targets",
        type=Path,
        required=True,
        help="FASTA file containing full target sequences. The script will auto-trim them.",
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=Path,
        default=Path("."),
        help="Input GenBank/FASTA file OR directory to scan (Default: current directory)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Optional: Output FASTA file path (Default: prints to terminal)",
    )
    return parser.parse_args()


def load_target_peptides(fasta_path: Path) -> dict[str, str]:
    """Loads targets and automatically runs the biochemical trimmer on them."""
    if not fasta_path.is_file():
        raise ValueError(f"Target FASTA file not found: {fasta_path}")

    targets = {}
    print(f"[*] Processing Target File: {fasta_path.name}", file=sys.stderr)

    with open(fasta_path, "r", encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_str = str(record.seq).strip()
            if seq_str:
                calculated_core = calculate_mature_core(seq_str)
                targets[record.id] = calculated_core

                print(
                    f"     -> {record.id}: Auto-trimmed to {len(calculated_core)}aa probe",
                    file=sys.stderr,
                )
                print(f"        Calculated core: {calculated_core}", file=sys.stderr)

    if not targets:
        raise ValueError(f"No valid sequences found in target file: {fasta_path.name}")

    return targets


def extract_orthologs(
    ref_path: Path, target_peptides: dict[str, str]
) -> Iterator[tuple[str, str, str, str, str]]:
    """Dynamically parses GenBank or FASTA format to yield ortholog data."""
    try:
        is_fasta = ref_path.suffix.lower() in (".fasta", ".fa", ".faa")
        fmt = "fasta" if is_fasta else "genbank"

        with open(ref_path, "r", encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, fmt):

                # BRANCH A: FASTA processing
                if is_fasta:
                    full_translation = str(record.seq).upper()
                    for target_name, core_peptide in target_peptides.items():
                        if not core_peptide:
                            continue

                        if core_peptide in full_translation:
                            full_header_desc = record.description.replace(
                                record.id, ""
                            ).strip()
                            product = (
                                full_header_desc
                                if full_header_desc
                                else "Unannotated FASTA sequence"
                            )

                            yield ref_path.stem, record.id, target_name, product, full_translation

                # BRANCH B: GenBank processing
                else:
                    for feature in record.features:
                        if feature.type == "CDS":
                            full_translation = feature.qualifiers.get(
                                "translation", [""]
                            )[0]
                            if not full_translation:
                                continue

                            for target_name, core_peptide in target_peptides.items():
                                if not core_peptide:
                                    continue

                                if core_peptide in full_translation:
                                    locus_tag = feature.qualifiers.get(
                                        "locus_tag", ["UNKNOWN"]
                                    )[0]
                                    product = feature.qualifiers.get(
                                        "product", ["Unknown product"]
                                    )[0]

                                    yield record.id, locus_tag, target_name, product, full_translation

    except Exception as e:
        raise ValueError(f"Failed to parse {ref_path.name}: {e}") from e


def main() -> None:
    args = get_args()

    try:
        targets = load_target_peptides(args.targets)

        print(f"\n[*] Scanning reference space: {args.reference}\n", file=sys.stderr)

        total_extracted_loci = 0

        with smart_open(args.output) as out_handle:
            for file_path in stream_reference_files(args.reference):
                print(f"  -> Scanning {file_path.name}...", file=sys.stderr)

                # file_hits aggregates multiple target probes hitting the same locus
                file_hits = {}

                for seq_id, locus, target_name, product, full_prot in extract_orthologs(
                    file_path, targets
                ):
                    if locus not in file_hits:
                        file_hits[locus] = {
                            "seq_id": seq_id,
                            "product": product,
                            "full_prot": full_prot,
                            "mapped_targets": [target_name],
                        }
                    else:
                        file_hits[locus]["mapped_targets"].append(target_name)

                # Process aggregated hits and write output
                for locus, data in file_hits.items():
                    total_extracted_loci += 1

                    # Clean the target names for the FASTA header
                    clean_names = [
                        t.replace(" ", "_").replace("/", "_").replace(",", "")
                        for t in data["mapped_targets"]
                    ]

                    # If it hit multiple probes, join them with a comma
                    targets_str = ",".join(clean_names)

                    fasta_header = f">{locus} | {data['seq_id']} | {data['product']} | [Ortholog_of_{targets_str}]"

                    out_handle.write(f"{fasta_header}\n{data['full_prot']}\n")
                    print(
                        f"      [Hit] {locus} contains core(s) from {len(data['mapped_targets'])} target probe(s)! ({len(data['full_prot'])} aa)",
                        file=sys.stderr,
                    )

        print("\n" + "=" * 50, file=sys.stderr)
        print(
            f"[*] SUCCESS: {total_extracted_loci} unique loci extracted.",
            file=sys.stderr,
        )

        if args.output:
            print(f"[*] Output saved to: {args.output.resolve()}", file=sys.stderr)
        print("=" * 50, file=sys.stderr)

    except ValueError as e:
        sys.exit(f"\n[!] Pipeline Error: {e}")
    except KeyboardInterrupt:
        sys.exit("\n[!] Pipeline gracefully interrupted by user.")


if __name__ == "__main__":
    main()
