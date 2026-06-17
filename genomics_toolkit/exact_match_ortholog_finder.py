"""
Ortholog Extractor

Calculates mature peptide cores and identifies orthologs across reference genomes
using EXACT SUBSTRING MATCHING (100% identity required).

This script accepts pre-peptide sequences and calculates their mature,
membrane-inserting cores based on biochemical cleavage sites. It then
uses these core sequences to scan target reference genomes, identifying
and extracting full orthologous protein sequences. If multiple targets map to
the same physical locus, it aggregates them into a single FASTA header
to prevent duplicate sequence outputs.

Important:
    This tool uses EXACT substring matching. A single amino acid substitution
    in the core region will cause a miss. For divergent homolog detection at
    lower identity (e.g., ≥35%), use gbk_ortholog_finder.py (Smith-Waterman
    alignment with BLOSUM62 scoring) instead.

    Use --raw to skip bacteriocin-specific core trimming and search using
    the full provided sequence. Required for non-bacteriocin targets.

License: MIT

Reproducibility:
    Associated with upcoming research (manuscript in preparation).
    Correct attribution is requested when used in derivative works.
    See LICENSE in the repository root for full details.

Example Usage:
    $ python3 ortholog_extractor.py -t target.faa -r references/ -o extracted_orthologs.faa

    # For non-bacteriocin proteins (housekeeping genes, TFs, kinases, etc.):
    $ python3 ortholog_extractor.py -t target.faa -r references/ --raw -o extracted_orthologs.faa
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.1.0"

import sys
import argparse
from pathlib import Path
from typing import Iterator

try:
    from Bio import SeqIO
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )
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
    parser.add_argument(
        "--raw",
        action="store_true",
        help=(
            "Skip calculate_mature_core() and use the full sequences exactly as provided. "
            "Required for non-bacteriocin targets (transcription factors, housekeeping genes, "
            "eukaryotic proteins, etc.) that don't use double-glycine cleavage logic."
        ),
    )
    return parser.parse_args()


def load_target_peptides(fasta_path: Path, use_raw: bool = False) -> dict[str, str]:
    """Loads targets and optionally runs the biochemical trimmer on them.

    Args:
        fasta_path: Path to the FASTA file containing target sequences.
        use_raw:    If True, skips calculate_mature_core() and uses the
                    full sequence as the search probe. Use for non-bacteriocin
                    targets that don't have double-glycine cleavage sites.
    """
    if not fasta_path.is_file():
        raise ValueError(f"Target FASTA file not found: {fasta_path}")

    targets = {}
    skipped = 0
    print(f"[*] Processing Target File: {fasta_path.name}", file=sys.stderr)
    if use_raw:
        print(
            "[*] Mode: --raw (full sequences used, no core trimming)", file=sys.stderr
        )
    else:
        print(
            "[*] Mode: bacteriocin core trimming (use --raw to disable)",
            file=sys.stderr,
        )

    with open(fasta_path, "r", encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_str = str(record.seq).strip()
            if not seq_str:
                continue

            if use_raw:
                probe = seq_str
                print(
                    f"     -> {record.id}: Using full {len(probe)}aa sequence as probe",
                    file=sys.stderr,
                )
            else:
                probe = calculate_mature_core(seq_str)

                if not probe:
                    # Empty core means calculate_mature_core found "GG" at the very
                    # end of the sequence (nothing left after the cleavage site).
                    # Rather than silently skipping, warn the user explicitly.
                    print(
                        f"     -> [!] WARNING: {record.id} — core trimming returned an "
                        f"empty sequence. This sequence will be SKIPPED.\n"
                        f"        Tip: If this is not a bacteriocin, re-run with --raw.",
                        file=sys.stderr,
                    )
                    skipped += 1
                    continue

                print(
                    f"     -> {record.id}: Trimmed to {len(probe)}aa probe",
                    file=sys.stderr,
                )
                print(f"        Core: {probe}", file=sys.stderr)

            targets[record.id] = probe

    if skipped > 0:
        print(
            f"\n[!] {skipped} sequence(s) skipped due to empty core. "
            f"Use --raw to include full sequences.",
            file=sys.stderr,
        )

    if not targets:
        raise ValueError(
            f"No valid probe sequences found in target file: {fasta_path.name}"
        )

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
        targets = load_target_peptides(args.targets, use_raw=args.raw)

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

                    # Truncate gracefully if many probes hit the same locus —
                    # oversized headers break some downstream tools (HMMER, aligners)
                    if len(clean_names) > 3:
                        targets_str = (
                            ",".join(clean_names[:3])
                            + f",_and_{len(clean_names) - 3}_more"
                        )
                    else:
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
