"""
GenBank Upstream Promoter Finder

Extracts regulatory promoter regions from multi-contig GenBank assemblies.

This tool locates a target gene by its locus tag, calculates the strand
orientation to extract the correct upstream sequence (applying reverse
complementation where necessary), and scans the region for specific motif hits.

Author: Jan Ephraim R. Vallente (ephrvallente@gmail.com)
Date: 2026-06-07
License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).

Example Usage:
    $ python3 gbk_promoter_finder.py -i C5_genome.gbk -l ctg1_50 -u 150 -m "TATAAT" -o ctg1_50_promoter.fasta
"""

__version__ = "1.0.3"

import re
import sys
import traceback
from typing import Iterator
from utils import base_parser, wrap_fasta, extract_upstream_sequence


def find_motif_regex_iterator(
    sequence: str, regex_pattern: str
) -> Iterator[tuple[int, str]]:
    if not sequence or not regex_pattern:
        return
    try:
        safe_pattern = re.compile(f"(?=({regex_pattern}))", re.IGNORECASE)
    except re.error as e:
        raise ValueError(f"Invalid regex pattern provided: {regex_pattern}") from e

    for match in safe_pattern.finditer(sequence):
        yield match.start() + 1, match.group(1)


def main() -> None:
    parser = base_parser("GenBank Targeted Upstream Motif Scanner")
    parser.add_argument("-l", "--locus", required=True, help="Target gene Locus Tag")
    parser.add_argument(
        "-u", "--upstream", type=int, default=100, help="Upstream bases to extract"
    )
    parser.add_argument(
        "-m", "--motif", required=False, help="Regex motif to search (optional)"
    )
    args = parser.parse_args()

    try:
        if args.upstream < 1:
            raise ValueError("Upstream bases must be a positive integer.")

        upstream_seq, start, end, strand = extract_upstream_sequence(
            args.input, args.locus, args.upstream
        )

        print(f"[*] Found {args.locus} at {start}-{end} (Strand: {strand})")
        print(f"[*] Extracting {args.upstream}bp upstream...")

        motifs = []
        if args.motif:
            print(f"[*] Searching for motif: {args.motif}")
            motifs = list(find_motif_regex_iterator(upstream_seq, args.motif))

        if args.output:
            with open(args.output, "w", encoding="utf-8") as out_file:
                # Strictly formatted FASTA Header
                out_file.write(
                    f">{args.locus}_upstream_{args.upstream}bp_strand_{strand}\n"
                )
                out_file.write(f"{wrap_fasta(upstream_seq)}\n")

                # Machine-readable TSV motif appendix
                if motifs:
                    out_file.write("\n# Motif_Hit_Position\tSequence\n")
                    for pos, seq in motifs:
                        out_file.write(f"{pos}\t{seq}\n")

            print(
                f"[*] Success! {len(motifs)} motifs found. Data written to {args.output.name}"
            )

        else:
            print("\n--- UPSTREAM SEQUENCE ---")
            if len(upstream_seq) > 500:
                print(
                    f"{upstream_seq[:100]} ... [snip {len(upstream_seq)-200}bp] ... {upstream_seq[-100:]}"
                )
            else:
                print(upstream_seq)
            print("-------------------------\n")

            if motifs:
                for pos, seq in motifs:
                    print(
                        f"    -> Motif Found! Relative Position: {pos} | Sequence: {seq}"
                    )
            elif args.motif is None:
                pass
            else:
                print("    -> No motifs found.")

    except FileNotFoundError:
        sys.exit(f"\n[!] File not found: {args.input}")
    except ValueError as e:
        sys.exit(f"\n[!] Error: {e}")
    except KeyboardInterrupt:
        sys.exit("\n[!] Pipeline gracefully interrupted by user.")
    except Exception:
        print("\n[!] UNEXPECTED BUG ENCOUNTERED:")
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
