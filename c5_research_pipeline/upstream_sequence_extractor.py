"""
Upstream Sequence Extractor

Extracts the raw upstream DNA sequence (promoter region) of a target gene
from a GenBank file using its locus tag.

This tool is used to empirically verify structural conservation in the
regulatory regions flanking orthologous genes across different genomes.

Author: Jan Ephraim R. Vallente (ephrvallente@gmail.com)
Date: 2026-06-04
License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).

Example Usage:
    # Extract 150bp upstream of LEUM_RS10400
    $ python3 upstream_sequence_extractor.py -i GCF_000014445.1_genomic.gbff -l LEUM_RS10400 -u 150 -o upstream_LEUM.fasta
"""

__version__ = "1.0.2"

import sys
import argparse
from pathlib import Path
from utils import extract_upstream_sequence, wrap_fasta


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extracts raw upstream promoter sequences."
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the GenBank file (.gbk or .gbff).",
    )
    parser.add_argument(
        "-l",
        "--locus",
        type=str,
        required=True,
        help="The exact locus tag of the gene (e.g., LEUM_RS10400).",
    )
    parser.add_argument(
        "-u",
        "--upstream",
        type=int,
        default=150,
        help="Number of base pairs to extract upstream of the start codon. Default: 150.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=False,
        help="Optional: Output FASTA file. If omitted, prints to terminal.",
    )
    return parser.parse_args()


def main() -> None:
    """Main execution block for upstream sequence extraction."""
    args = get_args()

    if args.upstream < 1:
        sys.exit("[!] --upstream must be a positive integer.")

    try:
        # Utilize the centralized engine from utils.py
        sequence, start, end, strand = extract_upstream_sequence(
            args.input, args.target, args.upstream
        )

        strand_symbol = "+" if strand == 1 else "-"

        # Adjusted to 1-based closed interval for biological standardization
        fasta_header = (
            f">{args.input.name} | {args.target} | "
            f"Upstream_{args.upstream}bp | Coord:{start + 1}-{end}({strand_symbol})"
        )

        wrapped_seq = wrap_fasta(sequence)
        final_output = f"{fasta_header}\n{wrapped_seq}\n"

        if args.output:
            with open(args.output, "w", encoding="utf-8") as f:
                f.write(final_output)
            print(
                f"[*] Success! Extracted {args.upstream}bp upstream of {args.target}.",
                file=sys.stderr,
            )
            print(f"[*] Saved to {args.output.name}", file=sys.stderr)
        else:
            # Print sequence to stdout for clean piping, logs go to stderr
            sys.stdout.write(final_output)

    except FileNotFoundError:
        sys.exit(f"\n[!] File not found: {args.input}")
    except ValueError as e:
        sys.exit(f"\n[!] Extraction Error: {e}")
    except KeyboardInterrupt:
        sys.exit("\n[!] Extraction interrupted by user.")


if __name__ == "__main__":
    main()
