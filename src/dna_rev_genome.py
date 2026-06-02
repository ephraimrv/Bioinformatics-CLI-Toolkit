"""
Reverse Complement Pipeline

Reads multiple sequences from a FASTA file using a lazy generator, calculates
the reverse complement (3'-to-5') of each sequence using C-optimized translation,
and streams the new sequences to an output TSV file.

Usage:
    $ python3 dna_rev_genome.py -i raw_dna.fasta -o reverse_complements.tsv
"""

__author__ = 'Jan Ephraim R. Vallente'
__email__ = 'ephrvallente@gmail.com'
__version__ = '1.0.0'

import sys
from utils import lazy_parse_fasta
from utils import base_parser

_DNA_MAP = str.maketrans("ACGNTacgnt", "TGCNAtgcna")


def reverse_complement(dna: str) -> str:
    """
    Calculate the reverse complement of a DNA string using a dictionary map.

    Biological Context:
    Simulates the anti-parallel nature of double-stranded DNA by generating
    the 3'-to-5' complementary strand from a 5'-to-3' sequence.

    Architecture & Performance:
    Executes a rapid string reversal (`[::-1]`) followed by a C-optimized 
    `str.translate()` character swap. 

    Args:
        dna: A DNA sequence containing only standard nucleotides.

    Returns:
        The sequence of the reverse complement strand.

    Raises:
        ValueError: If the string contains any anomalous characters outside the 
            standard IUPAC base set.
    """

    unique_characters = set(dna)
    allowed_characters = set("ACGNTacgnt")

    if not unique_characters.issubset(allowed_characters):
        # We find exactly which characters are rogue and report them
        rogue_bases = unique_characters - allowed_characters
        raise ValueError(f"Invalid bases found: {rogue_bases}")

    return dna[::-1].translate(_DNA_MAP)


def main():
    """
    Pipeline manager for reverse complement I/O operations.

    Reads raw DNA sequences, orchestrates the reversal engine, 
    and safely streams the results to an output TSV file.
    """
    args = base_parser("Genome Reverse Complement").parse_args()
    input_path = args.input
    output_path = args.output

    sequences_processed = 0
    try:
        with open(output_path, "w", encoding="utf-8") as out_file:
            out_file.write("Sequence ID\tDNA Reverse Complement\n")
            for fasta_id, sequence in lazy_parse_fasta(input_path):
                sequences_processed += 1

                dna_rev = reverse_complement(sequence)
                print(f"> {fasta_id} reverse complement generated!")

                row = f"{fasta_id}\t{dna_rev}\n"
                out_file.write(row)

        if sequences_processed == 0:
            raise ValueError(
                "Pipeline Halted: The input file contained no valid FASTA records."
            )

    # EAFP: Catching all our errors gracefully here at the execution layer.
    except (ValueError, FileNotFoundError, PermissionError) as e:
        output_path.unlink(missing_ok=True)
        sys.exit(f"Pipeline Halted: {e}")

    print(f"\nSuccess! {sequences_processed} sequences processed with O(1) memory.")
    print(f"Machine-readable matrix written to: {output_path.name}")


if __name__ == "__main__":
    main()
