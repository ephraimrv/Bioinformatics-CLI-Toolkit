"""
DNA to RNA Transcription Pipeline

Reads multi-FASTA DNA datasets using a lazy generator, transcribes the coding
strand into mRNA using C-optimized translation, and streams the output to a TSV.
"""

import sys
from utils import lazy_parse_fasta
from utils import base_parser

_RNA_MAP = str.maketrans("Tt", "Uu")


def transcribe_to_rna(dna: str) -> str:
    """
    Transcribes a DNA string into RNA.

    Biological Context:
    Transcription is the first step of gene expression, where a segment of
    DNA is copied into RNA. Computationally, this is modeled by a direct
    character substitution of 'T' to 'U'.

    Args:
        dna_sequence: A raw string representing the coding DNA strand.

    Returns:
        The transcribed RNA sequence string.
    """

    unique_characters = set(dna)
    allowed_characters = set("ACGNTacgnt")

    if not unique_characters.issubset(allowed_characters):
        rogue_bases = unique_characters - allowed_characters
        raise ValueError(f"Invalid bases found: {rogue_bases}")

    return dna.translate(_RNA_MAP)


def main():
    """
    Pipeline manager for I/O operations.

    Reads raw DNA sequences from an input file, computes their reverse
    complements, and writes the aggregated results to an output file.
    """
    args = base_parser("DNA to RNA Transcription Pipeline").parse_args()
    input_path = args.input
    output_path = args.output

    sequences_processed = 0
    try:
        with open(output_path, "w", encoding="utf-8") as out_file:
            out_file.write("Sequence_ID\tTranscribed_RNA\n")
            for fasta_id, sequence in lazy_parse_fasta(input_path):
                sequences_processed += 1

                rna = transcribe_to_rna(sequence)
                print(f"> {fasta_id} transcribed to RNA!")

                row = f"{fasta_id}\t{rna}\n"
                out_file.write(row)

        if sequences_processed == 0:
            raise ValueError("The input file contained no valid FASTA records.")

    # EAFP: Catching all our errors gracefully here at the execution layer.
    except (ValueError, FileNotFoundError, PermissionError) as e:
        output_path.unlink(missing_ok=True)
        sys.exit(f"Pipeline Halted: {e}")

    print(f"\nSuccess! {sequences_processed} sequences processed with O(1) memory.")
    print(f"Machine-readable matrix written to: {output_path.name}")


if __name__ == "__main__":
    main()
