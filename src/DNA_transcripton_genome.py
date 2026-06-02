"""
DNA to RNA Transcription Pipeline

Reads multi-FASTA DNA datasets using a lazy generator, transcribes the coding
strand into mRNA using C-optimized translation, and streams the output to a TSV.

Usage:
    $ python3 DNA_transcripton_genome.py -i coding_dna.fasta -o transcribed_rna.tsv
"""

__author__ = 'Jan Ephraim R. Vallente'
__email__ = 'ephrvallente@gmail.com'
__version__ = '1.0.1'

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

    Architecture & Performance:
    Utilizes Python's native `str.translate()` backed by a C-optimized 
    translation map for instantaneous character substitution.

    Args:
        dna: A raw string representing the coding DNA strand.

    Returns:
        The transcribed RNA sequence.

    Raises:
        ValueError: If the input sequence contains any invalid characters.
    """

    unique_characters = set(dna)
    allowed_characters = set("ACGNTacgnt")

    if not unique_characters.issubset(allowed_characters):
        rogue_bases = unique_characters - allowed_characters
        raise ValueError(f"Invalid bases found: {rogue_bases}")

    return dna.translate(_RNA_MAP)


def main() -> None:
    """
    Pipeline manager for transcription I/O operations.

    Reads raw DNA sequences, orchestrates the transcription engine, 
    and writes the aggregated results to an output TSV file.
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

    except KeyboardInterrupt:
        output_path.unlink(missing_ok=True)
        sys.exit(
            "\nPipeline Halted: Scan interrupted by user. Partial output safely removed."
        )
    except (ValueError, FileNotFoundError, PermissionError) as e:
        output_path.unlink(missing_ok=True)
        sys.exit(f"\nSystem Exit: {e}")

    print(f"\nSuccess! {sequences_processed} sequences processed with O(1) memory.")
    print(f"Machine-readable matrix written to: {output_path.name}")


if __name__ == "__main__":
    main()
