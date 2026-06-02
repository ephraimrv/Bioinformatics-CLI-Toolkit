"""
Bioinformatics Standard Utilities

Centralized functions for file routing and genomic data parsing.
Designed to handle both sandbox (Rosalind) datasets and industrial
(NCBI/UniProt) genomic files with strict validation.
"""

import argparse
from pathlib import Path
from collections.abc import Iterator

# fmt: off
DNA_NUC_LIST = ["A", "C", "G", "T"]

DNA_CODON_TABLE: dict[str, str] = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': 'Stop', 'TAG': 'Stop',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'Stop', 'TGG': 'W',
}

RNA_CODON_TABLE: dict[str, str] = {
    "UUU": "F", "CUU": "L", "AUU": "I", "GUU": "V",
    "UUC": "F", "CUC": "L", "AUC": "I", "GUC": "V",
    "UUA": "L", "CUA": "L", "AUA": "I", "GUA": "V",
    "UUG": "L", "CUG": "L", "AUG": "M", "GUG": "V",
    "UCU": "S", "CCU": "P", "ACU": "T", "GCU": "A",
    "UCC": "S", "CCC": "P", "ACC": "T", "GCC": "A",
    "UCA": "S", "CCA": "P", "ACA": "T", "GCA": "A",
    "UCG": "S", "CCG": "P", "ACG": "T", "GCG": "A",
    "UAU": "Y", "CAU": "H", "AAU": "N", "GAU": "D",
    "UAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
    "UAA": "Stop", "CAA": "Q", "AAA": "K", "GAA": "E",
    "UAG": "Stop", "CAG": "Q", "AAG": "K", "GAG": "E",
    "UGU": "C", "CGU": "R", "AGU": "S", "GGU": "G",
    "UGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
    "UGA": "Stop", "CGA": "R", "AGA": "R", "GGA": "G",
    "UGG": "W", "CGG": "R", "AGG": "R", "GGG": "G"
}

MONOISOTOPIC_MASS_TABLE: dict[str, float] = {
    'A': 71.03711, 'C': 103.00919, 'D': 115.02694,
    'E': 129.04259, 'F': 147.06841, 'G': 57.02146,
    'H': 137.05891, 'I': 113.08406, 'K': 128.09496,
    'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
    'P': 97.05276, 'Q': 128.05858, 'R': 156.10111,
    'S': 87.03203, 'T': 101.04768, 'V': 99.06841,
    'W': 186.07931, 'Y': 163.06333
}
# fmt: on


def base_parser(description_text: str) -> argparse.ArgumentParser:
    """
    Creates a standard argument parser with default input/output arguments.
    """
    parser = argparse.ArgumentParser(description=description_text)
    parser.add_argument(
        "-i", "--input", type=Path, required=True, help="Path to input file."
    )
    parser.add_argument(
        "-o", "--output", type=Path, required=True, help="Path to save the output file."
    )
    return parser


def parse_fasta(fasta_string: str) -> dict[str, str]:
    """Parses a raw FASTA text string into a dictionary of ID to sequence mappings."""
    fasta_dict: dict[str, str] = {}
    fasta_id: str = ""
    seq_buffer: list[str] = []

    for line in fasta_string.splitlines():
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if fasta_id:
                fasta_dict[fasta_id] = "".join(seq_buffer)

            header_parts = line[1:].strip().split(None, 1)
            if not header_parts:
                raise ValueError(
                    "Malformed FASTA: Empty identifier header ('>') found."
                )

            fasta_id = header_parts[0]
            seq_buffer = []
        else:
            if not fasta_id:
                raise ValueError(
                    "Malformed FASTA: Sequence data found before an identifier header."
                )
            seq_buffer.append(line)

    if fasta_id:
        fasta_dict[fasta_id] = "".join(seq_buffer)

    return fasta_dict


def lazy_parse_fasta(file_path: str | Path) -> Iterator[tuple[str, str]]:
    """Yields one (ID, sequence) tuple at a time from a FASTA file on disk."""
    fasta_id = ""
    seq_buffer = []

    with open(file_path, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if fasta_id:
                    yield fasta_id, "".join(seq_buffer)

                header_parts = line[1:].strip().split(None, 1)
                if not header_parts:
                    raise ValueError(
                        "Malformed FASTA: Empty identifier header ('>') found."
                    )

                fasta_id = header_parts[0]
                seq_buffer = []
            else:
                if not fasta_id:
                    raise ValueError(
                        "Malformed FASTA: Sequence data found before an identifier header."
                    )
                seq_buffer.append(line)

        if fasta_id:
            yield fasta_id, "".join(seq_buffer)


def get_rosalind_paths(input_filename: str, output_filename: str) -> tuple[Path, Path]:
    """Constructs absolute paths for standardized Rosalind input/output routing."""
    utils_dir = Path(__file__).resolve().parent
    input_path = utils_dir.parent / "dataset_bioinformatics_stronghold" / input_filename
    output_path = (
        utils_dir.parent / "outputs_bioinformatics_stronghold" / output_filename
    )

    if not input_path.exists():
        raise FileNotFoundError(f"Critical Error: could not find {input_path}")

    return input_path, output_path
