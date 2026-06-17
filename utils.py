#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

"""Bioinformatics Standard Utilities

Centralized utility functions for file routing, parsing, and sequence
manipulation. Provides standardized handling for GenBank and FASTA data
with integrated error validation and memory-safe processing.

Note:
    Associated with ongoing, unpublished research (manuscript in
    preparation). Correct attribution is requested when used in
    derivative works.
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.0.7"

import argparse
import contextlib
import sys
from pathlib import Path
from collections.abc import Iterator
from typing import IO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO

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


def extract_upstream_sequence(
    gbk_path: Path, locus_tag: str, upstream_bp: int
) -> tuple[str, int, int, int]:
    """Extracts the upstream promoter region of a specific gene.

    Args:
        gbk_path: Path to the GenBank file.
        locus_tag: The unique locus tag of the target gene.
        upstream_bp: Number of base pairs to extract upstream of the start codon.

    Returns:
        A tuple containing: (sequence, start, end, strand).

    Raises:
        ValueError: If the locus tag is not found or the file cannot be parsed.
    """
    try:
        for record in SeqIO.parse(gbk_path, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    if locus_tag in feature.qualifiers.get("locus_tag", []):
                        start = int(feature.location.start)
                        end = int(feature.location.end)
                        strand = feature.location.strand

                        if strand == 1:
                            slice_start = max(0, start - upstream_bp)
                            return (
                                str(record.seq[slice_start:start]),
                                start,
                                end,
                                strand,
                            )
                        else:
                            slice_end = min(len(record.seq), end + upstream_bp)
                            raw_seq = record.seq[end:slice_end]
                            return str(raw_seq.reverse_complement()), start, end, strand

    except (FileNotFoundError, OSError, UnicodeDecodeError, ValueError) as e:
        raise ValueError(f"Failed to parse GenBank file {gbk_path.name}: {e}") from e

    raise ValueError(f"Locus tag '{locus_tag}' not found in {gbk_path.name}.")


def stream_reference_files(target_path: Path) -> Iterator[Path]:
    """
    Yields valid GenBank or FASTA files from a file or directory.

    If given a directory, it recursively searches all sub-directories.

    Args:
        target_path: Path object pointing to a single file or a directory.

    Yields:
        Path objects for every valid genomic file found.

    Raises:
        ValueError: If the provided path does not exist.
    """
    valid_exts = (".gbk", ".gbff", ".fasta", ".fa", ".faa")
    import sys  # Imported locally for the warning stream

    valid_exts = (".gbk", ".gbff", ".fasta", ".fa", ".faa")

    if target_path.is_file():
        if target_path.suffix.lower() in valid_exts:
            yield target_path
        else:
            print(
                f"[!] Warning: {target_path.name} is not a valid GenBank/FASTA extension.",
                file=sys.stderr,
            )

    elif target_path.is_dir():
        # rglob recursively searches the folder and all sub-folders
        for ext in ("*.gbk", "*.gbff", "*.fasta", "*.fa", "*.faa"):
            yield from target_path.rglob(ext)

    else:
        raise ValueError(f"Provided path does not exist: {target_path}")


def calculate_mature_core(full_protein: str) -> str:
    """
    Calculates the mature peptide core based on double-glycine cleavage sites.

    Trims the leader sequence at the first 'GG' site and removes C-terminal
    hydrophilic tails based on a Kyte-Doolittle hydrophobicity gradient.

    Args:
        full_protein: The full protein sequence string.

    Returns:
        The trimmed mature protein sequence.
    """
    if "GG" not in full_protein:
        return full_protein  # If no cut site, use the whole sequence

    # Split at the FIRST 'GG' and keep everything after it
    parts = full_protein.split("GG", 1)
    mature_peptide = parts[1]

    MIN_LENGTH = 25
    if len(mature_peptide) < MIN_LENGTH:
        return mature_peptide

    # Scan the sequence for structural boundaries
    for i in range(MIN_LENGTH, len(mature_peptide)):
        current_residue = mature_peptide[i]

        # Rule A: Proline Helix-Breaker
        if current_residue == "P":
            return mature_peptide[: i + 1]

        # Rule B: Kyte-Doolittle Hydrophobicity Drop
        if i + 5 <= len(mature_peptide):
            window = mature_peptide[i : i + 5]
            analyzer = ProteinAnalysis(window)

            # FIX: Use the built-in GRAVY (Grand Average of Hydropathy) method
            # This perfectly replaces the broken manual kd dictionary math.
            avg_hydro = analyzer.gravy()

            # If the window becomes highly charged/hydrophilic, snip it!
            if avg_hydro < -0.5:
                return mature_peptide[:i]

    return mature_peptide


@contextlib.contextmanager
def smart_open(filename: Path | None) -> Iterator[IO]:
    """
    Routes output to a file or standard output.

    Args:
        filename: Target path. If None, routes to sys.stdout.

    Yields:
        An open file handle or sys.stdout.
    """
    if filename:
        handle = open(filename, "w", encoding="utf-8")
        try:
            yield handle
        finally:
            handle.close()
    else:
        yield sys.stdout


def wrap_fasta(sequence: str, width: int = 60) -> str:
    """
    Wraps a sequence string into multiple lines for FASTA formatting.

    Args:
        sequence: The raw nucleotide or amino acid string.
        width: Maximum number of characters per line.

    Returns:
        The line-wrapped sequence string.
    """
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def base_parser(
    description_text: str, include_input: bool = True, include_output: bool = True
) -> argparse.ArgumentParser:
    """
    Creates a standard CLI argument parser for pipeline scripts.

    Args:
        description_text: Description displayed in the help menu.
        include_input: If True, adds the -i/--input flag.
        include_output: If True, adds the -o/--output flag.

    Returns:
        A configured argparse.ArgumentParser instance.
    """
    parser = argparse.ArgumentParser(description=description_text)

    if include_input:
        parser.add_argument(
            "-i", "--input", type=Path, required=True, help="Path to input file."
        )

    if include_output:
        parser.add_argument(
            "-o",
            "--output",
            type=Path,
            required=False,
            default=None,
            help="Path to save the output file (optional — prints to terminal if omitted).",
        )

    return parser


def parse_fasta(fasta_string: str) -> dict[str, str]:
    """
    Parses a raw multi-FASTA text string into an ID-to-sequence dictionary.

    Args:
        fasta_string: Multi-FASTA formatted text.

    Returns:
        A dictionary mapping sequence IDs to sequence strings.

    Raises:
        ValueError: If the FASTA structure is malformed.
    """
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
    """
    Reads a FASTA file line-by-line to minimize memory usage.

    Args:
        file_path: Path to the FASTA file.

    Yields:
        A tuple of (sequence_id, sequence).

    Raises:
        ValueError: If the FASTA structure is malformed.
    """
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
    """
    Constructs absolute paths for Rosalind dataset management.

    Args:
        input_filename: Dataset file name.
        output_filename: Results file name.

    Returns:
        A tuple of (input_path, output_path).

    Raises:
        FileNotFoundError: If the input path does not exist.
    """
    utils_dir = Path(__file__).resolve().parent
    input_path = utils_dir.parent / "dataset_bioinformatics_stronghold" / input_filename
    output_path = (
        utils_dir.parent / "outputs_bioinformatics_stronghold" / output_filename
    )

    if not input_path.exists():
        raise FileNotFoundError(f"Critical Error: could not find {input_path}")

    return input_path, output_path
