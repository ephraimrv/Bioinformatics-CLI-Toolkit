#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

"""Bioinformatics Standard Utilities

Note:
    v1.3.0: Added ``extract_upstream_sequence_with_length()``. Found while
    reviewing comparative_kmer_analyzer.py: ``extract_upstream_sequence()``
    has always discarded ``actual_upstream`` (the real available length,
    shorter than requested only when truncated by a contig boundary) to
    preserve its 4-tuple return contract for existing callers
    (gbk_promoter_finder.py also unpacks exactly 4 values). This meant a
    silently-truncated upstream window — e.g. a gene 40bp from its
    contig's start returning only 40bp when 150bp was requested — was
    undetectable by any caller. The new function exposes the 5th value;
    ``extract_upstream_sequence()`` now delegates to it internally and is
    otherwise unchanged.
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.3.0"

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
}

RNA_CODON_TABLE: dict[str, str] = {
    "UUU": "F",
}

MONOISOTOPIC_MASS_TABLE: dict[str, float] = {
    'A': 71.03711,
}

_REVCOMP_SRC = "ACGTNRYSWKMBDHVacgtnryswkmbdhv"
_REVCOMP_DST = "TGCANYRSWMKVHDBtgcanyrswmkvhdb"
_REVCOMP_TABLE = str.maketrans(_REVCOMP_SRC, _REVCOMP_DST)
# fmt: on


def revcomp(seq: str) -> str:
    """Returns the reverse complement of a DNA string."""
    return seq.translate(_REVCOMP_TABLE)[::-1]


def resolve_upstream_window(
    seq_len: int,
    start: int,
    end: int,
    strand: int,
    upstream_bp: int,
) -> tuple[int, int, int]:
    """Resolves strand-aware upstream-window slice boundaries (arithmetic only)."""
    if strand == 1:
        slice_start = max(0, start - upstream_bp)
        actual_upstream = start - slice_start
        return slice_start, start, actual_upstream
    else:
        slice_end = min(seq_len, end + upstream_bp)
        actual_upstream = slice_end - end
        return end, slice_end, actual_upstream


def extract_upstream_window(
    record,
    start: int,
    end: int,
    strand: int,
    upstream_bp: int,
) -> tuple[str, int, int, int]:
    """Extracts the strand-corrected upstream DNA sequence for a feature."""
    slice_start, slice_end, actual_upstream = resolve_upstream_window(
        len(record.seq), start, end, strand, upstream_bp
    )
    if strand == 1:
        upstream_seq = str(record.seq[slice_start:slice_end])
    else:
        upstream_seq = str(record.seq[slice_start:slice_end].reverse_complement())
    return upstream_seq, actual_upstream, slice_start, slice_end


def extract_upstream_sequence_with_length(
    gbk_path: Path, locus_tag: str, upstream_bp: int
) -> tuple[str, int, int, int, int]:
    """Extracts the upstream promoter region of a gene, with truncation info.

    Identical to ``extract_upstream_sequence()`` except it additionally
    returns ``actual_upstream`` — the real number of upstream bases
    available, which is less than ``upstream_bp`` only when the feature
    sits within ``upstream_bp`` of its contig's edge.

    ``extract_upstream_sequence()`` discards this value for backward
    compatibility (its callers unpack a fixed 4-tuple); use this function
    instead when truncation matters to your analysis — e.g. comparing two
    genes' upstream windows, where one being silently shorter than
    requested would otherwise go unnoticed and skew any length-normalized
    statistic computed from it.

    Args:
        gbk_path: Path to the GenBank file.
        locus_tag: The unique locus tag of the target gene.
        upstream_bp: Number of base pairs to extract upstream of the start codon.

    Returns:
        A tuple containing: (sequence, start, end, strand, actual_upstream).

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

                        seq, actual_upstream, _slice_start, _slice_end = (
                            extract_upstream_window(
                                record, start, end, strand, upstream_bp
                            )
                        )
                        return seq, start, end, strand, actual_upstream

    except (FileNotFoundError, OSError, UnicodeDecodeError, ValueError) as e:
        raise ValueError(f"Failed to parse GenBank file {gbk_path.name}: {e}") from e

    raise ValueError(f"Locus tag '{locus_tag}' not found in {gbk_path.name}.")


def extract_upstream_sequence(
    gbk_path: Path, locus_tag: str, upstream_bp: int
) -> tuple[str, int, int, int]:
    """Extracts the upstream promoter region of a specific gene.

    Delegates to ``extract_upstream_sequence_with_length()`` and drops the
    ``actual_upstream`` element to preserve this function's original
    4-tuple contract for existing callers. Callers that need to detect
    contig-boundary truncation should use
    ``extract_upstream_sequence_with_length()`` instead.

    Args:
        gbk_path: Path to the GenBank file.
        locus_tag: The unique locus tag of the target gene.
        upstream_bp: Number of base pairs to extract upstream of the start codon.

    Returns:
        A tuple containing: (sequence, start, end, strand).

    Raises:
        ValueError: If the locus tag is not found or the file cannot be parsed.
    """
    seq, start, end, strand, _actual_upstream = extract_upstream_sequence_with_length(
        gbk_path, locus_tag, upstream_bp
    )
    return seq, start, end, strand


def stream_reference_files(target_path: Path) -> Iterator[Path]:
    """Yields valid GenBank or FASTA files from a file or directory."""
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
        for ext in ("*.gbk", "*.gbff", "*.fasta", "*.fa", "*.faa"):
            yield from target_path.rglob(ext)

    else:
        raise ValueError(f"Provided path does not exist: {target_path}")


def calculate_mature_core(full_protein: str) -> str:
    """Calculates the mature peptide core based on double-glycine cleavage sites."""
    if "GG" not in full_protein:
        return full_protein

    parts = full_protein.split("GG", 1)
    mature_peptide = parts[1]

    MIN_LENGTH = 25
    if len(mature_peptide) < MIN_LENGTH:
        return mature_peptide

    for i in range(MIN_LENGTH, len(mature_peptide)):
        current_residue = mature_peptide[i]

        if current_residue == "P":
            return mature_peptide[: i + 1]

        if i + 5 <= len(mature_peptide):
            window = mature_peptide[i : i + 5]
            analyzer = ProteinAnalysis(window)
            avg_hydro = analyzer.gravy()

            if avg_hydro < -0.5:
                return mature_peptide[:i]

    return mature_peptide


@contextlib.contextmanager
def smart_open(filename: Path | None) -> Iterator[IO]:
    """Routes output to a file or standard output."""
    if filename:
        handle = open(filename, "w", encoding="utf-8")
        try:
            yield handle
        finally:
            handle.close()
    else:
        yield sys.stdout


def wrap_fasta(sequence: str, width: int = 60) -> str:
    """Wraps a sequence string into multiple lines for FASTA formatting."""
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def base_parser(
    description_text: str, include_input: bool = True, include_output: bool = True
) -> argparse.ArgumentParser:
    """Creates a standard CLI argument parser for pipeline scripts."""
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
    """Parses a raw multi-FASTA text string into an ID-to-sequence dictionary."""
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
    """Reads a FASTA file line-by-line to minimize memory usage."""
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
