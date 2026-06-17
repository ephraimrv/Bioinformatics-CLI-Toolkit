"""
Multi-FASTA Quality Control Engine

Processes massive multi-FASTA datasets (including multi-gigabyte chromosomes).
Reads the file line-by-line, calculates nucleotide counts on the fly, and 
discards the sequence string immediately to guarantee strict O(1) memory usage.

Usage:
    $ python3 stream_nuc_count.py -i massive_genome.fasta -o counts.tsv
"""

__author__ = 'Jan Ephraim R. Vallente'
__email__ = 'ephrvallente@gmail.com'
__version__ = '1.0.1'

import sys
from pathlib import Path
from typing import Iterator
from utils import base_parser


def stream_nuc_count(file_path: str | Path) -> Iterator[tuple[str, dict[str, int]]]:
    """
    A memory-safe generator that calculates sequence metadata on the fly.

    Architecture & Performance:
    Unlike standard parsers that buffer a sequence until the next '>' is found,
    this engine processes the math line-by-line. The raw sequence string is
    destroyed immediately after the math is computed. Memory usage remains 
    near zero, even if a single sequence is 50 Gigabytes long.

    Args:
        file_path: The system path to the target FASTA file.

    Yields:
        A tuple containing the FASTA ID and a dictionary of its base counts.
        
    Raises:
        ValueError: If the input file cannot be found or the FASTA is malformed.
    """
    current_id = ""
    totals = {"A": 0, "C": 0, "G": 0, "T": 0, "Anomalies": 0}
    DNA_TRANSLATOR = str.maketrans("acgtn", "ACGTN", " \n\r")

    try:
        with open(file_path, "r", encoding="utf-8") as file:
            for line in file:
                line = line.strip()

                if not line:
                    continue

                if line.startswith(">"):
                    # If we were tracking a sequence, yield its final math before resetting
                    if current_id:
                        yield current_id, totals.copy()

                    header_parts: list = line[1:].strip().split(None, 1)

                    if not header_parts:
                        raise ValueError(
                            "Malformed FASTA: Empty identifier header ('>') found."
                        )
                    # Reset for the new sequence
                    current_id: str = header_parts[0]
                    totals = {"A": 0, "C": 0, "G": 0, "T": 0, "Anomalies": 0}

                else:
                    # The Malformed File Guard Clause
                    if not current_id:
                        raise ValueError(
                            "Malformed FASTA: Sequence data found before an identifier header."
                        )

                    # Clean the line and do the math immediately
                    clean_line = line.translate(DNA_TRANSLATOR)

                    a_count = clean_line.count("A")
                    c_count = clean_line.count("C")
                    g_count = clean_line.count("G")
                    t_count = clean_line.count("T")

                    totals["A"] += a_count
                    totals["C"] += c_count
                    totals["G"] += g_count
                    totals["T"] += t_count

                    # Ambiguous bases
                    totals["Anomalies"] += len(clean_line) - (
                        a_count + c_count + g_count + t_count
                    )

            # The Final Flush: Yield the math for the last sequence in the file
            if current_id:
                yield current_id, totals

    except FileNotFoundError as e:
        # Preserve traceback and system context using 'from e'
        raise ValueError(f"Input file not found at '{file_path}'") from e


def main() -> None:
    """
    Pipeline manager for I/O routing and TSV formatting.
    
    Reads standard CLI arguments, processes the streaming nucleotide counts, 
    and safely outputs a machine-readable tab-separated matrix.
    """

    args = base_parser("Nucleotide Count Pipeline for Genome Assemblies").parse_args()

    input_path = args.input
    output_path = args.output

    # Formatting
    print("--- TRUE INDUSTRIAL MULTI-FASTA QC REPORT ---")
    sequences_processed = 0

    try:
        with open(output_path, "w", encoding="utf-8") as out_file:
            out_file.write(
                "Sequence_ID\tA_count\tC_count\tG_count\tT_count\tAnomalies\n"
            )

            for fasta_id, counts in stream_nuc_count(input_path):
                sequences_processed += 1

                print(f"> {fasta_id}")
                print(
                    f"  A: {counts['A']:,} | C: {counts['C']:,} | G: {counts['G']:,} | T: {counts['T']:,} | Anomalies: {counts['Anomalies']:,}"
                )

                # Machine Readable Output (TSV)
                row = (
                    f"{fasta_id}\t{counts['A']}\t{counts['C']}\t"
                    f"{counts['G']}\t{counts['T']}\t{counts['Anomalies']}\n"
                )
                out_file.write(row)

        if sequences_processed == 0:
            raise ValueError(
                "Pipeline Halted: The input file contained no valid FASTA records."
            )

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
