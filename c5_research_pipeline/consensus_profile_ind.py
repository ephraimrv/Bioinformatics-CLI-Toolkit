"""
Consensus and Motif Conservation Profiler

Constructs a Position Probability Matrix (PPM) and calculates Shannon entropy.

This script ingests aligned FASTA sequences to perform column-wise frequency
counts. It computes a PPM and Shannon entropy to quantify positional
conservation, outputting a machine-readable TSV matrix for downstream analysis.

License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).

Example Usage:
    $ python3 consensus_profile_ind.py -i alignment.fasta -o profile_matrix.tsv
"""


__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.0.2"

import math
import sys
from typing import Iterator
from utils import base_parser, lazy_parse_fasta


def build_industrial_profile(
    fasta_iterator: Iterator[tuple[str, str]],
) -> tuple[str, dict[str, list[float]], list[float]]:
    try:
        first_id, first_seq = next(fasta_iterator)
    except StopIteration:
        raise ValueError("Pipeline Halted: The FASTA stream is empty.")

    if not first_seq:
        raise ValueError(
            f"Pipeline Halted: First sequence '{first_id}' contains no bases."
        )

    sequence_length = len(first_seq)
    total_sequences = 1
    valid_chars = "ACGT-"
    counts = {char: [0] * sequence_length for char in valid_chars}

    for index, char in enumerate(first_seq.upper()):
        if char not in counts:
            raise ValueError(
                f"Pipeline Halted: Invalid character '{char}' found in {first_id}."
            )
        counts[char][index] += 1

    for fasta_id, seq in fasta_iterator:
        clean_seq = seq.upper()
        if len(clean_seq) != sequence_length:
            raise ValueError(
                f"Pipeline Halted: Alignment broken at {fasta_id}. "
                f"Expected length {sequence_length}, found {len(clean_seq)}."
            )

        for index, char in enumerate(clean_seq):
            if char not in counts:
                raise ValueError(
                    f"Pipeline Halted: Invalid character '{char}' found in {fasta_id}."
                )
            counts[char][index] += 1

        total_sequences += 1

    consensus_list = []
    ppm = {char: [0.0] * sequence_length for char in valid_chars}
    information_content = [0.0] * sequence_length
    max_bits = math.log2(len(valid_chars.replace("-", "")))

    for i in range(sequence_length):
        max_count = -1
        winning_char = ""
        column_entropy_sum = 0.0

        for char in valid_chars:
            current_count = counts[char][i]
            probability = current_count / total_sequences
            ppm[char][i] = probability

            if char != "-" and current_count > max_count:
                max_count = current_count
                winning_char = char

            if probability > 0:
                column_entropy_sum += probability * math.log2(probability)

        consensus_list.append(winning_char)
        information_content[i] = max(0.0, max_bits + column_entropy_sum)

    consensus_string = "".join(consensus_list)
    return consensus_string, ppm, information_content


def main() -> None:
    parser = base_parser("Industrial Consensus and Motif Profiling Pipeline")
    args = parser.parse_args()

    print(
        f"[*] Building PPM and calculating Information Content for: {args.input.name}"
    )

    try:
        fasta_stream = lazy_parse_fasta(args.input)
        consensus, ppm, entropy = build_industrial_profile(fasta_stream)

    except ValueError as e:
        sys.exit(f"\n[!] Pipeline Halted: {e}")
    except FileNotFoundError:
        sys.exit(f"\n[!] Pipeline Halted: Could not find file {args.input}")
    except KeyboardInterrupt:
        sys.exit("\n[!] Pipeline gracefully interrupted by user.")

    seq_length = len(consensus)
    col_width = 8
    pos_header_visual = "".join(f"{i:<{col_width}}" for i in range(1, seq_length + 1))
    visual_rows = {
        char: "".join(f"{val:<{col_width}.3f}" for val in ppm[char]) for char in "ACGT-"
    }
    visual_entropy = "".join(f"{val:<{col_width}.3f}" for val in entropy)

    print("\n--- CONSERVATION PROFILE ---")
    print(f"Consensus: {consensus}\n")
    print("Position Probability Matrix (PPM):")
    print(f"{'pos':<{col_width}}{pos_header_visual}")
    for char in "ACGT-":
        print(f"{f'{char}:':<{col_width}}{visual_rows[char]}")

    print("\nInformation Content (Bits):")
    print(f"{'IC:':<{col_width}}{visual_entropy}")

    if args.output:
        pos_header_tsv = "\t".join(str(i) for i in range(1, seq_length + 1))
        tabbed_consensus = "\t".join(consensus)
        tsv_rows = {
            char: "\t".join(f"{val:.3f}" for val in ppm[char]) for char in "ACGT-"
        }
        tsv_entropy = "\t".join(f"{val:.3f}" for val in entropy)

        output_lines = [
            f"Position\t{pos_header_tsv}",
            f"Consensus\t{tabbed_consensus}",
        ]
        for char in "ACGT-":
            output_lines.append(f"{char}\t{tsv_rows[char]}")
        output_lines.append(f"IC\t{tsv_entropy}")

        final_output = "\n".join(output_lines)

        try:
            args.output.write_text(final_output, encoding="utf-8")
            print(f"\n[*] Success! Machine-readable TSV written to: {args.output.name}")
        except OSError as e:
            sys.exit(f"\n[!] Error: Could not write to {args.output.name}. Reason: {e}")
    else:
        print(
            "\n[*] Note: No output file specified (-o). Results printed to terminal only."
        )


if __name__ == "__main__":
    main()
