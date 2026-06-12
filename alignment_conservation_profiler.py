"""
Alignment Conservation Profiler

Builds a Position Probability Matrix (PPM) and Information Content profile
from a multiple sequence alignment (MSA) in FASTA format.

This tool accepts a gap-aligned FASTA file and computes position-by-position
conservation metrics for sequence logo generation and motif discovery. Outputs
include the consensus sequence, PPM (probabilities of each nucleotide at each
position), and gap-penalized Information Content using the WebLogo standard
(Crooks et al. 2004).

Handles real-world alignments gracefully: ambiguous IUPAC characters (N, Y, R, W, etc.)
are silently skipped without crashing the pipeline. Gaps in columns are penalized
in the IC calculation to correctly report conservation.

License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).

Example Usage:
    $ python3 alignment_conservation_profiler.py -i alignment.fasta -o profile_matrix.tsv
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.1.0"

import math
import sys

try:
    from Bio import SeqIO
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )
from utils import base_parser


def build_industrial_profile(
    sequence_records,
) -> tuple[str, dict[str, list[float]], list[float]]:
    """Builds PPM and Information Content from aligned sequences.

    Args:
        sequence_records: An iterator yielding SeqRecord objects
                          (e.g., from SeqIO.parse(..., "fasta")).

    Returns:
        A tuple of (consensus_string, ppm_dict, information_content_list).

    Raises:
        ValueError: If sequences are unaligned, contain invalid characters, or stream is empty.
    """
    try:
        first_record = next(sequence_records)
    except StopIteration:
        raise ValueError("Pipeline Halted: The FASTA stream is empty.")

    first_seq = str(first_record.seq).upper()
    first_id = first_record.id

    if not first_seq:
        raise ValueError(
            f"Pipeline Halted: First sequence '{first_id}' contains no bases."
        )

    sequence_length = len(first_seq)
    total_sequences = 1
    valid_chars = "ACGT-"
    counts = {char: [0] * sequence_length for char in valid_chars}

    for index, char in enumerate(first_seq.upper()):
        if char in counts:
            counts[char][index] += 1
        # Ambiguous characters (N, Y, R, etc.) are silently skipped

    for record in sequence_records:
        clean_seq = str(record.seq).upper()
        if len(clean_seq) != sequence_length:
            raise ValueError(
                f"Pipeline Halted: Alignment broken at {record.id}. "
                f"Expected length {sequence_length}, found {len(clean_seq)}."
            )

        for index, char in enumerate(clean_seq):
            if char in counts:
                # Only count known ACGT/gap characters.
                # Ambiguous IUPAC codes (N, Y, R, W, etc.) are silently skipped —
                # they are not added to any count, which means they are treated as
                # absent from that column. This prevents pipeline crashes on
                # real-world draft assemblies and eukaryotic alignments.
                counts[char][index] += 1

        total_sequences += 1

    consensus_list = []
    ppm = {char: [0.0] * sequence_length for char in valid_chars}
    information_content = [0.0] * sequence_length
    max_bits = math.log2(4.0)  # Always 2.0 bits for a 4-letter DNA alphabet

    for i in range(sequence_length):
        max_count = -1
        winning_char = "-"  # Default to gap; prevents wrong nucleotide appearing
        # as consensus for all-gap columns (which have all
        # nucleotide counts at 0, beating initial max_count=-1)
        column_entropy_sum = 0.0

        # Count only ACGT bases in this column — gaps and skipped ambiguous
        # characters do not contribute to the nucleotide entropy calculation.
        acgt_count = sum(counts[c][i] for c in "ACGT")

        for char in valid_chars:
            current_count = counts[char][i]

            # PPM: probability across ALL sequences (gaps included in denominator)
            ppm[char][i] = current_count / total_sequences

            # Consensus: allow gap to win if it dominates
            if current_count > max_count:
                max_count = current_count
                winning_char = char

            # Entropy: calculated ONLY over nucleotides, normalized against
            # the number of sequences that have an actual base at this position.
            # This ensures gap-heavy columns don't artificially suppress entropy.
            if char in "ACGT" and acgt_count > 0:
                base_prob = current_count / acgt_count
                if base_prob > 0:
                    column_entropy_sum += base_prob * math.log2(base_prob)

        consensus_list.append(winning_char)

        # Gap-penalized Information Content (Crooks et al. 2004, WebLogo standard):
        #   IC = (log2(4) - H_nucleotides) * fraction_present
        #
        # The fraction_present term scales IC down proportionally to how many
        # sequences have an actual base at this position. A 100% gap column
        # gets IC = 0.0 (not 2.0 bits as the unpenalized formula would give).
        # A 50% gap, 100% conserved column gets IC = 1.0, not 2.0.
        fraction_present = acgt_count / total_sequences
        raw_ic = (
            max_bits + column_entropy_sum
        )  # entropy_sum is negative; adding gives IC
        information_content[i] = max(0.0, raw_ic * fraction_present)

    consensus_string = "".join(consensus_list)
    return consensus_string, ppm, information_content


def main() -> None:
    parser = base_parser("Alignment Conservation Profiler")
    args = parser.parse_args()

    print(
        f"[*] Building PPM and calculating Information Content for: {args.input.name}",
        file=sys.stderr,
    )

    try:
        # Use BioPython's standard SeqIO.parse() for publication-grade reproducibility
        sequence_records = SeqIO.parse(args.input, "fasta")
        consensus, ppm, entropy = build_industrial_profile(sequence_records)

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
            print(
                f"\n[*] Success! Machine-readable TSV written to: {args.output.resolve()}",
                file=sys.stderr,
            )
        except OSError as e:
            sys.exit(f"\n[!] Error: Could not write to {args.output.name}. Reason: {e}")
    else:
        print(
            "\n[*] Note: No output file specified (-o). Results printed to terminal only.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
