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

ALGORITHM AND IMPLEMENTATION NOTES:

    PPM builder uses zip(*sequences) + Counter for column-wise processing.
    This is significantly faster than nested Python loops for large alignments
    because Counter's inner counting loop is implemented in C (CPython). For
    a 500-sequence × 1000bp alignment, the Python nested loop would evaluate
    500,000 Python-level iterations; zip+Counter collapses that to 1000
    C-backed Counter operations. The asymptotic complexity is the same
    O(N × L), but the constant factor is substantially smaller.

    Gaps in the alignment are included in the PPM (as '-' column probabilities)
    but excluded from entropy and IC calculations. IC is penalised by the
    fraction of non-gap bases at each position (WebLogo formula).

Note on scope:
    This tool is a PROFILE SCANNER, not a motif discovery engine. It quantifies
    conservation across a pre-built alignment. For discovering unknown regulatory
    motifs from unaligned upstream sequences, use motif_discovery.py, which
    implements the Expectation-Maximization algorithm (log-odds + background
    model, MEME-style OOPS model with bidirectional strand scanning and seed
    clustering).

License: MIT

Note:
    This module is part of ongoing research and is associated with an upcoming
    publication. Please cite appropriately when used in derivative works.
    See LICENSE file in the repository root for full license terms.

Example Usage:
    $ python3 alignment_conservation_profiler.py -i alignment.fasta -o profile_matrix.tsv

Notes on Excel import:
    TSV output uses UTF-8 encoding (plain, no BOM). When you open the file in
    Excel, it will automatically recognize the header row. If the Transform Data
    wizard appears, you can close it and Excel will display the data correctly.

"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.2.0"

import math
import sys
from collections import Counter

try:
    from Bio import SeqIO
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )
from utils import base_parser


def _build_profile(
    sequences: list[str],
) -> tuple[str, dict[str, list[float]], list[float]]:
    """Core PPM + Information Content builder.

    Accepts a pre-validated list of equal-length uppercase DNA strings.
    Ambiguous IUPAC characters (N, Y, R, W, etc.) are silently skipped.
    Gaps are legal and tracked in the PPM but excluded from entropy.

    Uses ``zip(*sequences)`` + ``Counter`` for column-wise counting rather
    than a nested Python loop. Counter's inner loop runs in C (CPython),
    making this substantially faster for large alignments (hundreds of
    sequences × hundreds of positions).

    Args:
        sequences: Non-empty list of equal-length uppercase strings.

    Returns:
        (consensus_string, ppm_dict, information_content_list)
        IC uses WebLogo gap-penalized formula (Crooks et al. 2004).
    """
    seq_len = len(sequences[0])
    n = len(sequences)
    valid_chars = "ACGT-"
    counts = {c: [0] * seq_len for c in valid_chars}
    max_bits = math.log2(4.0)

    # Column-wise counting via zip(*sequences) + Counter.
    # Each iteration yields one column (a tuple of n characters) and counts
    # all characters in a single C-backed call, replacing the inner Python loop.
    for i, column in enumerate(zip(*sequences)):
        col_counts = Counter(column)
        for c in valid_chars:
            counts[c][i] = col_counts[c]  # Counter returns 0 for missing keys

    consensus_list = []
    ppm = {c: [0.0] * seq_len for c in valid_chars}
    information_content = [0.0] * seq_len

    for i in range(seq_len):
        best_count = -1
        winning_char = "-"
        entropy_sum = 0.0
        acgt_count = sum(counts[c][i] for c in "ACGT")

        for char in valid_chars:
            cnt = counts[char][i]
            ppm[char][i] = cnt / n
            if cnt > best_count:
                best_count = cnt
                winning_char = char
            if char in "ACGT" and acgt_count > 0:
                p = cnt / acgt_count
                if p > 0:
                    entropy_sum += p * math.log2(p)

        consensus_list.append(winning_char)
        fraction_present = acgt_count / n
        information_content[i] = max(0.0, (max_bits + entropy_sum) * fraction_present)

    return "".join(consensus_list), ppm, information_content


def build_industrial_profile(
    sequence_records,
) -> tuple[str, dict[str, list[float]], list[float]]:
    """Builds PPM and IC from an aligned FASTA SeqRecord iterator.

    Args:
        sequence_records: Iterator of BioPython SeqRecord objects.

    Returns:
        (consensus_string, ppm_dict, information_content_list)

    Raises:
        ValueError: If stream is empty, a sequence is empty, or lengths differ.
    """
    sequences = []
    ref_len = None

    for record in sequence_records:
        seq = str(record.seq).upper()
        if not seq:
            raise ValueError(f"Sequence '{record.id}' contains no bases.")
        if ref_len is None:
            ref_len = len(seq)
        elif len(seq) != ref_len:
            raise ValueError(
                f"Alignment broken at '{record.id}'. "
                f"Expected {ref_len}, found {len(seq)}."
            )
        sequences.append(seq)

    if not sequences:
        raise ValueError("The FASTA stream is empty.")

    return _build_profile(sequences)


def main() -> None:
    parser = base_parser("Alignment Conservation Profiler (-i = FASTA file,-o = tsv)")
    args = parser.parse_args()

    print(
        f"[*] Building PPM and Information Content for: {args.input.name}",
        file=sys.stderr,
    )

    try:
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
        try:
            # Use open() + write() like regulon_scanner.py to ensure Excel
            # recognizes row 1 as headers. Write headers first, separately.
            with open(args.output, "w", encoding="utf-8") as tsv:
                # Header row
                pos_header_tsv = "\t".join(str(i) for i in range(1, seq_length + 1))
                tsv.write(f"Position\t{pos_header_tsv}\n")

                # Data rows
                tabbed_consensus = "\t".join(consensus)
                tsv.write(f"Consensus\t{tabbed_consensus}\n")

                for char in "ACGT-":
                    row = "\t".join(f"{ppm[char][i]:.3f}" for i in range(seq_length))
                    tsv.write(f"{char}\t{row}\n")

                ic_row = "\t".join(f"{val:.3f}" for val in entropy)
                tsv.write(f"IC\t{ic_row}\n")

            print(
                f"\n[*] Success! TSV written to: {args.output.resolve()}",
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
