"""
GenBank Upstream Promoter Finder

Extracts regulatory promoter regions from multi-contig GenBank assemblies.

This tool locates a target gene by its locus tag, calculates the strand
orientation to extract the correct upstream sequence (applying reverse
complementation where necessary), and scans the region for specific motif hits.

License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).

Example Usage:
    $ python3 gbk_promoter_finder.py -i C5_genome.gbk -l ctg1_50 -u 150 -m "TATAAT" -o ctg1_50_promoter.fasta
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.1.0"

import re
import sys
import traceback
from typing import Iterator

try:
    from Bio.Seq import Seq
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )
from utils import base_parser, wrap_fasta, extract_upstream_sequence


def find_motif_regex_iterator(
    sequence: str, regex_pattern: str, actual_len: int
) -> Iterator[tuple[int, str, str]]:
    """Scans a sequence for a motif on both strands.

    Positions are reported as negative integers relative to the Translation
    Start Site (TSS), following standard molecular biology convention
    (e.g., -35 means 35 bases upstream of the ATG).

    Both the coding (+) and template (-) strands are scanned so TF binding
    sites in either orientation are detected.

    Args:
        sequence:       The upstream DNA sequence string (coding strand, 5'→3').
        regex_pattern:  IUPAC/regex motif (IGNORECASE applied).
        actual_len:     Actual length of the upstream sequence (may be shorter
                        than requested if near a contig boundary).

    Yields:
        Tuple of (rel_pos, matched_sequence, strand_indicator).
        rel_pos is a negative integer relative to TSS.
        strand_indicator is '+' or '-'.
    """
    if not sequence or not regex_pattern:
        return
    try:
        safe_pattern = re.compile(f"(?=({regex_pattern}))", re.IGNORECASE)
    except re.error as e:
        raise ValueError(f"Invalid regex pattern provided: {regex_pattern}") from e

    # Forward (coding) strand scan
    for match in safe_pattern.finditer(sequence):
        rel_pos = -(actual_len - match.start())
        yield rel_pos, match.group(1), "+"

    # Reverse complement (template) strand scan
    rc_seq = str(Seq(sequence).reverse_complement())
    for match in safe_pattern.finditer(rc_seq):
        fwd_pos = len(sequence) - match.end()
        rel_pos = -(actual_len - fwd_pos)
        yield rel_pos, match.group(1), "-"


def main() -> None:
    parser = base_parser("GenBank Targeted Upstream Motif Scanner")
    parser.add_argument("-l", "--locus", required=True, help="Target gene locus tag")
    parser.add_argument(
        "-u",
        "--upstream",
        type=int,
        default=100,
        help=(
            "Upstream bases to extract. Default: 100. "
            "For eukaryotes, consider --upstream 2000 or higher."
        ),
    )
    parser.add_argument(
        "-m",
        "--motif",
        required=False,
        help="Regex/IUPAC motif to search for on both strands (optional)",
    )
    args = parser.parse_args()

    try:
        if args.upstream < 1:
            raise ValueError("Upstream bases must be a positive integer.")

        upstream_seq, start, end, strand = extract_upstream_sequence(
            args.input, args.locus, args.upstream
        )

        actual_len = len(upstream_seq)

        print(
            f"[*] Found {args.locus} at {start}-{end} (Gene strand: {strand})",
            file=sys.stderr,
        )
        print(
            f"[*] Requested: {args.upstream}bp upstream | Extracted: {actual_len}bp",
            file=sys.stderr,
        )

        if actual_len < args.upstream:
            print(
                f"[!] Warning: Upstream truncated to {actual_len}bp (contig boundary).",
                file=sys.stderr,
            )

        motifs = []
        if args.motif:
            print(
                f"[*] Searching for motif: {args.motif} (both strands)",
                file=sys.stderr,
            )
            # Sort results biologically: most negative (farthest from TSS) first
            motifs = sorted(
                find_motif_regex_iterator(upstream_seq, args.motif, actual_len),
                key=lambda x: x[0],
            )

        if args.output:
            with open(args.output, "w", encoding="utf-8") as out_file:
                # FASTA header reflects actual extracted length
                out_file.write(
                    f">{args.locus}_upstream_{actual_len}bp_gene_strand_{strand}\n"
                )
                out_file.write(f"{wrap_fasta(upstream_seq)}\n")

                # TSV motif appendix: position is TSS-relative (negative = upstream)
                if motifs:
                    out_file.write(
                        "\n# Position_Relative_to_TSS\tMotif_Strand\tSequence\n"
                    )
                    for pos, seq, motif_strand in motifs:
                        out_file.write(f"{pos}\t{motif_strand}\t{seq}\n")

            print(
                f"[*] Success! {len(motifs)} motif(s) found. Written to {args.output.resolve()}",
                file=sys.stderr,
            )

        else:
            print("\n--- UPSTREAM SEQUENCE ---", file=sys.stderr)
            if len(upstream_seq) > 500:
                print(
                    f"{upstream_seq[:100]} ... [snip {len(upstream_seq)-200}bp] ... {upstream_seq[-100:]}",
                    file=sys.stderr,
                )
            else:
                print(upstream_seq, file=sys.stderr)
            print("-------------------------\n", file=sys.stderr)

            if motifs:
                for pos, seq, motif_strand in motifs:
                    print(
                        f"    -> Motif Found! Position: {pos} ({motif_strand} strand) | Sequence: {seq}",
                        file=sys.stderr,
                    )
            elif args.motif is None:
                pass
            else:
                print("    -> No motifs found.", file=sys.stderr)

    except FileNotFoundError:
        sys.exit(f"\n[!] File not found: {args.input}")
    except ValueError as e:
        sys.exit(f"\n[!] Error: {e}")
    except KeyboardInterrupt:
        sys.exit("\n[!] Pipeline gracefully interrupted by user.")
    except Exception:
        print("\n[!] UNEXPECTED BUG ENCOUNTERED:")
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
