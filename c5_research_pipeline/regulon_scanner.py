"""
Genome-Wide Regulon Mapper

Maps transcriptional networks by identifying operator motifs in upstream regions.

The pipeline isolates the upstream sequence (default 150bp) of every CDS in
a GenBank assembly. It performs a regex-based motif search for a provided
IUPAC/Regex operator footprint and compiles matches into a genomic matrix
suitable for network analysis.

Author: Jan Ephraim R. Vallente (ephrvallente@gmail.com)
Date: 2026-06-07
License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).

Example Usage:
    $ python3 regulon_scanner.py -i C5_genome.gbk -u 200 -m "GCGCAG[CT]G[GT]T[TA]AAAT" -o regulon.tsv
"""

__version__ = "1.0.2"

import re
import sys
import traceback
from pathlib import Path
from typing import Iterator
from Bio import SeqIO
from utils import base_parser


def stream_regulon_hits(
    gbk_path: Path, regex_pattern: str, upstream_bp: int
) -> Iterator[dict]:
    try:
        safe_pattern = re.compile(f"(?=({regex_pattern}))", re.IGNORECASE)
    except re.error as e:
        raise ValueError(f"Invalid regex pattern: '{regex_pattern}'") from e

    try:
        for record in SeqIO.parse(gbk_path, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":

                    locus_tag = feature.qualifiers.get("locus_tag", ["UNKNOWN"])[0]
                    product = feature.qualifiers.get(
                        "product", ["hypothetical protein"]
                    )[0]
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    strand = feature.location.strand

                    # Extract Upstream with Boundary Tracking
                    if strand == 1:
                        slice_start = max(0, start - upstream_bp)
                        actual_upstream = start - slice_start
                        upstream_seq = str(record.seq[slice_start:start])
                    else:
                        slice_end = min(len(record.seq), end + upstream_bp)
                        actual_upstream = slice_end - end
                        raw_seq = record.seq[end:slice_end]
                        upstream_seq = str(raw_seq.reverse_complement())

                    if actual_upstream < upstream_bp:
                        print(
                            f"    [!] Warning: {locus_tag} upstream truncated to {actual_upstream}bp (contig boundary).",
                            file=sys.stderr,
                        )

                    matches = []
                    for match in safe_pattern.finditer(upstream_seq):
                        matches.append((match.start() + 1, match.group(1)))

                    if matches:
                        yield {
                            "locus_tag": locus_tag,
                            "product": product,
                            "contig": record.id,
                            "strand": strand,
                            "matches": matches,
                        }

    except (FileNotFoundError, OSError, UnicodeDecodeError, ValueError) as e:
        raise ValueError(f"GenBank Parsing Error in {gbk_path.name}: {e}") from e


def main() -> None:
    parser = base_parser("Genome-Wide Regulon Scanner")
    parser.add_argument(
        "-u", "--upstream", type=int, default=150, help="Upstream bases to extract"
    )
    parser.add_argument("-m", "--motif", required=True, help="Regex motif to search")
    args = parser.parse_args()

    print(f"[*] Scanning entire genome: {args.input.name}")
    print(f"[*] Upstream region: {args.upstream}bp")
    print(f"[*] Searching for Motif: {args.motif}")

    total_genes_hit = 0
    total_motifs_found = 0

    try:
        results_iterator = stream_regulon_hits(args.input, args.motif, args.upstream)

        if args.output:
            with open(args.output, "w", encoding="utf-8") as tsv:
                tsv.write(
                    "Locus_Tag\tContig\tStrand\tMotif_Count\tRelative_Positions\tMatched_Sequences\tProduct\n"
                )
                for hit in results_iterator:
                    total_genes_hit += 1
                    total_motifs_found += len(hit["matches"])
                    positions = ",".join([str(m[0]) for m in hit["matches"]])
                    sequences = ",".join([m[1] for m in hit["matches"]])
                    tsv.write(
                        f"{hit['locus_tag']}\t{hit['contig']}\t{hit['strand']}\t{len(hit['matches'])}\t{positions}\t{sequences}\t{hit['product']}\n"
                    )
                    print(
                        f"    -> Regulon member found: {hit['locus_tag']} ({hit['product'][:40]}...)"
                    )
        else:
            for hit in results_iterator:
                total_genes_hit += 1
                total_motifs_found += len(hit["matches"])
                print(
                    f"    -> Regulon member found: {hit['locus_tag']} ({hit['product'][:40]}...)"
                )
            print(
                "\n[*] Note: No output file specified (-o). Results printed to terminal only."
            )

        print("\n" + "=" * 40)
        print("          PIPELINE COMPLETE")
        print("=" * 40)
        print(f"Total Genes in Regulon: {total_genes_hit}")
        print(f"Total Motifs Bound:     {total_motifs_found}")
        if args.output:
            print(f"Results written to:     {args.output.name}")
        print("=" * 40)

    except (ValueError, FileNotFoundError, PermissionError) as e:
        sys.exit(f"\n[!] Pipeline Halted: {e}")
    except KeyboardInterrupt:
        sys.exit("\n[!] Pipeline gracefully interrupted by user.")
    except Exception:
        print("\n[!] UNEXPECTED BUG ENCOUNTERED:")
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
