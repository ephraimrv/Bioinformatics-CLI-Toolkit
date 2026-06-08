"""
Genomic Open Reading Frame (ORF) Extraction Pipeline

This pipeline scans multi-FASTA DNA assemblies across all 6 reading frames
to identify and extract valid Open Reading Frames. It utilizes a heuristic
"Longest ORF" algorithm to filter out nested and redundant protein fragments,
yielding a biologically accurate, non-bloated FASTA file of distinct predicted
proteins.

Usage:
    Standard execution (defaults to minimum 100 amino acids):
        $ python3 orf_pipeline.py -i raw_contigs.fasta -o predicted_proteins.fasta

    Targeted execution (e.g., hunting 30-AA viral smORFs):
        $ python3 orf_pipeline.py -i virus.fasta -o viral_smORFs.fasta -m 30
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.0.3"

import sys
from typing import Generator
from utils import DNA_CODON_TABLE, base_parser, lazy_parse_fasta, wrap_fasta
from dna_rev_genome import reverse_complement


def prot_frame_heuristic(
    dna: str, min_prot_len: int = 100
) -> Generator[tuple[str, int], None, None]:
    """
    Scans a DNA sequence and yields the longest distinct Open Reading Frames (ORFs).

    Biological Context:
    Mimics standard heuristic gene predictors by applying the "Longest ORF" rule.
    If multiple Start codons (ATG) share the exact same reading frame and Stop
    codon, only the longest continuous protein is yielded. Internal/nested
    Start codons are mathematically ignored to prevent data bloat.

    Architecture & Performance:
    Utilizes a Python set (`processed_stops`) to track spatial coordinates of
    utilized Stop codons in O(1) time. Stream processing via `yield` guarantees
    a flat, near-zero memory footprint even on multi-gigabyte chromosomes.

    Args:
        dna: A raw string of genomic DNA to be scanned.
        min_prot_len: The minimum amino acid length required to pass the
            quality control filter. Defaults to 100.

    Yields:
        A package containing:
            - The translated amino acid sequence.
            - The 1-based reading frame (+1, +2, or +3) where the gene resides.
    """
    seq_len = len(dna)
    start = dna.find("ATG")
    processed_stops = set()

    while start != -1:
        prot_seq = []

        for i in range(start, seq_len - 2, 3):
            codon = dna[i : i + 3]
            amino = DNA_CODON_TABLE.get(codon)

            if amino == "Stop":
                if i not in processed_stops:
                    if len(prot_seq) >= min_prot_len:
                        yield "".join(prot_seq), (start % 3) + 1
                        processed_stops.add(i)
                break

            elif amino is None:
                break

            else:
                prot_seq.append(amino)

        start = dna.find("ATG", start + 1)


def main() -> None:

    parser = base_parser("Open Reading Frames Pipeline")

    parser.add_argument(
        "-m",
        "--min_length",
        type=int,
        default=100,
        help="Minimum protein length in amino acids (default: 100)",
    )

    args = parser.parse_args()

    input_path = args.input
    output_path = args.output
    min_length = args.min_length

    if min_length < 1:
        sys.exit("Error: Minimum protein length must be at least 1 amino acid.")

    seq_processed = 0
    total_prot = 0

    try:
        with open(output_path, "w", encoding="utf-8") as out_file:

            for fasta_id, sequence in lazy_parse_fasta(input_path):
                seq_processed += 1
                orf_counter = 1

                for protein, frame in prot_frame_heuristic(
                    sequence, min_prot_len=min_length
                ):
                    header = (
                        f">{fasta_id}_ORF{orf_counter} Frame=+{frame} Strand=Forward"
                    )

                    out_file.write(f"{header}\n{wrap_fasta(protein)}\n")

                    orf_counter += 1
                    total_prot += 1

                dna_rev = reverse_complement(sequence)

                for protein, frame in prot_frame_heuristic(
                    dna_rev, min_prot_len=min_length
                ):
                    header = (
                        f">{fasta_id}_ORF{orf_counter} Frame=-{frame} Strand=Reverse"
                    )
                    out_file.write(f"{header}\n{wrap_fasta(protein)}\n")
                    orf_counter += 1
                    total_prot += 1

                print(f"Processed {fasta_id} | Found {orf_counter -1} ORFs")

        if seq_processed == 0:
            raise ValueError("Pipeline Halted: Input file contained no valid FASTA")

    except KeyboardInterrupt:
        output_path.unlink(missing_ok=True)
        sys.exit(
            "\nPipeline Halted: Scan interrupted by user. Partial output safely removed."
        )
    except (ValueError, FileNotFoundError, PermissionError) as e:
        output_path.unlink(missing_ok=True)
        sys.exit(f"\nSystem Exit: {e}")

    print("-" * 43)
    print(f"Success! {total_prot:,} total proteins extracted.")
    print(f"Results safely written to: {output_path.name}")


if __name__ == "__main__":
    main()
