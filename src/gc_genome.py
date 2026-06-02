"""
GC Content & Assembly QC Pipeline

Calculates per-contig GC percentage and the statistically true whole-genome
GC content of a multi-FASTA assembly. Outputs results to the terminal and
saves a TSV file suitable for downstream contamination checks (e.g., GC vs Length).

Usage:
    $ python3 gc_genome.py -i input_assembly.fasta -o gc_results.tsv
"""

__author__ = Jan Ephraim R. Vallente
__email__ = ephrvallente@gmail.com
__version__ = 1.0.0

import sys
from stream_nuc_count import stream_nuc_count
from utils import base_parser


def main() -> None:
    """
    Command-line interface for calculating GC content.

    Parses arguments, orchestrates the memory-safe stream_nuc_count generator,
    aggregates global genome statistics, and writes a formatted TSV report 
    along with a summary log file.
    """
    args = base_parser("GC Content Pipeline for Genome Assemblies").parse_args()
    input_path = args.input
    output_path = args.output

    # Format output
    highest_gc_id = ""
    highest_gc_value = -1.0
    sequences_processed = 0
    global_g = 0
    global_c = 0
    global_valid = 0

    print(("-" * 13) + (" GC CONTENT PIPELINE RUNNING ") + ("-" * 13))

    try:
        with open(output_path, "w", encoding="utf-8") as out_file:
            out_file.write(("Sequence_ID\tGC_Content_Percent\tValid_Length_bp\n"))

            for fasta_id, counts in stream_nuc_count(input_path):
                sequences_processed += 1
                local_valid = counts["A"] + counts["C"] + counts["G"] + counts["T"]

                if local_valid > 0:
                    local_gc = ((counts["G"] + counts["C"]) / local_valid) * 100.0
                else:
                    local_gc = 0.0

                if local_gc > highest_gc_value:
                    highest_gc_value = local_gc
                    highest_gc_id = fasta_id

                # Terminal progress
                print(
                    f"> {fasta_id} | GC: {local_gc:.6f}% | Length: {local_valid:,} bp"
                )

                # --- GLOBAL MATH  ---
                global_g += counts["G"]
                global_c += counts["C"]
                global_valid += local_valid

                out_file.write(f"{fasta_id}\t{local_gc:.6f}\t{local_valid}\n")

        if sequences_processed == 0:
            raise ValueError(
                "Pipeline Halted: The input file contained no valid FASTA records."
            )

        if global_valid == 0:
            raise ValueError(
                "Pipeline Halted: No valid ACGT bases found across all sequences."
            )

    except (ValueError, PermissionError) as e:
        output_path.unlink(missing_ok=True)
        sys.exit(f"Pipeline Halted: {e}")

    true_genome_gc = ((global_g + global_c) / global_valid) * 100.0

    log_path = output_path.with_name(f"{output_path.stem}_summary.txt")
    log_content = (
        "----------- FINAL RESULTS --------------\n"
        f"Highest GC Content:        {highest_gc_value:.6f}\n"
        f"(ID: {highest_gc_id})\n\n"
        "========================================\n"
        "      FINAL ASSEMBLY STATISTICS\n"
        "========================================\n"
        f"Total Contigs Processed:  {sequences_processed}\n"
        f"Total Usable Genome Size: {global_valid:,} bp\n"
        f"True Whole Genome GC:     {true_genome_gc:.6f}%\n"
        f"Highest GC Content:       {highest_gc_value:.6f}%\n"
        f"(ID: {highest_gc_id})\n"
        "========================================"
    )

    log_path.write_text(log_content, encoding="utf-8")

    print("\n" + ("-" * 11) + " FINAL RESULTS " + ("-") * 14)
    print(f"Highest GC Content:        {highest_gc_value:.6f}% \n(ID: {highest_gc_id})")

    print("\n" + "=" * 40)
    print("      FINAL ASSEMBLY STATISTICS")
    print("=" * 40)
    print(f"Total Contigs Processed:  {sequences_processed}")
    print(f"Total Usable Genome Size: {global_valid:,} bp")
    print(f"True Whole Genome GC:     {true_genome_gc:.6f}%")
    print("=" * 40)
    print(f"\nMachine-readable matrix safely written to: {output_path.name}")
    print(f"Log files safely written to: {output_path.stem}_summary.txt")


if __name__ == "__main__":
    main()
