"""
antiSMASH BGC Region Explorer

Parses antiSMASH .region.gbk output files for genomic characterization.

This pipeline reads GenBank files produced by antiSMASH. It extracts gene
annotations, summarizes biosynthetic roles, and exports functional manifests
for secondary metabolite clusters.

Author: Jan Ephraim R. Vallente (ephrvallente@gmail.com)
Date: 2026-06-07
License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).

Example Usage:
    $ python3 bgc_explorer.py -i C5_gnlProkkaNHJNNNGJ_1.region001.gbk
    $ python3 bgc_explorer.py -i C5_prokka_result.gbk -o C5_manifest.tsv
"""

__version__ = "1.0.1"
import sys
import re
from pathlib import Path
from collections import Counter
from Bio import SeqIO
from utils import base_parser

# Compile regex at the module level for performance.
# Targets the final closing parenthesis and captures everything after it.
# Example: "biosynthetic (rule-based-clusters) RiPP-like: Bacteriocin_IIc" -> "RiPP-like: Bacteriocin_IIc"
FUNC_PATTERN = re.compile(r"\)\s*([^)]+)$")


def parse_bgc_region(gbk_path: Path) -> tuple[str, list[tuple]]:
    """
    Parses an antiSMASH v8 region file to extract cluster metadata and CDS features.

    Returns:
        A tuple containing:
        - cluster_type (str): The product type of the BGC (from the 'region' feature).
        - cds_results (list): A list of tuples containing CDS data:
          (locus_tag, gene, start, end, strand, gene_kind, role_str, product)

    Raises:
        ValueError: If Biopython or file IO explicitly fails to parse the file.
    """
    cluster_type = "Unknown Cluster Type"
    cds_results = []

    try:
        for record in SeqIO.parse(gbk_path, "genbank"):
            for feature in record.features:

                # Extract Top-Level Cluster Metadata
                if feature.type == "region":
                    cluster_type = feature.qualifiers.get("product", [cluster_type])[0]

                # Extract CDS Details
                elif feature.type == "CDS":
                    locus_tag = feature.qualifiers.get("locus_tag", ["UNKNOWN"])[0]
                    gene = feature.qualifiers.get("gene", ["-"])[0]
                    product = feature.qualifiers.get(
                        "product", ["Hypothetical protein"]
                    )[0]

                    # 1-based start coordinates for publication standards
                    start = int(feature.location.start) + 1
                    end = int(feature.location.end)
                    strand = "+" if feature.location.strand == 1 else "-"

                    gene_kind = feature.qualifiers.get("gene_kind", ["unassigned"])[0]

                    extra_detail = ""
                    gene_funcs = feature.qualifiers.get("gene_functions", [])
                    if gene_funcs:
                        raw_func = gene_funcs[0]
                        match = FUNC_PATTERN.search(raw_func)
                        if match:
                            extra_detail = match.group(1).strip()
                        else:
                            extra_detail = raw_func.strip()

                    if extra_detail:
                        role_str = f"{gene_kind} [{extra_detail}]"
                    else:
                        role_str = gene_kind

                    cds_results.append(
                        (
                            locus_tag,
                            gene,
                            start,
                            end,
                            strand,
                            gene_kind,
                            role_str,
                            product,
                        )
                    )

        return cluster_type, cds_results

    # Catch ONLY file-system and strict parsing errors, not our own bugs
    except (FileNotFoundError, OSError, UnicodeDecodeError, ValueError) as e:
        raise ValueError(
            f"Failed to parse BGC GenBank file '{gbk_path.name}': {e}"
        ) from e


def main() -> None:
    """CLI Entry point for the BGC Explorer."""
    parser = base_parser(
        description_text="antiSMASH v8 BGC Explorer", include_output=False
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Optional: Path to save the extracted manifest as a TSV file.",
    )

    args = parser.parse_args()

    try:
        cluster_type, results = parse_bgc_region(args.input)

        if not results:
            raise ValueError(f"No coding sequences found in {args.input.name}.")

        print(f"\n--- BGC FACTORY MANIFEST: {args.input.name} ---")
        print(f"[*] Cluster Type: {cluster_type}")
        print("-" * 135)
        print(
            f"{'':<4}{'Locus Tag':<17} | {'Gene':<6} | {'Coordinates':<18} | {'Role in BGC':<45} | {'Product'}"
        )
        print("-" * 135)

        for locus, gene, start, end, strand, kind, role, product in results:

            # Visually highlight the core weapons and cannons
            if "biosynthetic" in role.lower():
                highlight = ">>> "
            elif "transport" in role.lower():
                highlight = " -> "
            else:
                highlight = "    "

            coords = f"{start}-{end} ({strand})"

            # Truncate role and product for terminal neatness
            short_role = role[:42] + "..." if len(role) > 42 else role
            short_product = product[:35] + "..." if len(product) > 35 else product

            print(
                f"{highlight:<4}{locus:<17} | {gene:<6} | {coords:<18} | {short_role:<45} | {short_product}"
            )

        print("-" * 135)

        # Summary Counter
        counts = Counter(cds[5] for cds in results)  # Index 5 is gene_kind
        summary_str = " | ".join(
            f"{kind.capitalize()}: {count}" for kind, count in counts.items()
        )

        print(f"[*] Total Genes : {len(results)}")
        print(f"[*] Category    : {summary_str}")

        # TSV Export Logic with strict IO error handling
        if args.output:
            try:
                with open(args.output, "w", encoding="utf-8-sig") as tsv_out:
                    tsv_out.write(
                        "Locus_Tag\tGene\tStart\tEnd\tStrand\tRole_in_BGC\tProduct\n"
                    )
                    for locus, gene, start, end, strand, kind, role, product in results:
                        tsv_out.write(
                            f"{locus}\t{gene}\t{start}\t{end}\t{strand}\t{role}\t{product}\n"
                        )
                print(f"[*] Successfully exported TSV file to: {args.output.resolve()}")
            except OSError as e:
                sys.exit(f"\n[!] Could not write TSV to '{args.output.resolve()}': {e}")

        print()

    except ValueError as e:
        sys.exit(f"\n[!] Pipeline Error: {e}")
    except KeyboardInterrupt:
        sys.exit("\n[!] Pipeline interrupted by user.")


if __name__ == "__main__":
    main()
