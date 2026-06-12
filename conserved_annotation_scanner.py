"""
Extracts and groups CDS product annotations across GenBank genomes to identify conserved genes.

This tool acts as a text-based core proteome profiler. It aggregates all 'product'
qualifiers, normalizes them, and reports only the gene products that meet a
specified genome frequency threshold. It automatically filters out uninformative
'hypothetical protein' annotations by default.

Output matrices are strictly sorted in this order (least to most conserved):
1. Number of genomes found (Ascending)
2. Total physical copies found across all genomes (Ascending)
3. Alphabetical by product name (Ascending)

It can output both TSV matrices and matching FASTA sequence files.

License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).

Example Usage:
    # Standard run: Find genes conserved in at least 2 genomes, output TSV
    $ python3 conserved_annotation_scanner.py -i references/ --min_genomes 2 -o core.tsv

    # Auto-FASTA run: Will output both 'core.tsv' and 'core.fasta' automatically
    $ python3 conserved_annotation_scanner.py -i references/ --min_genomes 2 -o core.tsv -f

    # Exact run: Find genes present in EXACTLY 2 genomes
    $ python3 conserved_annotation_scanner.py -i references/ --min_genomes 2 --exact -o strict_two.tsv
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.0.2"

import sys
import argparse
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from utils import stream_reference_files, wrap_fasta


def extract_all_cds_features(file_path: Path) -> list[dict]:
    """Parses a GenBank file and extracts metadata for every CDS feature.

    Note:
        Parsing errors (OSError, ValueError) are caught internally and logged to stderr.
        Returns an empty list on failure to allow batch scanning to continue.

    Args:
        file_path: Path to the .gbk or .gbff file.

    Returns:
        A list of dictionaries containing locus tags, products, and sequences.
    """
    features = []
    try:
        if file_path.suffix.lower() not in (".gbk", ".gbff"):
            print(
                f"  [i] Skipping {file_path.name} (Not a GenBank format).",
                file=sys.stderr,
            )
            return features

        for record in SeqIO.parse(file_path, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    product = feature.qualifiers.get("product", [""])[0].strip()

                    if not product:
                        continue

                    features.append(
                        {
                            "original_product": product,
                            "normalized_key": product.lower(),
                            "locus_tag": feature.qualifiers.get(
                                "locus_tag", ["UNKNOWN"]
                            )[0],
                            "translation": feature.qualifiers.get(
                                "translation", [None]
                            )[0],
                            "locus": record.id,
                        }
                    )
    except (OSError, UnicodeDecodeError, ValueError) as e:
        print(f"  [!] Error reading {file_path.name}: {e}", file=sys.stderr)

    return features


def get_args() -> argparse.Namespace:
    """Configures CLI arguments."""
    parser = argparse.ArgumentParser(description="Conserved Annotation Scanner")
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Directory containing the GenBank reference files.",
    )
    parser.add_argument(
        "--min_genomes",
        type=int,
        default=2,
        help="Minimum number of genomes the product must appear in. Default: 2.",
    )
    parser.add_argument(
        "--exact",
        action="store_true",
        help="If flagged, restricts output to products present in EXACTLY the min_genomes value.",
    )
    parser.add_argument(
        "--keep_hypothetical",
        action="store_true",
        help="If flagged, includes 'hypothetical protein' in the results.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output TSV file to save the data matrix.",
    )
    parser.add_argument(
        "-f",
        "--fasta",
        action="store_true",
        help="If flagged, automatically generates a matching FASTA file alongside the TSV output.",
    )
    return parser.parse_args()


def main() -> None:
    """Main execution block for conserved annotation scanning."""
    args = get_args()

    if args.min_genomes < 1:
        sys.exit("[!] --min_genomes must be at least 1.")

    print(f"[*] Target Directory: {args.input.name}")
    condition_str = "==" if args.exact else ">="
    print(
        f"[*] Condition: Must be present in {condition_str} {args.min_genomes} genomes"
    )

    if not args.keep_hypothetical:
        print("[*] Filter: Ignoring 'hypothetical protein' annotations\n")

    master_results = defaultdict(lambda: defaultdict(list))
    scanned_files = 0

    try:
        for file_path in stream_reference_files(args.input):
            print(f"  -> Scanning {file_path.name}...")
            scanned_files += 1
            features = extract_all_cds_features(file_path)

            for feat in features:
                norm_key = feat["normalized_key"]
                if not args.keep_hypothetical and "hypothetical" in norm_key:
                    continue
                master_results[norm_key][file_path.name].append(feat)

        print(f"\n[*] Scan complete. Parsed {scanned_files} files.")
        print("[*] Aggregating and filtering data...")
        print("-" * 60)

        tsv_lines = [
            "Conserved_Product_Group\tGenomes_Found\tGenome_File\tLocus\tLocus_Tag\tOriginal_Product\tProtein_Sequence"
        ]
        fasta_lines = []
        conserved_count = 0

        # 3-Tier Sorting Logic:
        # Negative numbers sort descending (highest conservation first), strings sort ascending (A-Z)
        sorted_results = sorted(
            master_results.items(),
            key=lambda item: (
                len(item[1]),  # Tier 1: Genome count (Ascending)
                sum(
                    len(hits) for hits in item[1].values()
                ),  # Tier 2: Hit count (Ascending)
                item[0],  # Tier 3: Alphabetical (Ascending)
            ),
        )

        for norm_key, genomes_dict in sorted_results:
            genome_count = len(genomes_dict)
            meets_condition = (
                (genome_count == args.min_genomes)
                if args.exact
                else (genome_count >= args.min_genomes)
            )

            if meets_condition:
                conserved_count += 1
                for genome_name, hits in genomes_dict.items():
                    for hit in hits:
                        tsv_lines.append(
                            f"{norm_key}\t{genome_count}\t{genome_name}\t{hit['locus']}\t"
                            f"{hit['locus_tag']}\t{hit['original_product']}\t{hit['translation']}"
                        )
                        if args.fasta and hit["translation"] is not None:
                            header = f">{genome_name}|{hit['locus']}|{hit['locus_tag']}|{hit['original_product']}"
                            seq_wrapped = wrap_fasta(hit["translation"])
                            fasta_lines.append(f"{header}\n{seq_wrapped}")

        if conserved_count == 0:
            print(f"[!] No functional annotations met the threshold criteria.")
            print(
                f"[-] Output file {args.output.name} was not created to prevent empty datasets."
            )
        else:
            # utf-8-sig adds a BOM, enabling auto-detection in Excel without a manual import step
            with open(args.output, "w", encoding="utf-8-sig") as out_tsv:
                out_tsv.write("\n".join(tsv_lines) + "\n")

            print(
                f"[*] Success! {conserved_count} distinct functional groups met the threshold."
            )
            print(f"[*] TSV matrix written to: {args.output.name}")

            if args.fasta and fasta_lines:
                fasta_path = args.output.with_suffix(".fasta")
                # FASTA must use standard utf-8, as downstream aligners will fail if a BOM is present
                with open(fasta_path, "w", encoding="utf-8") as out_fasta:
                    out_fasta.write("\n".join(fasta_lines) + "\n")
                print(f"[*] FASTA sequences written to: {fasta_path.name}")

    except KeyboardInterrupt:
        sys.exit("\n[!] Scan interrupted by user.")


if __name__ == "__main__":
    main()
