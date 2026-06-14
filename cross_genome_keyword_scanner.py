"""
Scans multiple GenBank reference genomes for specific annotation keywords.

Filters and returns only the genes (and their protein sequences) that are
conserved across a minimum number of different genomes. This is used to
empirically verify the presence of specific protein families (e.g., 'lactobin',
'cerein') across a comparative genomic dataset.

Output matrices are strictly sorted in this order (least to most conserved):
1. Number of genomes found (Ascending)
2. Total physical copies found across all genomes (Ascending)
3. Alphabetical by keyword/product name (Ascending)

It supports exporting results as both a tabular TSV matrix and a sequence-ready FASTA file.

License: MIT

Reproducibility:
    Associated with upcoming research (manuscript in preparation).
    Correct attribution is requested when used in derivative works.
    See LICENSE in the repository root for full details.

Example Usage:
    # Standard run: Search for 'bacteriocin', output TSV
    $ python3 cross_genome_keyword_scanner.py -i references/ -k bacteriocin --min_genomes 3 -o core_annotations.tsv

    # Auto-FASTA run: Search for multiple keywords, output TSV and matching FASTA
    $ python3 cross_genome_keyword_scanner.py -i references/ -k cerein lactobin -o conserved_hits.tsv -f

    # Exact run: Search for keywords present in EXACTLY 3 genomes
    $ python3 cross_genome_keyword_scanner.py -i references/ -k bacteriocin --min_genomes 3 --exact -o strict_three.tsv
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.1.0"

import sys
import argparse
from pathlib import Path
from collections import defaultdict

try:
    from Bio import SeqIO
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )
from utils import stream_reference_files, wrap_fasta

try:
    from tqdm import tqdm

    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    tqdm = lambda x, **kwargs: x


def scan_file_for_keywords(
    file_path: Path, keywords: list[str]
) -> dict[str, list[dict]]:
    """Parses a GenBank file to find CDS features matching target keywords.

    Matching is case-insensitive substring search across the product, gene,
    and note qualifiers combined. A single CDS can match multiple keywords
    independently.

    Note:
        Substring matching means short keywords may over-match. For example,
        searching 'cin' will match 'bacteriocin', 'nisin', 'colistin', etc.
        Use specific keywords (e.g., 'bacteriocin' not 'cin') for clean results.

        Parsing errors (OSError, ValueError) are caught internally and logged
        to stderr, allowing batch scanning to continue.

    Args:
        file_path: Path to the .gbk or .gbff file.
        keywords:  A list of string keywords to search for (case-insensitive).

    Returns:
        A dictionary mapping each matched keyword to a list of feature metadata.
    """
    lower_keywords = [k.lower() for k in keywords]
    hits = defaultdict(list)

    try:
        if file_path.suffix.lower() not in (".gbk", ".gbff"):
            print(
                f"  [i] Skipping {file_path.name} (Not a GenBank format).",
                file=sys.stderr,
            )
            return hits

        for record in SeqIO.parse(file_path, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    product = feature.qualifiers.get("product", [""])[0]
                    gene = feature.qualifiers.get("gene", [""])[0]
                    note = feature.qualifiers.get("note", [""])[0]

                    searchable_text = f"{product} {gene} {note}".lower()

                    for kw in lower_keywords:
                        if kw in searchable_text:
                            locus_tag = feature.qualifiers.get(
                                "locus_tag", ["UNKNOWN"]
                            )[0]
                            translation = feature.qualifiers.get("translation", [None])[
                                0
                            ]

                            hits[kw].append(
                                {
                                    "locus_tag": locus_tag,
                                    "product": product,
                                    "translation": translation,
                                    "locus": record.id,
                                }
                            )
                            # A gene matching multiple keywords is added to each keyword's result list independently

    except (OSError, UnicodeDecodeError, ValueError) as e:
        print(f"  [!] Error reading {file_path.name}: {e}", file=sys.stderr)

    return hits


def get_args() -> argparse.Namespace:
    """Configures the CLI and returns parsed arguments."""
    parser = argparse.ArgumentParser(
        description="Cross-Genome Keyword Conservation Scanner"
    )
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Directory containing the GenBank reference files.",
    )
    parser.add_argument(
        "-k",
        "--keywords",
        type=str,
        nargs="+",
        required=True,
        help="One or more keywords to search for (e.g., cerein lactobin).",
    )
    parser.add_argument(
        "--min_genomes",
        type=int,
        default=2,
        help="Minimum number of genomes a keyword must appear in to be reported. Default: 2.",
    )
    parser.add_argument(
        "--exact",
        action="store_true",
        help="If flagged, restricts output to products present in EXACTLY the min_genomes value.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=False,
        help="Output TSV file to save the results.",
    )
    parser.add_argument(
        "-f",
        "--fasta",
        action="store_true",
        help="If flagged, automatically generates a matching FASTA file alongside the TSV output.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help=(
            "Increase output verbosity. By default, only major milestones and "
            "a progress bar are shown. With --verbose, every file being scanned "
            "is printed. Useful for debugging."
        ),
    )
    return parser.parse_args()


def main() -> None:
    """Main execution block for cross-genome keyword scanning."""
    args = get_args()

    if args.min_genomes < 1:
        sys.exit("[!] --min_genomes must be at least 1.")

    print(f"[*] Target Directory: {args.input.name}", file=sys.stderr)
    print(f"[*] Keywords to hunt: {', '.join(args.keywords)}", file=sys.stderr)
    condition_str = "==" if args.exact else ">="
    print(
        f"[*] Condition: Must be present in {condition_str} {args.min_genomes} genomes\n",
        file=sys.stderr,
    )

    master_results = defaultdict(lambda: defaultdict(list))
    scanned_files = 0

    try:
        ref_files = list(stream_reference_files(args.input))
        ref_iter = tqdm(
            ref_files,
            desc="Scanning genomes",
            disable=not HAS_TQDM or args.verbose,
        )

        for file_path in ref_iter:
            if args.verbose:
                print(f"  -> Scanning {file_path.name}...", file=sys.stderr)
            scanned_files += 1
            file_hits = scan_file_for_keywords(file_path, args.keywords)

            for kw, hit_list in file_hits.items():
                master_results[kw][file_path.name].extend(hit_list)

        print(f"\n[*] Scan complete. Parsed {scanned_files} files.", file=sys.stderr)
        print("-" * 60, file=sys.stderr)

        output_lines = [
            "Keyword\tGenomes_Found\tGenome_File\tLocus\tLocus_Tag\tProduct\tProtein_Sequence",
        ]
        fasta_lines = []
        valid_results = {}
        failed_keywords = []

        for kw in args.keywords:
            kw_lower = kw.lower()
            genomes_with_kw = master_results.get(kw_lower, {})
            genome_count = len(genomes_with_kw)

            meets_condition = (
                (genome_count == args.min_genomes)
                if args.exact
                else (genome_count >= args.min_genomes)
            )

            if meets_condition:
                valid_results[kw] = genomes_with_kw
            else:
                failed_keywords.append((kw, genome_count))

        found_anything = len(valid_results) > 0

        if found_anything:
            # All tiers sort ascending (lowest to highest, A-Z)
            # Negative numbers sort descending (highest conservation first), strings sort ascending (A-Z)
            sorted_results = sorted(
                valid_results.items(),
                key=lambda item: (
                    len(item[1]),  # Tier 1: Genome count (Ascending)
                    sum(
                        len(hits) for hits in item[1].values()
                    ),  # Tier 2: Total hit count (Ascending)
                    item[0].lower(),  # Tier 3: Alphabetical keyword (Ascending)
                ),
            )

            for kw, genomes_with_kw in sorted_results:
                print(
                    f"[+] Keyword '{kw}' is CONSERVED across {len(genomes_with_kw)} genomes.",
                    file=sys.stderr,
                )
                for genome_name, hits in genomes_with_kw.items():
                    for hit in hits:
                        output_lines.append(
                            f"{kw}\t{len(genomes_with_kw)}\t{genome_name}\t{hit['locus']}\t{hit['locus_tag']}\t{hit['product']}\t{hit['translation']}"
                        )
                        if args.fasta and hit["translation"] is not None:
                            header = f">{genome_name}|{hit['locus']}|{hit['locus_tag']}|{hit['product']}"
                            seq_wrapped = wrap_fasta(hit["translation"])
                            fasta_lines.append(f"{header}\n{seq_wrapped}")

        for kw, genome_count in failed_keywords:
            print(
                f"[-] Keyword '{kw}' failed threshold (Found in {genome_count} genomes).",
                file=sys.stderr,
            )

        print("-" * 60, file=sys.stderr)

        if not found_anything:
            print(
                "[!] No keywords met the minimum genome threshold. No output written.",
                file=sys.stderr,
            )
        elif args.output:
            try:
                # utf-8-sig adds a BOM, enabling auto-detection in Excel
                with open(args.output, "w", encoding="utf-8-sig") as f:
                    f.write("\n".join(output_lines) + "\n")
                print(
                    f"[*] Success! Matrix written to {args.output.resolve()}",
                    file=sys.stderr,
                )

                if args.fasta and fasta_lines:
                    fasta_path = args.output.with_suffix(".fasta")
                    # FASTA must use standard utf-8 for compatibility with downstream tools
                    with open(fasta_path, "w", encoding="utf-8") as f_fasta:
                        f_fasta.write("\n".join(fasta_lines) + "\n")
                    print(
                        f"[*] FASTA sequences written to {fasta_path.name}",
                        file=sys.stderr,
                    )

            except OSError as e:
                sys.exit(f"\n[!] Error writing file: {e}")
        else:
            # If no output file is provided but data exists, dump TSV to stdout
            print(
                "\n[!] Note: No output file specified (-o). Printing to stdout:\n",
                file=sys.stderr,
            )
            sys.stdout.write("\n".join(output_lines) + "\n")

    except KeyboardInterrupt:
        sys.exit("\n[!] Scan interrupted by user.")


if __name__ == "__main__":
    main()
