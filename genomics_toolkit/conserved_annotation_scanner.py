#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

"""Conserved Annotation Scanner

Extracts and groups CDS product annotations across GenBank genomes to identify
conserved genes. Acts as a text-based core proteome profiler: aggregates all
/product qualifiers, normalizes them, and reports only gene products that meet
a specified genome frequency threshold. Filters uninformative 'hypothetical
protein' annotations by default.

ANNOTATION NORMALIZATION:
    Raw /product strings are normalized before grouping using
    _normalize_product():
    - Lowercasing
    - Stripping common qualifier noise words ("putative", "probable",
      "predicted", "possible", "potential", "uncharacterized", "candidate")
    - Normalizing punctuation (hyphens, underscores, slashes) to spaces
    - Removing trailing bracketed organism specs (e.g., "[E. coli K-12]")
    - Removing trailing parenthetical specs (e.g., "(plasmid)")

    This handles the most common cross-pipeline annotation inconsistencies.
    NOTE: Abbreviation-vs-full-name discrepancies ("atpA" vs "ATP synthase
    subunit alpha") cannot be resolved by text normalization. They require
    sequence-based clustering via gbk_ortholog_finder.py.

OUTPUT SORT ORDER (most to least conserved):
    1. Number of genomes found        (Descending — most conserved first)
    2. Total physical copies          (Descending — single-copy before multicopy)
    3. Alphabetical by product name   (Ascending)

    Descending Tier 1 ensures core genes appear at the TOP of the output
    matrix. Descending Tier 2 prevents copy-number inflation (a gene with
    12 copies in one genome ranks BELOW a perfectly distributed single-copy
    universal gene with the same genome count).

PSEUDOGENE HANDLING:
    CDS features without a /translation qualifier (pseudogenes, frameshifted
    genes) are stored with an empty string rather than None. The TSV
    Protein_Sequence column will be blank for these entries, preventing the
    literal string "None" from being written and misinterpreted as a peptide
    sequence by downstream tools.

FILE KEY SAFETY:
    Genomes are tracked by their path relative to the input directory, not
    by filename alone. This prevents two files named genome.gbk in different
    subdirectories (e.g., wild_type/genome.gbk and mutant/genome.gbk) from
    being merged into the same tracking entry.

Note:
    Associated with ongoing, unpublished research (manuscript in
    preparation). Correct attribution is requested when used in
    derivative works.

Examples:
    # Standard run: Find genes conserved in at least 2 genomes, output TSV
    $ python3 conserved_annotation_scanner.py -i references/ --min_genomes 2 -o core.tsv

    # Auto-FASTA run: Output both 'core.tsv' and 'core.fasta'
    $ python3 conserved_annotation_scanner.py -i references/ --min_genomes 2 -o core.tsv -f

    # Exact run: Find genes present in EXACTLY 2 genomes
    $ python3 conserved_annotation_scanner.py -i references/ --min_genomes 2 --exact -o strict_two.tsv
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.2.0"

import re
import sys
import argparse
from contextlib import ExitStack
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


# ── Product name normalization ─────────────────────────────────────────────────

_NOISE_WORDS = frozenset(
    [
        "putative",
        "probable",
        "predicted",
        "possible",
        "potential",
        "uncharacterized",
        "candidate",
    ]
)


def _normalize_product(product: str) -> str:
    """Normalize a /product annotation string for cross-pipeline grouping.

    Handles the most common inconsistencies introduced by different annotation
    pipelines (Prokka, Bakta, NCBI RefSeq, RAST):
    - Lowercasing
    - Stripping noise qualifier words ("putative", "probable", etc.)
    - Normalizing punctuation (hyphens, underscores, slashes) to spaces
    - Removing trailing bracketed organism specs (e.g., "[Lactobacillus sp.]")
    - Removing trailing parenthetical details (e.g., "(plasmid)")

    Limitation: Cannot resolve abbreviation-vs-full-name discrepancies
    (e.g., "atpA" vs "ATP synthase subunit alpha"). These require
    sequence-based clustering. See gbk_ortholog_finder.py.

    Args:
        product: Raw /product qualifier string from a GenBank CDS feature.

    Returns:
        Normalized, lowercase string suitable for grouping.
    """
    name = product.lower().strip()
    # Remove trailing bracketed organism specs
    name = re.sub(r"\s*\[.*?\]\s*$", "", name)
    # Remove trailing parenthetical specs
    name = re.sub(r"\s*\(.*?\)\s*$", "", name)
    # Normalize punctuation to spaces
    name = re.sub(r"[-_/,;]+", " ", name)
    # Collapse multiple spaces
    name = re.sub(r"\s+", " ", name).strip()
    # Strip noise qualifier words
    words = [w for w in name.split() if w not in _NOISE_WORDS]
    return " ".join(words).strip()


# ── GenBank parsing ────────────────────────────────────────────────────────────


def extract_all_cds_features(file_path: Path) -> list[dict]:
    """Parses a GenBank file and extracts metadata for every CDS feature.

    Note:
        Parsing errors (OSError, ValueError) are caught internally and logged
        to stderr. Returns an empty list on failure to allow batch scanning
        to continue.

        CDS features without a /translation qualifier (pseudogenes, truncated
        genes) are stored with an empty string in the 'translation' field
        rather than None, preventing the literal string "None" from appearing
        in TSV output and being misinterpreted by downstream tools.

    Args:
        file_path: Path to the .gbk or .gbff file.

    Returns:
        A list of feature metadata dicts.
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

                    # Store "" instead of None to prevent "None" appearing in TSV
                    raw_translation = feature.qualifiers.get("translation", [None])[0]
                    translation = raw_translation if raw_translation is not None else ""

                    features.append(
                        {
                            "original_product": product,
                            "normalized_key": _normalize_product(product),
                            "locus_tag": feature.qualifiers.get(
                                "locus_tag", ["UNKNOWN"]
                            )[0],
                            "translation": translation,
                            "locus": record.id,
                        }
                    )
    except (OSError, UnicodeDecodeError, ValueError) as e:
        print(f"  [!] Error reading {file_path.name}: {e}", file=sys.stderr)

    return features


# ── CLI ────────────────────────────────────────────────────────────────────────


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
        help="Restrict output to products present in EXACTLY the min_genomes value.",
    )
    parser.add_argument(
        "--keep_hypothetical",
        action="store_true",
        help="Include 'hypothetical protein' in the results.",
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
        help="Also generate a matching FASTA file alongside the TSV output.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help=(
            "Print every file being scanned. By default only major milestones "
            "and a progress bar are shown."
        ),
    )
    return parser.parse_args()


# ── Main ───────────────────────────────────────────────────────────────────────


def main() -> None:
    """Main execution block for conserved annotation scanning."""
    args = get_args()

    if args.min_genomes < 1:
        sys.exit("[!] --min_genomes must be at least 1.")

    condition_str = "==" if args.exact else ">="
    print(f"[*] Target Directory : {args.input.name}", file=sys.stderr)
    print(
        f"[*] Condition        : Present in {condition_str} {args.min_genomes} genomes",
        file=sys.stderr,
    )
    if not args.keep_hypothetical:
        print(
            "[*] Filter           : Ignoring 'hypothetical protein' annotations",
            file=sys.stderr,
        )
    print("", file=sys.stderr)

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

            # Use relative path as genome key to prevent namespace collisions
            # when two files share the same filename in different subdirectories.
            try:
                genome_key = str(file_path.relative_to(args.input))
            except ValueError:
                genome_key = file_path.name

            features = extract_all_cds_features(file_path)

            for feat in features:
                norm_key = feat["normalized_key"]
                if not args.keep_hypothetical and "hypothetical" in norm_key:
                    continue
                master_results[norm_key][genome_key].append(feat)

        print(f"\n[*] Scan complete. Parsed {scanned_files} files.", file=sys.stderr)
        print("[*] Aggregating and filtering data...", file=sys.stderr)
        print("-" * 60, file=sys.stderr)

        # Filter first, then sort — avoids sorting items that will be discarded
        filtered_results = {
            norm_key: genomes_dict
            for norm_key, genomes_dict in master_results.items()
            if (
                (len(genomes_dict) == args.min_genomes)
                if args.exact
                else (len(genomes_dict) >= args.min_genomes)
            )
        }

        # 3-Tier sort: most conserved first
        #   Tier 1: genome count          DESCENDING (most conserved at top)
        #   Tier 2: total physical copies DESCENDING (single-copy before multicopy)
        #   Tier 3: alphabetical          ASCENDING
        # Negating Tiers 1 and 2 achieves descending order without reverse=True,
        # while keeping Tier 3 in natural ascending order.
        sorted_results = sorted(
            filtered_results.items(),
            key=lambda item: (
                -len(item[1]),
                -sum(len(hits) for hits in item[1].values()),
                item[0],
            ),
        )

        if not sorted_results:
            print(
                "[!] No functional annotations met the threshold criteria.",
                file=sys.stderr,
            )
            print(
                f"[-] Output file {args.output.name} was not created to prevent empty datasets.",
                file=sys.stderr,
            )
            return

        conserved_count = 0
        TSV_HEADER = (
            "Conserved_Product_Group\tGenomes_Found\tGenome_File\t"
            "Locus\tLocus_Tag\tOriginal_Product\tProtein_Sequence"
        )
        fasta_path = args.output.with_suffix(".fasta") if args.fasta else None

        # Stream-write rows directly rather than accumulating strings in RAM.
        # For large eukaryotic or metagenomic datasets, list accumulation +
        # join() can spike memory significantly. Writing row-by-row keeps
        # memory usage flat regardless of dataset size.
        with ExitStack() as stack:
            out_tsv = stack.enter_context(open(args.output, "w", encoding="utf-8-sig"))
            out_fasta = (
                stack.enter_context(open(fasta_path, "w", encoding="utf-8"))
                if fasta_path
                else None
            )

            out_tsv.write(TSV_HEADER + "\n")

            for norm_key, genomes_dict in sorted_results:
                genome_count = len(genomes_dict)
                conserved_count += 1

                for genome_name, hits in genomes_dict.items():
                    for hit in hits:
                        out_tsv.write(
                            f"{norm_key}\t{genome_count}\t{genome_name}\t"
                            f"{hit['locus']}\t{hit['locus_tag']}\t"
                            f"{hit['original_product']}\t{hit['translation']}\n"
                        )
                        if out_fasta and hit["translation"]:
                            header = (
                                f">{genome_name}|{hit['locus']}|"
                                f"{hit['locus_tag']}|{hit['original_product']}"
                            )
                            seq_wrapped = wrap_fasta(hit["translation"])
                            out_fasta.write(f"{header}\n{seq_wrapped}\n")

        print(
            f"[*] Success! {conserved_count} distinct functional groups met the threshold.",
            file=sys.stderr,
        )
        print(
            f"[*] TSV matrix written to : {args.output.resolve()}",
            file=sys.stderr,
        )
        if fasta_path:
            print(
                f"[*] FASTA sequences written to : {fasta_path.resolve()}",
                file=sys.stderr,
            )

    except KeyboardInterrupt:
        sys.exit("\n[!] Scan interrupted by user.")


if __name__ == "__main__":
    main()
