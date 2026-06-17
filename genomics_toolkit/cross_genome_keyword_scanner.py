#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

"""Cross-Genome Keyword Conservation Scanner

Scans multiple GenBank reference genomes for specific annotation keywords.
Filters and returns only the genes (and their protein sequences) that are
conserved across a minimum number of different genomes. Used to empirically
verify the presence of specific protein families (e.g., 'lactobin', 'cerein')
across a comparative genomic dataset.

Supports exporting results as both a tabular TSV matrix and a FASTA file.

OUTPUT SORT ORDER (most to least conserved):
    1. Number of genomes found        (Descending — most conserved first)
    2. Total physical copies          (Descending — single-copy before multicopy)
    3. Alphabetical by keyword        (Ascending)

    Descending Tier 1 ensures highly conserved hits appear at the TOP of the
    output matrix. Descending Tier 2 prevents copy-number inflation artifacts
    (a gene duplicated 12 times in one genome ranks BELOW a truly universal
    single-copy gene with the same genome count).

FASTA HEADER FORMAT:
    >{keyword}|{genome_file}|{contig_id}|{locus_tag}|{product}

    The keyword is prepended to every FASTA header. This guarantees unique
    headers when the same gene matches multiple keywords (e.g., a "nisin-
    controlled class II bacteriocin" matching both "nisin" and "bacteriocin").
    Without keyword prefixing, downstream alignment tools (MAFFT, Clustal
    Omega, IQ-TREE) would receive identical non-unique record headers and
    either crash or produce corrupt output.

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
    being merged into the same tracking entry, which would produce silently
    wrong conservation counts.

Note:
    Associated with ongoing, unpublished research (manuscript in
    preparation). Correct attribution is requested when used in
    derivative works.

Examples:
    # Standard run: Search for 'bacteriocin', output TSV
    $ python3 cross_genome_keyword_scanner.py -i references/ -k bacteriocin --min_genomes 3 -o core_annotations.tsv

    # Auto-FASTA run: Search for multiple keywords, output TSV and matching FASTA
    $ python3 cross_genome_keyword_scanner.py -i references/ -k cerein lactobin -o conserved_hits.tsv -f

    # Exact run: Search for keywords present in EXACTLY 3 genomes
    $ python3 cross_genome_keyword_scanner.py -i references/ -k bacteriocin --min_genomes 3 --exact -o strict_three.tsv
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.2.0"

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


# ── GenBank parsing ────────────────────────────────────────────────────────────


def scan_file_for_keywords(
    file_path: Path, keywords: list[str]
) -> dict[str, list[dict]]:
    """Parses a GenBank file to find CDS features matching target keywords.

    Matching is case-insensitive substring search across the product, gene,
    and note qualifiers combined. A single CDS can match multiple keywords
    independently — it will appear in each matching keyword's result list.

    Note:
        Substring matching means short keywords may over-match. For example,
        searching 'cin' will match 'bacteriocin', 'nisin', 'colistin', etc.
        Use specific keywords (e.g., 'bacteriocin' not 'cin') for clean results.

        CDS features without a /translation qualifier (pseudogenes) are stored
        with an empty string rather than None to prevent "None" appearing as a
        literal value in TSV output.

        Parsing errors (OSError, ValueError) are caught internally and logged
        to stderr, allowing batch scanning to continue.

    Args:
        file_path: Path to the .gbk or .gbff file.
        keywords:  A list of string keywords to search for (case-insensitive).

    Returns:
        A dictionary mapping each matched keyword to a list of feature metadata.
    """
    lower_keywords = [k.lower() for k in keywords]
    hits: dict[str, list[dict]] = defaultdict(list)

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

                            # Store "" instead of None to prevent "None" in TSV
                            raw_translation = feature.qualifiers.get(
                                "translation", [None]
                            )[0]
                            translation = (
                                raw_translation if raw_translation is not None else ""
                            )

                            hits[kw].append(
                                {
                                    "locus_tag": locus_tag,
                                    "product": product,
                                    "translation": translation,
                                    "locus": record.id,
                                }
                            )

    except (OSError, UnicodeDecodeError, ValueError) as e:
        print(f"  [!] Error reading {file_path.name}: {e}", file=sys.stderr)

    return hits


# ── CLI ────────────────────────────────────────────────────────────────────────


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
        help="Restrict output to keywords present in EXACTLY the min_genomes value.",
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

            # Use relative path as genome key to prevent namespace collisions
            # when two files share the same filename in different subdirectories.
            try:
                genome_key = str(file_path.relative_to(args.input))
            except ValueError:
                genome_key = file_path.name

            file_hits = scan_file_for_keywords(file_path, args.keywords)
            for kw, hit_list in file_hits.items():
                master_results[kw][genome_key].extend(hit_list)

        print(f"\n[*] Scan complete. Parsed {scanned_files} files.", file=sys.stderr)
        print("-" * 60, file=sys.stderr)

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

        for kw, genome_count in failed_keywords:
            print(
                f"[-] Keyword '{kw}' failed threshold (Found in {genome_count} genomes).",
                file=sys.stderr,
            )

        if not valid_results:
            print(
                "\n[!] No keywords met the minimum genome threshold. No output written.",
                file=sys.stderr,
            )
            print("-" * 60, file=sys.stderr)
            return

        # 3-Tier sort: most conserved first
        #   Tier 1: genome count          DESCENDING (most conserved at top)
        #   Tier 2: total physical copies DESCENDING (single-copy before multicopy)
        #   Tier 3: alphabetical keyword  ASCENDING
        sorted_results = sorted(
            valid_results.items(),
            key=lambda item: (
                -len(item[1]),
                -sum(len(hits) for hits in item[1].values()),
                item[0].lower(),
            ),
        )

        for kw, genomes_with_kw in sorted_results:
            print(
                f"[+] Keyword '{kw}' is CONSERVED across {len(genomes_with_kw)} genomes.",
                file=sys.stderr,
            )

        print("-" * 60, file=sys.stderr)

        TSV_HEADER = "Keyword\tGenomes_Found\tGenome_File\tLocus\tLocus_Tag\tProduct\tProtein_Sequence"
        fasta_path = (
            args.output.with_suffix(".fasta") if (args.fasta and args.output) else None
        )

        if args.output:
            try:
                # Stream-write rows directly rather than accumulating strings in RAM.
                # For large eukaryotic or metagenomic datasets, list accumulation +
                # join() can spike memory significantly. Writing row-by-row keeps
                # memory usage flat regardless of dataset size.
                with ExitStack() as stack:
                    out_tsv = stack.enter_context(
                        open(args.output, "w", encoding="utf-8-sig")
                    )
                    out_fasta = (
                        stack.enter_context(open(fasta_path, "w", encoding="utf-8"))
                        if fasta_path
                        else None
                    )

                    out_tsv.write(TSV_HEADER + "\n")

                    for kw, genomes_with_kw in sorted_results:
                        genome_count = len(genomes_with_kw)
                        for genome_name, hits in genomes_with_kw.items():
                            for hit in hits:
                                out_tsv.write(
                                    f"{kw}\t{genome_count}\t{genome_name}\t"
                                    f"{hit['locus']}\t{hit['locus_tag']}\t"
                                    f"{hit['product']}\t{hit['translation']}\n"
                                )
                                if out_fasta and hit["translation"]:
                                    # Keyword is prepended to the FASTA header to guarantee
                                    # uniqueness when the same gene matches multiple keywords.
                                    # Without the keyword prefix, a gene matching both
                                    # "nisin" and "bacteriocin" would produce two identical
                                    # headers, causing crashes in MAFFT/IQ-TREE.
                                    header = (
                                        f">{kw}|{genome_name}|{hit['locus']}|"
                                        f"{hit['locus_tag']}|{hit['product']}"
                                    )
                                    seq_wrapped = wrap_fasta(hit["translation"])
                                    out_fasta.write(f"{header}\n{seq_wrapped}\n")

                print(
                    f"[*] Success! Matrix written to {args.output.resolve()}",
                    file=sys.stderr,
                )
                if fasta_path:
                    print(
                        f"[*] FASTA sequences written to {fasta_path.resolve()}",
                        file=sys.stderr,
                    )

            except OSError as e:
                sys.exit(f"\n[!] Error writing file: {e}")

        else:
            # No output file — dump TSV to stdout
            print(
                "\n[!] Note: No output file specified (-o). Printing to stdout:\n",
                file=sys.stderr,
            )
            sys.stdout.write(TSV_HEADER + "\n")
            for kw, genomes_with_kw in sorted_results:
                genome_count = len(genomes_with_kw)
                for genome_name, hits in genomes_with_kw.items():
                    for hit in hits:
                        sys.stdout.write(
                            f"{kw}\t{genome_count}\t{genome_name}\t"
                            f"{hit['locus']}\t{hit['locus_tag']}\t"
                            f"{hit['product']}\t{hit['translation']}\n"
                        )

    except KeyboardInterrupt:
        sys.exit("\n[!] Scan interrupted by user.")


if __name__ == "__main__":
    main()
