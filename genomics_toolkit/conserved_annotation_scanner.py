#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

"""Conserved Annotation Scanner

Extracts and groups CDS product annotations across GenBank genomes to identify
conserved genes. Acts as a text-based core proteome profiler: aggregates all
/product qualifiers, normalizes them, and reports only gene products that meet
a specified genome frequency threshold. Filters uninformative 'hypothetical
protein' annotations by default.

EUKARYOTIC COMPATIBILITY:
    This script relies entirely on the /product and /translation qualifiers
    already present on each CDS feature, as annotated by the upstream tool
    (Prokka, Bakta, NCBI RefSeq, RAST) — it never re-derives a sequence from
    genomic coordinates the way contig_gene_profiler.py does. This means it
    is already safe for multi-exon eukaryotic CDS features out of the box:
    /translation is the correctly spliced protein regardless of how many
    exons the underlying CDS has, since the annotation tool computed it, not
    this script.

    The one eukaryotic-specific issue that DID exist — isoforms inflating
    the Tier-2 "physical copies" sort statistic — is fixed as of v1.3.0;
    see the changelog below.

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
    sequence-based clustering via pairwise_homolog_finder.py.

OUTPUT SORT ORDER (most to least conserved):
    1. Number of genomes found        (Descending — most conserved first)
    2. Total physical copies          (Descending — single-copy before multicopy)
    3. Alphabetical by product name   (Ascending)

    Descending Tier 1 ensures core genes appear at the TOP of the output
    matrix. Descending Tier 2 prevents copy-number inflation (a gene with
    12 copies in one genome ranks BELOW a perfectly distributed single-copy
    universal gene with the same genome count).

    "Physical copies" (Tier 2) counts DISTINCT locus_tags per genome, not
    raw CDS feature count — see v1.3.0 changelog below for why.

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
    This module is part of ongoing research and is associated with an upcoming
    publication. Please cite appropriately when used in derivative works.
    See LICENSE file in the repository root for full license terms.

    v1.3.0: Fixed two bugs and one stale reference found during review.
    (1) ``_normalize_product()`` could reduce a raw /product string to an
    empty string — "Putative", "(hypothetical)",
    and "predicted" (three different, unrelated raw annotations) ALL
    normalize to "". Without a fallback, every such low-information
    annotation across every genome in a scan would silently merge into one
    fake "conserved gene group" under the shared empty key, regardless of
    how unrelated the underlying genes actually were. This also broke the
    existing hypothetical-protein filter for phrasings like
    "(hypothetical)": "hypothetical" is not a substring of "", so the
    filter (which checks the NORMALIZED key) never caught it, letting an
    uninformative annotation leak into the output as if it were a real
    conserved gene. Fixed by falling back to the lowercased-but-otherwise-
    unmodified original product string whenever normalization yields "" —
    this preserves enough of the original text to avoid the
    cross-contamination, and keeps "hypothetical" as a substring so the
    existing filter still works on this fallback value too.
    (2) The Tier-2 sort key summed raw CDS feature count
    (``len(hits)``) per genome as "physical copies" — a eukaryotic gene
    with 3 annotated splice isoforms sharing one
    locus_tag in one genome, plus one normal single-copy ortholog in a
    second genome, reported "4 total physical copies" when both genomes
    actually have exactly 1 physical locus each (isoforms are alternative
    transcripts of the SAME locus, not separate gene copies). This did not
    corrupt the TSV output rows (every hit is still written, with its
    locus_tag visible) or the Tier-1 genome count, only the Tier-2 sort
    statistic and therefore the report's ranking. Fixed by counting
    DISTINCT locus_tags per genome (``len({hit["locus_tag"] for hit in
    hits})``) rather than raw hit count. No effect on prokaryotic genomes,
    where every locus_tag already maps to exactly one CDS feature, so
    distinct-count and raw-count are always identical. Edge case: multiple
    genuinely unannotated CDS features (no /locus_tag at all, stored under
    the literal fallback string "UNKNOWN") within the same genome and
    product group will now collapse to a sort-count of 1 rather than their
    true count — a deliberate tradeoff, since this is rarer than the
    eukaryotic isoform case and does not reduce TSV output completeness
    (every such hit is still written as its own row).
    (3) Corrected a stale reference: this docstring referenced
    "gbk_ortholog_finder.py" twice for sequence-based clustering — this
    was the script's former name before it was renamed to
    pairwise_homolog_finder.py, a more appropriate name for what it does
    (confirmed via that file's own docstring, which independently
    describes itself as the Smith-Waterman/BLOSUM62 sequence-based
    clustering tool). The rename was never propagated to these two
    references. Corrected both.

    v1.3.1: Added the missing shebang/SPDX-license/copyright header that
    every other script in this toolkit carries at the very top of the
    file — this one had skipped straight to the module docstring. Also
    removed the now-redundant inline "License: MIT" line from the
    docstring body, since the SPDX header covers it. Documentation only;
    no behavior change.

Example Usage:
    # Standard run: Find genes conserved in at least 2 genomes, output TSV
    $ python3 conserved_annotation_scanner.py -i references/ --min_genomes 2 -o core.tsv

    # Auto-FASTA run: Output both 'core.tsv' and 'core.fasta'
    $ python3 conserved_annotation_scanner.py -i references/ --min_genomes 2 -o core.tsv -f

    # Exact run: Find genes present in EXACTLY 2 genomes
    $ python3 conserved_annotation_scanner.py -i references/ --min_genomes 2 --exact -o strict_two.tsv
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.3.1"

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

    EMPTY-STRING FALLBACK (v1.3.0):
        The steps above can reduce some raw /product strings to nothing at
        all — "Putative", "(hypothetical)", and
        "predicted" all normalize to "". Returning "" directly would merge
        every such low-information annotation, across every genome scanned,
        into one shared (and meaningless) group key — regardless of how
        unrelated the actual underlying genes are. It would also silently
        defeat the caller's "hypothetical protein" filter for phrasings
        like "(hypothetical)", since "hypothetical" is not a substring of
        "". When normalization would yield "", this function instead
        returns the lowercased-but-otherwise-unmodified original string,
        which both avoids the cross-contamination and keeps "hypothetical"
        (or any other identifying text) intact for that filter to act on.

    Limitation: Cannot resolve abbreviation-vs-full-name discrepancies
    (e.g., "atpA" vs "ATP synthase subunit alpha"). These require
    sequence-based clustering. See pairwise_homolog_finder.py.

    Args:
        product: Raw /product qualifier string from a GenBank CDS feature.

    Returns:
        Normalized, lowercase string suitable for grouping. Falls back to
        the lowercased raw string if normalization would otherwise yield
        an empty string.
    """
    raw_lower = product.lower().strip()
    name = raw_lower
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
    normalized = " ".join(words).strip()

    # See EMPTY-STRING FALLBACK note above.
    return normalized if normalized else raw_lower


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
        #
        # Tier 2 counts DISTINCT locus_tags per genome (v1.3.0), not raw CDS
        # feature count. Eukaryotic GenBank annotations represent alternative
        # splice isoforms as separate CDS features sharing one locus_tag —
        # a gene with 3 isoforms in one genome plus
        # 1 normal ortholog in another previously reported "4 total physical
        # copies" when both genomes actually have exactly 1 physical locus
        # each. Counting distinct locus_tags collapses isoforms of the same
        # locus back down to 1 per genome. No effect on prokaryotic genomes,
        # where every locus_tag already maps to exactly one CDS feature
        # (distinct-count == raw-count always). See the module changelog for
        # the "UNKNOWN" locus_tag edge-case tradeoff this introduces.
        sorted_results = sorted(
            filtered_results.items(),
            key=lambda item: (
                -len(item[1]),
                -sum(
                    len({hit["locus_tag"] for hit in hits}) for hits in item[1].values()
                ),
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
