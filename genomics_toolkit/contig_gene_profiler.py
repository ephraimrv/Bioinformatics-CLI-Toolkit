#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

"""Extract genes from a specific contig in a BAKTA GBK/GBFF file.

This module provides tools to extract and analyze genes from a specific contig
in GenBank format files produced by BAKTA genome annotation. It supports multiple
output formats (table, TSV) and can filter by feature type (CDS, rRNA, tRNA,
mRNA, ncRNA, or all).

Supports output with or without protein sequences for detailed downstream analysis.

EUKARYOTIC COMPATIBILITY:
    Unlike the upstream-promoter family of scripts in this toolkit
    (gbk_promoter_finder.py, regulon_scanner.py,
    comparative_kmer_analyzer.py — all CDS-anchored, not TSS-anchored,
    and explicitly prokaryote-only), this script's sequence/length
    extraction is genuinely correct for multi-exon eukaryotic genes as of
    v1.2.0: it uses Biopython's ``feature.extract()``, which correctly
    concatenates only the exonic parts of a compound (joined) location,
    strand-aware — rather than a manual ``start:end`` slice of genomic
    DNA, which would silently include intron sequence. ``length`` is
    computed from the extracted sequence itself, so it always matches
    what's actually in the ``sequence`` column, whether the feature has a
    single-part location (every prokaryotic gene, and single-exon
    eukaryotic genes) or a compound/joined one (multi-exon eukaryotic
    genes).

    Remaining eukaryotic-specific gaps:
    - ``--type`` now includes ``mRNA``/``ncRNA`` (v1.2.0) alongside
      CDS/rRNA/tRNA, since eukaryotic GenBank files represent genes via
      mRNA features that prokaryotic Prokka/Bakta output never emits.
    - ``--type all`` does not attempt to deduplicate overlapping feature
      types (e.g. a single gene's ``mRNA`` and ``CDS`` features cover
      overlapping spans, and a ``gene`` wrapper duplicates its own
      CDS/rRNA/tRNA coordinates) — coding-density reporting is disabled
      specifically for ``--type all`` because of this; see
      ``print_statistics()``.
    - GC% is computed over the whole contig only, not feature-by-feature;
      this is organism-agnostic and unaffected by exon/intron structure.

Note:
    Associated with ongoing, unpublished research (manuscript in
    preparation). Correct attribution is requested when used in
    derivative works.

    v1.2.0: Fixed three bugs and two eukaryotic-compatibility gaps found
    during review.
    (1) ``extract_contig_genes()`` matched contigs with
    ``if contig_name in record.id`` and never broke out of the loop on a
    match — confirmed empirically that searching "contig_2" on a file
    containing contig_2/contig_20/contig_21/contig_200 would silently
    merge genes from all four into one combined list, while the printed
    contig_id/length/GC% reflected only whichever record matched last.
    Now collects all matches first; proceeds only if exactly one match is
    found, and exits with an explicit list of the ambiguous matches
    otherwise. Partial-match convenience is preserved for the common case
    where it resolves to exactly one contig.
    (2) Feature sequences were extracted via manual
    ``record.seq[start:end]`` slicing (+ reverse_complement for minus
    strand) instead of ``feature.extract(record.seq)``. For any
    multi-exon feature (a CompoundLocation/``join()``, which only occurs
    in eukaryotic annotation), manual slicing returns the full genomic
    span INCLUDING introns, mislabeled as the gene's sequence; ``length``
    was likewise the genomic envelope (end-start) rather than the actual
    feature length. No behavior change for prokaryotic (single-part
    location) features, where both approaches are identical.
    (3) ``--type all``'s coding-density statistic summed every feature's
    length regardless of type, including the contig-spanning ``source``
    feature and duplicate ``gene``/``CDS`` pairs at identical
    coordinates — confirmed empirically that a toy 2-gene contig already
    reported 144.8% "coding density" this way. Rather than introduce a
    new feature-type allowlist that could still double-count (e.g. a
    eukaryotic gene's overlapping ``mRNA`` and ``CDS`` spans),
    coding-density reporting is now skipped entirely for ``--type all``
    with an explanatory message; behavior for any single ``--type``
    (CDS/rRNA/tRNA/mRNA/ncRNA) is completely unchanged, since a
    single-type ``genes`` list was never affected by this bug.
    (4) Added ``mRNA`` and ``ncRNA`` to ``--type``'s choices — eukaryotic
    GenBank files represent genes via these feature types, which were
    previously unselectable except through the (buggy) ``all`` path.
    (5) GC% calculation (``base in "GC"``) was case-sensitive and would
    undercount GC on soft-masked sequences (lowercase repeat regions,
    more common in eukaryotic assemblies) — now uppercases the sequence
    once before counting.
    (6) Fixed a stale filename in the CLI epilog and the
    no-output-file runtime hint: both referenced
    ``extract_contig_genes_v2.py``, an earlier name for this script,
    making every shown example command non-copy-pasteable as written.
    Corrected to ``contig_gene_profiler.py`` throughout (the module
    docstring's own Examples section already had the correct name and
    was not affected).

Examples:
    Extract genes from contig_2 and display as table::

        python3 contig_gene_profiler.py C5_genome.gbff contig_2

    Extract CDS genes and save to TSV file with sequences::

        python3 contig_gene_profiler.py C5_genome.gbff contig_2 \\
            --output plasmid.tsv --format tsv

    Extract rRNA genes with sequences included::

        python3 contig_gene_profiler.py C5_genome.gbff contig_1 \\
            --type rRNA --output rrna_genes.tsv

    Extract mRNA genes from a eukaryotic GenBank file::

        python3 contig_gene_profiler.py yeast_genome.gbff chrI \\
            --type mRNA --output chrI_mrna.tsv
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.2.0"

import sys
import argparse
import csv

try:
    from Bio import SeqIO
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )

# Feature types whose extracted sequence represents real gene-body
# content. Used only to decide whether coding-density reporting is safe
# (see print_statistics()) — not used to filter what gets extracted or
# displayed, which is governed entirely by --type.
GENE_BODY_TYPES = frozenset({"CDS", "rRNA", "tRNA", "ncRNA", "tmRNA", "mRNA"})


def extract_contig_genes(gbk_file, contig_name, feature_type="CDS"):
    """Extract genes from a specific contig in a GenBank file.

    Parses a GenBank format file and extracts all features from a specified
    contig. Optionally filters by feature type (CDS, rRNA, tRNA, mRNA,
    ncRNA, or all). Also calculates contig-level statistics (length, GC
    content).

    Args:
        gbk_file (str):
            Path to the GBFF/GBK file.
        contig_name (str):
            Name/identifier of the contig to extract (e.g., "contig_2").
            Partial matches are supported (e.g., "contig" matches
            "contig_2") PROVIDED the partial name resolves to exactly one
            contig in the file — if it matches more than one, this raises
            an error listing every match rather than silently combining
            them (see v1.2.0 changelog).
        feature_type (str, optional):
            Type of features to extract. Defaults to "CDS".
            Options: "CDS", "rRNA", "tRNA", "mRNA", "ncRNA", or "all".

    Returns:
        tuple:
            A tuple containing:
                - genes (list): List of gene dictionaries with keys:
                    - locus_tag (str): Unique gene identifier
                    - gene_name (str): Gene name (if available)
                    - product (str): Gene product description
                    - type (str): Feature type (CDS, tRNA, etc.)
                    - start (int): Genomic start position in contig
                    - end (int): Genomic end position in contig
                    - strand (str): "+" or "-"
                    - length (int): Length of the EXTRACTED sequence in bp
                      (equal to end-start for single-part locations; for a
                      multi-exon/compound location, this is the spliced
                      coding length, NOT the genomic envelope end-start)
                    - sequence (str): Nucleotide sequence (if CDS/rRNA/
                      tRNA/mRNA/ncRNA), correctly spliced for multi-exon
                      features via feature.extract()
                - contig_length (int): Total contig length in bp
                - contig_gc (float): GC content percentage
                - contig_id (str): Full contig identifier

    Raises:
        SystemExit: If the specified contig is not found, or if
            ``contig_name`` matches more than one contig in the file.
    """
    # Collect every record whose ID contains contig_name BEFORE processing
    # any of them. Previously this loop processed (and accumulated genes
    # from) every match it found with no break, while contig_id/length/GC
    # were silently overwritten on each successive match — so a partial
    # name matching more than one contig merged their genes under a
    # single mislabeled identity. Resolving to exactly one match first
    # makes that ambiguity an explicit error instead of a silent merge.
    matched_records = [
        record
        for record in SeqIO.parse(gbk_file, "genbank")
        if contig_name in record.id
    ]

    if not matched_records:
        print(f"ERROR: Contig '{contig_name}' not found in {gbk_file}")
        sys.exit(1)

    if len(matched_records) > 1:
        matched_ids = [r.id for r in matched_records]
        print(
            f"ERROR: '{contig_name}' matches {len(matched_records)} contigs "
            f"in {gbk_file}, not exactly one:"
        )
        for mid in matched_ids:
            print(f"    - {mid}")
        print(
            "Partial matching only proceeds when it resolves to exactly "
            "one contig — matching more than one would otherwise merge "
            "their genes under a single contig_id/length/GC% from "
            "whichever happened to be processed last. Provide a more "
            "specific name (e.g. the full contig ID shown above) to "
            "disambiguate."
        )
        sys.exit(1)

    record = matched_records[0]
    contig_id = record.id
    contig_length = len(record.seq)

    # Case-insensitive GC count: lowercase bases (e.g. soft-masked repeat
    # regions, more common in eukaryotic assemblies) were previously
    # excluded from the count entirely, undercounting GC%.
    upper_seq = str(record.seq).upper()
    gc_count = sum(1 for base in upper_seq if base in "GC")
    contig_gc = (gc_count / contig_length) * 100 if contig_length else 0.0

    genes = []
    for feature in record.features:
        # Filter by feature type
        if feature_type != "all" and feature.type != feature_type:
            continue

        # Extract metadata
        locus_tag = feature.qualifiers.get("locus_tag", ["N/A"])[0]
        product = feature.qualifiers.get("product", ["Unknown"])[0]
        gene_name = feature.qualifiers.get("gene", [""])[0]

        # Genomic envelope coordinates (context/display only — see
        # `length` below for the actual extracted-sequence length).
        start = int(feature.location.start)
        end = int(feature.location.end)
        strand = "+" if feature.location.strand == 1 else "-"

        # Extract the feature's sequence via Biopython's feature.extract(),
        # which handles strand and compound (joined) locations correctly
        # by concatenating only the exonic parts in genomic order before
        # reverse-complementing as a whole for minus-strand features. The
        # previous manual `record.seq[start:end]` (+ reverse_complement)
        # is only equivalent to this for a single-part location (every
        # prokaryotic gene, and single-exon eukaryotic genes) — for a
        # multi-exon eukaryotic CDS it would instead return the full
        # genomic span INCLUDING introns, silently mislabeled as the
        # gene's sequence.
        try:
            sequence = str(feature.extract(record.seq))
        except Exception:
            # Fall back to manual slicing only if extract() itself raises
            # (e.g. a malformed or out-of-bounds location) — better to
            # report something for this one feature than abort the whole
            # run over it.
            if feature.location.strand == 1:
                sequence = str(record.seq[start:end])
            else:
                sequence = str(record.seq[start:end].reverse_complement())

        # Length now reflects the ACTUAL extracted sequence, not
        # end-start. Identical to end-start for single-part locations;
        # for a multi-exon feature, end-start would include intron
        # length while len(sequence) correctly does not. Keeping `length`
        # derived from `sequence` guarantees the two fields can never
        # silently disagree.
        length = len(sequence)

        gene_dict = {
            "locus_tag": locus_tag,
            "gene_name": gene_name,
            "product": product,
            "type": feature.type,
            "start": start,
            "end": end,
            "strand": strand,
            "length": length,
            "sequence": sequence,
        }

        genes.append(gene_dict)

    return genes, contig_length, contig_gc, contig_id


def print_header_info(contig_id, contig_length, contig_gc, genes, feature_type):
    """Print header information about the contig.

    Args:
        contig_id (str):
            Full contig identifier.
        contig_length (int):
            Length of the contig in bp.
        contig_gc (float):
            GC content percentage of the contig.
        genes (list):
            List of extracted gene dictionaries.
        feature_type (str):
            Type of features extracted (e.g., "CDS", "tRNA").
    """
    print(f"\n{'='*80}")
    print(f"CONTIG: {contig_id}")
    print(f"Length: {contig_length:,} bp | GC: {contig_gc:.2f}%")
    print(f"Feature type: {feature_type} | Genes found: {len(genes)}")
    print(f"{'='*80}\n")


def print_genes_table(genes):
    """Print genes in a nicely formatted table.

    Displays gene information in human-readable columns without sequences.
    Useful for quick terminal viewing. Column widths are optimized for
    readability.

    Args:
        genes (list):
            List of gene dictionaries to display.
    """
    if not genes:
        print("No genes found.")
        return

    # Define column widths
    col_locus = 15
    col_gene = 8
    col_type = 10
    col_product = 40
    col_start = 10
    col_end = 10
    col_strand = 6
    col_length = 10

    # Print header
    header = f"{'Locus Tag':<{col_locus}} | {'Gene':<{col_gene}} | {'Type':<{col_type}} | {'Product':<{col_product}} | {'Start':<{col_start}} | {'End':<{col_end}} | {'Strand':<{col_strand}} | {'Length':<{col_length}}"
    print(header)

    # Calculate total width for separator line
    total_width = len(header)
    print("-" * total_width)

    # Print each gene
    for gene in genes:
        locus = gene["locus_tag"]
        gene_name = gene["gene_name"] if gene["gene_name"] else ""
        gtype = gene["type"]
        product = gene["product"][:35] + "..."
        start = str(gene["start"])
        end = str(gene["end"])
        strand = gene["strand"]
        length = str(gene["length"])

        row = f"{locus:<{col_locus}} | {gene_name:<{col_gene}} | {gtype:<{col_type}} | {product:<{col_product}} | {start:<{col_start}} | {end:<{col_end}} | {strand:<{col_strand}} | {length:<{col_length}}"
        print(row)

    print("-" * total_width)
    print(f"Total genes: {len(genes)}\n")


def save_genes_tsv(genes, output_file, include_sequence=True):
    """Save genes to a TSV (Tab-Separated Values) file.

    Exports gene data in tab-separated format suitable for spreadsheet
    applications (Excel, LibreOffice Calc) and data analysis tools.
    Optionally includes nucleotide sequences.

    Args:
        genes (list):
            List of gene dictionaries to save.
        output_file (str):
            Path to output TSV file.
        include_sequence (bool, optional):
            Whether to include nucleotide sequences. Defaults to True.

    Raises:
        SystemExit: If file write operation fails.
    """
    if not genes:
        print("No genes to save.")
        return

    try:
        with open(output_file, "w", newline="") as f:
            # Define column headers
            if include_sequence:
                fieldnames = [
                    "locus_tag",
                    "gene_name",
                    "type",
                    "product",
                    "start",
                    "end",
                    "strand",
                    "length",
                    "sequence",
                ]
            else:
                fieldnames = [
                    "locus_tag",
                    "gene_name",
                    "type",
                    "product",
                    "start",
                    "end",
                    "strand",
                    "length",
                ]

            # Create TSV writer (delimiter='\t' makes it tab-separated)
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")

            # Write header row
            writer.writeheader()

            # Write data rows
            for gene in genes:
                if include_sequence:
                    writer.writerow(gene)
                else:
                    # Write only non-sequence columns
                    row = {k: v for k, v in gene.items() if k != "sequence"}
                    writer.writerow(row)

        print(f"✓ TSV file saved: {output_file}")
        print(f"  Rows: {len(genes)} genes")
        print(f"  Columns: {len(fieldnames)}")
        if include_sequence:
            print(f"  Sequences: INCLUDED")

    except Exception as e:
        print(f"ERROR saving TSV file: {e}")
        sys.exit(1)


def save_genes_table(genes, output_file, include_sequence=True):
    """Save genes in formatted table format to a text file.

    Exports gene data in a human-readable table format for easy viewing
    in text editors. Optionally includes nucleotide sequences.

    Args:
        genes (list):
            List of gene dictionaries to save.
        output_file (str):
            Path to output text file.
        include_sequence (bool, optional):
            Whether to include nucleotide sequences. Defaults to True.

    Raises:
        SystemExit: If file write operation fails.
    """
    if not genes:
        print("No genes to save.")
        return

    try:
        with open(output_file, "w") as f:
            # Write header
            f.write(
                f"{'Locus Tag':<15} | {'Type':<8} | {'Product':<50} | {'Location':<15}\n"
            )
            f.write("-" * 130 + "\n")

            # Write each gene
            for gene in genes:
                locus = gene["locus_tag"]
                gtype = gene["type"]
                product = gene["product"][:48]
                location = f"{gene['start']}..{gene['end']} ({gene['strand']})"

                f.write(f"{locus:<15} | {gtype:<8} | {product:<50} | {location:<15}\n")

                # Add sequence if requested
                if include_sequence:
                    sequence = gene["sequence"]
                    f.write(f"  Sequence ({len(sequence)} bp): {sequence}\n")
                    f.write("\n")

            f.write(f"\nTotal genes: {len(genes)}\n")
            if include_sequence:
                f.write("[Sequences included in output]\n")

        print(f"✓ Table file saved: {output_file}")
        print(f"  Rows: {len(genes)} genes")
        if include_sequence:
            print(f"  Sequences: INCLUDED")

    except Exception as e:
        print(f"ERROR saving table file: {e}")
        sys.exit(1)


def categorize_genes(genes):
    """Categorize genes by type and print summary.

    Args:
        genes (list):
            List of gene dictionaries to categorize.
    """
    categories = {}

    for gene in genes:
        ftype = gene["type"]
        if ftype not in categories:
            categories[ftype] = []
        categories[ftype].append(gene)

    print(f"{'='*80}")
    print("GENE CATEGORY SUMMARY")
    print(f"{'='*80}\n")

    for ftype in sorted(categories.keys()):
        genes_of_type = categories[ftype]
        print(f"{ftype.upper()}: {len(genes_of_type)} genes")


def print_statistics(genes, contig_length, contig_gc, feature_type):
    """Print genome statistics.

    Displays summary statistics about the contig and extracted genes,
    including coding density and total coding bases.

    Coding-density reporting is skipped for ``--type all`` (see v1.2.0
    changelog note in the module docstring): mixing feature types means
    summing their lengths would double-count overlapping spans (e.g. a
    eukaryotic gene's ``mRNA`` and ``CDS`` features, or any gene's
    ``gene`` wrapper duplicating its own CDS/rRNA/tRNA coordinates) and
    would also include the whole-contig ``source`` feature if present.
    For any single ``--type`` (CDS/rRNA/tRNA/mRNA/ncRNA), ``genes``
    contains only that one type, so no overlap is possible and the
    calculation is unchanged from previous versions.

    Args:
        genes (list):
            List of extracted gene dictionaries.
        contig_length (int):
            Total contig length in bp.
        contig_gc (float):
            GC content percentage.
        feature_type (str):
            The --type value used for extraction (e.g. "CDS", "all").
    """
    print(f"\n{'='*80}")
    print("STATISTICS")
    print(f"{'='*80}")
    print(f"Total genes found: {len(genes)}")
    print(f"Contig length: {contig_length:,} bp")
    print(f"GC content: {contig_gc:.2f}%")

    if genes:
        if feature_type == "all":
            print(
                "Coding density: not computed for --type all — mixing "
                "feature types (overlapping gene/CDS/mRNA spans, plus the "
                "whole-contig 'source' feature if present) would "
                "double-count the same bases. Re-run with a single "
                "--type (e.g. CDS) for a meaningful density figure."
            )
        else:
            total_bp = sum(gene["length"] for gene in genes)
            coding_density = (
                (total_bp / contig_length) * 100 if contig_length else 0.0
            )
            print(f"Total coding bases: {total_bp:,} bp")
            print(f"Coding density: {coding_density:.1f}%")


def main():
    """Parse arguments and execute gene extraction workflow."""
    parser = argparse.ArgumentParser(
        description="Extract genes from a specific contig in a BAKTA GBK/GBFF file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Display genes as table on terminal (no sequences)
  python3 contig_gene_profiler.py C5_genome.gbff contig_2

  # Save to TSV file with sequences included
  python3 contig_gene_profiler.py C5_genome.gbff contig_2 \\
      --output plasmid.tsv --format tsv

  # Save to text file with sequences included
  python3 contig_gene_profiler.py C5_genome.gbff contig_2 \\
      --output plasmid.txt --format table

  # Extract only rRNA genes with sequences
  python3 contig_gene_profiler.py C5_genome.gbff contig_1 \\
      --type rRNA --output rrna_genes.tsv

  # Extract mRNA genes from a eukaryotic GenBank file
  python3 contig_gene_profiler.py yeast_genome.gbff chrI \\
      --type mRNA --output chrI_mrna.tsv

  # Redirect to text file via shell (includes sequences)
  python3 contig_gene_profiler.py C5_genome.gbff contig_2 > output.txt
        """,
    )
    parser.add_argument("gbk_file", help="Path to GBFF or GBK file")
    parser.add_argument(
        "contig_name", help="Name of contig to extract (e.g., 'contig_2')"
    )
    parser.add_argument(
        "--type",
        default="CDS",
        choices=["CDS", "rRNA", "tRNA", "mRNA", "ncRNA", "all"],
        help=(
            "Feature type to extract (default: CDS). 'mRNA' and 'ncRNA' "
            "are for eukaryotic GenBank files, which represent genes via "
            "these feature types rather than (or in addition to) CDS."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output filename (if not specified, prints to screen only)",
    )
    parser.add_argument(
        "--format",
        default="tsv",
        choices=["table", "tsv"],
        help="Output format: 'table' for human-readable, 'tsv' for spreadsheet (default: table)",
    )

    args = parser.parse_args()

    # Extract genes
    genes, contig_length, contig_gc, contig_id = extract_contig_genes(
        args.gbk_file, args.contig_name, args.type
    )

    # Print to screen (always show header and summary)
    print_header_info(contig_id, contig_length, contig_gc, genes, args.type)

    # Decide: show genes on screen or save to file
    if args.output:
        # USER REQUESTED TO SAVE → Don't display genes, just save and confirm
        print(f"{'='*80}")
        print(f"Saving {len(genes)} genes to file (with sequences)...")
        print(f"{'='*80}\n")

        if args.format == "tsv":
            save_genes_tsv(genes, args.output, include_sequence=True)
        else:  # table format
            save_genes_table(genes, args.output, include_sequence=True)

        print(f"\n{'='*80}")
        print(f"✓ Successfully saved!")
        print(f"  File: {args.output}")
        print(f"  Format: {args.format.upper()}")
        print(f"  Genes: {len(genes)}")
        print(f"  Sequences: INCLUDED")
        print(f"{'='*80}\n")
    else:
        # NO OUTPUT FILE REQUESTED → Display genes on screen (clean, no sequences)
        print("[*] NOTE: Nucleotide sequences will be included if you save to a file:")
        print(
            "    python3 contig_gene_profiler.py C5_genome.gbff contig_2 -o output.tsv"
        )
        print(
            "    or redirect output: python3 contig_gene_profiler.py C5_genome.gbff contig_2 > output.txt\n"
        )

        if args.format == "table":
            print_genes_table(genes)
        else:  # TSV format
            print("\nTSV Format (tab-separated, no sequences):")
            print(
                "\t".join(
                    [
                        "locus_tag",
                        "gene_name",
                        "type",
                        "product",
                        "start",
                        "end",
                        "strand",
                        "length",
                    ]
                )
            )
            for gene in genes:
                print(
                    f"{gene['locus_tag']}\t{gene['gene_name']}\t{gene['type']}\t{gene['product']}\t{gene['start']}\t{gene['end']}\t{gene['strand']}\t{gene['length']}"
                )

    # Always print summary
    categorize_genes(genes)
    print_statistics(genes, contig_length, contig_gc, args.type)


if __name__ == "__main__":
    main()