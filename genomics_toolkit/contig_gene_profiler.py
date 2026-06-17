#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

"""Extract genes from a specific contig in a BAKTA GBK/GBFF file.

This module provides tools to extract and analyze genes from a specific contig
in GenBank format files produced by BAKTA genome annotation. It supports multiple
output formats (table, TSV) and can filter by feature type (CDS, rRNA, tRNA).

Supports output with or without protein sequences for detailed downstream analysis.

Note:
    Associated with ongoing, unpublished research (manuscript in
    preparation). Correct attribution is requested when used in
    derivative works.

Examples:
    Extract genes from contig_2 and display as table::

        python3 contig_gene_profiler.py C5_genome.gbff contig_2

    Extract CDS genes and save to TSV file with sequences::

        python3 contig_gene_profiler.py C5_genome.gbff contig_2 \\
            --output plasmid.tsv --format tsv

    Extract rRNA genes with sequences included::

        python3 contig_gene_profiler.py C5_genome.gbff contig_1 \\
            --type rRNA --output rrna_genes.tsv
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.1.0"

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


def extract_contig_genes(gbk_file, contig_name, feature_type="CDS"):
    """Extract genes from a specific contig in a GenBank file.

    Parses a GenBank format file and extracts all features from a specified
    contig. Optionally filters by feature type (CDS, rRNA, tRNA, or all).
    Also calculates contig-level statistics (length, GC content).

    Args:
        gbk_file (str):
            Path to the GBFF/GBK file.
        contig_name (str):
            Name/identifier of the contig to extract (e.g., "contig_2").
            Partial matches are supported (e.g., "contig" matches "contig_2").
        feature_type (str, optional):
            Type of features to extract. Defaults to "CDS".
            Options: "CDS", "rRNA", "tRNA", or "all".

    Returns:
        tuple:
            A tuple containing:
                - genes (list): List of gene dictionaries with keys:
                    - locus_tag (str): Unique gene identifier
                    - gene_name (str): Gene name (if available)
                    - product (str): Gene product description
                    - type (str): Feature type (CDS, tRNA, etc.)
                    - start (int): Start position in contig
                    - end (int): End position in contig
                    - strand (str): "+" or "-"
                    - length (int): Feature length in bp
                    - sequence (str): Nucleotide sequence (if CDS/rRNA/tRNA)
                - contig_length (int): Total contig length in bp
                - contig_gc (float): GC content percentage
                - contig_id (str): Full contig identifier

    Raises:
        SystemExit: If the specified contig is not found in the file.
    """

    genes = []
    contig_found = False
    contig_length = 0
    contig_gc = 0.0
    contig_id = ""

    # Parse the GenBank file
    for record in SeqIO.parse(gbk_file, "genbank"):
        # Check if this is the contig we're looking for
        if contig_name in record.id:
            contig_found = True
            contig_id = record.id
            contig_length = len(record.seq)

            # Calculate GC content
            gc_count = sum(1 for base in record.seq if base in "GC")
            contig_gc = (gc_count / len(record.seq)) * 100

            # Extract features
            for feature in record.features:
                # Filter by feature type
                if feature_type != "all" and feature.type != feature_type:
                    continue

                # Extract metadata
                locus_tag = feature.qualifiers.get("locus_tag", ["N/A"])[0]
                product = feature.qualifiers.get("product", ["Unknown"])[0]
                gene_name = feature.qualifiers.get("gene", [""])[0]

                # Get location
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = "+" if feature.location.strand == 1 else "-"

                # Get length
                length = end - start

                # Extract sequence from the contig
                if feature.location.strand == 1:
                    sequence = str(record.seq[start:end])
                else:
                    # Reverse complement for minus strand
                    sequence = str(record.seq[start:end].reverse_complement())

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

    if not contig_found:
        print(f"ERROR: Contig '{contig_name}' not found in {gbk_file}")
        sys.exit(1)

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


def print_statistics(genes, contig_length, contig_gc):
    """Print genome statistics.

    Displays summary statistics about the contig and extracted genes,
    including coding density and total coding bases.

    Args:
        genes (list):
            List of extracted gene dictionaries.
        contig_length (int):
            Total contig length in bp.
        contig_gc (float):
            GC content percentage.
    """
    print(f"\n{'='*80}")
    print("STATISTICS")
    print(f"{'='*80}")
    print(f"Total genes found: {len(genes)}")
    print(f"Contig length: {contig_length:,} bp")
    print(f"GC content: {contig_gc:.2f}%")

    if genes:
        total_bp = sum(gene["length"] for gene in genes)
        coding_density = (total_bp / contig_length) * 100
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
  python extract_contig_genes_v2.py C5_genome.gbff contig_2

  # Save to TSV file with sequences included
  python extract_contig_genes_v2.py C5_genome.gbff contig_2 \\
      --output plasmid.tsv --format tsv

  # Save to text file with sequences included
  python extract_contig_genes_v2.py C5_genome.gbff contig_2 \\
      --output plasmid.txt --format table

  # Extract only rRNA genes with sequences
  python extract_contig_genes_v2.py C5_genome.gbff contig_1 \\
      --type rRNA --output rrna_genes.tsv

  # Redirect to text file via shell (includes sequences)
  python extract_contig_genes_v2.py C5_genome.gbff contig_2 > output.txt
        """,
    )
    parser.add_argument("gbk_file", help="Path to GBFF or GBK file")
    parser.add_argument(
        "contig_name", help="Name of contig to extract (e.g., 'contig_2')"
    )
    parser.add_argument(
        "--type",
        default="CDS",
        choices=["CDS", "rRNA", "tRNA", "all"],
        help="Feature type to extract (default: CDS)",
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
            "    python extract_contig_genes_v2.py C5_genome.gbff contig_2 -o output.tsv"
        )
        print(
            "    or redirect output: python extract_contig_genes_v2.py C5_genome.gbff contig_2 > output.txt\n"
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
    print_statistics(genes, contig_length, contig_gc)


if __name__ == "__main__":
    main()
