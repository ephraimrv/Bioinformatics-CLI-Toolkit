#!/usr/bin/env python3
"""
Extract genes from a specific contig in a BAKTA GBK/GBFF file

Usage:
    python extract_contig_genes.py <gbk_file> <contig_name> [--type CDS]

Example:
    python extract_contig_genes.py C5_genome.gbff contig_2
    python extract_contig_genes.py C5_genome.gbff contig_2 --type CDS
    python extract_contig_genes.py C5_genome.gbff contig_2 --type all
"""

from Bio import SeqIO
import sys
import argparse


def extract_contig_genes(gbk_file, contig_name, feature_type="CDS"):
    """
    Extract genes from a specific contig in a GenBank file.

    Args:
        gbk_file (str): Path to the GBFF/GBK file
        contig_name (str): Name of the contig to extract (e.g., "contig_2")
        feature_type (str): Type of features to extract ("CDS", "rRNA", "tRNA", "all")

    Returns:
        list: List of gene dictionaries with metadata
    """

    genes = []
    contig_found = False
    contig_length = 0
    contig_gc = 0.0

    print(f"\nSearching for contig: {contig_name}")
    print(f"Feature type: {feature_type}\n")

    # Parse the GenBank file
    for record in SeqIO.parse(gbk_file, "genbank"):
        # Check if this is the contig we're looking for
        if contig_name in record.id:
            contig_found = True
            contig_length = len(record.seq)

            # Calculate GC content
            gc_count = sum(1 for base in record.seq if base in "GC")
            contig_gc = (gc_count / len(record.seq)) * 100

            print(f"{'='*80}")
            print(f"CONTIG: {record.id}")
            print(f"Length: {contig_length:,} bp | GC: {contig_gc:.2f}%")
            print(f"{'='*80}\n")

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

                gene_dict = {
                    "locus_tag": locus_tag,
                    "gene_name": gene_name,
                    "product": product,
                    "type": feature.type,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "length": length,
                }

                genes.append(gene_dict)

    if not contig_found:
        print(f"ERROR: Contig '{contig_name}' not found in {gbk_file}")
        sys.exit(1)

    return genes, contig_length, contig_gc


def print_genes_table(genes):
    """
    Print genes in a nice formatted table.
    """
    if not genes:
        print("No genes found.")
        return

    # Print header
    print(f"{'Locus Tag':<15} | {'Type':<8} | {'Product':<50} | {'Location':<15}")
    print("-" * 130)

    # Print each gene
    for gene in genes:
        locus = gene["locus_tag"]
        gtype = gene["type"]
        product = gene["product"][:48]  # Truncate long product names
        location = f"{gene['start']}..{gene['end']} ({gene['strand']})"

        print(f"{locus:<15} | {gtype:<8} | {product:<50} | {location:<15}")

    print(f"\nTotal genes: {len(genes)}")


def print_genes_detailed(genes):
    """
    Print genes with more detailed information.
    """
    if not genes:
        print("No genes found.")
        return

    for i, gene in enumerate(genes, 1):
        print(f"\n{i}. {gene['locus_tag']}")
        print(f"   Type: {gene['type']}")
        print(f"   Product: {gene['product']}")
        if gene["gene_name"]:
            print(f"   Gene: {gene['gene_name']}")
        print(f"   Location: {gene['start']}..{gene['end']} ({gene['strand']})")
        print(f"   Length: {gene['length']} bp")


def categorize_genes(genes):
    """
    Categorize genes by type and print summary.
    """
    categories = {}

    for gene in genes:
        ftype = gene["type"]
        if ftype not in categories:
            categories[ftype] = []
        categories[ftype].append(gene)

    print(f"\n{'='*80}")
    print("GENE CATEGORY SUMMARY")
    print(f"{'='*80}\n")

    for ftype in sorted(categories.keys()):
        genes_of_type = categories[ftype]
        print(f"{ftype.upper()}: {len(genes_of_type)} genes")


def main():
    parser = argparse.ArgumentParser(
        description="Extract genes from a specific contig in a BAKTA GBK/GBFF file"
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
        "--format",
        default="table",
        choices=["table", "detailed"],
        help="Output format (default: table)",
    )

    args = parser.parse_args()

    # Extract genes
    genes, contig_length, contig_gc = extract_contig_genes(
        args.gbk_file, args.contig_name, args.type
    )

    # Print results
    if args.format == "table":
        print_genes_table(genes)
    else:
        print_genes_detailed(genes)

    # Print category summary
    categorize_genes(genes)

    # Print statistics
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


if __name__ == "__main__":
    main()
