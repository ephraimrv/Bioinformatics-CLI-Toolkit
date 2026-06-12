#!/usr/bin/env python3
"""
Extract plasmid genes from a BAKTA GBK file
"""

from Bio import SeqIO
import sys


def extract_plasmid_info(gbk_file, plasmid_name="NHJNNNGJ_2"):
    """
    Extract genes from the secondary contig (plasmid) in a GBK file
    """
    plasmid_genes = []

    for record in SeqIO.parse(gbk_file, "genbank"):
        # Check if this is the plasmid contig
        if plasmid_name in record.id or "contig" in record.id.lower():
            print(f"\n=== PLASMID CONTIG: {record.id} ===")
            print(f"Length: {len(record.seq)} bp")
            print(
                f"GC content: {sum(1 for base in record.seq if base in 'GC')/len(record.seq)*100:.2f}%"
            )
            print(
                f"\nGenes found: {len([f for f in record.features if f.type == 'CDS'])}\n"
            )

            # Extract all CDS (genes)
            for feature in record.features:
                if feature.type == "CDS":
                    gene_name = feature.qualifiers.get("product", ["Unknown"])[0]
                    locus_tag = feature.qualifiers.get("locus_tag", ["N/A"])[0]
                    location = f"{feature.location.start}..{feature.location.end}"

                    plasmid_genes.append(
                        {
                            "locus_tag": locus_tag,
                            "product": gene_name,
                            "location": location,
                            "strand": "+" if feature.location.strand == 1 else "-",
                        }
                    )

                    print(f"{locus_tag:15} | {gene_name:50} | {location}")

    return plasmid_genes


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python extract_plasmid.py <gbk_file>")
        sys.exit(1)

    gbk_file = sys.argv[1]
    extract_plasmid_info(gbk_file)
