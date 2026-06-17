import sys
from pathlib import Path

try:
    from Bio import SeqIO  # type: ignore[import]  # pylint: disable=import-error
except Exception:
    sys.stderr.write(
        "Error: Biopython is not installed or cannot be imported. Install it with: pip install biopython\n"
    )
    sys.exit(1)

# Define your file path and criteria
script_dir = Path(__file__).resolve().parent
gbff_file = (
    script_dir / "references" / "probiomin_ani" / "GCA_029823215.1_genomic.gbff"
)  # <--- Change this to your actual file name
target_locus = "PUV52_02875"
upstream_bp = 150


def extract_upstream(filename, locus_tag, num_bases):
    # Iterate through all genome records/scaffolds in the gbff file
    for record in SeqIO.parse(filename, "genbank"):
        for feature in record.features:
            # Check features that hold locus tags (CDS or gene)
            if "locus_tag" in feature.qualifiers:
                if locus_tag in feature.qualifiers["locus_tag"]:

                    strand = feature.location.strand
                    start = int(feature.location.start)
                    end = int(feature.location.end)

                    print(f"Found {locus_tag} on record: {record.id}")
                    print(f"Coordinates: {start}-{end} | Strand: {strand}")

                    # Positive (Forward) Strand: Upstream is to the left (before start)
                    if strand == 1:
                        upstream_start = max(0, start - num_bases)
                        upstream_end = start
                        upstream_seq = record.seq[upstream_start:upstream_end]

                    # Negative (Reverse) Strand: Upstream is to the right (after end)
                    # We must take the reverse complement to keep 5' -> 3' orientation
                    elif strand == -1:
                        upstream_start = end
                        upstream_end = min(len(record.seq), end + num_bases)
                        upstream_seq = record.seq[
                            upstream_start:upstream_end
                        ].reverse_complement()

                    else:
                        print("Error: Unknown strand orientation.")
                        return

                    # Output the result in FASTA format
                    print(
                        f"\n>Upstream_{num_bases}bp_of_{locus_tag} | Record: {record.id} | Strand: {strand}"
                    )
                    print(upstream_seq)
                    return

    print(f"Error: Locus tag '{locus_tag}' not found in the file.")


# Run the function
extract_upstream(gbff_file, target_locus, upstream_bp)
