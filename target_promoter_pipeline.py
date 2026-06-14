"""
Homology-Anchored Regulatory Extractor

Bypasses text annotation errors by hunting for target genes using their
exact amino acid sequences. Once found, extracts the 150bp upstream regulatory region.
"""

import sys
from pathlib import Path
from typing import Iterator
from Bio import SeqIO


def stream_genbank_files(base_dir: Path) -> Iterator[Path]:
    """Yields all GenBank files (.gbk, .gbff)."""
    for ext in ("*.gbk", "*.gbff"):
        yield from base_dir.rglob(ext)


def extract_homology_upstream(
    gbk_path: Path, target_peptides: dict[str, str], upstream_bp: int = 150
) -> Iterator[tuple[str, str, str, str]]:
    """
    Hunts for CDS translations containing specific core peptide strings.
    Yields (Contig ID, Locus Tag, Target Name, Upstream Sequence).
    """
    try:
        with open(gbk_path, "r", encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":

                        # EAFP: Safely grab the translation
                        translation = feature.qualifiers.get("translation", [""])[0]
                        if not translation:
                            continue

                        # Check if any of our target mature peptides are inside this translation
                        for target_name, core_peptide in target_peptides.items():
                            if core_peptide in translation:
                                locus_tag = feature.qualifiers.get(
                                    "locus_tag", ["UNKNOWN"]
                                )[0]
                                start = int(feature.location.start)
                                end = int(feature.location.end)
                                strand = feature.location.strand

                                # Strand-aware slicing
                                if strand == 1:
                                    slice_start = max(0, start - upstream_bp)
                                    upstream_seq = str(record.seq[slice_start:start])
                                else:
                                    slice_end = min(len(record.seq), end + upstream_bp)
                                    raw_seq = record.seq[end:slice_end]
                                    upstream_seq = str(raw_seq.reverse_complement())

                                yield record.id, locus_tag, target_name, upstream_seq

    except Exception as e:
        raise ValueError(f"Failed to parse {gbk_path.name}") from e


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    output_fasta = script_dir / "bacteriocin_upstream_MEME.fasta"

    # THE UPGRADE: We define our targets by their un-mutatable sequences
    # You can add the Lactobin sequence here too if you want to include it in MEME
    targets = {
        "Novel_IIc_Weapon": "VMTGRGLARAIGGGIVGGVIRGIPGGP",  # C5 specific
        # MLNNFEKMNIEQLETVVGGVMTGRGLARAIGGGIVGGVIRGIPGGPAGMFVGAHLGAAAGAATYAVTHY
        "House_Bacteriocin": "FPLLPIVVPIIAGGATYVAKDAWNHLDQIR",  # The sequence you just found (post-GG cleavage)
        # Add the Lactobin A/Cerein 7B mature core here if needed!
    }  # MEKLSEQELAKVSGGFPLLPIVVPIIAGGATYVAKDAWNHLDQIRSGWRKAGNSKW

    upstream_length = 150

    print(f"[*] Scanning project directory: {script_dir}")
    print(f"[*] Hunting by homology to bypass annotation errors...")

    hits_found = 0

    try:
        with open(output_fasta, "w", encoding="utf-8") as out_file:
            for file_path in stream_genbank_files(script_dir):
                print(f"  -> Parsing {file_path.name}...")

                for seq_id, locus, target_name, seq in extract_homology_upstream(
                    file_path, targets, upstream_length
                ):
                    hits_found += 1

                    # Clean FASTA header
                    fasta_header = (
                        f">{seq_id}_{locus}_{target_name}_up{upstream_length}"
                    )
                    out_file.write(f"{fasta_header}\n{seq}\n")
                    print(f"      [Hit] {locus} matches {target_name}")

        print("\n" + "=" * 50)
        print(f"[*] SUCCESS: {hits_found} regulatory regions extracted.")
        print(f"[*] Output saved to: {output_fasta.name}")
        print("=" * 50)

    except ValueError as e:
        import traceback

        print("\n[!] RAW TRACEBACK DUMP:")
        traceback.print_exc()
        sys.exit(1)
    except KeyboardInterrupt:
        sys.exit("\n[!] Pipeline gracefully interrupted by user.")


if __name__ == "__main__":
    main()
