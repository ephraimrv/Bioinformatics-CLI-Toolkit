"""
Genome-Wide Regulon Mapper

Maps transcriptional networks by identifying operator motifs in upstream regions.

The pipeline isolates the upstream sequence (default 150bp) of every CDS in
a GenBank assembly. It performs a regex-based motif search for a provided
IUPAC/Regex operator footprint and compiles matches into a genomic matrix
suitable for network analysis.

Both strands are scanned independently. Motif positions are reported as
negative integers relative to the Translation Start Site (TSS), following
standard molecular biology convention (e.g., -10 and -35 boxes). The motif
strand column (+/-) indicates which DNA strand the binding site was found on.

Eukaryotic note:
    Prokaryotic promoters are tightly packed within 150-300bp upstream.
    In eukaryotes, enhancers and distal regulatory elements can reside
    1,000-50,000bp upstream. Use --upstream 5000 or higher for eukaryotic
    genomes to avoid missing distal regulatory sites.

License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).

Example Usage:
    # Prokaryotic (default upstream is 150bp)
    $ python3 regulon_scanner.py -i C5_genome.gbk -u 200 -m "GCGCAG[CT]G[GT]T[TA]AAAT" -o regulon.tsv

    # Eukaryotic (increase upstream window significantly)
    $ python3 regulon_scanner.py -i yeast_genome.gbff -u 2000 -m "TATAAA" -o regulon.tsv
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.1.2"

import re
import sys
import traceback
from pathlib import Path
from typing import Iterator

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )

from utils import base_parser

# ── IUPAC ambiguity code → regex character class mapping ──────────────────────
# Standard IUPAC-IUB nucleotide ambiguity codes (1970).
# Python's re module treats every letter literally — "W" matches only the
# character W, not [AT]. This table enables transparent auto-translation
# so users can supply either IUPAC notation or regex directly.
_IUPAC_TO_REGEX: dict[str, str] = {
    "R": "[AG]",  # puRine
    "Y": "[CT]",  # pYrimidine
    "S": "[GC]",  # Strong (3 H-bonds)
    "W": "[AT]",  # Weak (2 H-bonds)
    "K": "[GT]",  # Keto
    "M": "[AC]",  # aMino
    "B": "[CGT]",  # not A  (B follows A in alphabet)
    "D": "[AGT]",  # not C  (D follows C)
    "H": "[ACT]",  # not G  (H follows G)
    "V": "[ACG]",  # not T/U (V follows U)
    "N": "[ACGT]",  # aNy base
    "U": "T",  # Uracil → Thymine (RNA ↔ DNA equivalence)
}


def _translate_iupac(pattern: str) -> str:
    """Translate IUPAC ambiguity codes in a motif string to regex equivalents.

    Characters already inside regex brackets ``[...]`` are left untouched
    so that user-defined character classes are preserved. Standard ACGT
    bases and all regex metacharacters pass through unchanged.

    Args:
        pattern: Motif string — pure IUPAC, pure regex, or mixed.

    Returns:
        Regex-compatible pattern with all IUPAC codes expanded.

    Examples::

        "TATAWAW"          →  "TATA[AT]A[AT]"
        "GCGCAG[CT]GWTAAAT" →  "GCGCAG[CT]G[AT]TAAAT"
        "TATAAA"           →  "TATAAA"    (plain ACGT — unchanged)
        "GCGCAG[CT]G[GT]T" →  "GCGCAG[CT]G[GT]T"  (brackets protected)
    """
    result = []
    in_brackets = 0

    for i, char in enumerate(pattern):
        prev = pattern[i - 1] if i > 0 else ""

        if char == "[" and prev != "\\":
            in_brackets += 1
            result.append(char)
        elif char == "]" and in_brackets > 0 and prev != "\\":
            in_brackets -= 1
            result.append(char)
        elif in_brackets == 0 and char.upper() in _IUPAC_TO_REGEX:
            result.append(_IUPAC_TO_REGEX[char.upper()])
        else:
            result.append(char)

    return "".join(result)


def stream_regulon_hits(
    gbk_path: Path, regex_pattern: str, upstream_bp: int
) -> Iterator[dict]:
    """Scans every CDS upstream region for a motif on both DNA strands.

    Motif positions are returned as negative integers relative to the
    Translation Start Site (TSS), following standard molecular biology
    convention (e.g., the -10 and -35 boxes in prokaryotic promoters).

    Both the coding strand (+) and the template strand (-) are scanned.
    This ensures Transcription Factor binding sites in either orientation
    are detected, including palindromic and non-palindromic motifs.

    Args:
        gbk_path:       Path to the GenBank file.
        regex_pattern:  IUPAC/regex motif string (IGNORECASE applied).
        upstream_bp:    Number of bases upstream of each CDS start to extract.

    Yields:
        A dict per CDS with at least one motif hit, containing locus_tag,
        product, contig, gene strand, and a sorted list of
        (rel_pos, matched_seq, motif_strand) tuples.
    """
    # Translate IUPAC ambiguity codes before compiling.
    # Python's re module treats letters literally — "W" matches only the
    # character W, not A or T. _translate_iupac() converts IUPAC codes to
    # their regex equivalents while leaving regex metacharacters untouched.
    translated_pattern = _translate_iupac(regex_pattern)
    try:
        safe_pattern = re.compile(f"(?=({translated_pattern}))", re.IGNORECASE)
    except re.error as e:
        raise ValueError(
            f"Invalid regex/IUPAC pattern: '{regex_pattern}' "
            f"(translated: '{translated_pattern}')"
        ) from e

    try:
        for record in SeqIO.parse(gbk_path, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":

                    # Use coordinate-based fallback for missing locus_tag.
                    # "UNKNOWN" would be ambiguous in TSV output if multiple
                    # CDS features lack annotation.
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    strand = feature.location.strand
                    locus_tag = feature.qualifiers.get(
                        "locus_tag", [f"UNKNOWN_CDS_{start}"]
                    )[0]
                    product = feature.qualifiers.get(
                        "product", ["hypothetical protein"]
                    )[0]

                    # Extract upstream with boundary tracking
                    if strand == 1:
                        slice_start = max(0, start - upstream_bp)
                        actual_upstream = start - slice_start
                        upstream_seq = str(record.seq[slice_start:start])
                    else:
                        slice_end = min(len(record.seq), end + upstream_bp)
                        actual_upstream = slice_end - end
                        raw_seq = record.seq[end:slice_end]
                        upstream_seq = str(raw_seq.reverse_complement())

                    if actual_upstream < upstream_bp:
                        print(
                            f"    [!] Warning: {locus_tag} upstream truncated to "
                            f"{actual_upstream}bp (contig boundary).",
                            file=sys.stderr,
                        )

                    matches = []

                    # Forward (coding) strand scan
                    # Position reported as negative distance from TSS:
                    #   match.start() 0  → -(actual_upstream)  (farthest from ATG)
                    #   match.start() L-1 → -1                 (one base before ATG)
                    for match in safe_pattern.finditer(upstream_seq):
                        rel_pos = -(actual_upstream - match.start())
                        matches.append((rel_pos, match.group(1), "+"))

                    # Reverse complement (template) strand scan
                    # The RC of the upstream sequence is scanned with the same pattern.
                    # The coordinate is mapped back to the forward-strand TSS origin:
                    #
                    # IMPORTANT: Because we use a zero-width lookahead assertion (?=(...)),
                    # match.end() always equals match.start() — the outer match consumes
                    # zero characters. Using match.end() directly would place every RC hit
                    # at the wrong position (off by the motif length). We must calculate
                    # the true end from the captured group's length instead.
                    #
                    # Math note: len(upstream_seq) == actual_upstream always (the sequence
                    # IS the extracted upstream window), so the formula simplifies:
                    #   fwd_pos = len(upstream_seq) - true_match_end
                    #   rel_pos = -(actual_upstream - fwd_pos)
                    #           = -(actual_upstream - actual_upstream + true_match_end)
                    #           = -true_match_end
                    rc_seq = str(Seq(upstream_seq).reverse_complement())
                    for match in safe_pattern.finditer(rc_seq):
                        true_match_end = match.start() + len(match.group(1))
                        rel_pos = -true_match_end
                        matches.append((rel_pos, match.group(1), "-"))

                    if matches:
                        # Sort biologically: 5' → 3' relative to TSS (most negative first)
                        matches.sort(key=lambda x: x[0])
                        yield {
                            "locus_tag": locus_tag,
                            "product": product,
                            "contig": record.id,
                            "strand": strand,
                            "matches": matches,
                        }

    except (FileNotFoundError, OSError, UnicodeDecodeError, ValueError) as e:
        raise ValueError(f"GenBank Parsing Error in {gbk_path.name}: {e}") from e


def main() -> None:
    parser = base_parser("Genome-Wide Regulon Scanner")
    parser.add_argument(
        "-u",
        "--upstream",
        type=int,
        default=150,
        help=(
            "Bases upstream of each CDS to extract and scan. "
            "Default: 150 (appropriate for prokaryotes). "
            "For eukaryotes, enhancers can reside 1,000-50,000bp upstream — "
            "use --upstream 5000 or higher for eukaryotic genomes."
        ),
    )
    parser.add_argument(
        "-m",
        "--motif",
        required=True,
        help=(
            "Motif to search for on both strands. Accepts three formats: "
            "(1) Pure regex: 'TATA[AT][AT]'. "
            "(2) IUPAC ambiguity codes: 'TATAWW' — automatically translated to "
            "regex before scanning (W=[AT], R=[AG], Y=[CT], S=[GC], K=[GT], "
            "M=[AC], B=[CGT], D=[AGT], H=[ACT], V=[ACG], N=[ACGT]). "
            "(3) Mixed: 'GCGCAG[CT]GWTTAAAT' (brackets protect existing regex). "
            "The translated pattern is printed to stderr for verification."
        ),
    )
    args = parser.parse_args()

    print(f"[*] Scanning genome    : {args.input.name}", file=sys.stderr)
    print(f"[*] Upstream region    : {args.upstream}bp", file=sys.stderr)
    print(f"[*] Motif (input)      : {args.motif}", file=sys.stderr)
    translated_motif = _translate_iupac(args.motif)
    if translated_motif != args.motif:
        print(
            f"[*] Motif (translated) : {translated_motif}  ← IUPAC codes expanded",
            file=sys.stderr,
        )
    print(f"[*] Strands scanned    : Both (+) coding and (-) template", file=sys.stderr)
    print(f"[*] Position reference : TSS (negative = upstream of ATG)", file=sys.stderr)

    total_genes_hit = 0
    total_motifs_found = 0

    try:
        results_iterator = stream_regulon_hits(args.input, args.motif, args.upstream)

        if args.output:
            with open(args.output, "w", encoding="utf-8") as tsv:
                tsv.write(
                    "Locus_Tag\tContig\tGene_Strand\tMotif_Count\t"
                    "Positions_Relative_to_TSS\tMatched_Sequences\tProduct\n"
                )
                for hit in results_iterator:
                    total_genes_hit += 1
                    total_motifs_found += len(hit["matches"])
                    # Format: -35(+), -10(-) — position(motif_strand)
                    positions = ",".join([f"{m[0]}({m[2]})" for m in hit["matches"]])
                    sequences = ",".join([m[1] for m in hit["matches"]])
                    tsv.write(
                        f"{hit['locus_tag']}\t{hit['contig']}\t{hit['strand']}\t"
                        f"{len(hit['matches'])}\t{positions}\t{sequences}\t{hit['product']}\n"
                    )
                    print(
                        f"    -> Regulon member found: {hit['locus_tag']} "
                        f"({hit['product'][:40]}...)",
                        file=sys.stderr,
                    )
        else:
            for hit in results_iterator:
                total_genes_hit += 1
                total_motifs_found += len(hit["matches"])
                print(
                    f"    -> Regulon member found: {hit['locus_tag']} "
                    f"({hit['product'][:40]}...)",
                    file=sys.stderr,
                )
            print(
                "\n[*] Note: No output file specified (-o). Results printed to terminal only.",
                file=sys.stderr,
            )

        print("\n" + "=" * 40, file=sys.stderr)
        print("          PIPELINE COMPLETE", file=sys.stderr)
        print("=" * 40, file=sys.stderr)
        print(f"Total Genes in Regulon : {total_genes_hit}", file=sys.stderr)
        print(f"Total Motifs Bound     : {total_motifs_found}", file=sys.stderr)
        if args.output:
            print(f"Results written to     : {args.output.resolve()}", file=sys.stderr)
        print("=" * 40, file=sys.stderr)

    except (ValueError, FileNotFoundError, PermissionError) as e:
        sys.exit(f"\n[!] Pipeline Halted: {e}")
    except KeyboardInterrupt:
        sys.exit("\n[!] Pipeline gracefully interrupted by user.")
    except Exception:
        print("\n[!] UNEXPECTED BUG ENCOUNTERED:", file=sys.stderr)
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
