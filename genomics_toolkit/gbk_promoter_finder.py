#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

"""GenBank Upstream Promoter Finder

Extracts regulatory promoter regions from multi-contig GenBank assemblies.

This tool locates a target gene by its locus tag, calculates the strand
orientation to extract the correct upstream sequence (applying reverse
complementation where necessary), and scans the region for specific motif hits.

PROKARYOTE-ONLY ANCHOR — NOT A WINDOW-SIZE PROBLEM:
    This script anchors the upstream window on the CDS start (the
    translation start / ATG) of whichever feature utils.py resolves for
    the given locus tag — not on the Transcription Start Site (TSS). In
    prokaryotes these coincide, since there is no 5' UTR separating them.
    In eukaryotes they do not: the TSS sits upstream of the CDS start,
    often separated by a 5' UTR that itself contains introns.

    Increasing --upstream on a eukaryotic genome does NOT fix this — it
    just extracts a longer stretch of 5' UTR/intron sequence anchored at
    the wrong coordinate, not the actual promoter. This script previously
    suggested "--upstream 2000 or higher" for eukaryotic use with no such
    caveat; that guidance was incomplete and has been corrected as of
    v1.2.1 (see changelog below). There is no eukaryote mode here, the
    same way there is none in the sibling script regulon_scanner.py,
    which uses the identical CDS-anchored mechanism and already documents
    this limitation. For eukaryotic promoter/regulatory work, extract
    upstream regions with universal_promoter_extractor.py or
    target_promoter_pipeline.py instead (both TSS-anchored, resolving the
    TSS from mRNA features across isoforms), then search that output with
    MEME/FIMO directly.

Note:
    Associated with ongoing, unpublished research (manuscript in
    preparation). Correct attribution is requested when used in
    derivative works.

    v1.2.0: Fixed three issues found during code review.
    (1) ``-m/--motif`` advertised IUPAC support but raw ambiguity codes
    (W, R, Y, etc.) were passed directly to ``re.compile()``, which has no
    concept of them — "TATAWAW" matched zero times against a sequence
    containing the valid instance "TATAAAA".
    Now uses ``utils.translate_iupac_to_regex()`` (shared with
    ``regulon_scanner.py``, which had the identical bug) before compiling.
    (2) The motif TSV table was previously appended to the bottom of the
    output ``.fasta`` file, corrupting it for any standard FASTA parser
    (BLAST, MEME, FIMO) — the appended lines contain characters
    outside any valid sequence alphabet. Motifs now go to a sibling
    ``.tsv`` file instead, matching ``target_promoter_pipeline.py``.
    (3) ``_get_genome_label()`` only matched ``CDS`` features when looking
    up the organism/strain label, so any mRNA-only or non-coding RNA
    locus tag fell back to the generic "unknown" label. The actual
    coordinate extraction (``utils.extract_upstream_sequence()``) had the
    same restriction and was fixed the same way in utils.py v1.4.0, which
    this script now relies on transitively — without that fix, this
    script would have raised "Locus tag not found" before ever reaching
    the label lookup, regardless of this function's own fix.

    v1.2.1: Corrected misleading eukaryotic guidance, after comparing
    this script against its sibling regulon_scanner.py. Both
    scripts ultimately anchor their upstream window on a CDS-start
    coordinate (this script via utils.extract_upstream_sequence(), which
    delegates to utils.extract_upstream_sequence_with_length()'s CDS-first
    feature resolution) — but regulon_scanner.py already documented
    itself as PROKARYOTE-ONLY for exactly that reason, while this
    script's own --upstream help text read "For eukaryotes, consider
    --upstream 2000 or higher," implying a bigger window alone makes the
    result eukaryote-correct. It does not: CDS start and the TSS are
    different coordinates in eukaryotes, separated by a 5' UTR (itself
    possibly containing introns), and no amount of extra upstream bp
    recovers the correct anchor point. Added the PROKARYOTE-ONLY ANCHOR
    docstring section above (mirroring regulon_scanner.py's wording,
    since it is the same root cause), corrected the --upstream help text,
    and added a one-time runtime warning — via the new shared
    utils.looks_eukaryotic() heuristic (mRNA-feature detection) — printed
    before extraction proceeds, pointing to universal_promoter_extractor.py
    / target_promoter_pipeline.py for actual eukaryotic TSS-anchored
    extraction. No change to prokaryote behavior or output format.

Example Usage:
    $ python3 gbk_promoter_finder.py -i C5_genome.gbk -l ctg1_50 -u 150 -m "TATAAT" -o ctg1_50_promoter.fasta
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.2.1"

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
from utils import (
    base_parser,
    wrap_fasta,
    extract_upstream_sequence,
    translate_iupac_to_regex,
    looks_eukaryotic,
)


def _get_genome_label(gbk_path: Path, locus_tag: str) -> tuple[str, str]:
    """Extract record ID and organism/strain label for a given locus tag.

    Scans the GenBank file to find the record containing the target locus,
    then pulls the organism name and strain qualifier from the source feature
    or record annotations. Strain is appended only when not already embedded
    in the organism string (e.g. NCBI already includes it in many entries).

    Args:
        gbk_path:  Path to the GenBank file.
        locus_tag: Target locus tag to locate the correct record.

    Returns:
        Tuple of (record_id, genome_label).
        genome_label is "Organism strain" or just "Organism" if strain is
        absent or already included. Falls back gracefully if not found.
    """
    try:
        for record in SeqIO.parse(gbk_path, "genbank"):
            for feature in record.features:
                # No feature-type restriction here deliberately: a
                # /locus_tag is unique within a genome regardless of which
                # feature type carries it. The previous CDS-only check
                # caused this function to fall back to "unknown" for any
                # mRNA-only or non-coding RNA (tRNA/rRNA/ncRNA) locus tag,
                # even though utils.extract_upstream_sequence() (called
                # before this function, and fixed the same way) already
                # successfully resolves those same locus tags.
                if locus_tag not in feature.qualifiers.get("locus_tag", []):
                    continue

                # Found the correct record — extract organism and strain
                organism = ""
                strain = ""

                # The /source feature is the most reliable location
                for src in record.features:
                    if src.type == "source":
                        organism = src.qualifiers.get("organism", [""])[0]
                        strain = src.qualifiers.get("strain", [""])[0]
                        if not strain:
                            strain = src.qualifiers.get("isolate", [""])[0]
                        break

                # Fall back to record-level annotations
                if not organism:
                    organism = record.annotations.get("organism", "")
                if not strain:
                    strain = record.annotations.get("strain", "")

                # Skip Prokka placeholder values ("." or blank)
                if organism in ("", "."):
                    organism = ""
                if strain in ("", "."):
                    strain = ""

                if organism and strain and strain not in organism:
                    label = f"{organism} {strain}"
                elif organism:
                    label = organism
                else:
                    # Prokka/local assemblies: use the file stem as fallback
                    label = gbk_path.stem

                return record.id, label.strip()

    except Exception:
        pass

    return "unknown", gbk_path.stem


def find_motif_regex_iterator(
    sequence: str, regex_pattern: str, actual_len: int
) -> Iterator[tuple[int, str, str]]:
    """Scans a sequence for a motif on both strands.

    Positions are reported as negative integers relative to the Translation
    Start Site (TSS), following standard molecular biology convention
    (e.g., -35 means 35 bases upstream of the ATG).

    Both the coding (+) and template (-) strands are scanned so TF binding
    sites in either orientation are detected.

    Args:
        sequence:       The upstream DNA sequence string (coding strand, 5'→3').
        regex_pattern:  IUPAC/regex motif (IGNORECASE applied).
        actual_len:     Actual length of the upstream sequence (may be shorter
                        than requested if near a contig boundary).

    Yields:
        Tuple of (rel_pos, matched_sequence, strand_indicator).
        rel_pos is a negative integer relative to TSS.
        strand_indicator is '+' or '-'.
    """
    if not sequence or not regex_pattern:
        return
    # Translate IUPAC ambiguity codes (W, R, Y, S, K, M, B, D, H, V, N) into
    # regex character classes before compiling. Python's re module has no
    # concept of these codes — compiling the raw pattern would search for
    # the literal letter (e.g. "W"), which never appears in a real DNA
    # sequence, so any motif using ambiguity codes would silently match
    # nothing — "TATAWAW" found zero hits against a
    # sequence containing the valid TATA-box instance "TATAAAA" before this
    # fix. Hand-written regex syntax in the motif (brackets, quantifiers,
    # groups) is left untouched — see translate_iupac_to_regex()'s docstring.
    translated_pattern = translate_iupac_to_regex(regex_pattern)
    try:
        safe_pattern = re.compile(f"(?=({translated_pattern}))", re.IGNORECASE)
    except re.error as e:
        raise ValueError(f"Invalid regex pattern provided: {regex_pattern}") from e

    # Forward (coding) strand scan
    for match in safe_pattern.finditer(sequence):
        rel_pos = -(actual_len - match.start())
        yield rel_pos, match.group(1), "+"

    # Reverse complement (template) strand scan
    rc_seq = str(Seq(sequence).reverse_complement())
    for match in safe_pattern.finditer(rc_seq):
        # IMPORTANT: Because we use a zero-width lookahead assertion (?=(...)),
        # match.end() always equals match.start() — the outer match consumes
        # zero characters. Using match.end() directly would place every RC hit
        # at the wrong position (off by the motif length). We must calculate
        # the true end from the captured group's length instead.
        #
        # Math note: len(sequence) == actual_len always (sequence IS the extracted
        # upstream), so the two-step formula simplifies:
        #   fwd_pos = len(sequence) - true_match_end
        #   rel_pos = -(actual_len - fwd_pos)
        #           = -(actual_len - actual_len + true_match_end)
        #           = -true_match_end
        true_match_end = match.start() + len(match.group(1))
        rel_pos = -true_match_end
        yield rel_pos, match.group(1), "-"


def main() -> None:
    parser = base_parser("GenBank Targeted Upstream Motif Scanner")
    parser.add_argument("-l", "--locus", required=True, help="Target gene locus tag")
    parser.add_argument(
        "-u",
        "--upstream",
        type=int,
        default=100,
        help=(
            "Upstream bases to extract. Default: 100. NOTE: this script "
            "anchors on CDS start, not the transcription start site (TSS) "
            "— increasing this value does NOT adapt it for eukaryotic use "
            "(see the PROKARYOTE-ONLY ANCHOR note in this script's "
            "docstring). Use universal_promoter_extractor.py for "
            "eukaryotic, TSS-anchored upstream extraction."
        ),
    )
    parser.add_argument(
        "-m",
        "--motif",
        required=False,
        help="Regex/IUPAC motif to search for on both strands (optional)",
    )
    args = parser.parse_args()

    try:
        if args.upstream < 1:
            raise ValueError("Upstream bases must be a positive integer.")

        # One-time heuristic warning: mRNA features are a strong signal of
        # eukaryotic annotation (Prokka/Bakta prokaryote output never emits
        # them). This script is CDS-anchored, not TSS-anchored (see the
        # PROKARYOTE-ONLY ANCHOR docstring section) — regulon_scanner.py
        # already warns on this exact signal; utils.looks_eukaryotic() is
        # the shared helper extracted so both scripts (and any future one)
        # stay consistent without duplicating the scan logic.
        if looks_eukaryotic(args.input):
            print(
                "[!] Warning: mRNA features detected — this looks like a "
                "eukaryotic genome. This script anchors the upstream "
                "window on CDS start, not the transcription start site "
                "(TSS), so the true promoter will likely be missed. See "
                "the PROKARYOTE-ONLY ANCHOR note in this script's "
                "docstring. For eukaryotic promoter extraction, use "
                "universal_promoter_extractor.py or "
                "target_promoter_pipeline.py instead.",
                file=sys.stderr,
            )

        upstream_seq, start, end, strand = extract_upstream_sequence(
            args.input, args.locus, args.upstream
        )

        actual_len = len(upstream_seq)
        record_id, genome_label = _get_genome_label(args.input, args.locus)
        strand_symbol = "+" if strand == 1 else "-"

        print(
            f"[*] Found {args.locus} at {start}-{end} (Gene strand: {strand})",
            file=sys.stderr,
        )
        print(
            f"[*] Requested: {args.upstream}bp upstream | Extracted: {actual_len}bp",
            file=sys.stderr,
        )

        if actual_len < args.upstream:
            print(
                f"[!] Warning: Upstream truncated to {actual_len}bp (contig boundary).",
                file=sys.stderr,
            )

        motifs = []
        if args.motif:
            print(
                f"[*] Searching for motif: {args.motif} (both strands)",
                file=sys.stderr,
            )
            # Sort results biologically: most negative (farthest from TSS) first
            motifs = sorted(
                find_motif_regex_iterator(upstream_seq, args.motif, actual_len),
                key=lambda x: x[0],
            )

        if args.output:
            # Write strict FASTA: header + sequence ONLY. Previously a TSV
            # motif table was appended below the sequence in the same
            # .fasta file. FASTA is a strict two-line-type format (header
            # lines starting with '>', sequence lines); appending tab-
            # delimited text with a '#' comment line corrupts the file for
            # any standard parser (BLAST, MEME, FIMO) — the appended lines
            # contain characters outside any valid sequence alphabet. The
            # motif table now goes to a sibling .tsv file instead, matching
            # the pattern already used in target_promoter_pipeline.py.
            with open(args.output, "w", encoding="utf-8") as out_file:
                # NCBI-style FASTA header with | separators
                fasta_header = (
                    f">{args.locus}"
                    f" | {record_id}"
                    f" | {genome_label}"
                    f" | {actual_len}bp upstream"
                    f" | strand {strand_symbol}"
                )
                out_file.write(f"{fasta_header}\n")
                out_file.write(f"{wrap_fasta(upstream_seq)}\n")

            print(
                f"[*] Success! FASTA written to {args.output.resolve()}",
                file=sys.stderr,
            )

            if motifs:
                tsv_path = args.output.with_suffix(".tsv")
                with open(tsv_path, "w", encoding="utf-8") as tsv_file:
                    tsv_file.write("Position_Relative_to_TSS\tMotif_Strand\tSequence\n")
                    for pos, seq, motif_strand in motifs:
                        tsv_file.write(f"{pos}\t{motif_strand}\t{seq}\n")
                print(
                    f"[*] {len(motifs)} motif(s) found. Written to {tsv_path.resolve()}",
                    file=sys.stderr,
                )
            else:
                print("[*] No motifs found.", file=sys.stderr)

        else:
            print("\n--- UPSTREAM SEQUENCE ---", file=sys.stderr)
            if len(upstream_seq) > 500:
                print(
                    f"{upstream_seq[:100]} ... [snip {len(upstream_seq)-200}bp] ... {upstream_seq[-100:]}",
                    file=sys.stderr,
                )
            else:
                print(upstream_seq, file=sys.stderr)
            print("-------------------------\n", file=sys.stderr)

            if motifs:
                for pos, seq, motif_strand in motifs:
                    print(
                        f"    -> Motif Found! Position: {pos} ({motif_strand} strand) | Sequence: {seq}",
                        file=sys.stderr,
                    )
            elif args.motif is None:
                pass
            else:
                print("    -> No motifs found.", file=sys.stderr)

    except FileNotFoundError:
        sys.exit(f"\n[!] File not found: {args.input}")
    except ValueError as e:
        sys.exit(f"\n[!] Error: {e}")
    except KeyboardInterrupt:
        sys.exit("\n[!] Pipeline gracefully interrupted by user.")
    except Exception:
        print("\n[!] UNEXPECTED BUG ENCOUNTERED:")
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
