#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

"""Interactive Protein Presence/Absence Scanner

An interactive REPL tool for scanning reference genomes for the presence
or absence of a protein of interest.

Paste a full protein sequence and press Enter. The tool automatically
calculates the mature peptide core (stripping signal peptides via
biochemical cleavage logic), then performs an exact substring search
across all reference genomes. Results are reported as a PRESENT/ABSENT
matrix per genome, with locus tag and product for each hit.

Important:
    This tool uses EXACT substring matching. A single amino acid substitution
    will cause a miss. For divergent homolog detection, use pairwise_homolog_finder.py.

    Use --raw to skip bacteriocin core trimming and search with the full
    pasted sequence. Required for non-bacteriocin targets.

    SCOPE — bacterial cleavage logic only: calculate_mature_core() models
    bacterial Sec/Tat/RiPP-leader cleavage rules. It has no knowledge of
    eukaryotic secretory-pathway signal peptides (ER/Golgi-targeted),
    which use different cleavage motifs entirely. There is no reliable
    way to auto-detect this from a pasted sequence alone (it carries no
    organism metadata), so always pass --raw when the query peptide is
    not a bacterial bacteriocin/RiPP — otherwise the trimmed "core" may
    be wrong and the exact-match search will silently miss real hits.

    Only GenBank files (.gbk, .gbff) and protein FASTA files (.faa, .mpfa) are
    supported as references. Raw nucleotide FASTA files (.fa, .fasta) contain
    DNA strings — protein probes will never match them. These files are
    automatically detected and skipped with a warning.

    For eukaryotic genomes: use annotated GenBank files (NCBI .gbff) rather
    than raw genomic FASTA. Unannotated DNA files cannot be protein-searched
    without six-frame translation.

Note:
    Associated with ongoing, unpublished research (manuscript in
    preparation). Correct attribution is requested when used in
    derivative works.

    v1.3.0: Fixed two bugs found during review.
    (1) Pasting a standard FASTA record (">header" line + wrapped
    sequence lines — the format any sequence comes in when copied from
    NCBI/UniProt) into the plain input()-based prompt had no defense at
    all: depending on terminal/readline paste handling, it either split
    into multiple nonsensical line-by-line queries, or merged into one
    string containing the header text and an embedded newline that can
    never match any real translation — a silent, guaranteed false
    ABSENT result with no indication that the paste FORMAT, not the
    target protein, was the problem. Added ``_read_query_sequence()``,
    which detects either case (an embedded newline in a single input, or
    a bare ">"-prefixed line arriving on its own) and strips the header/
    joins the wrapped lines before treating the result as the query.
    (2) ``scan_for_peptide()`` re-raised parsing failures as ``ValueError``,
    which propagated straight out of the per-file loop in main() — one
    unreadable or malformed reference file silently aborted the scan for
    every file after it in the same directory, with no indication those
    later files were never attempted. This broke from this toolkit's own
    established defensive pattern (see conserved_annotation_scanner.py's
    docstring: "Returns an empty list on failure to allow batch scanning
    to continue"). Fixed by extracting the duplicated with-output/
    without-output scanning loops into one shared ``_scan_and_report()``
    helper that catches ``ValueError`` per file and continues to the
    next one, writing an "ERROR" status row to the TSV (if any) so the
    failure is visible rather than the file just silently vanishing from
    the matrix. Neither fix changes the SCOPE caveat, the .mpfa/.faa
    reference handling, or any other behavior documented above.

Examples:
    # Interactive scan against a directory of reference genomes
    $ python3 protein_presence_scanner.py -i references/ -o presence_matrix.tsv

    # For non-bacteriocin proteins (skip core trimming)
    $ python3 protein_presence_scanner.py -i references/ --raw -o presence_matrix.tsv
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.3.0"

import sys
import argparse
from pathlib import Path
from typing import TextIO

try:
    from Bio import SeqIO
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )
from utils import stream_reference_files, calculate_mature_core


def scan_for_peptide(gbk_path: Path, target_peptide: str) -> list[tuple[str, str]]:
    hits = []
    try:
        suffix = gbk_path.suffix.lower()
        is_protein_fasta = suffix in (".faa", ".mpfa")
        is_nucleotide_fasta = suffix in (".fasta", ".fa", ".fna")

        # Nucleotide FASTA files contain DNA strings (ATCG...).
        # A protein probe (MKKTLV...) will never match them — silently
        # returning ABSENT for every genome. Catch this and warn explicitly.
        if is_nucleotide_fasta:
            print(
                f"  [!] Skipping {gbk_path.name}: nucleotide FASTA files cannot be "
                f"protein-searched. Provide a GenBank (.gbff/.gbk) or protein FASTA (.faa) instead.",
                file=sys.stderr,
            )
            return hits

        if is_protein_fasta:
            with open(gbk_path, "r", encoding="utf-8") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    translation = str(record.seq).upper()
                    if target_peptide in translation:
                        full_header_desc = record.description.replace(
                            record.id, ""
                        ).strip()
                        product = (
                            full_header_desc
                            if full_header_desc
                            else "Unannotated FASTA sequence"
                        )
                        hits.append((record.id, product))
        else:
            # GenBank format
            with open(gbk_path, "r", encoding="utf-8") as handle:
                for record in SeqIO.parse(handle, "genbank"):
                    for feature in record.features:
                        if feature.type == "CDS":
                            translation = feature.qualifiers.get("translation", [""])[
                                0
                            ].upper()
                            if target_peptide in translation:
                                locus_tag = feature.qualifiers.get(
                                    "locus_tag", ["UNKNOWN"]
                                )[0]
                                product = feature.qualifiers.get(
                                    "product", ["Unknown product"]
                                )[0]
                                hits.append((locus_tag, product))
        return hits

    except (FileNotFoundError, OSError, UnicodeDecodeError, ValueError) as e:
        raise ValueError(f"Failed to parse or extract from {gbk_path.name}: {e}") from e


def _read_query_sequence() -> str | None:
    """Reads one full pasted protein sequence from the user, tolerating
    both a bare amino-acid string and a FASTA-formatted paste (a
    ">header" line followed by sequence, often wrapped across multiple
    lines — the standard format any sequence comes in when copied from
    NCBI/UniProt).

    Pasting a FASTA record into a plain input()-based prompt can fail in
    one of two ways depending on the terminal/readline's paste handling:
    either each line arrives as its own separate input() call (so the
    header line and each wrapped sequence fragment would previously have
    been treated as independent, nonsensical queries), or the whole
    block arrives merged into one string containing the header text and
    an embedded newline (so the previous code would search for that
    literal, header-contaminated string, which can never match any real
    translation — a silent, guaranteed false ABSENT result with no
    indication that the paste format, not the target protein, was the
    problem). This function handles both cases by detecting a leading
    ">" — either as the start of an embedded multi-line block, or as a
    bare line of its own — and stripping it before joining the remaining
    sequence lines into one continuous string.

    Note: a sequence pasted WITHOUT a ">" header but still wrapped across
    multiple lines is only handled correctly if the terminal delivers it
    as one merged multi-line string (Case A below). If lines arrive
    separately with no header to signal "more lines follow," there is no
    reliable way to distinguish "this is the whole query" from "more is
    coming" — paste it as a single continuous line in that case, or
    include the ">" header so this function knows to keep reading.

    Returns:
        The cleaned, header-stripped, single-line amino acid sequence, or
        None if the user typed a quit command or hit EOF.
    """
    try:
        first_line = input("Paste protein sequence: ").strip()
    except EOFError:
        return None

    if first_line.lower() in ("quit", "exit", "q"):
        return None
    if not first_line:
        return ""

    # Case A: the whole pasted block (header and all) arrived as ONE
    # string with embedded newlines.
    if "\n" in first_line:
        lines = [ln.strip() for ln in first_line.splitlines() if ln.strip()]
        seq_lines = [ln for ln in lines if not ln.startswith(">")]
        return "".join(seq_lines)

    # Case B: a bare FASTA header arrived as its own line — the terminal
    # fed the paste one line at a time. Keep reading subsequent lines as
    # sequence fragments until a blank line, quit command, or EOF.
    if first_line.startswith(">"):
        print(
            "  [i] FASTA header detected — reading wrapped sequence lines "
            "below it (blank line or 'quit' to stop)...",
            file=sys.stderr,
        )
        seq_lines = []
        while True:
            try:
                line = input().strip()
            except EOFError:
                break
            if not line or line.lower() in ("quit", "exit", "q"):
                break
            seq_lines.append(line)
        return "".join(seq_lines)

    # Case C: a single bare amino-acid line — the original, always-
    # supported usage pattern, unchanged.
    return first_line


def _scan_and_report(
    file_paths,
    core_target: str,
    tsv: TextIO | None = None,
) -> int:
    """Scans every reference file for core_target, reporting PRESENT/ABSENT
    status per file and (if ``tsv`` is given) writing matrix rows.

    Each file's scan is wrapped in its own try/except so that one
    unreadable or malformed reference file does not abort the scan for
    every other file in the batch — before
    this fix, scan_for_peptide()'s re-raised ValueError on a single bad
    file propagated out of the per-file loop entirely, silently
    truncating the matrix for every file after the one that failed, with
    no indication to the user that those later files were never actually
    scanned. This mirrors the defensive pattern
    conserved_annotation_scanner.py already uses for the same reason
    (its own docstring: "Returns an empty list on failure to allow batch
    scanning to continue"). Also de-duplicates what were previously two
    near-identical with-output/without-output loops in main() into one
    shared implementation.

    Args:
        file_paths: Iterable of reference file paths to scan.
        core_target: The (already core-trimmed or raw) peptide to search for.
        tsv: Open file handle to write matrix rows to, or None to skip
            writing (terminal-only mode).

    Returns:
        Total number of hits found across all files.
    """
    total_hits = 0
    for file_path in file_paths:
        print(f"  [*] Scanning {file_path.name}...", file=sys.stderr)
        try:
            hits = scan_for_peptide(file_path, core_target)
        except ValueError as e:
            print(f"      [X] Error scanning {file_path.name}: {e}", file=sys.stderr)
            print(
                "      [i] Skipping this file and continuing with the "
                "remaining reference files.",
                file=sys.stderr,
            )
            if tsv:
                tsv.write(f"{core_target}\t{file_path.name}\t-\t-\tERROR\n")
            continue

        if hits:
            print(
                f"      [!] ALERT: Found {len(hits)} match(es) in {file_path.name}",
                file=sys.stderr,
            )
            total_hits += len(hits)
            for locus, product in hits:
                short_prod = product[:45] + "..." if len(product) > 45 else product
                print(
                    f"          -> Locus: {locus:<15} | Product: {short_prod}",
                    file=sys.stderr,
                )
                if tsv:
                    tsv.write(
                        f"{core_target}\t{file_path.name}\t{locus}\t{product}\tPRESENT\n"
                    )
        else:
            print(f"      [✓] ABSENT in {file_path.name}", file=sys.stderr)
            if tsv:
                tsv.write(f"{core_target}\t{file_path.name}\t-\t-\tABSENT\n")

    return total_hits


def main() -> None:
    parser = argparse.ArgumentParser(description="Bacteriocin Presence/Absence Hunter")
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        default=Path("."),
        help="Input GenBank file OR a directory to scan (Default: current directory)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=False,
        help="Output TSV file for the presence/absence matrix",
    )
    parser.add_argument(
        "--raw",
        action="store_true",
        help=(
            "Skip bacteriocin core trimming and use the full pasted sequence as the probe. "
            "Required for non-bacteriocin targets (housekeeping genes, TFs, kinases, etc.)."
        ),
    )
    args = parser.parse_args()

    print("=========================================", file=sys.stderr)
    print(f"  Target Scope: {args.input.name}", file=sys.stderr)
    if args.output:
        print(f"  Output File:  {args.output.resolve()}", file=sys.stderr)
    print("=========================================", file=sys.stderr)
    print("Paste a full protein sequence and press Enter.", file=sys.stderr)
    print("Type 'quit' or press Ctrl+C to exit.\n", file=sys.stderr)

    if args.raw:
        print(
            "[*] Mode: --raw (full sequences used, no core trimming)", file=sys.stderr
        )
    else:
        print(
            "[*] Mode: bacteriocin core trimming (bacterial Sec/Tat/RiPP "
            "cleavage rules; NOT valid for eukaryotic secretome proteins "
            "— use --raw for those)",
            file=sys.stderr,
        )

    # Track whether this is the first successful query in this session.
    # We open the output file in "w" mode for the first query (creating/overwriting)
    # and "a" mode for all subsequent queries (appending rows to same file).
    # This prevents the session eraser bug where "w" inside the loop would
    # wipe the previous query's results each time a new sequence is pasted.
    first_query = True

    while True:
        try:
            user_input = _read_query_sequence()

            if user_input is None:
                break
            if not user_input:
                continue

            if args.raw:
                core_target = user_input.upper()
                print(
                    f"\n  [+] Using full sequence: {len(core_target)}aa\n",
                    file=sys.stderr,
                )
            else:
                print("\n  [*] Calculating structural core...", file=sys.stderr)
                core_target = calculate_mature_core(user_input.upper())

                if not core_target:
                    print(
                        "  [!] Error: Mature core calculation returned an empty sequence.\n"
                        "      This can happen if the signal peptide spans the entire protein,\n"
                        "      or if 'GG' appears at the very end of the sequence.\n"
                        "      Try a longer sequence, check the input, or use --raw.",
                        file=sys.stderr,
                    )
                    continue

                print(f"  [+] Core Extracted : {core_target}", file=sys.stderr)
                print(
                    f"  [i] Core Length    : {len(core_target)} amino acids\n",
                    file=sys.stderr,
                )

            if args.output:
                file_mode = "w" if first_query else "a"
                with open(args.output, file_mode, encoding="utf-8-sig") as tsv:
                    if first_query:
                        # Write header only on the first query of this session
                        tsv.write(
                            "Query_Core\tGenome_File\tLocus_Tag\tProduct\tStatus\n"
                        )
                        first_query = False
                    total_input_hits = _scan_and_report(
                        stream_reference_files(args.input), core_target, tsv=tsv
                    )
            else:
                total_input_hits = _scan_and_report(
                    stream_reference_files(args.input), core_target, tsv=None
                )

            if args.output:
                print(
                    f"\n  [=] BATCH SUMMARY: {total_input_hits} matches found. "
                    f"Matrix saved to {args.output.name}.\n",
                    file=sys.stderr,
                )
            print("-" * 60, file=sys.stderr)

        except KeyboardInterrupt:
            print("\n[!] Force quitting. Goodbye!", file=sys.stderr)
            break
        except ValueError as e:
            print(f"\n  [X] Error: {e}\n", file=sys.stderr)


if __name__ == "__main__":
    main()
