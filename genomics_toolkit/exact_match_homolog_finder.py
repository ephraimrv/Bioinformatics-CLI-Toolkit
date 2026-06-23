#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

"""Exact Match Homolog Finder — exact-substring peptide search across reference genomes.

Calculates mature peptide cores and finds candidate homologs across reference
genomes using EXACT SUBSTRING MATCHING (100%% identity required).

This script accepts pre-peptide sequences and calculates their mature,
membrane-inserting cores based on biochemical cleavage sites. It then
uses these core sequences to scan target reference genomes, identifying
and extracting the full matching protein sequences. If multiple targets map
to the same physical locus, it aggregates them into a single FASTA header
to prevent duplicate sequence outputs.

Important:
    This tool uses EXACT substring matching. A single amino acid substitution
    in the core region will cause a miss. For divergent homolog detection at
    lower identity (e.g., ≥35%%), use pairwise_homolog_finder.py
    (Smith-Waterman alignment with BLOSUM62 scoring) instead.

    You must choose explicitly between --raw (search full sequences as
    provided) and --mature (apply bacterial Sec/Tat/RiPP double-glycine
    leader-cleavage trimming first) — there is no default. See the v1.3.0
    changelog below for why.

    Reference files must be protein sequences (GenBank with /translation
    qualifiers, or protein FASTA: .faa/.mpfa). Nucleotide FASTA files
    (.fasta/.fa/.fna) are skipped with a warning — see "Reference format
    handling" below.

Terminology note:
    "Exact match" here means byte-identical sequence, not validated
    orthology. An exact match can equally be a true ortholog, a recent
    paralog, or simply a short, highly conserved motif shared by unrelated
    proteins. Treat hits as candidate homologs requiring further
    confirmation (e.g. reciprocal-best-hit or phylogenetic analysis) before
    describing them as orthologs in a manuscript or other formal report.

Note:
    Associated with ongoing, unpublished research (manuscript in
    preparation). Correct attribution is requested when used in derivative
    works.

    v1.2.1: Corrected two stale references found while auditing
    conserved_annotation_scanner.py for an identical issue. (1) The
    "Important" section pointed to "gbk_ortholog_finder.py" for divergent
    homolog detection — a filename that has never existed in this
    project. Confirmed pairwise_homolog_finder.py is the renamed target
    (its own docstring independently describes itself as the
    Smith-Waterman/BLOSUM62 sequence-based clustering tool). (2) This
    file's own Example section instructed users to run
    "exact_match_ortholog_finder.py" — the deprecated predecessor this
    file superseded — instead of its own actual filename,
    exact_match_homolog_finder.py. No behavior change; documentation only.

    v1.3.0: BREAKING CHANGE — found while reasoning through whether this
    script's sibling, pairwise_homolog_finder.py, implies this one is
    eukaryote-incompatible (it doesn't, structurally, but the two scripts'
    defaults diverged in a way that mattered). pairwise_homolog_finder.py
    defaults to NOT applying calculate_mature_core() (its --mature flag
    is opt-in); this script defaulted to APPLYING it unless --raw was
    explicitly passed — the riskier direction, since calculate_mature_core()
    models bacterial Sec/Tat/RiPP double-glycine cleavage specifically.
    The danger is narrower than "every eukaryotic target breaks": most
    eukaryotic proteins contain no "GG" dipeptide acting as a meaningful
    signal, so the function is a no-op for them even without --raw. The
    real failure mode is a non-bacteriocin protein that happens to
    contain a COINCIDENTAL "GG" substring (not rare in a 200+ residue
    protein) — in that case the function silently chops the sequence at
    that coincidental site and runs bacterial-specific hydrophobicity-
    window logic on the remainder, producing a biologically meaningless
    fragment with no indication anything went wrong. There is no
    reliable way to detect "is this protein bacterial" from a bare FASTA
    target (no organism metadata), so rather than pick a new default in
    either direction (which only swaps which silent-failure mode is
    possible), --raw and --mature are now a REQUIRED mutually exclusive
    pair — running this script without explicitly choosing one now
    raises an argparse error instead of silently assuming either
    behavior. Any existing call site that previously relied on the
    implicit default (running with neither flag) must add one
    explicitly; there is no backward-compatible default to fall back on
    by design.

Example:
    Bacteriocin/RiPP core search (explicit, no longer the default)::

        python3 exact_match_homolog_finder.py \\
            -t target.faa -r references/ --mature -o extracted_homologs.faa

    Non-bacteriocin proteins (housekeeping genes, TFs, kinases, eukaryotic
    proteins in general)::

        python3 exact_match_homolog_finder.py \\
            -t target.faa -r references/ --raw -o extracted_homologs.faa
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.3.0"

import sys
import argparse
from pathlib import Path
from typing import Iterator

try:
    from Bio import SeqIO
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )
from utils import stream_reference_files, calculate_mature_core, smart_open


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Searches for exact-match peptide cores against reference genomes."
    )
    parser.add_argument(
        "-t",
        "--targets",
        type=Path,
        required=True,
        help="FASTA file containing full target sequences. The script will auto-trim them.",
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=Path,
        default=Path("."),
        help="Input GenBank/FASTA file OR directory to scan (Default: current directory)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Optional: Output FASTA file path (Default: prints to terminal)",
    )
    core_group = parser.add_mutually_exclusive_group(required=True)
    core_group.add_argument(
        "--raw",
        action="store_true",
        help=(
            "Skip calculate_mature_core() and use the full sequences exactly "
            "as provided. Required for non-bacteriocin targets — "
            "transcription factors, housekeeping genes, kinases, and "
            "eukaryotic proteins in general — since none of these use the "
            "bacterial double-glycine leader-cleavage logic this trimmer "
            "models."
        ),
    )
    core_group.add_argument(
        "--mature",
        action="store_true",
        help=(
            "Apply calculate_mature_core() before searching: trims each "
            "target to its predicted mature core using bacterial "
            "Sec/Tat/RiPP double-glycine cleavage rules. Use this ONLY for "
            "bacteriocin/RiPP-type targets. This does NOT model eukaryotic "
            "secretory-pathway signal peptides (ER/Golgi-targeted), which "
            "use entirely different cleavage motifs — there is no reliable "
            "way to auto-detect this from a bare sequence (it carries no "
            "organism metadata), so the choice between --raw and --mature "
            "is mandatory rather than defaulted, to avoid silently applying "
            "bacterial-specific trimming to a non-bacterial target."
        ),
    )
    return parser.parse_args()


def load_target_peptides(fasta_path: Path, use_raw: bool = False) -> dict[str, str]:
    """Loads targets and optionally runs the biochemical trimmer on them.

    Args:
        fasta_path: Path to the FASTA file containing target sequences.
        use_raw:    If True, skips calculate_mature_core() and uses the
                    full sequence as the search probe. Use for non-bacteriocin
                    targets that don't have double-glycine cleavage sites.
    """
    if not fasta_path.is_file():
        raise ValueError(f"Target FASTA file not found: {fasta_path}")

    targets = {}
    skipped = 0
    print(f"[*] Processing Target File: {fasta_path.name}", file=sys.stderr)
    if use_raw:
        print(
            "[*] Mode: --raw (full sequences used, no core trimming)", file=sys.stderr
        )
    else:
        print(
            "[*] Mode: --mature (bacterial Sec/Tat/RiPP double-glycine "
            "cleavage rules; NOT valid for eukaryotic secretory-pathway "
            "targets)",
            file=sys.stderr,
        )

    with open(fasta_path, "r", encoding="utf-8") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_str = str(record.seq).strip()
            if not seq_str:
                continue

            if use_raw:
                probe = seq_str
                print(
                    f"     -> {record.id}: Using full {len(probe)}aa sequence as probe",
                    file=sys.stderr,
                )
            else:
                probe = calculate_mature_core(seq_str)

                if not probe:
                    # Empty core means calculate_mature_core found "GG" at the very
                    # end of the sequence (nothing left after the cleavage site).
                    # Rather than silently skipping, warn the user explicitly.
                    print(
                        f"     -> [!] WARNING: {record.id} — core trimming returned an "
                        f"empty sequence. This sequence will be SKIPPED.\n"
                        f"        Tip: If this is not a bacteriocin, re-run with --raw.",
                        file=sys.stderr,
                    )
                    skipped += 1
                    continue

                print(
                    f"     -> {record.id}: Trimmed to {len(probe)}aa probe",
                    file=sys.stderr,
                )
                print(f"        Core: {probe}", file=sys.stderr)

            targets[record.id] = probe

    if skipped > 0:
        print(
            f"\n[!] {skipped} sequence(s) skipped due to empty core. "
            f"Use --raw to include full sequences.",
            file=sys.stderr,
        )

    if not targets:
        raise ValueError(
            f"No valid probe sequences found in target file: {fasta_path.name}"
        )

    return targets


def extract_homologs(
    ref_path: Path, target_peptides: dict[str, str]
) -> Iterator[tuple[str, str, str, str, str, str]]:
    """Dynamically parses GenBank or FASTA format to yield exact-match hits.

    Reference format handling:
        GenBank files always work — protein sequences come from the
        /translation qualifier regardless of the file's nucleotide backbone.
        For bare FASTA files, extension determines how the content is
        treated (same convention as protein_presence_scanner.py): ``.faa``/
        ``.mpfa`` are protein FASTA and are searched directly. ``.fasta``/
        ``.fa``/``.fna`` are treated as nucleotide FASTA and are skipped
        with a warning rather than silently exact-matching a protein probe
        against a DNA string (which always fails with zero hits and no
        explanation otherwise). If you do have protein sequences saved with
        a ``.fasta``/``.fa`` extension, rename to ``.faa`` to search them.

    Yields:
        Tuples of ``(seq_id, locus, target_name, product, full_prot,
        protein_id)``. ``protein_id`` is the RefSeq-style ``/protein_id``
        qualifier when present (distinguishes splice isoforms that share
        the same ``locus_tag``); empty string when absent or not applicable
        (always empty for the FASTA branch, where ``locus`` is already the
        FASTA record's own unique ID).
    """
    try:
        suffix = ref_path.suffix.lower()
        is_protein_fasta = suffix in (".faa", ".mpfa")
        is_nucleotide_fasta = suffix in (".fasta", ".fa", ".fna")
        is_fasta = is_protein_fasta or is_nucleotide_fasta

        if is_nucleotide_fasta:
            # Guard against nucleotide FASTA silently scanned as protein
            # (same bug class already fixed in protein_presence_scanner.py;
            # same extension-based convention applied here for consistency).
            # Skip-and-warn rather than hard-exit, since a directory scan
            # shouldn't abort over one bad file among many valid ones.
            print(
                f"  [!] Skipping {ref_path.name}: nucleotide FASTA "
                f"(.{suffix.lstrip('.')}) cannot be searched with a protein "
                f"probe. Provide a protein FASTA (.faa) or an annotated "
                f"GenBank file instead.",
                file=sys.stderr,
            )
            return

        fmt = "fasta" if is_fasta else "genbank"

        with open(ref_path, "r", encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, fmt):

                # BRANCH A: FASTA processing
                if is_fasta:
                    full_translation = str(record.seq).upper()

                    for target_name, core_peptide in target_peptides.items():
                        if not core_peptide:
                            continue

                        if core_peptide in full_translation:
                            full_header_desc = record.description.replace(
                                record.id, ""
                            ).strip()
                            product = (
                                full_header_desc
                                if full_header_desc
                                else "Unannotated FASTA sequence"
                            )

                            yield (
                                ref_path.stem,
                                record.id,
                                target_name,
                                product,
                                full_translation,
                                "",  # protein_id: not applicable, record.id is unique here
                            )

                # BRANCH B: GenBank processing
                else:
                    for feature in record.features:
                        if feature.type == "CDS":
                            full_translation = feature.qualifiers.get(
                                "translation", [""]
                            )[0]
                            if not full_translation:
                                continue

                            for target_name, core_peptide in target_peptides.items():
                                if not core_peptide:
                                    continue

                                if core_peptide in full_translation:
                                    locus_tag = feature.qualifiers.get(
                                        "locus_tag", ["UNKNOWN"]
                                    )[0]
                                    product = feature.qualifiers.get(
                                        "product", ["Unknown product"]
                                    )[0]
                                    protein_id = feature.qualifiers.get(
                                        "protein_id", [""]
                                    )[0]

                                    yield (
                                        record.id,
                                        locus_tag,
                                        target_name,
                                        product,
                                        full_translation,
                                        protein_id,
                                    )

    except Exception as e:
        raise ValueError(f"Failed to parse {ref_path.name}: {e}") from e


def main() -> None:
    args = get_args()

    try:
        targets = load_target_peptides(args.targets, use_raw=args.raw)

        print(f"\n[*] Scanning reference space: {args.reference}\n", file=sys.stderr)

        total_extracted_entries = 0

        with smart_open(args.output) as out_handle:
            for file_path in stream_reference_files(args.reference):
                print(f"  -> Scanning {file_path.name}...", file=sys.stderr)

                # file_hits aggregates multiple target probes hitting the same
                # exact sequence. Keyed by (locus, full_prot) rather than bare
                # locus: a locus_tag is NOT guaranteed unique to one sequence —
                # eukaryotic splice isoforms routinely share a locus_tag while
                # having different translations. Keying on locus alone would
                # let a second isoform's match silently overwrite/be discarded
                # under the first isoform's sequence (a biologically false
                # result). Keying on the actual sequence guarantees distinct
                # isoforms get distinct entries, while truly identical
                # sequences correctly still collapse into one entry.
                file_hits = {}

                for (
                    seq_id,
                    locus,
                    target_name,
                    product,
                    full_prot,
                    protein_id,
                ) in extract_homologs(file_path, targets):
                    key = (locus, full_prot)
                    if key not in file_hits:
                        file_hits[key] = {
                            "seq_id": seq_id,
                            "product": product,
                            "full_prot": full_prot,
                            "protein_id": protein_id,
                            "mapped_targets": [target_name],
                        }
                    else:
                        file_hits[key]["mapped_targets"].append(target_name)

                # Count how many distinct sequences share each locus, so headers
                # can disambiguate isoforms (e.g. LOCUS_XP_001 vs LOCUS_XP_002)
                # instead of all displaying the same bare locus tag.
                locus_variant_count: dict[str, int] = {}
                for locus, _seq in file_hits:
                    locus_variant_count[locus] = locus_variant_count.get(locus, 0) + 1
                locus_running_index: dict[str, int] = {}

                # Process aggregated hits and write output
                for (locus, _seq), data in file_hits.items():
                    total_extracted_entries += 1

                    # Clean the target names for the FASTA header
                    clean_names = [
                        t.replace(" ", "_").replace("/", "_").replace(",", "")
                        for t in data["mapped_targets"]
                    ]

                    # Truncate gracefully if many probes hit the same locus —
                    # oversized headers break some downstream tools (HMMER, aligners)
                    if len(clean_names) > 3:
                        targets_str = (
                            ",".join(clean_names[:3])
                            + f",_and_{len(clean_names) - 3}_more"
                        )
                    else:
                        targets_str = ",".join(clean_names)

                    # Disambiguate the displayed identifier when multiple
                    # distinct sequences (isoforms) share this locus_tag.
                    if locus_variant_count[locus] > 1:
                        locus_running_index[locus] = (
                            locus_running_index.get(locus, 0) + 1
                        )
                        if data["protein_id"]:
                            display_locus = f"{locus}_{data['protein_id']}"
                        else:
                            display_locus = (
                                f"{locus}_isoform{locus_running_index[locus]}"
                            )
                    else:
                        display_locus = locus

                    fasta_header = (
                        f">{display_locus} | {data['seq_id']} | {data['product']} | "
                        f"[Homolog_of_{targets_str}]"
                    )

                    out_handle.write(f"{fasta_header}\n{data['full_prot']}\n")
                    print(
                        f"      [Hit] {display_locus} contains core(s) from "
                        f"{len(data['mapped_targets'])} target probe(s)! "
                        f"({len(data['full_prot'])} aa)",
                        file=sys.stderr,
                    )

        print("\n" + "=" * 50, file=sys.stderr)
        print(
            f"[*] SUCCESS: {total_extracted_entries} unique sequence(s) extracted.",
            file=sys.stderr,
        )

        if args.output:
            print(f"[*] Output saved to: {args.output.resolve()}", file=sys.stderr)
        print("=" * 50, file=sys.stderr)

    except ValueError as e:
        sys.exit(f"\n[!] Pipeline Error: {e}")
    except KeyboardInterrupt:
        sys.exit("\n[!] Pipeline gracefully interrupted by user.")


if __name__ == "__main__":
    main()
