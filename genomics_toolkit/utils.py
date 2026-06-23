#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

"""Bioinformatics Standard Utilities

Note:
    v1.3.0: Added ``extract_upstream_sequence_with_length()``. While
    working on comparative_kmer_analyzer.py I noticed
    ``extract_upstream_sequence()`` has always discarded
    ``actual_upstream`` (the real available length, shorter than
    requested only when truncated by a contig boundary) to preserve its
    4-tuple return contract for existing callers (gbk_promoter_finder.py
    also unpacks exactly 4 values). That meant a silently-truncated
    upstream window — e.g. a gene 40bp from its contig's start returning
    only 40bp when 150bp was requested — was undetectable by any caller.
    The new function exposes the 5th value; ``extract_upstream_sequence()``
    now delegates to it internally and is otherwise unchanged.

    v1.4.0: Two additions, both arising from a closer look at
    gbk_promoter_finder.py's downstream consumers.

    (1) Added ``IUPAC_TO_REGEX`` and ``translate_iupac_to_regex()``. Raw
    IUPAC ambiguity codes (W, R, Y, S, K, M, B, D, H, V, N) passed
    directly to Python's ``re.compile()`` are matched as literal
    characters, not ambiguity classes — searching "TATAWAW" finds zero
    matches even against a sequence containing the valid instance
    "TATAAAA". Both ``gbk_promoter_finder.py`` and ``regulon_scanner.py``
    advertise "IUPAC/regex motif" support and had this exact bug; I put
    the translator here so both (and any future motif-scanning script)
    stay consistent.

    (2) Fixed ``extract_upstream_sequence_with_length()``'s CDS-only
    restriction. Previously any locus_tag that exists only on a non-CDS
    feature (mRNA, tRNA, rRNA, ncRNA — i.e. any non-coding RNA gene,
    which by definition has no CDS) raised "Locus tag not found" even
    though the file plainly contains it. Now falls back through
    CDS -> mRNA -> tRNA -> rRNA -> ncRNA, trying each only if the
    previous type isn't present for this locus_tag in the record. CDS is
    tried first and used whenever present, so prokaryote file behavior
    (the only case before this fix) is completely unchanged.

    v1.4.1: Fixed a second, independent gap in
    ``extract_upstream_sequence_with_length()``, in the same root-cause
    class I'd already fixed in find_gbk_features.py v1.7.5: a eukaryotic
    gene can have MULTIPLE features of the SAME type (typically CDS)
    sharing one locus_tag — alternate transcript isoforms with different
    first coding exons, legitimate per NCBI's own submission conventions
    (a locus_tag is shared by every CDS/mRNA/exon belonging to one
    gene). The v1.4.0 fix correctly resolved which FEATURE TYPE to
    prefer (CDS before mRNA before tRNA...) via
    ``candidates.setdefault(feature.type, feature)``, but ``setdefault``
    silently keeps only the FIRST feature seen for a given type — so if
    two CDS features carried the same locus_tag at different
    coordinates, the second was discarded with no warning and no record
    of the collision.

    UNLIKE the find_gbk_features.py fix: that script was solving a
    display/context problem, where taking the full genomic envelope
    (min start, max end across all isoforms) is the correct answer — it
    just needs the gene's whole footprint. This function solves a
    PROMOTER-extraction problem, where isoforms can have genuinely
    different, biologically real transcription start sites. Silently
    unioning their coordinates would manufacture a single answer that
    might be biologically wrong for either isoform, not just imprecise.
    There is no universally "correct" choice here without knowing which
    transcript the caller actually wants.

    So this fix does NOT pick a "correct" isoform — it detects the
    collision, prints a loud warning naming every discarded coordinate
    set, and documents (here and in every caller's --help text) that the
    current behavior is "first feature encountered in file order" and
    that there is no transcript-ID disambiguation option yet. This keeps
    prokaryote behavior (where no locus_tag has more than one CDS) fully
    silent and unchanged, while making the eukaryotic multi-isoform case
    loud instead of silently wrong.

    v1.4.2: Added ``looks_eukaryotic()``, after noticing
    gbk_promoter_finder.py and comparative_kmer_analyzer.py both anchor
    their upstream window, like regulon_scanner.py, on whatever feature
    this module's CDS-first resolution picks (CDS start, i.e. the
    translation start / ATG) — never the transcription start site
    (TSS). regulon_scanner.py already documents itself as
    PROKARYOTE-ONLY for exactly this reason and emits a one-time runtime
    warning when it detects mRNA features mid-scan (a strong eukaryotic
    signal — Prokka/Bakta prokaryote output never emits mRNA features).
    The other two scripts shared the identical underlying anchor problem
    but neither documented nor warned about it; gbk_promoter_finder.py's
    own --upstream help text even read "For eukaryotes, consider
    --upstream 2000 or higher," implying a bigger window alone fixes it,
    which it does not. Rather than duplicate the mRNA-detection scan in
    both files, I extracted it here once and now use it in both (see
    their own changelogs for what changed on their end).
    regulon_scanner.py's own inline check is left as-is: it's folded
    into a scan it already performs for other reasons, so calling this
    separately there would add a redundant full-genome pass for no
    benefit.

    v1.4.3: Added ".mpfa" to ``stream_reference_files()``'s recognized
    extensions. Both protein_presence_scanner.py and
    exact_match_homolog_finder.py document ".mpfa" as a supported
    protein-FASTA reference format and both rely on this function for
    directory-mode scanning (``-i references_dir/``) — but neither the
    single-file check nor the directory glob loop here ever included
    ".mpfa", so a reference file with that extension was silently never
    even handed to either script's own per-file logic. Purely additive:
    every previously-recognized extension is unchanged, so no existing
    caller's behavior is affected for any file that isn't ".mpfa".

    v1.4.4: Added ``disambiguate_isoform_id()``. The same disambiguated-
    identifier logic (``{locus_tag}#{per-isoform-id}``, or
    ``{locus_tag}#isoform{N}`` as a counter-based fallback) had ended up
    independently duplicated inline, twice each, in
    pairwise_homolog_finder.py's ``--keep-all-isoforms`` and
    universal_promoter_extractor.py's ``--all-isoforms`` — four copies
    of essentially the same logic across two files, differing only in
    which qualifier (``protein_id`` vs ``transcript_id``) each happened
    to look up. All four call sites now use this one shared function
    instead. No behavior change at any of the four call sites — I
    re-ran each script's existing isoform-mode tests before and after
    the swap to be sure.
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.4.4"

import argparse
import contextlib
import sys
from pathlib import Path
from collections.abc import Iterator
from typing import IO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO

# fmt: off
DNA_NUC_LIST = ["A", "C", "G", "T"]

DNA_CODON_TABLE: dict[str, str] = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': 'Stop', 'TAG': 'Stop',
    'TGC': 'C', 'TGT': 'C', 'TGA': 'Stop', 'TGG': 'W',
}

RNA_CODON_TABLE: dict[str, str] = {
    "UUU": "F", "CUU": "L", "AUU": "I", "GUU": "V",
    "UUC": "F", "CUC": "L", "AUC": "I", "GUC": "V",
    "UUA": "L", "CUA": "L", "AUA": "I", "GUA": "V",
    "UUG": "L", "CUG": "L", "AUG": "M", "GUG": "V",
    "UCU": "S", "CCU": "P", "ACU": "T", "GCU": "A",
    "UCC": "S", "CCC": "P", "ACC": "T", "GCC": "A",
    "UCA": "S", "CCA": "P", "ACA": "T", "GCA": "A",
    "UCG": "S", "CCG": "P", "ACG": "T", "GCG": "A",
    "UAU": "Y", "CAU": "H", "AAU": "N", "GAU": "D",
    "UAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
    "UAA": "Stop", "CAA": "Q", "AAA": "K", "GAA": "E",
    "UAG": "Stop", "CAG": "Q", "AAG": "K", "GAG": "E",
    "UGU": "C", "CGU": "R", "AGU": "S", "GGU": "G",
    "UGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
    "UGA": "Stop", "CGA": "R", "AGA": "R", "GGA": "G",
    "UGG": "W", "CGG": "R", "AGG": "R", "GGG": "G"
}

MONOISOTOPIC_MASS_TABLE: dict[str, float] = {
    'A': 71.03711, 'C': 103.00919, 'D': 115.02694,
    'E': 129.04259, 'F': 147.06841, 'G': 57.02146,
    'H': 137.05891, 'I': 113.08406, 'K': 128.09496,
    'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
    'P': 97.05276, 'Q': 128.05858, 'R': 156.10111,
    'S': 87.03203, 'T': 101.04768, 'V': 99.06841,
    'W': 186.07931, 'Y': 163.06333
}

_REVCOMP_SRC = "ACGTNRYSWKMBDHVacgtnryswkmbdhv"
_REVCOMP_DST = "TGCANYRSWMKVHDBtgcanyrswmkvhdb"
_REVCOMP_TABLE = str.maketrans(_REVCOMP_SRC, _REVCOMP_DST)
# fmt: on


def revcomp(seq: str) -> str:
    """Returns the reverse complement of a DNA string."""
    return seq.translate(_REVCOMP_TABLE)[::-1]


# fmt: off
IUPAC_TO_REGEX: dict[str, str] = {
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ATCG]",
}
# fmt: on


def translate_iupac_to_regex(motif: str) -> str:
    """Translates IUPAC nucleotide ambiguity codes into regex character classes.

    Python's ``re`` module has no built-in concept of IUPAC ambiguity codes —
    compiling a raw pattern like ``"TATAWAW"`` searches for the literal
    letter "W", which never appears in a DNA sequence (only A/C/G/T do), so
    every motif containing an ambiguity code silently matches nothing:
    ``re.compile("(?=(TATAWAW))")`` against a sequence containing the
    valid TATA-box instance ``TATAAAA`` returns zero matches.

    Translation happens character-by-character, so genuine regex syntax the
    caller has written by hand — brackets, quantifiers, groups, anchors —
    passes through completely unchanged; only the 11 single-letter IUPAC
    ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N) are rewritten into
    their equivalent character class. This lets callers freely mix shorthand
    IUPAC codes with hand-written regex in the same motif string, e.g.
    ``"TATA[AT]W"`` translates to ``"TATA[AT][AT]"`` — the explicit class is
    left alone and only the bare ``W`` is expanded.

    This function does not itself compile or anchor the pattern; callers
    are responsible for ``re.compile()`` and any flags (e.g. IGNORECASE).

    Args:
        motif: A motif string that may mix literal bases (A/C/G/T), IUPAC
               ambiguity codes, and/or raw regex syntax. Case-insensitive —
               internally uppercased before translation.

    Returns:
        The motif with every IUPAC ambiguity code expanded to a regex
        character class. Characters not in the ambiguity-code table
        (plain bases, regex metacharacters, digits, etc.) are returned
        unchanged.
    """
    return "".join(IUPAC_TO_REGEX.get(ch, ch) for ch in motif.upper())


def resolve_upstream_window(
    seq_len: int,
    start: int,
    end: int,
    strand: int,
    upstream_bp: int,
) -> tuple[int, int, int]:
    """Resolves strand-aware upstream-window slice boundaries (arithmetic only)."""
    if strand == 1:
        slice_start = max(0, start - upstream_bp)
        actual_upstream = start - slice_start
        return slice_start, start, actual_upstream
    else:
        slice_end = min(seq_len, end + upstream_bp)
        actual_upstream = slice_end - end
        return end, slice_end, actual_upstream


def extract_upstream_window(
    record,
    start: int,
    end: int,
    strand: int,
    upstream_bp: int,
) -> tuple[str, int, int, int]:
    """Extracts the strand-corrected upstream DNA sequence for a feature."""
    slice_start, slice_end, actual_upstream = resolve_upstream_window(
        len(record.seq), start, end, strand, upstream_bp
    )
    if strand == 1:
        upstream_seq = str(record.seq[slice_start:slice_end])
    else:
        upstream_seq = str(record.seq[slice_start:slice_end].reverse_complement())
    return upstream_seq, actual_upstream, slice_start, slice_end


# Priority order for resolving a locus_tag to a feature when more than one
# feature type in a record could carry it. CDS first preserves prokaryote
# behavior exactly (every prokaryote file has CDS); RNA types are only
# consulted when no CDS with this locus_tag exists in the record at all.
# See extract_upstream_sequence_with_length()'s "FEATURE TYPE FALLBACK" note.
_UPSTREAM_FEATURE_PRIORITY = ("CDS", "mRNA", "tRNA", "rRNA", "ncRNA")


def extract_upstream_sequence_with_length(
    gbk_path: Path, locus_tag: str, upstream_bp: int
) -> tuple[str, int, int, int, int]:
    """Extracts the upstream promoter region of a gene, with truncation info.

    Identical to ``extract_upstream_sequence()`` except it additionally
    returns ``actual_upstream`` — the real number of upstream bases
    available, which is less than ``upstream_bp`` only when the feature
    sits within ``upstream_bp`` of its contig's edge.

    ``extract_upstream_sequence()`` discards this value for backward
    compatibility (its callers unpack a fixed 4-tuple); use this function
    instead when truncation matters to your analysis — e.g. comparing two
    genes' upstream windows, where one being silently shorter than
    requested would otherwise go unnoticed and skew any length-normalized
    statistic computed from it.

    FEATURE TYPE FALLBACK (v1.4.0):
        Previously this function only matched ``CDS`` features, so any
        locus tag that exists only on an ``mRNA``, ``tRNA``, ``rRNA``, or
        ``ncRNA`` feature — e.g. a non-coding RNA gene, which by definition
        has no CDS at all — raised ``ValueError: Locus tag not found``,
        even though the file plainly contains that locus tag — I ran into
        this directly against ``gbk_promoter_finder.py``, the primary caller.

        A single record's CDS/mRNA/tRNA/etc. features for one gene
        typically share the same locus_tag but have DIFFERENT coordinates
        (e.g. an mRNA spans 5' UTR + exons + 3' UTR; its CDS spans only the
        coding portion) — so which feature type is used changes which
        coordinate gets anchored as "upstream of this gene" starts from.
        To keep prokaryote behavior byte-for-byte unchanged, CDS is tried
        first and used if present at all (every prokaryote file has it, so
        no prokaryote caller is affected by this change in any way); the
        RNA feature types below are tried, in order, only when no CDS with
        this locus_tag exists in the record at all.

        ``gene`` is deliberately NOT in this fallback list: its span
        typically includes introns/UTRs and would be a different "upstream
        anchor" than any of the types below, with no clear caller need for
        it today — better to fail loudly than silently extract upstream of
        a likely-wrong coordinate.

    MULTI-ISOFORM CAVEAT (v1.4.1):
        The fallback above resolves which feature TYPE to use. It does NOT
        resolve which feature to use when MULTIPLE features of that same
        type share this locus_tag at different coordinates — e.g. two CDS
        features for two transcript isoforms of one gene, each with a
        different first coding exon. This is legitimate GenBank structure
        per NCBI's own submission conventions (a locus_tag is shared by
        every CDS/mRNA/exon belonging to one gene), and it has no single
        "correct" resolution here: each isoform can have a genuinely
        different, biologically real transcription start site, so unioning
        their coordinates (as find_gbk_features.py's full-envelope fix does
        for its own, different, display-only use case) would manufacture a
        promoter region that may be wrong for either isoform — not just
        imprecise.

        Current behavior: the FIRST matching feature encountered in file
        order is used. If more than one feature of the chosen type shares
        this locus_tag at genuinely different coordinates, a warning is
        printed to stderr naming every discarded coordinate set, so the
        collision is visible instead of silent. There is currently no way
        to request a specific transcript/isoform by ID — if you need a
        specific isoform's promoter, locate its coordinates directly
        (e.g. via its /protein_id or /transcript_id) and extract by
        genomic coordinate instead of by locus_tag. Prokaryote behavior is
        unaffected: no locus_tag in a prokaryote file maps to more than one
        CDS, so the warning path is never reached on prokaryote input.

    Args:
        gbk_path: Path to the GenBank file.
        locus_tag: The unique locus tag of the target gene.
        upstream_bp: Number of base pairs to extract upstream of the start codon.

    Returns:
        A tuple containing: (sequence, start, end, strand, actual_upstream).

    Raises:
        ValueError: If the locus tag is not found or the file cannot be parsed.
    """
    try:
        for record in SeqIO.parse(gbk_path, "genbank"):
            # Collect every candidate feature matching this locus_tag, keyed
            # by feature type, so a deterministic priority order can be
            # applied rather than just using whichever type happens to
            # appear first in the file's feature table. Unlike v1.4.0, this
            # keeps ALL matching features per type (not just the first) so
            # same-type collisions (isoforms) can be detected below instead
            # of silently discarded.
            candidates: dict[str, list] = {}
            for feature in record.features:
                if feature.type not in _UPSTREAM_FEATURE_PRIORITY:
                    continue
                if locus_tag in feature.qualifiers.get("locus_tag", []):
                    candidates.setdefault(feature.type, []).append(feature)

            if not candidates:
                continue

            feature = None
            chosen_type = None
            feature_list: list = []
            for ftype in _UPSTREAM_FEATURE_PRIORITY:
                if ftype in candidates:
                    chosen_type = ftype
                    feature_list = candidates[ftype]
                    feature = feature_list[0]
                    break

            # MULTI-ISOFORM CAVEAT (v1.4.1): warn if the chosen type has more
            # than one feature for this locus_tag at genuinely different
            # coordinates. Identical coordinates (e.g. a duplicate annotation)
            # are not a collision and stay silent — only distinct coordinate
            # sets indicate isoforms with potentially different TSS.
            if len(feature_list) > 1:
                coord_set = {
                    (int(f.location.start), int(f.location.end), f.location.strand)
                    for f in feature_list
                }
                if len(coord_set) > 1:
                    chosen_coord = (
                        int(feature.location.start),
                        int(feature.location.end),
                        feature.location.strand,
                    )
                    other_coords = sorted(coord_set - {chosen_coord})
                    print(
                        f"[!] Warning: locus_tag '{locus_tag}' matches "
                        f"{len(feature_list)} '{chosen_type}' features in "
                        f"{gbk_path.name} at different coordinates — likely "
                        f"alternate transcript isoforms sharing one "
                        f"locus_tag. Using the FIRST one encountered "
                        f"(start={chosen_coord[0]}, end={chosen_coord[1]}, "
                        f"strand={chosen_coord[2]}); other isoform(s) at "
                        f"{other_coords} were ignored. No transcript-ID "
                        f"disambiguation is available yet — see "
                        f"extract_upstream_sequence_with_length()'s "
                        f"MULTI-ISOFORM CAVEAT docstring.",
                        file=sys.stderr,
                    )

            start = int(feature.location.start)
            end = int(feature.location.end)
            strand = feature.location.strand

            seq, actual_upstream, _slice_start, _slice_end = extract_upstream_window(
                record, start, end, strand, upstream_bp
            )
            return seq, start, end, strand, actual_upstream

    except (FileNotFoundError, OSError, UnicodeDecodeError, ValueError) as e:
        raise ValueError(f"Failed to parse GenBank file {gbk_path.name}: {e}") from e

    raise ValueError(f"Locus tag '{locus_tag}' not found in {gbk_path.name}.")


def extract_upstream_sequence(
    gbk_path: Path, locus_tag: str, upstream_bp: int
) -> tuple[str, int, int, int]:
    """Extracts the upstream promoter region of a specific gene.

    Delegates to ``extract_upstream_sequence_with_length()`` and drops the
    ``actual_upstream`` element to preserve this function's original
    4-tuple contract for existing callers. Callers that need to detect
    contig-boundary truncation should use
    ``extract_upstream_sequence_with_length()`` instead.

    See ``extract_upstream_sequence_with_length()``'s MULTI-ISOFORM CAVEAT
    (v1.4.1) — this function inherits that same behavior and warning.

    Args:
        gbk_path: Path to the GenBank file.
        locus_tag: The unique locus tag of the target gene.
        upstream_bp: Number of base pairs to extract upstream of the start codon.

    Returns:
        A tuple containing: (sequence, start, end, strand).

    Raises:
        ValueError: If the locus tag is not found or the file cannot be parsed.
    """
    seq, start, end, strand, _actual_upstream = extract_upstream_sequence_with_length(
        gbk_path, locus_tag, upstream_bp
    )
    return seq, start, end, strand


def disambiguate_isoform_id(
    locus_tag: str,
    feature,
    isoform_counters: dict[str, int],
    id_qualifier: str = "protein_id",
) -> str:
    """Builds a disambiguated identifier for one isoform of a multi-isoform locus.

    The same logic was independently duplicated, inline, in two places
    each across two different scripts (``pairwise_homolog_finder.py``'s
    ``--keep-all-isoforms`` and ``universal_promoter_extractor.py``'s
    ``--all-isoforms`` — 4 copies of essentially the same logic total)
    before I consolidated it here. Both flags exist for the same reason:
    a "keep only the longest/5'-most representative isoform" default is
    a reasonable, deterministic choice, but not necessarily the
    biologically dominant one — when a caller wants every isoform
    individually instead, each one needs its own addressable identifier
    rather than colliding on the shared ``locus_tag``.

    Resolution order:
      1. ``feature.qualifiers[id_qualifier]`` — typically unique per
         isoform (``/protein_id`` on a CDS feature;
         ``/transcript_id`` on an mRNA feature — callers pass whichever
         qualifier their feature type actually carries). Returns
         ``f"{locus_tag}#{isoform_id}"``.
      2. If that qualifier is absent, falls back to a per-locus running
         counter, mutating ``isoform_counters`` in place (the caller owns
         this dict and must reuse the SAME dict across every call within
         one extraction pass, so the counter persists per locus rather
         than resetting). Returns ``f"{locus_tag}#isoform{N}"``.

    Args:
        locus_tag:        The resolved base identifier for this isoform's
                           gene (shared across all its isoforms).
        feature:           A Biopython SeqFeature for this specific isoform.
        isoform_counters: A ``{locus_tag: count}`` dict the caller owns and
                           reuses across every call in one extraction pass
                           — mutated in place to track how many
                           qualifier-less isoforms of this locus have
                           already been assigned a counter-based ID.
        id_qualifier:      Which qualifier to look up for a per-isoform ID.
                           ``"protein_id"`` for CDS features, ``"transcript_id"``
                           for mRNA features. Default: ``"protein_id"``.

    Returns:
        A disambiguated string identifier, never colliding with another
        isoform of the same locus within one extraction pass.
    """
    isoform_id = feature.qualifiers.get(id_qualifier, [""])[0]
    if isoform_id:
        return f"{locus_tag}#{isoform_id}"
    isoform_counters[locus_tag] = isoform_counters.get(locus_tag, 0) + 1
    return f"{locus_tag}#isoform{isoform_counters[locus_tag]}"


def looks_eukaryotic(gbk_path: Path) -> bool:
    """Cheap heuristic check for eukaryotic content: any mRNA feature present.

    Used by CDS-anchored upstream-extraction tools (gbk_promoter_finder.py,
    comparative_kmer_analyzer.py) to emit a one-time warning that their
    upstream window is anchored on CDS start, not the TSS, before the user
    silently gets a "promoter" that is actually 5' UTR/intron sequence
    anchored at the wrong coordinate. See those scripts' PROKARYOTE-ONLY
    ANCHOR docstring sections, and regulon_scanner.py, which performs the
    equivalent check inline during a single-pass genome scan it already
    does for other reasons.

    This is a heuristic, not a guarantee: presence of an mRNA feature is a
    strong signal of eukaryotic annotation (Prokka/Bakta prokaryote output
    does not emit mRNA features), but its ABSENCE does not prove a genome
    is prokaryotic — it only means this particular signal wasn't found.
    Stops at the first mRNA feature found rather than scanning the entire
    file, since one occurrence is sufficient to warrant the warning.

    Args:
        gbk_path: Path to the GenBank file.

    Returns:
        True if any mRNA feature is found in any record; False otherwise,
        including on any parse error — this check is advisory-only, so a
        failed heuristic should never block real extraction work, which
        will raise its own clear error separately if the file is actually
        unreadable.
    """
    try:
        for record in SeqIO.parse(gbk_path, "genbank"):
            for feature in record.features:
                if feature.type == "mRNA":
                    return True
    except (OSError, UnicodeDecodeError, ValueError):
        return False
    return False


def stream_reference_files(target_path: Path) -> Iterator[Path]:
    """Yields valid GenBank or FASTA files from a file or directory."""
    valid_exts = (".gbk", ".gbff", ".fasta", ".fa", ".faa", ".mpfa")

    if target_path.is_file():
        if target_path.suffix.lower() in valid_exts:
            yield target_path
        else:
            print(
                f"[!] Warning: {target_path.name} is not a valid GenBank/FASTA extension.",
                file=sys.stderr,
            )

    elif target_path.is_dir():
        for ext in ("*.gbk", "*.gbff", "*.fasta", "*.fa", "*.faa", "*.mpfa"):
            yield from target_path.rglob(ext)

    else:
        raise ValueError(f"Provided path does not exist: {target_path}")


def calculate_mature_core(full_protein: str) -> str:
    """Calculates the mature peptide core based on double-glycine cleavage sites."""
    if "GG" not in full_protein:
        return full_protein

    parts = full_protein.split("GG", 1)
    mature_peptide = parts[1]

    MIN_LENGTH = 25
    if len(mature_peptide) < MIN_LENGTH:
        return mature_peptide

    for i in range(MIN_LENGTH, len(mature_peptide)):
        current_residue = mature_peptide[i]

        if current_residue == "P":
            return mature_peptide[: i + 1]

        if i + 5 <= len(mature_peptide):
            window = mature_peptide[i : i + 5]
            analyzer = ProteinAnalysis(window)
            avg_hydro = analyzer.gravy()

            if avg_hydro < -0.5:
                return mature_peptide[:i]

    return mature_peptide


@contextlib.contextmanager
def smart_open(filename: Path | None) -> Iterator[IO]:
    """Routes output to a file or standard output."""
    if filename:
        handle = open(filename, "w", encoding="utf-8")
        try:
            yield handle
        finally:
            handle.close()
    else:
        yield sys.stdout


def wrap_fasta(sequence: str, width: int = 60) -> str:
    """Wraps a sequence string into multiple lines for FASTA formatting."""
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def base_parser(
    description_text: str, include_input: bool = True, include_output: bool = True
) -> argparse.ArgumentParser:
    """Creates a standard CLI argument parser for pipeline scripts."""
    parser = argparse.ArgumentParser(description=description_text)

    if include_input:
        parser.add_argument(
            "-i", "--input", type=Path, required=True, help="Path to input file."
        )

    if include_output:
        parser.add_argument(
            "-o",
            "--output",
            type=Path,
            required=False,
            default=None,
            help="Path to save the output file (optional — prints to terminal if omitted).",
        )

    return parser


def parse_fasta(fasta_string: str) -> dict[str, str]:
    """Parses a raw multi-FASTA text string into an ID-to-sequence dictionary."""
    fasta_dict: dict[str, str] = {}
    fasta_id: str = ""
    seq_buffer: list[str] = []

    for line in fasta_string.splitlines():
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if fasta_id:
                fasta_dict[fasta_id] = "".join(seq_buffer)

            header_parts = line[1:].strip().split(None, 1)
            if not header_parts:
                raise ValueError(
                    "Malformed FASTA: Empty identifier header ('>') found."
                )

            fasta_id = header_parts[0]
            seq_buffer = []
        else:
            if not fasta_id:
                raise ValueError(
                    "Malformed FASTA: Sequence data found before an identifier header."
                )
            seq_buffer.append(line)

    if fasta_id:
        fasta_dict[fasta_id] = "".join(seq_buffer)

    return fasta_dict


def lazy_parse_fasta(file_path: str | Path) -> Iterator[tuple[str, str]]:
    """Reads a FASTA file line-by-line to minimize memory usage."""
    fasta_id = ""
    seq_buffer = []

    with open(file_path, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if fasta_id:
                    yield fasta_id, "".join(seq_buffer)

                header_parts = line[1:].strip().split(None, 1)
                if not header_parts:
                    raise ValueError(
                        "Malformed FASTA: Empty identifier header ('>') found."
                    )

                fasta_id = header_parts[0]
                seq_buffer = []
            else:
                if not fasta_id:
                    raise ValueError(
                        "Malformed FASTA: Sequence data found before an identifier header."
                    )
                seq_buffer.append(line)

        if fasta_id:
            yield fasta_id, "".join(seq_buffer)
