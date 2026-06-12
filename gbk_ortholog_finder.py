# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

"""GBK Ortholog Finder — pairwise protein ortholog detection from GenBank files.

Extracts protein-coding sequences from a query GBK/GBFF file and compares
them against one or more reference genomes using local pairwise alignment
(Smith-Waterman with BLOSUM62) to identify orthologs by percent identity.

Works on prokaryotic and eukaryotic GenBank files. Accepts a single
reference file or a directory of reference files.

Performance pipeline (each stage cheaper than the last):
    1. Length ratio filter  — skips pairs outside 0.5–2.0 length ratio (O(1)).
    2. k-mer Jaccard filter — skips pairs sharing too few 4-mers to meet the
       identity threshold (O(n+m) Python, no alignment needed).
    3. Score pre-filter     — calls ``aligner.score()`` which runs entirely in
       optimised C with no traceback. Skips pairs scoring ≤ 0 or below a
       length-scaled minimum. This prevents the Biopython traceback combinatorial
       explosion that hangs the process on unrelated sequences.
    4. Full alignment       — ``aligner.align()`` with traceback only for the
       small fraction of pairs that survive stages 1–3.

    All reference protein mature sequences and k-mers are pre-computed once
    per reference file before the inner query loop runs, so these operations
    are never repeated inside the nested comparison.

Flag interaction — ``--mature`` and ``--min-coverage``:
    ``--mature`` trims signal peptides from both query and reference proteins
    before alignment, so only the bioactive mature core is compared.
    ``--min-coverage`` then measures how much of the shorter *trimmed* sequence
    is covered by the alignment. Used together, they require that the entire
    mature core — not just a short conserved motif — is genuinely present in
    the reference. This combination is strongly recommended for bacteriocin
    and small-peptide searches.

Note:
    This script is part of ongoing research and is associated with an upcoming
    publication. Correct attribution is requested when used in derivative works.
    Released under the MIT License. See LICENSE in the repository root.

Example:
    Compare a region GBK against a single reference genome::

        python3 gbk_ortholog_finder.py \\
            -q region001.gbk -r ATCC8293.gbff -o results.tsv

    Bacteriocin screen with signal peptide trimming::

        python3 gbk_ortholog_finder.py \\
            -q region001.gbk -r references/ \\
            --mature --max-length 150 --identity 0.35 --min-coverage 0.65 \\
            -o bacteriocin_hits.tsv

    Full genome vs directory of references::

        python3 gbk_ortholog_finder.py \\
            -q genome.gbff -r references/ --identity 0.40 -o results.tsv
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.4"

import sys
import argparse
import csv
from pathlib import Path
from dataclasses import dataclass

try:
    from Bio import SeqIO
    from Bio.Align import substitution_matrices, PairwiseAligner
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )

from utils import stream_reference_files, calculate_mature_core, smart_open

# ── Module-level aligner and matrix (created ONCE, reused for all comparisons) ─
# Moving PairwiseAligner() out of the comparison function avoids re-instantiating
# it millions of times, which was the primary constant-factor bottleneck.
_BLOSUM62 = substitution_matrices.load("BLOSUM62")

_ALIGNER = PairwiseAligner()
_ALIGNER.mode = "local"  # Smith-Waterman; no forced end-to-end gaps
_ALIGNER.substitution_matrix = _BLOSUM62
_ALIGNER.open_gap_score = -11  # Standard protein local alignment penalties
_ALIGNER.extend_gap_score = -1

_K: int = 4  # k-mer length; 4-mers are more selective than 3-mers for proteins


# ─────────────────────────────────────────────────────────────────────────────
# DATA STRUCTURES
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class Protein:
    """A protein extracted from the query GBK file."""

    locus_tag: str
    product: str
    sequence: str  # Full /translation= sequence
    mature_sequence: (
        str  # After calculate_mature_core() — equals sequence if --mature not used
    )
    source_file: str
    length: int


@dataclass
class _RefProtein:
    """Internal reference protein with all comparison fields pre-computed.

    ``cmp_seq`` and ``kmers`` are derived from the mature sequence when
    ``use_mature=True``, otherwise from the full sequence. They are computed
    ONCE at load time so the inner comparison loop never recomputes them.
    """

    locus_tag: str
    product: str
    sequence: str  # Full original sequence (stored for TSV output)
    cmp_seq: str  # Sequence used for alignment (full or mature)
    kmers: frozenset  # Pre-computed k-mers of cmp_seq
    length: int  # len(cmp_seq)


@dataclass
class OrthoHit:
    """A single ortholog comparison result."""

    query_locus: str
    query_product: str
    query_seq: str
    ref_locus: str
    ref_product: str
    ref_seq: str
    ref_file: str
    identity: float  # 0.0–1.0
    alignment_length: int
    mismatches: int  # Number of non-matching positions in the alignment
    query_length: int
    ref_length: int


# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: ARGUMENT PARSER
# ─────────────────────────────────────────────────────────────────────────────


def get_args() -> argparse.Namespace:
    """Configures and returns the CLI argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Find orthologs by comparing proteins from a query GenBank file "
            "against one or more reference GenBank files. "
            "Uses local pairwise alignment (Smith-Waterman / BLOSUM62)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-q",
        "--query",
        type=Path,
        required=True,
        help="Query GBK or GBFF file. Proteins are extracted from every CDS feature.",
    )
    parser.add_argument(
        "-r",
        "--reference",
        type=Path,
        required=True,
        help=(
            "Reference GBK/GBFF file, OR a folder containing multiple such files. "
            "All .gbk and .gbff files in the folder are scanned automatically."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help=(
            "Save results as a TSV file (includes full protein sequences). "
            "If not used, results are printed as a table in the terminal "
            "(without sequences)."
        ),
    )
    parser.add_argument(
        "--identity",
        type=float,
        default=0.30,
        help=(
            "Minimum percent identity to report a hit, expressed as a decimal "
            "(e.g. 0.30 = 30%%, 0.50 = 50%%). "
            "Lower values find more distant relatives; higher values find closer ones. "
            "Default: 0.30."
        ),
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.50,
        help=(
            "How much of the shorter protein must be covered by the alignment "
            "(0.0–1.0, i.e. 0.50 = 50%%). "
            "This prevents hits where only a small conserved domain (e.g. a signal "
            "peptide or a zinc-finger motif) matches, while the rest of the protein "
            "is unrelated. When used with --mature, coverage is measured against the "
            "trimmed mature core, so 0.70 means '70%% of the active peptide aligns'. "
            "Recommended range: 0.50 (permissive) to 0.80 (strict). Default: 0.50."
        ),
    )
    parser.add_argument(
        "--mature",
        action="store_true",
        default=False,
        help=(
            "Strip signal peptides and pro-sequences before comparing. "
            "Bacteriocins and many secreted peptides have a signal peptide region "
            "at the N-terminus that is cleaved off before the protein becomes active. "
            "This region is often poorly conserved even between close relatives, "
            "which can lower identity scores or cause missed hits. "
            "Enabling --mature trims both the query and reference proteins to their "
            "predicted mature cores before alignment, so you are comparing only the "
            "biologically active portion. Strongly recommended for bacteriocin, "
            "lantibiotic, and small secreted peptide comparisons. "
            "Combine with --min-coverage to ensure the entire mature core aligns, "
            "not just a fragment of it."
        ),
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=None,
        help=(
            "Only include query proteins at or below this length (in amino acids). "
            "Useful to focus a search on small peptides: bacteriocins are typically "
            "20–150 aa, so --max-length 150 screens out larger unrelated proteins "
            "and speeds up the run considerably."
        ),
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=10,
        help="Skip query proteins shorter than this many amino acids. Default: 10.",
    )

    return parser.parse_args()


# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: EXTRACT PROTEINS FROM QUERY GBK
# ─────────────────────────────────────────────────────────────────────────────


def extract_proteins_from_gbk(
    gbk_path: Path,
    apply_mature: bool = False,
    min_length: int = 10,
    max_length: int | None = None,
) -> list[Protein]:
    """Reads a GBK/GBFF file and returns all filtered CDS protein sequences.

    Args:
        gbk_path:     Path to the query GBK/GBFF file.
        apply_mature: If ``True``, applies ``calculate_mature_core()`` to each protein.
        min_length:   Skip proteins shorter than this.
        max_length:   Skip proteins longer than this (``None`` = no limit).

    Returns:
        List of ``Protein`` objects.
    """
    proteins: list[Protein] = []
    print(f"\n[*] Extracting proteins from query: {gbk_path.name}", file=sys.stderr)

    skipped_zero_length = 0

    for record in SeqIO.parse(gbk_path, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue
            translation = feature.qualifiers.get("translation", [""])[0]
            if not translation:
                continue

            locus_tag = feature.qualifiers.get("locus_tag", ["UNKNOWN"])[0]
            product = feature.qualifiers.get("product", ["Unknown product"])[0]
            full_length = len(translation)

            if full_length < min_length:
                continue
            if max_length is not None and full_length > max_length:
                continue

            mature_seq = (
                calculate_mature_core(translation) if apply_mature else translation
            )

            # Skip proteins whose mature core is zero-length
            if apply_mature and (not mature_seq or len(mature_seq) == 0):
                print(
                    f"   [!] Skipping {locus_tag}: mature core is zero-length.",
                    file=sys.stderr,
                )
                skipped_zero_length += 1
                continue

            proteins.append(
                Protein(
                    locus_tag=locus_tag,
                    product=product,
                    sequence=translation,
                    mature_sequence=mature_seq,
                    source_file=gbk_path.name,
                    length=full_length,
                )
            )

            if apply_mature and mature_seq != translation:
                print(
                    f"   {locus_tag} ({full_length:>3} aa)"
                    f" \u2192 mature: {len(mature_seq):>3} aa | {product[:55]}",
                    file=sys.stderr,
                )
            else:
                print(
                    f"   {locus_tag} ({full_length:>3} aa) | {product[:60]}",
                    file=sys.stderr,
                )

    if skipped_zero_length > 0:
        print(
            f"  [{skipped_zero_length} query protein(s) skipped due to zero-length mature core]",
            file=sys.stderr,
        )

    print(f"\n[*] Extracted {len(proteins)} protein(s) from query.\n", file=sys.stderr)
    return proteins


# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: PRE-FILTER HELPERS
# ─────────────────────────────────────────────────────────────────────────────


def _build_kmers(seq: str, k: int = _K) -> frozenset[str]:
    """Returns the frozen set of all k-mers in ``seq``.

    Args:
        seq: Amino acid sequence.
        k:   k-mer length.

    Returns:
        Frozen set of k-mer strings. Empty if ``seq`` is shorter than ``k``.
    """
    if len(seq) < k:
        return frozenset()
    return frozenset(seq[i : i + k] for i in range(len(seq) - k + 1))


def _passes_kmer_filter(
    kmers_a: frozenset,
    kmers_b: frozenset,
    min_identity: float,
    k: int = _K,
) -> bool:
    """Fast k-mer pre-screen before running the alignment.

    Uses the inequality ``Jaccard(k-mers) ≥ identity^k`` with a 50 %%
    safety slack. Pairs whose k-mer overlap is mathematically too low to
    ever meet ``min_identity`` are rejected here without any alignment.

    Args:
        kmers_a:      Pre-computed k-mers for sequence A.
        kmers_b:      Pre-computed k-mers for sequence B.
        min_identity: Identity threshold (0.0–1.0).
        k:            k-mer length used when building the sets.

    Returns:
        ``True`` to proceed to alignment; ``False`` to skip the pair.
    """
    if not kmers_a or not kmers_b:
        return True  # Short sequences pass; alignment will handle them
    union = len(kmers_a | kmers_b)
    if union == 0:
        return True
    jaccard = len(kmers_a & kmers_b) / union
    threshold = (min_identity**k) * 0.50
    return jaccard >= threshold


# ─────────────────────────────────────────────────────────────────────────────
# STEP 4: SCORE PRE-FILTER AND FULL ALIGNMENT
# ─────────────────────────────────────────────────────────────────────────────


def _passes_score_filter(
    seq_a: str,
    seq_b: str,
    min_identity: float,
) -> bool:
    """Score pre-screen using aligner.score() and S_max normalisation.

    Runs two ``aligner.score()`` calls, both entirely in Biopython's C engine
    with no traceback.  Together they are still many times faster than a single
    full ``aligner.align()`` call.

    Theory (BLOSUM62):
        For two proteins with fraction ``p`` of identical positions in a local
        alignment of length L, the expected alignment score is approximately::

            actual_score ≈ L × (6p − 1)

        derived from avg. BLOSUM62 diagonal ≈ +5 (matches) and avg. off-diagonal
        ≈ −1 (substitutions).

        The self-alignment score of the shorter sequence (S_max) is the absolute
        ceiling for any alignment between the two sequences::

            S_max ≈ L × avg_diagonal ≈ L × 5

        The correct safety factor must satisfy::

            factor < (6 × min_identity − 1) / (5 × min_identity)

        At min_identity = 0.30 this ceiling is 0.53.  Using **0.40** gives a
        25 %% margin below the ceiling and works correctly for any
        ``--identity`` value from 0.25 upwards.  A factor of 0.90 would cause
        false negatives at any threshold below ~65 %% identity.

    Args:
        seq_a:        First protein sequence.
        seq_b:        Second protein sequence.
        min_identity: Minimum identity threshold (0.0–1.0).

    Returns:
        ``True`` to proceed to full alignment; ``False`` to skip the pair.
    """
    # Step 1: actual local alignment score (C-level, no traceback)
    actual_score = _ALIGNER.score(seq_a, seq_b)

    # Fast-fail: score ≤ 0 means no positive local alignment exists at all
    if actual_score <= 0:
        return False

    # Step 2: theoretical maximum score (self-alignment of the shorter seq)
    # In local alignment, the shorter sequence caps the maximum match length.
    shorter_seq = seq_a if len(seq_a) <= len(seq_b) else seq_b
    max_possible_score = _ALIGNER.score(shorter_seq, shorter_seq)

    if max_possible_score <= 0:
        return False

    # Step 3: safety factor 0.40 — derived from BLOSUM62 expected score
    # equations. This is mathematically correct; 0.90 would cause false
    # negatives for any --identity setting below ~0.65.
    scaled_threshold = max_possible_score * min_identity * 0.40
    return actual_score >= scaled_threshold


def calculate_identity(seq_a: str, seq_b: str) -> tuple[float, int, int]:
    """Computes percent identity via Smith-Waterman alignment (traceback).

    This function is only called on pairs that have already passed the
    length ratio filter, k-mer filter, and score pre-filter.  The small
    number of surviving pairs makes the traceback safe from combinatorial
    explosion.

    Args:
        seq_a: First protein sequence.
        seq_b: Second protein sequence.

    Returns:
        A tuple of ``(percent_identity, alignment_length, mismatches)`` where
        ``percent_identity`` is 0.0–1.0,
        ``alignment_length`` is the number of positions in the local alignment
        including any gaps, and
        ``mismatches`` is the number of non-matching aligned positions
        (not counting gaps).
    """
    if not seq_a or not seq_b:
        return 0.0, 0, 0

    alignments = _ALIGNER.align(seq_a, seq_b)
    best = next(iter(alignments), None)
    if best is None:
        return 0.0, 0, 0

    aligned_a = best[0]
    aligned_b = best[1]
    alignment_length = len(aligned_a)

    if alignment_length == 0:
        return 0.0, 0, 0

    identical = sum(1 for a, b in zip(aligned_a, aligned_b) if a == b and a != "-")
    mismatches = sum(
        1 for a, b in zip(aligned_a, aligned_b) if a != b and a != "-" and b != "-"
    )
    return identical / alignment_length, alignment_length, mismatches


# ─────────────────────────────────────────────────────────────────────────────
# STEP 5: LOAD REFERENCE PROTEINS (pre-compute cmp_seq and k-mers ONCE)
# ─────────────────────────────────────────────────────────────────────────────


def _load_reference_proteins(
    ref_path: Path,
    use_mature: bool = False,
) -> list[_RefProtein]:
    """Loads all CDS proteins from a reference GenBank file.

    The comparison sequence (``cmp_seq``) and its k-mers are computed here,
    once per protein, BEFORE the inner query loop runs.  This means mature
    trimming and k-mer building are never executed inside the nested loop,
    regardless of how many query proteins are being compared.

    If ``use_mature=True`` and a protein's mature core becomes zero-length
    (e.g. the signal peptide is longer than the entire protein), the protein
    is silently skipped with a warning printed to stderr.

    Args:
        ref_path:   Path to the reference GBK/GBFF file.
        use_mature: If ``True``, applies ``calculate_mature_core()`` to each
                    protein so that ``cmp_seq`` is the trimmed mature core.

    Returns:
        List of ``_RefProtein`` objects with pre-computed comparison data.
    """
    proteins: list[_RefProtein] = []
    skipped_zero_length = 0

    for record in SeqIO.parse(ref_path, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue
            translation = feature.qualifiers.get("translation", [""])[0]
            if not translation:
                continue

            # Pre-compute the comparison sequence and its k-mers right here.
            # This is the core fix for Flaw 1: these operations now run exactly
            # once per reference protein, not once per (query × reference) pair.
            cmp_seq = calculate_mature_core(translation) if use_mature else translation

            # Skip proteins whose mature core is zero-length or None
            if not cmp_seq or len(cmp_seq) == 0:
                locus_tag = feature.qualifiers.get("locus_tag", ["UNKNOWN"])[0]
                print(
                    f"  [!] Skipping {locus_tag}: mature core is zero-length. "
                    f"(Signal peptide may be longer than the full protein.)",
                    file=sys.stderr,
                )
                skipped_zero_length += 1
                continue

            proteins.append(
                _RefProtein(
                    locus_tag=feature.qualifiers.get("locus_tag", ["UNKNOWN"])[0],
                    product=feature.qualifiers.get("product", ["Unknown product"])[0],
                    sequence=translation,
                    cmp_seq=cmp_seq,
                    kmers=_build_kmers(cmp_seq),
                    length=len(cmp_seq),
                )
            )

    if skipped_zero_length > 0:
        print(
            f"  [{skipped_zero_length} protein(s) skipped due to zero-length mature core]",
            file=sys.stderr,
        )

    return proteins


# ─────────────────────────────────────────────────────────────────────────────
# STEP 6: COMPARE QUERY PROTEINS AGAINST A REFERENCE FILE
# ─────────────────────────────────────────────────────────────────────────────


def find_orthologs(
    query_proteins: list[Protein],
    ref_path: Path,
    min_identity: float,
    use_mature: bool,
    min_coverage: float = 0.50,
) -> list[OrthoHit]:
    """Compares query proteins against all CDS in a reference file.

    Applies four filters in order of increasing computational cost. The
    majority of pairs (typically >99 %%) are rejected before the full
    alignment ever runs.

    Args:
        query_proteins: ``Protein`` objects from the query GBK.
        ref_path:       Path to the reference file.
        min_identity:   Minimum identity (0.0–1.0) to report a hit.
        use_mature:     If ``True``, uses mature sequences for comparison.
        min_coverage:   Minimum alignment coverage of the shorter sequence.

    Returns:
        List of ``OrthoHit`` objects for all pairs passing all filters.
    """
    hits: list[OrthoHit] = []

    # Load reference proteins with cmp_seq and k-mers pre-computed.
    # calculate_mature_core() and _build_kmers() are called here ONCE per
    # reference protein — never again inside the nested query loop below.
    ref_proteins = _load_reference_proteins(ref_path, use_mature=use_mature)
    if not ref_proteins:
        return hits

    for qprotein in query_proteins:
        query_seq = qprotein.mature_sequence if use_mature else qprotein.sequence
        query_kmers = _build_kmers(query_seq)
        query_len = len(query_seq)

        for ref in ref_proteins:

            # Filter 1 — length ratio (O(1)) ──────────────────────────────────
            ratio = query_len / max(ref.length, 1)
            if ratio < 0.5 or ratio > 2.0:
                continue

            # Filter 2 — k-mer Jaccard (O(n+m) Python) ────────────────────────
            if not _passes_kmer_filter(query_kmers, ref.kmers, min_identity):
                continue

            # Filter 3 — score pre-filter (O(n×m) C, NO traceback) ────────────
            # This is the critical guard against Biopython's traceback
            # combinatorial explosion. aligner.score() runs entirely in C and
            # returns immediately without enumerating alignment paths.
            if not _passes_score_filter(query_seq, ref.cmp_seq, min_identity):
                continue

            # Filter 4 — full alignment with traceback ─────────────────────────
            # Only a tiny fraction of pairs reach this point.
            identity, aln_length, mismatches = calculate_identity(
                query_seq, ref.cmp_seq
            )
            if identity < min_identity:
                continue

            # Filter 5 — coverage check ────────────────────────────────────────
            shorter = min(query_len, ref.length)
            coverage = aln_length / shorter if shorter > 0 else 0.0
            if coverage < min_coverage:
                continue

            hits.append(
                OrthoHit(
                    query_locus=qprotein.locus_tag,
                    query_product=qprotein.product,
                    query_seq=query_seq,
                    ref_locus=ref.locus_tag,
                    ref_product=ref.product,
                    ref_seq=ref.cmp_seq,  # Consistent: same seq as used for alignment
                    ref_file=ref_path.stem,
                    identity=identity,
                    alignment_length=aln_length,
                    mismatches=mismatches,
                    query_length=query_len,
                    ref_length=ref.length,
                )
            )

    return hits


# ─────────────────────────────────────────────────────────────────────────────
# STEP 7: PRETTY PRINT TABLE
# ─────────────────────────────────────────────────────────────────────────────


def print_hits_table(hits: list[OrthoHit]) -> None:
    """Prints ortholog hits as a formatted, aligned-column terminal table.

    Sorted by query locus tag, then by identity descending.

    Args:
        hits: List of ``OrthoHit`` objects to display.
    """
    if not hits:
        print("No hits found above the identity and coverage thresholds.")
        return

    cQL, cQP = 15, 23
    cRL, cRP = 14, 25
    cRF = 23
    cID, cAL, cQL2, cRL2 = 10, 7, 9, 7

    header = (
        f"{'Query Locus':<{cQL}} | {'Query Product':<{cQP}} | "
        f"{'Ref Locus':<{cRL}} | {'Ref Product':<{cRP}} | "
        f"{'Ref File':<{cRF}} | {'Identity %':<{cID}} | "
        f"{'Aln Len':<{cAL}} | {'Query Len':<{cQL2}} | {'Ref Len':<{cRL2}}"
    )
    sep = "-" * len(header)
    print(header)
    print(sep)

    for hit in sorted(hits, key=lambda h: (h.query_locus, -h.identity)):
        qp = (
            (hit.query_product[:20] + "...")
            if len(hit.query_product) > 23
            else hit.query_product
        )
        rp = (
            (hit.ref_product[:22] + "...")
            if len(hit.ref_product) > 25
            else hit.ref_product
        )
        print(
            f"{hit.query_locus:<{cQL}} | {qp:<{cQP}} | "
            f"{hit.ref_locus:<{cRL}} | {rp:<{cRP}} | "
            f"{hit.ref_file:<{cRF}} | {hit.identity*100:>{cID-1}.2f}% | "
            f"{hit.alignment_length:>{cAL}} | {hit.query_length:>{cQL2}} | "
            f"{hit.ref_length:>{cRL2}}"
        )

    print(sep)
    print(f"Total hits: {len(hits)}\n")


# ─────────────────────────────────────────────────────────────────────────────
# STEP 8: TSV OUTPUT
# ─────────────────────────────────────────────────────────────────────────────

_TSV_HEADERS = [
    "query_locus",
    "query_product",
    "ref_locus",
    "ref_product",
    "ref_file",
    "query_sequence",
    "ref_sequence",
    "mismatches",
    "identity_pct",
    "alignment_length",
    "query_length",
    "ref_length",
]


def write_tsv(hits: list[OrthoHit], out_handle) -> None:
    """Writes all ortholog hits to a TSV file with protein sequences.

    Args:
        hits:       ``OrthoHit`` results to write.
        out_handle: Open file handle or ``sys.stdout``.
    """
    writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
    writer.writerow(_TSV_HEADERS)
    for hit in sorted(hits, key=lambda h: (h.query_locus, -h.identity)):
        writer.writerow(
            [
                hit.query_locus,
                hit.query_product,
                hit.ref_locus,
                hit.ref_product,
                hit.ref_file,
                hit.query_seq,
                hit.ref_seq,
                hit.mismatches,
                f"{hit.identity * 100:.2f}",
                hit.alignment_length,
                hit.query_length,
                hit.ref_length,
            ]
        )


# ─────────────────────────────────────────────────────────────────────────────
# STEP 9: MAIN PIPELINE
# ─────────────────────────────────────────────────────────────────────────────


def main() -> None:
    """Parses arguments and runs the full ortholog-finding pipeline."""
    args = get_args()

    print("=" * 100, file=sys.stderr)
    print(f"GBK ORTHOLOG FINDER v{__version__}", file=sys.stderr)
    print("=" * 100, file=sys.stderr)
    print(f"  Query        : {args.query}", file=sys.stderr)
    print(f"  Reference    : {args.reference}", file=sys.stderr)
    print(f"  Min identity : {args.identity * 100:.0f}%", file=sys.stderr)
    print(f"  Min coverage : {args.min_coverage * 100:.0f}%", file=sys.stderr)
    print(
        f"  Mature core  : {'YES — signal peptides trimmed' if args.mature else 'NO'}",
        file=sys.stderr,
    )
    if args.max_length:
        print(f"  Max length   : {args.max_length} aa", file=sys.stderr)
    print("=" * 100, file=sys.stderr)

    try:
        query_proteins = extract_proteins_from_gbk(
            gbk_path=args.query,
            apply_mature=args.mature,
            min_length=args.min_length,
            max_length=args.max_length,
        )
        if not query_proteins:
            sys.exit("[!] No proteins extracted from query. Check the file.")

        ref_files = list(stream_reference_files(args.reference))
        if not ref_files:
            sys.exit("[!] No valid reference files found.")

        all_hits: list[OrthoHit] = []
        for i, ref_file in enumerate(ref_files, 1):
            print(
                f"[*] [{i}/{len(ref_files)}] Scanning {ref_file.name} ...",
                file=sys.stderr,
            )
            file_hits = find_orthologs(
                query_proteins=query_proteins,
                ref_path=ref_file,
                min_identity=args.identity,
                use_mature=args.mature,
                min_coverage=args.min_coverage,
            )
            print(
                f"      \u2192 {len(file_hits)} hit(s) above "
                f"{args.identity*100:.0f}% identity / "
                f"{args.min_coverage*100:.0f}% coverage",
                file=sys.stderr,
            )
            all_hits.extend(file_hits)

        print(f"\n[*] Total hits: {len(all_hits)}\n", file=sys.stderr)

        if args.output:
            print(
                f"[*] Saving TSV \u2192 {args.output.resolve()}",
                file=sys.stderr,
            )
            with smart_open(args.output) as out_handle:
                write_tsv(all_hits, out_handle)
        else:
            print(
                "[*] Tip: add -o results.tsv to save results with protein sequences.",
                file=sys.stderr,
            )
            print()
            print_hits_table(all_hits)

        print("=" * 100, file=sys.stderr)
        print("[*] Done.", file=sys.stderr)
        print("=" * 100, file=sys.stderr)

    except ValueError as exc:
        sys.exit(f"\n[!] Error: {exc}")
    except KeyboardInterrupt:
        sys.exit("\n[!] Interrupted by user.")


if __name__ == "__main__":
    main()
