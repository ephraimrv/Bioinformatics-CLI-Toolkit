#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

"""GBK Homolog Finder — pairwise protein homology detection from GenBank files.

Extracts protein-coding sequences from a query GBK/GBFF file and compares
them against one or more reference genomes using local pairwise alignment
(Smith-Waterman with BLOSUM62) to identify homologs by percent identity.

Works on prokaryotic and eukaryotic GenBank files. Accepts a single
reference file or a directory of reference files.

ISOFORM HANDLING (eukaryotic genomes):
    Eukaryotic GenBank annotations represent alternative splice isoforms as
    separate CDS features that all share the same ``/locus_tag`` (the gene
    identifier) — only ``/protein_id`` differs per isoform. Without
    deduplication, a gene with N query isoforms and M reference isoforms
    produces N x M near-duplicate comparisons for what is biologically one
    gene, flooding output with redundant rows and risking a combinatorial
    blow-up in runtime at genome scale. Both ``extract_proteins_from_gbk()``
    and ``_load_reference_proteins()`` keep only the LONGEST isoform per
    ``locus_tag`` — the same fallback heuristic tools like OrthoFinder use
    when no canonical-transcript annotation is available. Prokaryotic
    genomes are unaffected: each locus_tag already maps to exactly one CDS.

Performance pipeline (each stage cheaper than the last):
    1. Length ratio filter  — skips pairs outside 0.5–2.0 length ratio (O(1)).
       Only active when ``--coverage-mode max``. Skipped in ``min`` mode to
       allow small peptides (bacteriocins) to match within larger proteins.
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

Flag interaction — ``--mature``, ``--min-coverage``, and ``--coverage-mode``:
    **Coverage mode choice determines your search goal:**

    ``--coverage-mode min`` (default; bacteriocin/domain-centric search):
        Coverage is measured against the shorter sequence. A 40 aa bacteriocin
        core may be fully present (100%% coverage) inside a 1000 aa multi-domain
        reference protein. The length ratio filter is skipped to allow this.
        Use this when searching for conserved peptides, domains, or bacteriocins.

    ``--coverage-mode max`` (whole-protein homolog search):
        Coverage is measured against the longer sequence. Both proteins must be
        substantially similar in size (0.5–2.0 ratio enforced). Use this for
        strict, size-matched 1:1 candidate identification in full-genome
        comparisons. (Confirming true 1:1 *orthology*, as opposed to homology,
        additionally requires reciprocal-best-hit or phylogenetic analysis —
        see "Terminology note" below.)

    ``--mature`` trims signal peptides from both query and reference proteins
    before alignment, so only the bioactive mature core is compared. Recommended
    for bacteriocins, lantibiotics, and secreted peptides.

    ``--min-coverage`` (0.0–1.0) gates hits on alignment coverage, computed
    from non-gap aligned residues only (gap characters in the alignment do
    not count toward either sequence's coverage). At 0.65 with
    ``--coverage-mode min``, you require at least 65%% of the shorter sequence
    to be covered by non-gap aligned residues. With ``--coverage-mode max``,
    the same applies to the longer sequence (much stricter for size-mismatched
    pairs).

    **Recommended combinations:**
        - Bacteriocin search: ``--mature --coverage-mode min --min-coverage 0.65``
        - Whole-genome homologs: ``--coverage-mode max --min-coverage 0.75``
        - General comparative: ``--coverage-mode min --min-coverage 0.50``

Terminology note:
    This tool performs pairwise local alignment, which establishes
    **homology** (sequences sharing common ancestry) — it does not by
    itself distinguish true **orthologs** (separated by a speciation event)
    from **paralogs** (separated by a gene duplication event). In
    gene-family-rich genomes (e.g. many eukaryotic kinases, cytochromes,
    transporters), a query will often match several paralogs above the
    identity/coverage threshold, all reported here as hits — not
    specifically as orthologs. Establishing true orthology requires
    Reciprocal Best Hit (RBH) analysis or phylogenetic tree reconciliation,
    both out of scope for this tool. Treat hits here as candidate homologs
    requiring that additional confirmation before describing them as
    orthologs in a manuscript or other formal report.

Note:
    This script is part of ongoing research and is associated with an upcoming
    publication. Correct attribution is requested when used in derivative works.
    Released under the MIT License. See LICENSE in the repository root.

    v1.7: Fixed a critical data-loss bug in ``extract_proteins_from_gbk()``
    and ``_load_reference_proteins()``. Both previously defaulted every CDS
    lacking a ``/locus_tag`` qualifier to the literal string ``"UNKNOWN"``.
    On files where ``/locus_tag`` is absent entirely — common for
    eukaryotic assemblies, GFF3-to-GenBank conversions, and some draft
    genomes — every CDS in the file collided on that one dictionary key,
    and the longest-isoform-wins deduplication silently discarded all but
    the single longest CDS in the whole file. Confirmed empirically: a
    simulated 15,000-CDS file with no locus_tag annotations dropped to
    exactly 1 surviving protein. Both functions now use
    ``_resolve_identifier()``, a fallback hierarchy (locus_tag ->
    protein_id -> gene -> contig-ID + coordinates) that always yields a
    safe, sufficiently-unique key — the final coordinate fallback includes
    the parent record's ID specifically because coordinates alone are only
    unique within a single contig, and fragmented draft assemblies
    routinely have many small contigs starting near position 0.

Example:
    Bacteriocin screen with signal peptide trimming and domain-centric search::

        python3 pairwise_homolog_finder.py \\
            -q region001.gbk -r references/ \\
            --mature --coverage-mode min --max-length 150 \\
            --identity 0.35 --min-coverage 0.65 \\
            -o bacteriocin_hits.tsv

    Whole-protein homolog search across genomes::

        python3 pairwise_homolog_finder.py \\
            -q genome.gbff -r references/ \\
            --coverage-mode max --identity 0.40 --min-coverage 0.75 \\
            -o homolog_hits.tsv

    Region vs single reference (default mode)::

        python3 pairwise_homolog_finder.py \\
            -q region001.gbk -r ATCC8293.gbff -o results.tsv
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.7"

import sys
import argparse
import csv
from pathlib import Path
from dataclasses import dataclass
from datetime import datetime
import time

try:
    from Bio import SeqIO
    from Bio.Align import substitution_matrices, PairwiseAligner
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )

# Optional: tqdm for progress bars. Works without it, but less elegant.
try:
    from tqdm import tqdm

    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False
    tqdm = lambda x, **kwargs: x  # Fallback: just iterate normally

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
    length: int  # len(cmp_seq) — used for ratio/coverage math (comparison-space length)
    full_length: (
        int  # len(sequence) — used only to pick the longest isoform per locus_tag
    )


@dataclass
class HomologHit:
    """A single homolog comparison result.

    Terminology note: this is the result of pairwise local alignment, which
    establishes homology (shared ancestry), not orthology specifically. See
    the module docstring's "Terminology note" section.
    """

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
    query_coverage: float  # Non-gap query residues / query_length (0.0-1.0)
    ref_coverage: float  # Non-gap ref residues / ref_length (0.0-1.0)
    coverage: float  # The coverage value actually used for filtering (mode-dependent)


# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: ARGUMENT PARSER
# ─────────────────────────────────────────────────────────────────────────────


def get_args() -> argparse.Namespace:
    """Configures and returns the CLI argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Find homologs by comparing proteins from a query GenBank file "
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
        "--output-fasta",
        type=Path,
        default=None,
        help=(
            "Save all homolog protein sequences to a FASTA file (.faa format). "
            "Useful for downstream analysis: multiple sequence alignment (MAFFT), "
            "phylogenetics (IQ-TREE), motif detection (HMMER), etc. "
            "If omitted, no FASTA file is generated. "
            "Example: --output-fasta homologs.faa"
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
            "How much of the protein must be covered by the alignment "
            "(0.0–1.0, i.e. 0.50 = 50%%). "
            "This prevents hits where only a small conserved domain (e.g. a signal "
            "peptide or a zinc-finger motif) matches, while the rest of the protein "
            "is unrelated. When used with --mature, coverage is measured against the "
            "trimmed mature core, so 0.70 means '70%% of the active peptide aligns'. "
            "Recommended range: 0.50 (permissive) to 0.80 (strict). Default: 0.50."
        ),
    )
    parser.add_argument(
        "--coverage-mode",
        choices=["min", "max"],
        default="min",
        help=(
            "Which sequence length to use as the denominator for coverage calculation. "
            "'min' uses the shorter sequence length — best for bacteriocin cores and "
            "conserved domain searches, where a small peptide (50 aa) may be fully "
            "present inside a much larger protein (1000 aa) and you want 100%% coverage "
            "reported in that case. "
            "'max' uses the longer sequence length — best for strict whole-protein "
            "homology searches, where you want both proteins to be substantially "
            "similar in length. Using 'max' also enforces a 0.5–2.0 length ratio "
            "pre-filter so extremely size-mismatched pairs are skipped early. "
            "Default: min."
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
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help=(
            "Increase output verbosity. By default, only major milestones and "
            "a progress bar are shown. With --verbose, every file scanned and "
            "every protein extracted is printed. Useful for debugging."
        ),
    )

    return parser.parse_args()


# ─────────────────────────────────────────────────────────────────────────────
# IDENTIFIER FALLBACK (fixes the "UNKNOWN" dictionary-collision data-loss bug)
# ─────────────────────────────────────────────────────────────────────────────


def _resolve_identifier(feature, record_id: str) -> str:
    """Resolves a unique grouping key for a CDS feature, with safe fallbacks.

    RefSeq bacterial GenBank files always carry ``/locus_tag``, but many
    eukaryotic assemblies, GFF3-to-GenBank conversions (MAKER, Augustus),
    and draft genomes omit it entirely. Previously this function's caller
    defaulted every such feature to the literal string ``"UNKNOWN"`` — on
    a file where NO feature has a locus_tag, every one of potentially
    thousands of CDS features collided on that single dictionary key, and
    the longest-isoform-wins logic in ``extract_proteins_from_gbk()`` /
    ``_load_reference_proteins()`` silently discarded all but the single
    longest CDS in the entire file. Confirmed empirically: a simulated
    15,000-CDS file with no locus_tag annotations dropped to exactly 1
    surviving protein.

    Fallback order (each step only used if the previous one is empty):
      1. ``/locus_tag``    — standard, shared across splice isoforms of one
                             gene (the correct grouping key when present).
      2. ``/protein_id``   — typically unique PER ISOFORM rather than per
                             gene, so falling back to it means isoform
                             deduplication won't collapse multiple isoforms
                             together for this feature. That's an accepted,
                             far smaller cost than the alternative (total
                             data loss): keeping a few extra near-duplicate
                             isoform entries is harmless; losing 99.99% of
                             the file's proteins is not.
      3. ``/gene``         — also typically shared across isoforms of one
                             gene, so this restores correct isoform grouping
                             for files that have gene symbols but neither of
                             the two qualifiers above.
      4. ``record.id`` + genomic coordinates — absolute last resort,
                             mathematically guaranteed unique. The contig/
                             record ID is included deliberately: coordinates
                             ALONE are only unique within a single contig,
                             and fragmented draft assemblies routinely have
                             many small contigs each starting near position
                             0 — omitting the contig ID would silently
                             reintroduce the same collision bug this
                             function exists to prevent, just at a smaller,
                             rarer scale. Confirmed empirically: two
                             different contigs each with a CDS at the same
                             coordinates collide under a coordinates-only
                             fallback.

    Args:
        feature:   A Biopython SeqFeature (CDS) to resolve an identifier for.
        record_id: The parent record's ``.id`` (contig/chromosome name),
                   used only in the final coordinate-based fallback.

    Returns:
        A non-empty string suitable for use as a deduplication dictionary
        key. Never ``"UNKNOWN"`` or any other shared constant.
    """
    identifier = feature.qualifiers.get("locus_tag", [""])[0]
    if identifier:
        return identifier

    identifier = feature.qualifiers.get("protein_id", [""])[0]
    if identifier:
        return identifier

    identifier = feature.qualifiers.get("gene", [""])[0]
    if identifier:
        return identifier

    start = int(feature.location.start)
    end = int(feature.location.end)
    return f"UNANNOTATED_{record_id}_{start}_{end}"


# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: EXTRACT PROTEINS FROM QUERY GBK
# ─────────────────────────────────────────────────────────────────────────────


def extract_proteins_from_gbk(
    gbk_path: Path,
    apply_mature: bool = False,
    min_length: int = 10,
    max_length: int | None = None,
    verbose: bool = False,
) -> list[Protein]:
    """Reads a GBK/GBFF file and returns the longest CDS per locus_tag.

    On a eukaryotic file, multiple CDS features (alternative splice
    isoforms) commonly share the same /locus_tag. Only the longest
    full-length translation per locus_tag is kept — see the module
    docstring's "ISOFORM HANDLING" section. On a prokaryotic file this is
    a no-op, since each locus_tag already has exactly one CDS.

    Args:
        gbk_path:     Path to the query GBK/GBFF file.
        apply_mature: If ``True``, applies ``calculate_mature_core()`` to each protein.
        min_length:   Skip proteins shorter than this.
        max_length:   Skip proteins longer than this (``None`` = no limit).
        verbose:      If ``True``, print details for each protein retained.

    Returns:
        List of ``Protein`` objects, one per distinct resolved identifier
        (locus_tag, or a safe fallback — see ``_resolve_identifier()``).
    """
    print(f"\n[*] Extracting proteins from query: {gbk_path.name}", file=sys.stderr)

    locus_best: dict[str, Protein] = {}
    skipped_zero_length = 0
    total_candidates = 0

    for record in SeqIO.parse(gbk_path, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue
            translation = feature.qualifiers.get("translation", [""])[0]
            if not translation:
                continue

            locus_tag = _resolve_identifier(feature, record.id)
            product = feature.qualifiers.get("product", ["Unknown product"])[0]
            full_length = len(translation)

            if full_length < min_length:
                continue
            if max_length is not None and full_length > max_length:
                continue

            total_candidates += 1

            mature_seq = (
                calculate_mature_core(translation) if apply_mature else translation
            )

            # Skip candidates whose mature core is zero-length
            if apply_mature and (not mature_seq or len(mature_seq) == 0):
                skipped_zero_length += 1
                continue

            existing = locus_best.get(locus_tag)
            if existing is not None and existing.length >= full_length:
                continue  # Already holding an equal-or-longer isoform for this locus

            locus_best[locus_tag] = Protein(
                locus_tag=locus_tag,
                product=product,
                sequence=translation,
                mature_sequence=mature_seq,
                source_file=gbk_path.name,
                length=full_length,
            )

    proteins = list(locus_best.values())

    if total_candidates > len(proteins):
        print(
            f"  [{total_candidates - len(proteins)} isoform(s) collapsed — "
            f"kept the longest CDS per locus_tag]",
            file=sys.stderr,
        )
    if skipped_zero_length > 0:
        print(
            f"  [{skipped_zero_length} candidate isoform(s) skipped due to "
            f"zero-length mature core]",
            file=sys.stderr,
        )

    if verbose:
        for p in proteins:
            if apply_mature and p.mature_sequence != p.sequence:
                print(
                    f"   {p.locus_tag} ({p.length:>3} aa)"
                    f" \u2192 mature: {len(p.mature_sequence):>3} aa | {p.product[:55]}",
                    file=sys.stderr,
                )
            else:
                print(
                    f"   {p.locus_tag} ({p.length:>3} aa) | {p.product[:60]}",
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


def calculate_identity(seq_a: str, seq_b: str) -> tuple[float, int, int, int, int]:
    """Computes percent identity via Smith-Waterman alignment (traceback).

    This function is only called on pairs that have already passed the
    length ratio filter, k-mer filter, and score pre-filter.  The small
    number of surviving pairs makes the traceback safe from combinatorial
    explosion.

    Args:
        seq_a: First protein sequence.
        seq_b: Second protein sequence.

    Returns:
        A tuple of ``(percent_identity, alignment_length, mismatches,
        a_nongap, b_nongap)`` where:
            ``percent_identity`` is 0.0–1.0 (identical positions over
            ``alignment_length`` — the standard, gap-inclusive definition,
            matching how BLAST/EMBOSS report identity).
            ``alignment_length`` is the number of columns in the local
            alignment INCLUDING gap columns on either side.
            ``mismatches`` is the number of non-matching aligned positions
            (not counting gaps).
            ``a_nongap``/``b_nongap`` are the count of non-gap residues each
            sequence actually contributes to the alignment. These are NOT
            the same as ``alignment_length`` — if one sequence has a large
            insertion relative to the other, ``alignment_length`` can exceed
            either original sequence's own length. Coverage must be computed
            from ``a_nongap``/``b_nongap`` against each sequence's own full
            length, never from ``alignment_length`` directly (doing so can
            produce coverage values exceeding 100%).
    """
    if not seq_a or not seq_b:
        return 0.0, 0, 0, 0, 0

    alignments = _ALIGNER.align(seq_a, seq_b)
    best = next(iter(alignments), None)
    if best is None:
        return 0.0, 0, 0, 0, 0

    aligned_a = best[0]
    aligned_b = best[1]
    alignment_length = len(aligned_a)

    if alignment_length == 0:
        return 0.0, 0, 0, 0, 0

    identical = sum(1 for a, b in zip(aligned_a, aligned_b) if a == b and a != "-")
    mismatches = sum(
        1 for a, b in zip(aligned_a, aligned_b) if a != b and a != "-" and b != "-"
    )
    a_nongap = sum(1 for c in aligned_a if c != "-")
    b_nongap = sum(1 for c in aligned_b if c != "-")
    return (
        identical / alignment_length,
        alignment_length,
        mismatches,
        a_nongap,
        b_nongap,
    )


# ─────────────────────────────────────────────────────────────────────────────
# STEP 5: LOAD REFERENCE PROTEINS (pre-compute cmp_seq and k-mers ONCE)
# ─────────────────────────────────────────────────────────────────────────────


def _load_reference_proteins(
    ref_path: Path,
    use_mature: bool = False,
) -> list[_RefProtein]:
    """Loads the longest-isoform CDS protein per locus_tag from a reference file.

    The comparison sequence (``cmp_seq``) and its k-mers are computed here,
    once per protein, BEFORE the inner query loop runs.  This means mature
    trimming and k-mer building are never executed inside the nested loop,
    regardless of how many query proteins are being compared.

    On a eukaryotic reference, multiple CDS features (splice isoforms)
    commonly share the same /locus_tag; only the longest full-length
    translation per locus_tag is kept (see module docstring, "ISOFORM
    HANDLING"). On a prokaryotic reference this is a no-op.

    If ``use_mature=True`` and a candidate's mature core becomes zero-length
    (e.g. the signal peptide is longer than the entire protein), the
    candidate is skipped and a warning is printed to stderr. This provides
    transparency about what the script did with your data.

    Args:
        ref_path:   Path to the reference GBK/GBFF file.
        use_mature: If ``True``, applies ``calculate_mature_core()`` to each
                    protein so that ``cmp_seq`` is the trimmed mature core.

    Returns:
        List of ``_RefProtein`` objects, one per distinct resolved
        identifier (locus_tag, or a safe fallback — see
        ``_resolve_identifier()``), with pre-computed comparison data.
    """
    locus_best: dict[str, _RefProtein] = {}
    skipped_zero_length = 0
    total_candidates = 0

    for record in SeqIO.parse(ref_path, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue
            translation = feature.qualifiers.get("translation", [""])[0]
            if not translation:
                continue

            total_candidates += 1
            locus_tag = _resolve_identifier(feature, record.id)
            full_length = len(translation)

            # Pre-compute the comparison sequence and its k-mers right here.
            # This is the core fix for Flaw 1: these operations now run exactly
            # once per reference protein, not once per (query × reference) pair.
            cmp_seq = calculate_mature_core(translation) if use_mature else translation

            # Skip candidates whose mature core is zero-length or None
            if not cmp_seq or len(cmp_seq) == 0:
                msg = (
                    f"  [!] Skipping {locus_tag}: mature core is zero-length. "
                    f"(Signal peptide may be longer than the full protein.)"
                )
                if HAS_TQDM:
                    tqdm.write(msg, file=sys.stderr)
                else:
                    print(msg, file=sys.stderr)
                skipped_zero_length += 1
                continue

            existing = locus_best.get(locus_tag)
            if existing is not None and existing.full_length >= full_length:
                continue  # Already holding an equal-or-longer isoform

            locus_best[locus_tag] = _RefProtein(
                locus_tag=locus_tag,
                product=feature.qualifiers.get("product", ["Unknown product"])[0],
                sequence=translation,
                cmp_seq=cmp_seq,
                kmers=_build_kmers(cmp_seq),
                length=len(cmp_seq),
                full_length=full_length,
            )

    proteins = list(locus_best.values())

    if total_candidates > len(proteins):
        collapse_msg = (
            f"  [{total_candidates - len(proteins)} isoform(s) collapsed — "
            f"kept the longest CDS per locus_tag]"
        )
        if HAS_TQDM:
            tqdm.write(collapse_msg, file=sys.stderr)
        else:
            print(collapse_msg, file=sys.stderr)

    if skipped_zero_length > 0:
        summary_msg = f"  [{skipped_zero_length} protein(s) skipped due to zero-length mature core]"
        # tqdm.write() handles output properly when progress bar is active
        if HAS_TQDM:
            tqdm.write(summary_msg, file=sys.stderr)
        else:
            print(summary_msg, file=sys.stderr)

    return proteins


# ─────────────────────────────────────────────────────────────────────────────
# STEP 6: COMPARE QUERY PROTEINS AGAINST A REFERENCE FILE
# ─────────────────────────────────────────────────────────────────────────────
def find_homologs(
    query_proteins: list[Protein],
    ref_path: Path,
    min_identity: float,
    use_mature: bool,
    min_coverage: float = 0.50,
    coverage_mode: str = "min",
) -> list[HomologHit]:
    """Compares query proteins against all CDS in a reference file.

    Applies filters in order of increasing computational cost. The
    majority of pairs (typically >99 %%) are rejected before the full
    alignment ever runs.

    Args:
        query_proteins: ``Protein`` objects from the query GBK.
        ref_path:       Path to the reference file.
        min_identity:   Minimum identity (0.0–1.0) to report a hit.
        use_mature:     If ``True``, uses mature sequences for comparison.
        min_coverage:   Minimum alignment coverage fraction to report a hit.
        coverage_mode:  ``'min'`` reports coverage relative to the shorter
                        sequence (correct for bacteriocin/domain searches).
                        ``'max'`` reports coverage relative to the longer
                        sequence (correct for whole-protein homolog search).
                        When ``'max'``, the length ratio pre-filter is also
                        enforced to reject size-mismatched pairs early.
        verbose:        If ``True``, print details during processing.

    Returns:
        List of ``HomologHit`` objects for all pairs passing all filters.
    """
    hits: list[HomologHit] = []

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
            # Only enforced in 'max' mode (whole-protein homolog search).
            # In 'min' mode (domain/bacteriocin search) a small peptide (40 aa)
            # may legitimately match one domain of a large protein (1000 aa);
            # a strict ratio filter would incorrectly reject those pairs.
            if coverage_mode == "max":
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
            identity, aln_length, mismatches, q_nongap, r_nongap = calculate_identity(
                query_seq, ref.cmp_seq
            )
            if identity < min_identity:
                continue

            # Filter 5 — coverage check (gap-aware) ─────────────────────────────
            # BUG THIS FIXES: previously used aln_length (gap-INCLUSIVE column
            # count) over min/max(query_len, ref.length) as the coverage
            # fraction. If the reference has a large insertion relative to
            # the query WITHIN the locally-aligned span, aln_length can exceed
            # either sequence's own real length, producing coverage > 100%
            # and corrupting downstream TSV filtering. Coverage must instead
            # be computed from each sequence's own NON-GAP residue count,
            # normalized against that same sequence's own full length.
            query_coverage = q_nongap / query_len if query_len > 0 else 0.0
            ref_coverage = r_nongap / ref.length if ref.length > 0 else 0.0
            coverage = (
                min(query_coverage, ref_coverage)
                if coverage_mode == "max"
                else max(query_coverage, ref_coverage)
            )
            if coverage < min_coverage:
                continue

            hits.append(
                HomologHit(
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
                    query_coverage=query_coverage,
                    ref_coverage=ref_coverage,
                    coverage=coverage,
                )
            )

    return hits


# ─────────────────────────────────────────────────────────────────────────────
# STEP 7: PRETTY PRINT TABLE
# ─────────────────────────────────────────────────────────────────────────────


def print_hits_table(hits: list[HomologHit]) -> None:
    """Prints homolog hits as a formatted, aligned-column terminal table.

    Sorted by query locus tag, then by identity descending.

    Args:
        hits: List of ``HomologHit`` objects to display.
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
    "coverage_pct",
    "alignment_length",
    "query_length",
    "ref_length",
]


def write_tsv(hits: list[HomologHit], out_handle) -> None:
    """Writes all homolog hits to a TSV file with protein sequences.

    Args:
        hits:       ``HomologHit`` results to write.
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
                f"{hit.coverage * 100:.2f}",
                hit.alignment_length,
                hit.query_length,
                hit.ref_length,
            ]
        )


def write_fasta(hits: list[HomologHit], fasta_path: Path) -> None:
    """Writes all homolog protein sequences to a FASTA file.

    Sequences are written in order: query first, then all references grouped
    by reference file. Each FASTA header includes locus tag, organism/file,
    product name, and identity percentage for easy identification.

    This output is suitable for downstream analysis:
        - Multiple sequence alignment (MAFFT, Clustal Omega)
        - Phylogenetic inference (IQ-TREE, RAxML)
        - Motif/domain detection (HMMER, InterProScan)
        - Sequence logos (WebLogo, ggseqlogo)

    Args:
        hits:       ``HomologHit`` results to write.
        fasta_path: Path to the output FASTA file.
    """
    if not hits:
        print("[!] No hits to write to FASTA.", file=sys.stderr)
        return

    with open(fasta_path, "w", encoding="utf-8") as fh:
        # Track which sequences we've written to avoid duplicates
        written = set()

        # Write query sequence(s) first
        # Group hits by query locus to handle cases where one query matched multiple refs
        query_hits = {}
        for hit in hits:
            if hit.query_locus not in query_hits:
                query_hits[hit.query_locus] = hit

        for hit in sorted(query_hits.values(), key=lambda hit: hit.query_locus):
            header = (
                f">{hit.query_locus} | {hit.query_product} | "
                f"Query | {hit.query_length}aa"
            )
            fh.write(f"{header}\n")
            fh.write(f"{hit.query_seq}\n")
            written.add((hit.query_locus, "query"))

        # Write reference sequences grouped by reference file
        # Sort by ref_file, then by ref_locus for consistent output
        sorted_hits = sorted(hits, key=lambda h: (h.ref_file, h.ref_locus, -h.identity))

        for hit in sorted_hits:
            seq_id = (hit.ref_locus, hit.ref_file)
            if seq_id in written:
                continue  # Skip if we've already written this ref sequence

            header = (
                f">{hit.ref_locus} | {hit.ref_product} | "
                f"{hit.ref_file} | {hit.identity*100:.2f}% identity | {hit.ref_length}aa"
            )
            fh.write(f"{header}\n")
            fh.write(f"{hit.ref_seq}\n")
            written.add(seq_id)


# ─────────────────────────────────────────────────────────────────────────────
# STEP 10: LOG FILE GENERATION
# ─────────────────────────────────────────────────────────────────────────────


def write_log_file(
    args: argparse.Namespace,
    query_proteins: list[Protein],
    ref_files: list[Path],
    all_hits: list[HomologHit],
    runtime_seconds: float,
    output_prefix: str,
) -> None:
    """Writes a comprehensive log/parameters file documenting the homolog search.

    This file mirrors the format used by BLAST, samtools, and other standard
    bioinformatics tools. It captures:
      - Command line used
      - All parameters and thresholds
      - Input file information
      - Summary statistics
      - Runtime and completion time
      - Output file locations

    The log file is saved alongside the main TSV output with a .log extension.

    Args:
        args:           Parsed command-line arguments.
        query_proteins: List of query proteins extracted from the input.
        ref_files:      List of reference GenBank file paths processed.
        all_hits:       List of homolog hits found.
        runtime_seconds: Elapsed time for the entire search (float).
        output_prefix:  Base name for output files (used to name the .log file).
    """
    # Construct log file path: same directory as TSV output, with .log extension
    if args.output:
        log_path = args.output.with_suffix(".log")
    else:
        log_path = Path(output_prefix).with_suffix(".log")

    # Reconstruct the command line for documentation
    import shlex

    cmd_parts = [
        "python3 pairwise_homolog_finder.py",
        f"-q {shlex.quote(str(args.query))}",
        f"-r {shlex.quote(str(args.reference))}",
        f"--identity {args.identity}",
        f"--min-coverage {args.min_coverage}",
        f"--coverage-mode {args.coverage_mode}",
    ]
    if args.mature:
        cmd_parts.append("--mature")
    if args.max_length:
        cmd_parts.append(f"--max-length {args.max_length}")
    if args.min_length != 10:
        cmd_parts.append(f"--min-length {args.min_length}")
    if args.output:
        cmd_parts.append(f"-o {shlex.quote(str(args.output))}")
    if args.output_fasta:
        cmd_parts.append(f"--output-fasta {shlex.quote(str(args.output_fasta))}")
    if args.verbose:
        cmd_parts.append("--verbose")

    command_line = " \\\n    ".join(cmd_parts)

    # Completion timestamp
    completion_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S UTC")

    # Gather reference file summaries
    ref_file_summaries = []
    for ref_file in ref_files:
        ref_file_summaries.append(f"  - {ref_file.name}")

    # Count hits per query
    hits_by_query = {}
    for hit in all_hits:
        if hit.query_locus not in hits_by_query:
            hits_by_query[hit.query_locus] = 0
        hits_by_query[hit.query_locus] += 1

    # Write the log file
    with open(log_path, "w", encoding="utf-8") as log:
        # Header
        log.write("=" * 100 + "\n")
        log.write(f"Homolog Search Results: {output_prefix}\n".center(100))
        log.write("=" * 100 + "\n\n")

        # Command line
        log.write("Command Line:\n")
        log.write(f"  {command_line}\n\n")

        # Parameters
        log.write("Parameters:\n")
        log.write(f"  Minimum Identity     : {args.identity * 100:.2f}%\n")
        log.write(f"  Minimum Coverage     : {args.min_coverage * 100:.2f}%\n")
        log.write(f"  Coverage Mode        : {args.coverage_mode}\n")
        log.write(
            f"  Mature Core Only     : {'YES (signal peptides trimmed)' if args.mature else 'NO'}\n"
        )
        log.write(f"  Min Query Length     : {args.min_length} aa\n")
        if args.max_length:
            log.write(f"  Max Query Length     : {args.max_length} aa\n")
        log.write("\n")

        # Query information
        log.write("Query Input:\n")
        log.write(f"  Query File           : {args.query.name}\n")
        log.write(f"  Proteins Extracted   : {len(query_proteins)}\n")
        if query_proteins:
            lengths = [p.length for p in query_proteins]
            log.write(f"  Query Length Range   : {min(lengths)}–{max(lengths)} aa\n")
        log.write("\n")

        # Reference information
        log.write("Reference Input:\n")
        log.write(f"  Reference Files      : {len(ref_files)} file(s)\n")
        for summary in ref_file_summaries:
            log.write(f"{summary}\n")
        log.write("\n")

        # Results summary
        log.write("Results Summary:\n")
        log.write(f"  Total Hits Found     : {len(all_hits)}\n")
        if hits_by_query:
            log.write(f"  Queries with Hits    : {len(hits_by_query)}\n")
            if query_proteins:
                pct_with_hits = len(hits_by_query) / len(query_proteins) * 100
                log.write(
                    f"  Coverage of Queries  : {pct_with_hits:.2f}% have at least one hit\n"
                )
        if all_hits:
            identities = [h.identity for h in all_hits]
            log.write(
                f"  Identity Range       : {min(identities)*100:.2f}–{max(identities)*100:.2f}%\n"
            )
            coverages = [h.coverage for h in all_hits]
            log.write(
                f"  Coverage Range       : {min(coverages)*100:.2f}–{max(coverages)*100:.2f}%\n"
            )
        log.write("\n")

        # Runtime information
        log.write("Execution:\n")
        log.write(f"  Runtime              : {runtime_seconds:.2f} seconds\n")
        log.write(f"  Completion Time      : {completion_time}\n")
        log.write("\n")

        # Output information
        log.write("Output Files:\n")
        if args.output:
            log.write(f"  TSV Results          : {args.output.resolve()}\n")
            log.write(
                "  TSV Columns          : query_locus, query_product, ref_locus,\n"
            )
            log.write(
                "                         ref_product, ref_file, query_sequence,\n"
            )
            log.write(
                "                         ref_sequence, mismatches, identity_pct,\n"
            )
            log.write(
                "                         alignment_length, query_length, ref_length\n"
            )
        if args.output_fasta:
            log.write(f"  FASTA Sequences      : {args.output_fasta.resolve()}\n")
            log.write("  FASTA Use Cases      : MAFFT, IQ-TREE, HMMER, InterProScan\n")
        log.write(f"  Log File             : {log_path.resolve()}\n")
        log.write("\n")

        # Footer
        log.write("=" * 100 + "\n")
        log.write(f"GBK Homolog Finder v{__version__}\n".center(100))
        log.write("=" * 100 + "\n")

    print(f"[*] Log saved to         : {log_path.resolve()}", file=sys.stderr)


def main() -> None:
    """Parses arguments and runs the full homolog-finding pipeline."""
    start_time = time.time()
    args = get_args()

    print("=" * 100, file=sys.stderr)
    print(f"GBK HOMOLOG FINDER v{__version__}", file=sys.stderr)
    print("=" * 100, file=sys.stderr)
    print(f"  Query        : {args.query}", file=sys.stderr)
    print(f"  Reference    : {args.reference}", file=sys.stderr)
    print(f"  Min identity : {args.identity * 100:.0f}%", file=sys.stderr)
    print(
        f"  Min coverage : {args.min_coverage * 100:.0f}% ({args.coverage_mode} sequence)",
        file=sys.stderr,
    )
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
            verbose=args.verbose,
        )
        if not query_proteins:
            sys.exit("[!] No proteins extracted from query. Check the file.")

        ref_files = list(stream_reference_files(args.reference))
        if not ref_files:
            sys.exit("[!] No valid reference files found.")

        all_hits: list[HomologHit] = []

        # Use progress bar if tqdm is available, otherwise just iterate
        ref_iter = tqdm(
            ref_files,
            desc="Scanning references",
            disable=not HAS_TQDM or args.verbose,
        )

        for ref_file in ref_iter:
            # Only print file name if verbose (tqdm shows progress in default mode)
            if args.verbose:
                print(
                    f"[*] Scanning {ref_file.name} ...",
                    file=sys.stderr,
                )
            file_hits = find_homologs(
                query_proteins=query_proteins,
                ref_path=ref_file,
                min_identity=args.identity,
                use_mature=args.mature,
                min_coverage=args.min_coverage,
                coverage_mode=args.coverage_mode,
            )
            # Only print results summary if verbose
            if args.verbose:
                print(
                    f"      \u2192 {len(file_hits)} hit(s) above "
                    f"{args.identity*100:.0f}% identity / "
                    f"{args.min_coverage*100:.0f}% coverage ({args.coverage_mode})",
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
            # Smart output: if there are many hits, print a preview instead of all
            if len(all_hits) > 50:
                print(
                    f"[!] Found {len(all_hits)} hits — too many to print to terminal.",
                    file=sys.stderr,
                )
                print(
                    "[*] Use: -o results.tsv to save results to a file (recommended).",
                    file=sys.stderr,
                )
                print(
                    "[*] Or use: --output-fasta results.faa for protein sequences.",
                    file=sys.stderr,
                )
                print(
                    "\n[*] Showing first 20 hits as preview:\n",
                    file=sys.stderr,
                )
                print_hits_table(all_hits[:20])
                print(
                    f"\n... ({len(all_hits) - 20} more hits omitted) ...\n",
                    file=sys.stderr,
                )
            else:
                print(
                    "[*] Tip: add -o results.tsv to save results with protein sequences.",
                    file=sys.stderr,
                )
                print()
                print_hits_table(all_hits)

        if args.output_fasta:
            print(
                f"[*] Saving FASTA \u2192 {args.output_fasta.resolve()}",
                file=sys.stderr,
            )
            write_fasta(all_hits, args.output_fasta)
            print(
                "      (For downstream: MAFFT, IQ-TREE, HMMER, InterProScan, etc.)",
                file=sys.stderr,
            )

        # Calculate runtime and write log file
        runtime_seconds = time.time() - start_time
        output_prefix = args.output.stem if args.output else args.query.stem
        write_log_file(
            args=args,
            query_proteins=query_proteins,
            ref_files=ref_files,
            all_hits=all_hits,
            runtime_seconds=runtime_seconds,
            output_prefix=output_prefix,
        )

        print("=" * 100, file=sys.stderr)
        print("[*] Done.", file=sys.stderr)
        print("=" * 100, file=sys.stderr)

    except ValueError as exc:
        sys.exit(f"\n[!] Error: {exc}")
    except KeyboardInterrupt:
        sys.exit("\n[!] Interrupted by user.")


if __name__ == "__main__":
    main()
