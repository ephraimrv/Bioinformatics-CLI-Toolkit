#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

r"""Remote BLAST Runner

Automates NCBI remote BLAST searches over sequences in a FASTA file.
Handles three sequence selection modes, writes TSV output with column
headers prepended automatically — no more typing the format string.

PROGRAMS:
    blastp  — protein query vs protein database (default)
    blastn  — nucleotide query vs nucleotide database
    blastx  — translated nucleotide query vs protein database

DEFAULT OUTPUT COLUMNS (outfmt 6 custom):
    qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore stitle

    Columns include:
    - Coordinates (qstart qend sstart send) for powerful pandas-based filtering
    - stitle: full subject description with organism in brackets
      Example: "WP_014324148.1 Blp family class II bacteriocin [Leuconostoc mesenteroides]"

CRITICAL NOTE — Taxonomy Fields (sscinames, scomnames) with -remote:
    Taxonomic columns (sscinames, scomnames) require a local Taxonomy Database
    (taxdb) to work. When using -remote, NCBI servers return only alignment
    data, not taxonomy lookups. Without local taxdb files, these columns return N/A.

    WORKAROUND: Extract organism name from stitle using regex in Python:

        import pandas as pd

        df = pd.read_csv('output.blast.tsv', sep='\t')
        df['organism'] = df['stitle'].str.extract(r'\[(.*?)\]', expand=False)
        print(df[['qseqid', 'sseqid', 'organism', 'evalue']])

    This extracts the organism name from the stitle brackets reliably.

SEQUENCE SELECTION MODES:
    (none)                         Blast every sequence in the file (default)
    --range LOCUS_START LOCUS_END  Blast a contiguous slice in file order
    --pick LOCUS_TAG [...]         Blast specific sequences by ID
    --list                         Preview sequence IDs and exit — does NOT
                                    run BLAST or touch NCBI. Without -o,
                                    shows the first 20 sequences only (a
                                    truncation notice covers the rest —
                                    important for eukaryotic-scale FASTA
                                    files with thousands of entries). With
                                    -o, writes the full list as a TSV
                                    instead of flooding the terminal with it.

MULTI-GENOME INPUT:
    -i accepts one or more FASTA files and/or directories, in any mix:
        -i genomeA.faa genomeB.faa genomeC.faa
        -i extracted_regions/
        -i genomeA.faa extracted_regions/

    A directory is expanded to every .fasta/.fa/.faa/.fna/.ffn file inside
    it (sorted alphabetically; non-FASTA files like a stray .txt are
    skipped automatically). This replaces manually running
    `cat file1 file2 file3 > combined.fasta` before this script.

    With more than one source file, every sequence ID is namespaced as
    "genome_label::original_id" (genome_label comes from the filename,
    e.g. "SCA72564" from SCA72564.faa). This is UNCONDITIONAL — it applies
    whether or not your specific IDs happen to collide today. The reason:
    Prokka/Bakta/SPAdes-style tools commonly restart numbering per genome
    (every genome can have its own "ctg1_1"), so a plain `cat` merge is one
    coincidence away from corrupting --pick/--range lookups and silently
    breaking --resume's per-ID progress tracking (genome A's ctg1_1
    finishing would wrongly mark genome B's ctg1_1 as done too, on a
    fully different sequence). Namespacing makes that structurally
    impossible rather than depending on you noticing it before it bites.
    A startup message reports exactly how many raw IDs actually collided,
    so you know how close you'd have come to that without it.

    With multiple sources, -o becomes REQUIRED — there's no single
    sensible default filename to guess once results are merged across
    genomes, so this is intentionally explicit rather than silently
    picking something like "combined.blast.tsv" on your behalf.

    With exactly one source file, behavior is identical to every prior
    version: no namespacing, no -o requirement, IDs untouched.

NCBI USAGE POLICY:
    Remote BLAST uses NCBI's public servers. A 10-second delay is enforced
    between searches by default (NCBI's stated guideline: no more than one
    contact every 10 seconds). By default, --batch-size bundles multiple
    sequences into each search — see BATCHING below — which also reduces
    how many separate searches you submit in a sitting. For very large
    inputs, consider downloading a local BLAST database instead:
    https://ftp.ncbi.nlm.nih.gov/blast/db/

    If results come back empty unexpectedly:
    1. NCBI rate limiting — you may have sent too many requests. Wait 30-60
       minutes before retrying.
    2. Timeout — refseq_protein and swissprot are slower than nr via remote.
       Increase timeout with --timeout (default: 600s, multiplied by
       --batch-size for the effective per-search ceiling).
    3. No hits — the sequence may genuinely have no hits in that database
       above your e-value threshold. Try a less stringent -e value.
    4. Use --debug to print the exact BLAST command and all server messages.

KNOWN LIMITATION — `-max_target_seqs` is NOT "top N best hits":
    This script passes --max-hits straight through to BLAST+'s
    -max_target_seqs, which is widely (and incorrectly) assumed to return
    the N highest-scoring database hits. Shah, Nute, Warnow & Pop (2019,
    Bioinformatics 35(9):1613-1614, doi:10.1093/bioinformatics/bty833)
    documented that the parameter instead limits how many hits are kept
    from the *search order* the algorithm encounters them in — not
    necessarily the N best by score. NCBI attributed part of the original
    severity to an over-aggressive internal optimization and patched it in
    BLAST+ 2.8.1 (2018); on current BLAST+ versions the worst outcome (the
    single best hit being silently dropped) is far less likely, but the
    parameter still does not guarantee a true top-N-by-score ranking.

    Practical takeaway for this toolkit: treat --max-hits as a ceiling on
    how many candidate hits to retrieve for downstream filtering, not as a
    statement that hit #1 in the output is necessarily the best possible
    match. If a manuscript's methods section reports "the top hit" from
    this script's output, that claim rests on the e-value/bitscore ranking
    of the RETURNED hits, not on an NCBI-side guarantee that no better hit
    was excluded before reaching you. For homology claims where this
    distinction matters, set --max-hits well above what you need (e.g. 250
    instead of 10) so the ranking step happens on your side, in pandas,
    after the fact.

RESUMING INTERRUPTED RUNS (--resume):
    A companion file <output>.progress is written alongside the TSV output,
    one line per sequence ID that has been DEFINITIVELY attempted (either it
    returned hits, or it returned zero hits cleanly — returncode 0 either
    way). Sequences that timed out or failed (non-zero returncode) are NOT
    marked as attempted, so a re-run will retry exactly those.

    With --batch-size > 1 (the default), the retry granularity is the
    BATCH, not the individual sequence: if a batch's search fails or times
    out, every sequence in that batch is left unmarked and all of them
    retry together on the next --resume run, even if some of them would
    have succeeded individually. Use --batch-size 1 for fully per-sequence
    retry granularity at the cost of one NCBI search per sequence.

    --resume reads that file, skips every sequence ID already in it, and
    APPENDS new results to the existing output TSV instead of overwriting
    it. Without --resume, output is always overwritten fresh, same as
    before.

    Only supported with the default outfmt 6 (tabular) — other formats
    don't have a parseable column to confirm a hit-free success vs. a
    silent failure, so --resume is ignored (with a warning) if -f is not 6.

    Safe to combine with --range/--pick/--list: resume only affects which
    of the SELECTED sequences are skipped, not which ones are selected.

BATCHING (--batch-size):
    By default, --batch-size sequences are bundled into ONE NCBI search
    instead of one search per sequence. This is NOT a risky workaround —
    NCBI's own developer documentation states plainly that "BLAST often
    runs more efficiently if multiple queries are sent as one search
    rather than if each query is sent as an individual search." Batching
    also directly reduces how many separate "searches" you submit in a
    sitting, which is the actual metric behind NCBI's stated 100-
    searches-per-24-hours throttling threshold (see NCBI USAGE POLICY
    above) — 200 sequences at --batch-size 10 is 20 searches, not 200.

    How it works: all sequences in a batch are written to one multi-FASTA
    query file and submitted as a single BLAST+ -remote invocation.
    outfmt 6's qseqid column already tags every hit row with which
    sequence it belongs to, so a multi-query submission's results need no
    extra parsing — they fall out of the existing TSV columns for free.

    --timeout is automatically multiplied by --batch-size for the
    effective per-search timeout (a 10-sequence batch needs more
    server-side processing time than one sequence). This is a ceiling,
    not an estimate of typical wait — most batches finish well under it.

    NCBI also limits total query length PER SEARCH (not per sequence):
    100,000 residues for protein queries (blastp), 1,000,000 bp for
    nucleotide queries (blastn/blastx). An oversized batch triggers a
    non-blocking warning before submission rather than a silent failure
    several minutes later. If you hit this, lower --batch-size.

    Use --batch-size 1 to restore the original one-sequence-per-search
    behavior exactly (e.g. if you need fully granular --resume retries,
    or are diagnosing a single problematic sequence with --debug).

    The output TSV format is completely unaffected by --batch-size — same
    columns, same content, regardless of how many sequences shared a
    submission. Only the request grouping, terminal logging granularity,
    and --resume retry granularity change.

SEQUENCE-TYPE CHECK:
    Before starting, a sample of the input is checked against the chosen
    -p program (protein vs nucleotide) using a cheap heuristic: amino acids
    E, F, I, L, P, Q, Z, J never appear in IUPAC nucleotide codes, so their
    presence/absence is a fast (not perfect) signal of which alphabet the
    file actually contains. This is advisory only — it prints a warning if
    the input looks mismatched (e.g. a .fna file run with -p blastp) but
    does not block the run, since short sequences can occasionally trip a
    false positive.

Reproducibility:
    Associated with upcoming research (manuscript in preparation).
    Correct attribution is requested when used in derivative works.
    See LICENSE in the repository root for full details.

Example Usage:
    # Blast all sequences with defaults (blastp, nr, e-value 1e-10)
    $ python3 remote_blast.py -i proteins.faa

    # Blast a range of sequences with less stringent e-value
    $ python3 remote_blast.py -i proteins.faa --range ctg1_1 ctg1_20 -e 1e-5

    # Blast three specific sequences
    $ python3 remote_blast.py -i proteins.faa --pick ctg1_1 ctg1_36 ctg1_45 -o picked.tsv

    # Nucleotide blast against nt database
    $ python3 remote_blast.py -i genome.fna -p blastn

    # Translated blast, custom max hits
    $ python3 remote_blast.py -i contigs.fna -p blastx --max-hits 50

    # Standard nr (equivalent to web BLAST "Clustered nr" option)
    $ python3 remote_blast.py -i proteins.faa -e 1e-10

    # RefSeq only (curated, fewer sequences, faster)
    $ python3 remote_blast.py -i proteins.faa --db refseq_protein -e 1e-10

    # Swiss-Prot (manually reviewed, highest quality)
    $ python3 remote_blast.py -i proteins.faa --db swissprot -e 1e-10

    # Protein Data Bank
    $ python3 remote_blast.py -i proteins.faa --db pdb -e 1e-10

    # Debug mode — prints exact command, NCBI messages, raw output length
    $ python3 remote_blast.py -i proteins.faa --debug

    # Increase timeout for slow databases (default: 600s)
    $ python3 remote_blast.py -i proteins.faa --db refseq_protein --timeout 900

    # Resume an interrupted run — skips sequences already completed,
    # appends new results to the existing output instead of overwriting it
    $ python3 remote_blast.py -i proteins.faa --resume

    # Preview sequence IDs without running BLAST. Use -o on a large
    # (e.g. eukaryotic) FASTA to save the full list instead of truncating
    # the terminal preview to the first 20.
    $ python3 remote_blast.py -i eukaryote_proteome.faa --list
    $ python3 remote_blast.py -i eukaryote_proteome.faa --list -o seq_list.tsv

    # 200 sequences, default batching (10/search = 20 searches, not 200)
    $ python3 remote_blast.py -i many_proteins.faa --resume

    # Larger batches for an even faster, NCBI-recommended bulk run
    $ python3 remote_blast.py -i many_proteins.faa --batch-size 50 --resume

    # Restore the original one-sequence-per-search behavior exactly
    $ python3 remote_blast.py -i proteins.faa --batch-size 1

    # Multiple genomes' .faa files in one BLAST run — no more manual
    # `cat file1 file2 file3 > combined.fasta` first. IDs are automatically
    # namespaced as "genome_label::id" to prevent cross-genome collisions.
    $ python3 remote_blast.py -i genomeA.faa genomeB.faa genomeC.faa \
        -o combined_results.tsv --resume

    # Or point at a whole directory of extracted regions/genomes
    $ python3 remote_blast.py -i extracted_regions/ -o combined_results.tsv

v1.4.2 changes (bugfix): A sequence whose remote search was silently
    rejected by NCBI server-side (e.g. CPU usage limit, overload,
    maintenance) used to be indistinguishable from a genuine zero-hit
    result — both produce empty stdout with returncode 0. That meant a
    server-side failure got permanently recorded as "done" in .progress
    and silently skipped on every future --resume, losing the sequence's
    data with no further warning. _blast_one() now checks stderr for known
    NCBI failure signatures whenever stdout is empty, and treats a match as
    a retryable failure (returns None) rather than a completed result.

v1.4.3 changes: (1) NCBI_MIN_DELAY corrected from 5s to 10s to match
    NCBI's current, actual BLAST URL API guideline ("do not contact the
    server more often than once every 10 seconds") — the old 5s figure was
    based on a rate limit that belongs to a different NCBI service
    (E-utilities), not the one -remote uses; (2) a run of more than 50
    sequences now prints an advisory citing NCBI's stated 100-searches-
    per-24-hours throttling/blocking threshold and their own
    recommendation to run large batches off-peak. Neither change alters
    output format or blocks a run — both are corrections/warnings only.
    NOTE: this script still submits one sequence per NCBI search. For
    100-200+ sequence runs, batching multiple sequences into one search is
    the real fix still pending (NCBI's own guidance explicitly recommends
    this over many single-query searches) — not yet implemented as of
    this version.

v1.5.0 changes (the batching fix promised in v1.4.3's note above):
    _blast_one() generalized to _blast_batch(), which submits --batch-size
    sequences (default 10) as ONE multi-FASTA NCBI search instead of one
    search per sequence. This directly addresses both the wall-clock pain
    of 100-200+ sequence runs AND the actual NCBI policy risk of submitting
    that many separate searches in one sitting (see BATCHING and NCBI
    USAGE POLICY above). --timeout now scales with batch size; an
    oversized batch's total query length is checked against NCBI's stated
    per-search limits before submission. --resume's retry granularity is
    now the batch, not the individual sequence, when --batch-size > 1 —
    documented explicitly since it's a real, intentional trade-off, not an
    oversight. --batch-size 1 reproduces the exact pre-v1.5.0 behavior.
    The output TSV format itself is unchanged regardless of batch size.

v1.6.0 changes: -i now accepts multiple files and/or directories instead
    of exactly one file, replacing the previous manual `cat file1 file2
    > combined.fasta` workflow. With more than one resolved source file,
    every sequence ID is unconditionally namespaced as
    "genome_label::original_id" to prevent silent cross-genome ID
    collisions (common with Prokka/Bakta/SPAdes-style per-genome
    numbering) from corrupting --pick/--range lookups or --resume's
    per-ID progress tracking — see MULTI-GENOME INPUT above. -o is now
    required when multiple sources are given. Single-file input behavior
    is completely unchanged.

v1.6.1 changes: Added the missing shebang/SPDX-license/copyright header
    that every other script in this toolkit carries at the very top of
    the file — this one had skipped straight to the module docstring.
    Also removed the now-redundant inline "License: MIT" line from the
    docstring body, since the SPDX header covers it. Documentation only;
    no behavior change.
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.6.1"

import argparse
import contextlib
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )


# ── Constants ─────────────────────────────────────────────────────────────────

# BLAST output columns. NOTE: taxonomic fields (sscinames, scomnames) require
# a local Taxonomy Database (taxdb) to populate. With -remote, NCBI servers
# return alignment data only, not taxonomy lookups, so those fields return N/A.
# Use stitle instead and extract organism names via regex in pandas.
COLUMNS = (
    "qseqid sseqid pident qcovs length mismatch gapopen "
    "qstart qend sstart send evalue bitscore stitle"
)
HEADER_ROW = COLUMNS.replace(" ", "\t")

DEFAULT_DB: dict[str, str] = {
    "blastp": "nr",
    "blastn": "nt",
    "blastx": "nr",
}

NCBI_MIN_DELAY = 10  # seconds between requests. NCBI's current developer
# guidelines (blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html)
# state plainly: "Do not contact the server more often than once every 10
# seconds." (Earlier versions of this script used 5s, based on a "<=3
# requests/sec" figure that belongs to NCBI's E-utilities service, a
# different API with a different, looser limit — not the BLAST URL API
# that -remote actually uses.)
DEFAULT_TIMEOUT = 600  # seconds — 10 minutes (remote BLAST can be slow)

DEFAULT_BATCH_SIZE = 10  # sequences bundled into one NCBI search by default

# NCBI's stated per-search total query length limits (combined across every
# sequence in one submission, not per individual sequence). Exceeding these
# means, per NCBI's own support documentation, "the search will simply not
# run" — checked before submitting a batch so an oversized --batch-size
# fails loud and early instead of silently wasting a submission.
MAX_QUERY_LENGTH: dict[str, int] = {
    "blastp": 100_000,  # protein queries (blastp, tblastn)
    "blastn": 1_000_000,  # nucleotide queries
    "blastx": 1_000_000,  # blastx query is nucleotide, translated server-side
}

# Programs whose QUERY is nucleotide. blastp's query is protein; blastn's and
# blastx's query is nucleotide (blastx translates it on the fly before
# comparing against a protein database). Used by the sequence-type check.
EXPECTS_NUCLEOTIDE_QUERY: dict[str, bool] = {
    "blastp": False,
    "blastn": True,
    "blastx": True,
}

# Amino-acid letters that are never valid IUPAC nucleotide codes (which only
# use A C G T U N R Y S W K M B D H V). Their presence is a strong signal
# that a sequence is protein, not nucleotide.
_PROTEIN_SIGNATURE_CHARS = frozenset("EFILPQZJ")
_NUCLEOTIDE_CHARS = frozenset("ACGTUN")
_SEQUENCE_TYPE_SAMPLE_SIZE = 20  # cheap check — no need to scan the whole file


def _looks_like_nucleotide(seq: str) -> bool:
    """Heuristic guess at whether a sequence is nucleotide or protein.

    Not a real alphabet parser — just a cheap, fast filter to catch the
    "ran blastp on a .fna file" class of mistake before it burns NCBI
    rate-limit allowance and a --timeout wait on garbage input.

    Args:
        seq: The sequence string to check.

    Returns:
        True if the sequence looks like nucleotide data, False if it
        looks like protein (or is too ambiguous to call confidently —
        defaults to False, i.e. "not clearly nucleotide", so the
        nucleotide-vs-protein vote stays conservative).
    """
    seq_upper = seq.upper()
    if not seq_upper:
        return False
    if any(c in _PROTEIN_SIGNATURE_CHARS for c in seq_upper):
        return False
    alpha_chars = [c for c in seq_upper if c.isalpha()]
    if not alpha_chars:
        return False
    non_nucleotide = sum(1 for c in alpha_chars if c not in _NUCLEOTIDE_CHARS)
    return (non_nucleotide / len(alpha_chars)) < 0.05


def _check_sequence_type(records: list[SeqRecord], program: str) -> None:
    """Prints a non-blocking warning if the input alphabet looks mismatched
    with the chosen BLAST program (e.g. protein FASTA run with -p blastn).

    Advisory only — never aborts the run. The heuristic in
    _looks_like_nucleotide() can be wrong on short or unusual sequences,
    so this is insurance against an obvious mix-up, not a hard gate.

    Args:
        records: The full list of parsed input sequences (sampled, not
            all scanned, to keep this cheap).
        program: The chosen BLAST program (blastp/blastn/blastx).
    """
    expects_nucleotide = EXPECTS_NUCLEOTIDE_QUERY[program]
    sample = records[:_SEQUENCE_TYPE_SAMPLE_SIZE]
    nuc_votes = sum(1 for r in sample if _looks_like_nucleotide(str(r.seq)))
    looks_nucleotide = (nuc_votes / len(sample)) > 0.5 if sample else False

    if expects_nucleotide and not looks_nucleotide:
        print(
            f"\n[!] WARNING: -p {program} expects a NUCLEOTIDE query, but the "
            f"sampled input looks like PROTEIN data.\n"
            f"    If this is a .faa file, you likely want -p blastp instead.\n"
            f"    Continuing anyway — this is a heuristic, not a hard stop.",
            file=sys.stderr,
        )
    elif not expects_nucleotide and looks_nucleotide:
        print(
            f"\n[!] WARNING: -p {program} expects a PROTEIN query, but the "
            f"sampled input looks like NUCLEOTIDE data.\n"
            f"    If this is a .fna/.fasta file, you likely want -p blastn "
            f"or -p blastx instead.\n"
            f"    Continuing anyway — this is a heuristic, not a hard stop.",
            file=sys.stderr,
        )


# ── Multi-source input handling ────────────────────────────────────────────────

# Extensions treated as FASTA when expanding a directory passed to -i.
# Covers protein (.faa), nucleotide contigs (.fna), CDS nucleotide (.ffn),
# and generic (.fasta/.fa) outputs — the formats Prokka/Bakta-style
# annotation pipelines actually produce.
FASTA_EXTENSIONS = frozenset({".fasta", ".fa", ".faa", ".fna", ".ffn"})

# Separator between a genome label and its original sequence ID when
# multiple sources are merged (e.g. "SCA72564::ctg1_49"). Deliberately NOT
# "|" — NCBI's classic FASTA ID format uses "|" as a field separator
# (gi|12345|ref|XYZ), so a literal pipe in a custom ID risks being
# misinterpreted by -remote's ID parsing. "::" is a common bioinformatics
# namespace separator (seqkit, samtools region syntax) and never appears in
# ordinary locus tags/contig names, so splitting back out is unambiguous.
GENOME_LABEL_SEPARATOR = "::"


def _resolve_input_files(input_paths: list[Path]) -> list[Path]:
    """Expands a mix of files and directories into a flat, deduplicated,
    deterministically-ordered list of FASTA files.

    Directories are expanded to every file inside them matching
    FASTA_EXTENSIONS (case-insensitive), sorted alphabetically for
    reproducibility. Individual files are used as given. The relative
    order of multiple -i arguments is preserved; only within-directory
    contents get sorted. Resolves to absolute paths for deduplication, so
    the same file reached two different ways (e.g. listed directly AND
    via a directory containing it) is only included once.

    Args:
        input_paths: Raw paths from -i, in the order given on the CLI —
            each may be a file or a directory.

    Returns:
        Flat list of FASTA file paths, no duplicates.
    """
    resolved: list[Path] = []
    seen: set[Path] = set()

    for p in input_paths:
        if not p.exists():
            sys.exit(f"[!] Input path not found: {p}")

        if p.is_dir():
            matches = sorted(
                f
                for f in p.iterdir()
                if f.is_file() and f.suffix.lower() in FASTA_EXTENSIONS
            )
            if not matches:
                sys.exit(
                    f"[!] No FASTA files found in directory: {p}\n"
                    f"    Looked for extensions: {sorted(FASTA_EXTENSIONS)}"
                )
            candidates = matches
        else:
            candidates = [p]

        for f in candidates:
            real = f.resolve()
            if real not in seen:
                seen.add(real)
                resolved.append(f)

    return resolved


def _load_multi_source(files: list[Path]) -> list[SeqRecord]:
    """Loads and merges sequences from one or more FASTA files.

    With a single file, behavior is identical to plain SeqIO.parse — IDs
    are untouched, fully backward compatible with every single-genome
    workflow used throughout this toolkit.

    With multiple files, every record's ID is rewritten to
    "<genome_label>::<original_id>" UNCONDITIONALLY — not only when a
    collision is detected. This is deliberate: behavior that only changes
    when IDs happen to collide today is exactly the kind of "works until
    it doesn't" inconsistency this toolkit has been burned by before
    (e.g. find_gbk_features.py's isoform-anchor bug). Always namespacing
    means --pick/--range/--resume/output qseqid are unambiguous regardless
    of what any particular run's IDs happen to look like.

    A collision count IS still reported (even though namespacing always
    applies) so you know how much silent corruption a plain `cat` would
    have caused on this specific input — useful context, not a
    conditional safety net.

    Args:
        files: FASTA files to load, in the order they should be merged.

    Returns:
        Merged list of SeqRecord objects, IDs namespaced if len(files) > 1.
    """
    if len(files) == 1:
        return list(SeqIO.parse(files[0], "fasta"))

    # Build genome labels from filenames, disambiguating any label clashes
    # (e.g. two different directories each containing a "genomeA.faa").
    used_labels: dict[str, int] = {}
    merged: list[SeqRecord] = []
    id_to_labels: dict[str, set[str]] = {}

    for f in files:
        base_label = f.stem
        used_labels[base_label] = used_labels.get(base_label, 0) + 1
        label = (
            base_label
            if used_labels[base_label] == 1
            else f"{base_label}_{used_labels[base_label]}"
        )

        file_records = list(SeqIO.parse(f, "fasta"))
        if not file_records:
            print(
                f"[!] WARNING: no sequences found in {f} — skipping.", file=sys.stderr
            )
            continue

        for record in file_records:
            id_to_labels.setdefault(record.id, set()).add(label)
            original_id = record.id
            # Preserve the human-readable part of the description (e.g.
            # Prokka/Bakta product annotations) — only the ID token at the
            # front needs namespacing, not the whole header. Losing this
            # would make --list noticeably less useful for picking targets
            # across genomes, which defeats the point of merging them.
            extra_desc = record.description[len(original_id) :].strip()
            new_id = f"{label}{GENOME_LABEL_SEPARATOR}{original_id}"
            record.id = new_id
            record.description = f"{new_id} {extra_desc}" if extra_desc else new_id
            merged.append(record)

    n_collisions = sum(1 for labels in id_to_labels.values() if len(labels) > 1)
    print(
        f"[*] Merged {len(files)} source file(s) into {len(merged):,} "
        f"sequence(s), namespaced as 'genome_label{GENOME_LABEL_SEPARATOR}id'.",
        file=sys.stderr,
    )
    if n_collisions:
        print(
            f"[*] {n_collisions} raw ID(s) appeared in more than one source "
            f"file — without namespacing, a plain `cat` merge would have "
            f"silently corrupted --pick/--range lookups and --resume "
            f"tracking for all {n_collisions} of them. Namespacing avoided that.",
            file=sys.stderr,
        )

    return merged


# ── Argument parsing ───────────────────────────────────────────────────────────


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="remote_blast.py",
        description=(
            "Remote NCBI BLAST with automatic TSV headers and flexible "
            "sequence selection."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
sequence selection (mutually exclusive — pick one or use none for all):
  --range LOCUS_START LOCUS_END  e.g. --range ctg1_1 ctg1_20
  --pick  LOCUS_TAG [...]        e.g. --pick ctg1_1 ctg1_36 ctg1_45

examples:
  python3 remote_blast.py -i proteins.faa
  python3 remote_blast.py -i proteins.faa -p blastp --range ctg1_1 ctg1_20 -e 1e-5
  python3 remote_blast.py -i proteins.faa -p blastp --pick ctg1_1 ctg1_36 ctg1_45
  python3 remote_blast.py -i genome.fna -p blastn
  python3 remote_blast.py -i contigs.fna -p blastx --db nr --max-hits 100
  python3 remote_blast.py -i proteins.faa --debug
        """,
    )

    # ── Required ──────────────────────────────────────────────────────────────
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        nargs="+",
        type=Path,
        metavar="FASTA",
        help=(
            "Input FASTA file(s) and/or directory/directories "
            "(.faa for protein, .fna/.fasta for nucleotide). "
            "Accepts multiple paths: files, directories, or a mix — e.g. "
            "-i genomeA.faa genomeB.faa, or -i extracted_regions/. "
            "A directory is expanded to every .fasta/.fa/.faa/.fna/.ffn "
            "file inside it. With more than one source file, sequence IDs "
            "are automatically namespaced as 'genome_label::id' to prevent "
            "cross-genome ID collisions (e.g. two genomes both having a "
            "'ctg1_1') from silently corrupting --pick/--resume."
        ),
    )

    # ── Optional ──────────────────────────────────────────────────────────────
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        metavar="TSV",
        help=(
            "Output TSV file. "
            "If not specified, defaults to <input_filename>.blast.tsv "
            "in the same directory as the input."
        ),
    )
    parser.add_argument(
        "-p",
        "--program",
        choices=["blastp", "blastn", "blastx"],
        default="blastp",
        help="BLAST program. Default: blastp",
    )
    parser.add_argument(
        "-e",
        "--evalue",
        type=float,
        default=1e-10,
        metavar="FLOAT",
        help="E-value threshold. Default: 1e-10",
    )
    parser.add_argument(
        "-f",
        "--outfmt",
        type=int,
        default=6,
        metavar="INT",
        help=(
            "BLAST output format number. Default: 6 (tabular with custom columns). "
            "If a non-6 format is selected, column headers are not written."
        ),
    )
    parser.add_argument(
        "--db",
        type=str,
        default=None,
        metavar="DB",
        help=(
            "BLAST database name (NCBI remote). "
            "Defaults: nr (blastp/blastx), nt (blastn). "
            "Other options: refseq_protein, swissprot, pdb, refseq_rna, etc. "
            "Note: 'Clustered nr' in the web BLAST UI is a presentation option on the "
            "standard 'nr' database — use --db nr for the equivalent."
        ),
    )
    parser.add_argument(
        "--max-hits",
        type=int,
        default=100,
        metavar="N",
        help="Maximum number of hits to return per sequence. Default: 100",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=DEFAULT_BATCH_SIZE,
        metavar="N",
        help=(
            f"Number of sequences bundled into ONE NCBI search. "
            f"Default: {DEFAULT_BATCH_SIZE}. This is NCBI's own recommended "
            f"approach for many queries (fewer, larger searches over many "
            f"small ones) — see the module docstring's BATCHING section. "
            f"--timeout is automatically scaled by batch size. Use "
            f"--batch-size 1 to restore the original one-sequence-per-search "
            f"behavior. With --resume, a failed/timed-out batch retries as "
            f"a whole batch, not per individual sequence."
        ),
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=DEFAULT_TIMEOUT,
        metavar="SEC",
        help=(
            f"Seconds to wait per sequence-equivalent before giving up on a "
            f"search. Default: {DEFAULT_TIMEOUT}s. Automatically multiplied "
            f"by --batch-size for the actual per-search timeout (e.g. "
            f"default timeout x default batch-size = "
            f"{DEFAULT_TIMEOUT * DEFAULT_BATCH_SIZE}s ceiling per search). "
            f"Increase for slow databases (refseq_protein, swissprot)."
        ),
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=NCBI_MIN_DELAY,
        metavar="SEC",
        help=(
            f"Seconds to wait between searches (i.e. between BATCHES, not "
            f"individual sequences within a batch). "
            f"Default: {NCBI_MIN_DELAY} (NCBI usage policy)."
        ),
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help=(
            "Print the exact BLAST command, all NCBI server messages (stderr), "
            "and raw output length for each sequence. Use when results are "
            "unexpectedly empty to diagnose rate limiting or server errors."
        ),
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help=(
            "Preview sequence IDs and exit immediately — no BLAST, no NCBI "
            "contact. Without -o, shows only the first 20 (truncated for "
            "large files); with -o, writes the full list as a TSV instead "
            "of printing it all. Use the ID with --range or --pick."
        ),
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help=(
            "Skip sequences already completed in a prior run (tracked in "
            "<output>.progress) and append new results to the existing "
            "output TSV instead of overwriting it. Only sequences that "
            "timed out or failed are retried. Requires outfmt 6 (ignored "
            "with a warning otherwise)."
        ),
    )

    # ── Sequence selection (mutually exclusive) ───────────────────────────────
    sel = parser.add_mutually_exclusive_group()
    sel.add_argument(
        "--range",
        nargs=2,
        type=str,
        metavar=("LOCUS_START", "LOCUS_END"),
        help=(
            "Blast all sequences from LOCUS_START to LOCUS_END (inclusive), "
            "in the order they appear in the file. "
            "Use --list to see available locus tags. "
            "Example: --range ctg1_1 ctg1_80"
        ),
    )
    sel.add_argument(
        "--pick",
        nargs="+",
        type=str,
        metavar="LOCUS_TAG",
        help=(
            "Blast specific sequences by locus tag. "
            "Use --list to see available locus tags. "
            "Example: --pick ctg1_1 ctg1_36 ctg1_45"
        ),
    )

    return parser


# ── Sequence selection ────────────────────────────────────────────────────────

# Maximum sequence IDs shown on the terminal for --list without -o. Beyond
# this, a truncation notice points to -o instead — important for eukaryotic
# proteomes/transcriptomes, which can run into the tens of thousands of
# sequences. Mirrors the --list-sequences convention used elsewhere in this
# toolkit (find_gbk_features.py, extract_genome_region.py).
_LIST_MAX_DISPLAY = 20


def _run_list(
    records: list[SeqRecord], args: argparse.Namespace, input_files: list[Path]
) -> None:
    """Preview sequence IDs in the input FASTA, or save the full list to a TSV.

    Without -o: prints the first _LIST_MAX_DISPLAY sequences to the
    terminal for quick scouting, then a truncation notice — never the
    entire file, which could be tens of thousands of entries for a
    eukaryotic genome.

    With -o: the TSV file IS the output. The terminal only gets a one-line
    confirmation, not a second full copy of the list — same "file is the
    output, terminal isn't a second output" rule used by
    find_gbk_features.py's run_list_sequences().

    This never touches BLAST or NCBI — it's pure local FASTA inspection,
    so no '[*] Output: <path>.blast.tsv' line should ever appear here; that
    line implies a BLAST run happened, which --list never does.

    Args:
        records: All sequences parsed/merged from the input source(s).
        args: Parsed argument namespace (uses args.output).
        input_files: Resolved source file list, for an accurate header
            when multiple files/a directory were given to -i.
    """
    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            f.write("seq_num\tsequence_id\tlength\tdescription\n")
            for i, rec in enumerate(records, 1):
                desc = rec.description[len(rec.id) :].strip()
                f.write(f"{i}\t{rec.id}\t{len(rec.seq)}\t{desc}\n")

        print(
            f"[*] {len(records):,} sequence(s) written to '{args.output}'.\n"
            f"    Use the ID (sequence_id column) with --range or --pick.",
            file=sys.stderr,
        )
        return

    source_desc = (
        input_files[0].name
        if len(input_files) == 1
        else f"{len(input_files)} source files"
    )
    print(f"\nSequences in {source_desc} ({len(records):,} total):\n")
    for i, rec in enumerate(records[:_LIST_MAX_DISPLAY], 1):
        desc = rec.description[len(rec.id) :].strip()
        desc_col = f"  {desc[:60]}" if desc else ""
        print(f"  [{i:>4}]  {rec.id:<20}{desc_col}")

    if len(records) > _LIST_MAX_DISPLAY:
        hidden = len(records) - _LIST_MAX_DISPLAY
        print(f"  ... and {hidden:,} more. Use -o to save the full list to a TSV.")

    print(
        f"\n  Use the ID (second column) with --range or --pick.\n"
        f"  Example: --range {records[0].id} {records[-1].id}"
    )


def _select_sequences(
    records: list[SeqRecord], args: argparse.Namespace
) -> list[tuple[int, SeqRecord]]:
    """Return (1-indexed position, record) pairs based on the selection mode."""
    id_to_idx: dict[str, int] = {rec.id: i for i, rec in enumerate(records)}

    if args.range:
        start_tag, end_tag = args.range
        missing = [t for t in (start_tag, end_tag) if t not in id_to_idx]
        if missing:
            sys.exit(
                f"[!] Locus tag(s) not found in FASTA: {missing}\n"
                f"    Run with --list to see all available IDs."
            )
        i = id_to_idx[start_tag]
        j = id_to_idx[end_tag]
        if i > j:
            sys.exit(
                f"[!] '{start_tag}' (position {i + 1}) comes AFTER "
                f"'{end_tag}' (position {j + 1}) in the file.\n"
                f"    Swap the arguments or check --list for the correct order."
            )
        return [(idx + 1, records[idx]) for idx in range(i, j + 1)]

    if args.pick:
        missing = [t for t in args.pick if t not in id_to_idx]
        if missing:
            sys.exit(
                f"[!] Locus tag(s) not found in FASTA: {missing}\n"
                f"    Run with --list to see all available IDs."
            )
        seen: set[str] = set()
        selected = []
        for tag in args.pick:
            if tag not in seen:
                idx = id_to_idx[tag]
                selected.append((idx + 1, records[idx]))
                seen.add(tag)
        return selected

    return [(i + 1, rec) for i, rec in enumerate(records)]


# ── BLAST execution ───────────────────────────────────────────────────────────


def _batch_label(records: list[SeqRecord]) -> str:
    """Short human-readable label for a batch, used in log/error messages.

    Args:
        records: The records in this batch (at least one).

    Returns:
        The single ID if there's one record, or "first..last (N seqs)" for
        a multi-record batch.
    """
    if len(records) == 1:
        return records[0].id
    return f"{records[0].id}..{records[-1].id} ({len(records)} seqs)"


def _check_batch_query_length(records: list[SeqRecord], program: str) -> None:
    """Warns (non-blocking) if a batch's combined query length exceeds
    NCBI's stated per-search limit for this program.

    Per NCBI's own support documentation, exceeding the limit means the
    search simply won't run server-side — this is advisory so the failure
    is explained immediately rather than discovered as a mystery empty
    result several minutes later.

    Args:
        records: The records about to be submitted as one batch.
        program: The chosen BLAST program (blastp/blastn/blastx).
    """
    limit = MAX_QUERY_LENGTH.get(program)
    if limit is None:
        return
    total_len = sum(len(r.seq) for r in records)
    if total_len > limit:
        print(
            f"\n[!] WARNING: This batch's combined query length is "
            f"{total_len:,} (residues/bp for {program}), exceeding NCBI's "
            f"stated limit of {limit:,} per search. Per NCBI's own "
            f"documentation, an oversized search simply will not run. "
            f"Consider lowering --batch-size.",
            file=sys.stderr,
        )


def _blast_batch(
    records: list[SeqRecord],
    program: str,
    db: str,
    evalue: float,
    outfmt: int,
    max_hits: int,
    timeout: int,
    debug: bool,
) -> str | None:
    """Runs one NCBI remote BLAST search over a batch of sequences and
    returns the raw output string.

    A "batch" of one sequence is the original, fully backward-compatible
    behavior (--batch-size 1). With more than one record, all of them are
    written to a single multi-FASTA query file and submitted as ONE NCBI
    search — this is NCBI's own recommended approach for many queries (see
    the module docstring's BATCHING section), not a workaround. BLAST+'s
    outfmt 6 already tags every hit row with `qseqid`, so a multi-query
    submission's results are inherently disambiguated per sequence with no
    extra parsing needed here.

    Writes the sequence(s) to a temporary file and passes it to BLAST via
    -query. The stdin approach (-query -) is documented but unreliable
    with -remote: BLAST+ does not consistently pipe stdin through to
    NCBI's servers across all versions. Temp file is the safe approach.

    Args:
        records:  The SeqRecord(s) to query, as one batch/submission.
        program:  BLAST program name (blastp/blastn/blastx).
        db:       Database name.
        evalue:   E-value cutoff.
        outfmt:   Output format number.
        max_hits: Maximum hits to return PER SEQUENCE (-max_target_seqs is
                  a per-query parameter in BLAST+, unaffected by batching).
        timeout:  Seconds before giving up on this request. Callers are
                  expected to scale this with batch size themselves (a
                  10-sequence batch needs more server-side time than one
                  sequence) — this function just enforces whatever value
                  it's given.
        debug:    If True, print command and all NCBI messages to stderr.

    Returns:
        On a completed request (returncode 0) with a genuine result: the
        BLAST output string, which may be "" if every sequence in the
        batch truly found zero hits — this counts as DEFINITIVELY DONE
        for --resume purposes, for EVERY sequence in the batch.
        On timeout, a non-zero returncode, OR a detected silent server-side
        failure (empty stdout + a known failure signature on stderr, e.g.
        NCBI's CPU usage limit rejection): None — this counts as NOT done
        for the WHOLE batch, so --resume will retry every sequence in it
        on the next run (batch is the retry granularity, not individual
        sequence, when batch-size > 1). Callers must check `is None`, not
        truthiness, to tell these apart (an empty string and None are
        different outcomes here).
    """
    label = _batch_label(records)
    tmp_file = tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False)
    tmp_path = tmp_file.name

    try:
        SeqIO.write(records, tmp_file, "fasta")
        tmp_file.close()  # must close before BLAST opens it (critical on WSL2)

        fmt_arg = f"6 {COLUMNS}" if outfmt == 6 else str(outfmt)
        cmd = [
            program,
            "-query",
            tmp_path,
            "-db",
            db,
            "-remote",
            "-outfmt",
            fmt_arg,
            "-evalue",
            str(evalue),
            "-max_target_seqs",
            str(max_hits),
        ]

        if debug:
            print(
                f"\n    [DEBUG] Command: {' '.join(cmd)}",
                file=sys.stderr,
            )

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
        )

        # Always print NCBI's stderr messages in debug mode.
        # In normal mode, print only if non-empty — NCBI uses stderr for
        # rate limiting warnings and other server messages that explain
        # why results may be empty even with returncode 0.
        if result.stderr.strip():
            if debug:
                print(
                    f"    [DEBUG] NCBI messages:\n"
                    + "\n".join(
                        f"            {line}"
                        for line in result.stderr.strip().splitlines()
                    ),
                    file=sys.stderr,
                )
            else:
                # Always show server messages — they explain empty results
                print(
                    f"\n    [!] NCBI server message for {label}:\n"
                    + "\n".join(
                        f"        {line}" for line in result.stderr.strip().splitlines()
                    ),
                    file=sys.stderr,
                )

        if debug:
            print(
                f"    [DEBUG] Return code: {result.returncode}",
                file=sys.stderr,
            )
            print(
                f"    [DEBUG] stdout length: {len(result.stdout)} chars",
                file=sys.stderr,
            )

        if result.returncode != 0:
            print(
                f"\n    [!] BLAST failed (returncode {result.returncode}) "
                f"for {label}.",
                file=sys.stderr,
            )
            return None

        # Guard against a silent server-side failure masquerading as a
        # genuine zero-hit result. BLAST+ -remote can complete the HTTP
        # round trip and exit 0 even when NCBI rejected or aborted the
        # search server-side (e.g. CPU usage limit hit, overload, scheduled
        # maintenance) — the only trace left is a message on stderr; stdout
        # is empty exactly like a real "no hits above threshold" result.
        # Without this check, that data loss is permanent: it gets written
        # to .progress as "done" and silently skipped on every future
        # --resume. This list is a heuristic, not exhaustive — known NCBI
        # failure signatures observed in the wild (e.g. "[blastsrv4.REAL]:
        # Error: CPU usage limit was exceeded"). Use --debug to see the raw
        # message if a batch is being retried unexpectedly so the list
        # can be extended.
        if not result.stdout.strip() and result.stderr.strip():
            lower_err = result.stderr.lower()
            failure_signatures = (
                "error",
                "failed",
                "timeout",
                "timed out",
                "cpu usage limit",
            )
            if any(sig in lower_err for sig in failure_signatures):
                print(
                    f"\n    [!] Treating this as a server-side failure, not "
                    f"a genuine zero-hit result, for {label} — will retry "
                    f"on --resume rather than marking it done.",
                    file=sys.stderr,
                )
                return None

        return result.stdout

    except subprocess.TimeoutExpired:
        print(
            f"\n    [!] Timeout ({timeout}s) for {label} — skipping.\n"
            f"        Increase with --timeout (current: {timeout}s"
            f"{', scaled by batch size' if len(records) > 1 else ''}).\n"
            f"        If using refseq_protein or swissprot, try --db nr first.\n"
            f"        If using nr and still timing out, NCBI may be rate-limiting\n"
            f"        your IP, or --batch-size may be too large. Wait 30-60\n"
            f"        minutes and retry, or lower --batch-size.",
            file=sys.stderr,
        )
        return None

    except FileNotFoundError:
        sys.exit(
            f"\n[!] '{program}' not found in PATH.\n"
            "    Install BLAST+ from:\n"
            "    https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
        )

    finally:
        Path(tmp_path).unlink(missing_ok=True)


# ── Main ──────────────────────────────────────────────────────────────────────


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    # ── Resolve and load sequences ─────────────────────────────────────────────
    # Loaded before anything else, since --list only needs this much and
    # shouldn't have to wait on a BLAST+ install check or print BLAST-run
    # settings (database, e-value, timeout...) that don't apply to it.
    input_files = _resolve_input_files(args.input)
    records = _load_multi_source(input_files)
    if not records:
        sys.exit(f"[!] No sequences found in: {input_files}")

    # ── --list: preview and exit, before any BLAST-specific setup ────────────
    if args.list:
        _run_list(records, args, input_files)
        sys.exit(0)

    # ── Fail-fast: check BLAST is installed before doing anything else ─────────
    if shutil.which(args.program) is None:
        sys.exit(
            f"\n[!] ERROR: '{args.program}' is not installed or not in your PATH.\n"
            "    Install NCBI BLAST+ from:\n"
            "    https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/\n"
        )

    # ── Resolve defaults ──────────────────────────────────────────────────────
    if args.output:
        output: Path = args.output
    elif len(input_files) == 1:
        output = input_files[0].with_suffix(".blast.tsv")
    else:
        sys.exit(
            "[!] -o/--output is required when -i resolves to more than one "
            f"source file (got {len(input_files)} files). With multiple "
            "genomes merged together, there's no single sensible default "
            "filename to guess — name it explicitly, e.g. "
            "-o combined_blast_results.tsv"
        )
    db: str = args.db or DEFAULT_DB[args.program]
    progress_path: Path = output.with_name(output.stem + ".progress")

    if args.batch_size < 1:
        sys.exit(f"[!] --batch-size must be >= 1 (got {args.batch_size}).")

    # --resume requires outfmt 6: that's the only format where a parseable
    # qseqid-bearing row tells us a sequence's hits were genuinely written
    # (vs. the run crashing mid-sequence). Other formats can still run fine,
    # they just don't support --resume.
    resume_active = args.resume
    if args.resume and args.outfmt != 6:
        print(
            "\n[!] WARNING: --resume requires outfmt 6 (tabular). "
            f"Ignoring --resume since -f {args.outfmt} was given.",
            file=sys.stderr,
        )
        resume_active = False

    effective_timeout = args.timeout * args.batch_size

    print(f"[*] Version       : {__version__}", file=sys.stderr)
    print(f"[*] Program       : {args.program}", file=sys.stderr)
    print(f"[*] Database      : {db} (remote)", file=sys.stderr)
    print(f"[*] E-value       : {args.evalue}", file=sys.stderr)
    print(f"[*] Output format : {args.outfmt}", file=sys.stderr)
    print(f"[*] Max hits/seq  : {args.max_hits}", file=sys.stderr)
    print(f"[*] Batch size    : {args.batch_size} sequence(s)/search", file=sys.stderr)
    print(
        f"[*] Timeout       : {args.timeout}s/seq -> {effective_timeout}s/search "
        f"(scaled by batch size)",
        file=sys.stderr,
    )
    print(f"[*] Debug mode    : {'ON' if args.debug else 'off'}", file=sys.stderr)
    print(f"[*] Resume        : {'ON' if resume_active else 'off'}", file=sys.stderr)
    if len(input_files) == 1:
        print(f"[*] Input         : {input_files[0]}", file=sys.stderr)
    else:
        print(
            f"[*] Input         : {len(input_files)} source file(s):", file=sys.stderr
        )
        for f in input_files[:10]:
            print(f"      - {f}", file=sys.stderr)
        if len(input_files) > 10:
            print(f"      ... and {len(input_files) - 10} more.", file=sys.stderr)
    print(f"[*] Output        : {output}", file=sys.stderr)
    print(f"    Total in file   : {len(records)}", file=sys.stderr)

    # ── Sequence-type sanity check (advisory only) ────────────────────────────
    _check_sequence_type(records, args.program)

    # ── Apply selection ───────────────────────────────────────────────────────
    selected = _select_sequences(records, args)

    if args.range:
        mode = f"range [{args.range[0]} -> {args.range[1]}]"
    elif args.pick:
        mode = f"pick {args.pick}"
    else:
        mode = "all"

    print(f"    Selection mode  : {mode}", file=sys.stderr)
    print(f"    To blast        : {len(selected)} sequence(s)", file=sys.stderr)

    # Each BATCH (not each sequence) is one NCBI search. NCBI's own developer
    # guidelines state they move users who submit more than 100 searches in
    # a 24-hour period to a slower queue, or in extreme cases block them —
    # and recommend running large batches off-peak (weekends, or 9pm-5am
    # Eastern on weekdays). This is advisory only; it does not block the run.
    n_batches_preview = -(-len(selected) // args.batch_size)  # ceil division
    if n_batches_preview > 50:
        print(
            f"\n[!] {len(selected)} sequences in {n_batches_preview} batches "
            f"= {n_batches_preview} separate NCBI searches in this run. "
            f"NCBI's own guidelines warn that submitting more than 100 "
            f"searches in 24 hours can get an IP moved to a slower queue, "
            f"or blocked in extreme cases, and recommend running large "
            f"batches on weekends or between 9pm-5am Eastern on weekdays. "
            f"Consider --resume if this gets interrupted, and avoid "
            f"stacking multiple large runs back-to-back today.",
            file=sys.stderr,
        )

    # ── Resume: filter out already-completed sequences ───────────────────────
    completed_ids: set[str] = set()
    if resume_active and progress_path.exists():
        completed_ids = {
            line.strip()
            for line in progress_path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        }

    if completed_ids:
        before = len(selected)
        selected = [(n, r) for n, r in selected if r.id not in completed_ids]
        skipped = before - len(selected)
        print(
            f"    Resuming        : {skipped} already-completed sequence(s) "
            f"skipped (from {progress_path.name})",
            file=sys.stderr,
        )
        if not selected:
            print(
                "\n[*] All selected sequences were already completed in a "
                "previous run. Nothing to do.",
                file=sys.stderr,
            )
            sys.exit(0)

    # Append (not overwrite) only when there's genuine prior progress to
    # continue. A --resume run against a file with no progress yet behaves
    # exactly like a fresh run, just starts the bookkeeping for next time.
    append_mode = bool(completed_ids)

    if args.delay < NCBI_MIN_DELAY:
        print(
            f"\n[!] Warning: --delay {args.delay}s is below the NCBI minimum "
            f"({NCBI_MIN_DELAY}s). Your IP may be rate-limited.",
            file=sys.stderr,
        )

    # ── Build batches ──────────────────────────────────────────────────────────
    # Chunk the (already resume-filtered) selection into groups of
    # --batch-size. A batch of 1 reproduces the original one-sequence-per-
    # search behavior exactly.
    batches: list[list[tuple[int, SeqRecord]]] = [
        selected[i : i + args.batch_size]
        for i in range(0, len(selected), args.batch_size)
    ]

    # ── BLAST loop ────────────────────────────────────────────────────────────
    print(
        f"\n[*] Starting remote BLAST: {len(selected)} sequence(s) in "
        f"{len(batches)} batch(es) ({args.delay}s delay between batches)...",
        file=sys.stderr,
    )

    hits_total = 0
    out_mode = "a" if append_mode else "w"

    # Progress tracking is written whenever outfmt 6 is used, regardless of
    # whether --resume was passed THIS run — so a future --resume always has
    # something to work from, even if this particular run never needed it.
    write_progress = args.outfmt == 6
    progress_ctx = (
        open(progress_path, out_mode, encoding="utf-8")
        if write_progress
        else contextlib.nullcontext()
    )

    with open(output, out_mode, encoding="utf-8") as out, progress_ctx as progress:
        if args.outfmt == 6 and not append_mode:
            out.write(HEADER_ROW + "\n")

        for batch_idx, batch in enumerate(batches):
            batch_records = [record for _, record in batch]
            label = _batch_label(batch_records)

            print(
                f"\n  [Batch {batch_idx + 1}/{len(batches)}] "
                f"{len(batch_records)} sequence(s): {label}...",
                file=sys.stderr,
            )

            _check_batch_query_length(batch_records, args.program)

            blast_output = _blast_batch(
                batch_records,
                args.program,
                db,
                args.evalue,
                args.outfmt,
                args.max_hits,
                effective_timeout,
                args.debug,
            )

            if blast_output is None:
                # Failed/timed out — NOT marked as progress. With
                # batch-size > 1, EVERY sequence in this batch retries
                # together on the next --resume run (batch is the retry
                # granularity, not the individual sequence).
                print(
                    f"      -> Batch failed/timed out "
                    f"({len(batch_records)} sequence(s) affected)."
                    + (
                        " Will retry this whole batch on the next --resume run."
                        if write_progress
                        else ""
                    ),
                    file=sys.stderr,
                )
            else:
                if blast_output.strip():
                    out.write(blast_output)
                    out.flush()
                    n_hits = len([l for l in blast_output.strip().splitlines() if l])
                    hits_total += n_hits
                    print(
                        f"      -> {n_hits} hit(s) across this batch",
                        file=sys.stderr,
                    )
                else:
                    print(
                        f"      -> No hits for any sequence in this batch. "
                        f"(Run with --debug to see what NCBI returned.)",
                        file=sys.stderr,
                    )
                # Definitively completed (hit or clean zero-hit) — safe to
                # skip on a future --resume run. Marked for EVERY sequence
                # in the batch, since the batch as a whole succeeded.
                if progress is not None:
                    for record in batch_records:
                        progress.write(record.id + "\n")
                    progress.flush()

            if batch_idx < len(batches) - 1:
                print(f"      Waiting {args.delay}s...", file=sys.stderr)
                time.sleep(args.delay)

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"\n{'=' * 44}", file=sys.stderr)
    print(f"  BLAST COMPLETE", file=sys.stderr)
    print(f"{'=' * 44}", file=sys.stderr)
    if completed_ids:
        print(f"  Skipped (resumed) : {len(completed_ids)}", file=sys.stderr)
    print(f"  Sequences blasted : {len(selected)}", file=sys.stderr)
    print(f"  NCBI searches sent: {len(batches)}", file=sys.stderr)
    print(f"  Total hits        : {hits_total}", file=sys.stderr)
    print(f"  Output written to : {output.resolve()}", file=sys.stderr)
    if write_progress:
        print(f"  Progress file     : {progress_path.resolve()}", file=sys.stderr)
    print(f"{'=' * 44}", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print(
            "\n[!] Pipeline interrupted by user (Ctrl+C).\n"
            "    Partial results may have been written to the output file.",
            file=sys.stderr,
        )
        sys.exit(130)
