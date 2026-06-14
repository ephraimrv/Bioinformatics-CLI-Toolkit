"""
Remote BLAST Runner

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

CRITICAL NOTE — Taxonomy Fields (sscinames, scomnames) with `-remote`:
    Taxonomic columns (sscinames, scomnames) require a local Taxonomy Database
    (taxdb) to work. When using `-remote`, NCBI servers return only alignment
    data, not taxonomy lookups. Without local taxdb files, these columns return N/A.

    WORKAROUND: Extract organism name from stitle using regex in Python:

        import pandas as pd
        import re

        df = pd.read_csv('output.blast.tsv', sep='\t')
        df['organism'] = df['stitle'].str.extract(r'\[(.*?)\]', expand=False)
        print(df[['qseqid', 'sseqid', 'organism', 'evalue']])

    This extracts the organism name from the stitle brackets reliably.

SEQUENCE SELECTION MODES:
    (none)                      Blast every sequence in the file (default)
    --range LOCUS_START LOCUS_END  Blast a contiguous slice in file order
    --pick LOCUS_TAG [...]         Blast specific sequences by ID
    --list                         Preview all IDs in the file and exit

NCBI USAGE POLICY:
    Remote BLAST uses NCBI's public servers. A 5-second delay is enforced
    between requests by default (NCBI guideline: no more than 3 requests/sec).
    For large batches (>100 sequences), consider downloading a local BLAST
    database: https://ftp.ncbi.nlm.nih.gov/blast/db/

License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).

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
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.2.1"

import argparse
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

NCBI_MIN_DELAY = 5  # seconds between requests (NCBI etiquette: ≤3 requests/sec)


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
  --range START END   e.g. --range 1 20  → sequences 1 through 20
  --pick  N [N ...]   e.g. --pick 1 36 45 → sequences 1, 36, and 45

examples:
  python3 remote_blast.py -i proteins.faa
  python3 remote_blast.py -i proteins.faa -p blastp --range 1 20 -e 1e-10
  python3 remote_blast.py -i proteins.faa -p blastp --pick 1 36 45
  python3 remote_blast.py -i genome.fna -p blastn
  python3 remote_blast.py -i contigs.fna -p blastx --db nr --max-hits 100
        """,
    )

    # ── Required ──────────────────────────────────────────────────────────────
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=Path,
        metavar="FASTA",
        help="Input FASTA file (.faa for protein, .fna/.fasta for nucleotide).",
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
        "--delay",
        type=float,
        default=NCBI_MIN_DELAY,
        metavar="SEC",
        help=(
            f"Seconds to wait between requests. "
            f"Default: {NCBI_MIN_DELAY} (NCBI usage policy)."
        ),
    )

    parser.add_argument(
        "--list",
        action="store_true",
        help=(
            "Print all sequence IDs in the input file with their position "
            "numbers, then exit. Use this to find locus tags for --range and --pick."
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


def _select_sequences(
    records: list[SeqRecord], args: argparse.Namespace
) -> list[tuple[int, SeqRecord]]:
    """Return (1-indexed position, record) pairs based on the selection mode.

    Sequences are identified by locus tag (the sequence ID from the FASTA
    header — the first word after '>'), not by positional index. This means
    the user never has to count sequences in the file manually.

    For --range LOCUS_START LOCUS_END: finds both tags in the file, then
    selects everything between them (inclusive) in file order.

    For --pick LOCUS_TAG [...]: selects each named sequence in the order
    the tags were given on the command line.

    Args:
        records: All SeqRecord objects parsed from the input file.
        args:    Parsed arguments (may contain .range or .pick).

    Returns:
        Ordered list of (seq_number, SeqRecord) tuples to blast.
    """
    # Build lookup: locus_tag → 0-indexed position in file
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

    # Default: all sequences
    return [(i + 1, rec) for i, rec in enumerate(records)]


# ── BLAST execution ───────────────────────────────────────────────────────────


def _blast_one(
    record: SeqRecord,
    program: str,
    db: str,
    evalue: float,
    outfmt: int,
    max_hits: int,
) -> str:
    """Run BLAST on a single sequence and return the raw output string.

    Writes the sequence to a temporary file and passes it to BLAST via
    -query. The stdin approach (-query -) is documented but unreliable
    with -remote: BLAST+ does not consistently pipe stdin through to
    NCBI's servers across all versions. Temp file is the safe approach.

    Args:
        record:   The SeqRecord to query.
        program:  BLAST program name (blastp/blastn/blastx).
        db:       Database name.
        evalue:   E-value cutoff.
        outfmt:   Output format number.
        max_hits: Maximum hits to return.

    Returns:
        BLAST output as a string, or "" on error/timeout.
    """
    tmp_path = Path(tempfile.mktemp(suffix=".fasta"))
    try:
        with open(tmp_path, "w") as tmp:
            SeqIO.write(record, tmp, "fasta")

        fmt_arg = f"6 {COLUMNS}" if outfmt == 6 else str(outfmt)
        cmd = [
            program,
            "-query",
            str(tmp_path),
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

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,  # 5 minutes — remote BLAST can be slow
        )

        if result.returncode != 0:
            print(
                f"\n    [!] BLAST error for {record.id}:\n"
                f"        {result.stderr.strip()}",
                file=sys.stderr,
            )
            return ""

        return result.stdout

    except subprocess.TimeoutExpired:
        print(
            f"\n    [!] Timeout (5 min) for {record.id} — skipping.",
            file=sys.stderr,
        )
        return ""

    except FileNotFoundError:
        sys.exit(
            f"\n[!] '{program}' not found in PATH.\n"
            "    Install BLAST+ from:\n"
            "    https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
        )

    finally:
        tmp_path.unlink(missing_ok=True)


# ── Main ──────────────────────────────────────────────────────────────────────


def main() -> None:
    parser = _build_parser()
    args = parser.parse_args()

    # ── Fail-fast: check BLAST is installed ───────────────────────────────────
    if shutil.which(args.program) is None:
        sys.exit(
            f"\n[!] ERROR: '{args.program}' is not installed or not in your PATH.\n"
            "    Install NCBI BLAST+ from:\n"
            "    https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/\n"
        )

    # ── Resolve defaults ──────────────────────────────────────────────────────
    if not args.input.exists():
        sys.exit(f"[!] Input file not found: {args.input}")

    output: Path = args.output or args.input.with_suffix(".blast.tsv")
    db: str = args.db or DEFAULT_DB[args.program]

    print(f"[*] Program       : {args.program}", file=sys.stderr)
    print(f"[*] Database      : {db} (remote)", file=sys.stderr)
    print(f"[*] E-value       : {args.evalue}", file=sys.stderr)
    print(f"[*] Output format : {args.outfmt}", file=sys.stderr)
    print(f"[*] Max hits/seq  : {args.max_hits}", file=sys.stderr)
    print(f"[*] Input         : {args.input}", file=sys.stderr)
    print(f"[*] Output        : {output}", file=sys.stderr)

    # ── Load sequences ────────────────────────────────────────────────────────
    print(f"\n[*] Loading sequences...", file=sys.stderr)
    records = list(SeqIO.parse(args.input, "fasta"))
    if not records:
        sys.exit(f"[!] No sequences found in {args.input}")
    print(f"    Total in file   : {len(records)}", file=sys.stderr)

    # ── --list: print all IDs and exit ────────────────────────────────────────
    if args.list:
        print(f"\nSequences in {args.input.name} ({len(records)} total):\n")
        # Strip the ID prefix from description (BioPython includes it)
        for i, rec in enumerate(records, 1):
            desc = rec.description[len(rec.id) :].strip()
            desc_col = f"  {desc[:60]}" if desc else ""
            print(f"  [{i:>4}]  {rec.id:<20}{desc_col}")
        print(
            f"\n  Use the ID (second column) with --range or --pick.\n"
            f"  Example: --range {records[0].id} {records[-1].id}"
        )
        sys.exit(0)

    # ── Apply selection ───────────────────────────────────────────────────────
    selected = _select_sequences(records, args)

    if args.range:
        mode = f"range [{args.range[0]} → {args.range[1]}]"
    elif args.pick:
        mode = f"pick {args.pick}"
    else:
        mode = "all"

    print(f"    Selection mode  : {mode}", file=sys.stderr)
    print(f"    To blast        : {len(selected)} sequence(s)", file=sys.stderr)

    if args.delay < NCBI_MIN_DELAY:
        print(
            f"\n[!] Warning: --delay {args.delay}s is below the NCBI minimum "
            f"({NCBI_MIN_DELAY}s). Your IP may be rate-limited.",
            file=sys.stderr,
        )

    # ── BLAST loop ────────────────────────────────────────────────────────────
    print(
        f"\n[*] Starting remote BLAST " f"({args.delay}s delay between requests)...",
        file=sys.stderr,
    )

    hits_total = 0

    with open(output, "w", encoding="utf-8") as out:
        # Write header only for tabular outfmt 6
        if args.outfmt == 6:
            out.write(HEADER_ROW + "\n")

        for idx, (seq_num, record) in enumerate(selected):
            print(
                f"\n  [{idx + 1}/{len(selected)}] "
                f"Seq {seq_num}: {record.id} "
                f"({len(record.seq)} residues)...",
                file=sys.stderr,
            )

            blast_output = _blast_one(
                record,
                args.program,
                db,
                args.evalue,
                args.outfmt,
                args.max_hits,
            )

            if blast_output.strip():
                out.write(blast_output)
                out.flush()  # persist each result as it arrives
                n_hits = len([l for l in blast_output.strip().splitlines() if l])
                hits_total += n_hits
                print(f"      → {n_hits} hit(s)", file=sys.stderr)
            else:
                print(f"      → No hits above threshold.", file=sys.stderr)

            # NCBI rate limiting — no delay needed after the last sequence
            if idx < len(selected) - 1:
                print(
                    f"      Waiting {args.delay}s...",
                    file=sys.stderr,
                )
                time.sleep(args.delay)

    # ── Summary ───────────────────────────────────────────────────────────────
    print(f"\n{'=' * 44}", file=sys.stderr)
    print(f"  BLAST COMPLETE", file=sys.stderr)
    print(f"{'=' * 44}", file=sys.stderr)
    print(f"  Sequences blasted : {len(selected)}", file=sys.stderr)
    print(f"  Total hits        : {hits_total}", file=sys.stderr)
    print(f"  Output written to : {output.resolve()}", file=sys.stderr)
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
        sys.exit(130)  # standard exit code for SIGINT
