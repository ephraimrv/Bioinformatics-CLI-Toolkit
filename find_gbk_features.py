# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

"""A human-readable viewer and search tool for GenBank genome files.

Provides six exploration modes for .gbk and .gbff files without requiring
manual inspection of raw GenBank format. Designed for use alongside
extract_genome_region.py as part of a no-antiSMASH genome exploration workflow.

antiSMASH region boundaries (e.g. 53317-78823) are computed by antiSMASH and
are not stored in the GBFF file. The ``--context`` mode provides a practical
manual equivalent for defining those boundaries.

Typical workflow for an unfamiliar genome:
    1. ``--list-sequences``    Identify available contigs and locus tag ranges.
    2. ``--list-products``     Enumerate searchable functional categories.
    3. ``-q "keyword"``        Find anchor genes of interest.
    4. ``--context LOCUS_TAG`` Explore the genomic neighbourhood.
    5. Pass coordinates to ``extract_genome_region.py``.

The ``-q`` flag accepts one or more keywords in a single run. When multiple
keywords are given, results are shown grouped by keyword in the terminal and
sorted by keyword order in the TSV file.

The ``--search-type`` flag controls where ``-q`` searches:
    - ``product`` (default): searches /product=, /gene_kind=, /gene_functions=,
      /sec_met_domain=, and /note= — use this when you know what a gene *does*.
    - ``locus``: substring match on /locus_tag= — use this when you know part
      of a locus tag, e.g. ``-q "RHP56_RS003"`` to find all genes in that range.
    - ``locus-exact``: exact match on /locus_tag= — use this to look up one
      specific gene before passing it to ``--context``.

When ``-o`` is specified, results are written as a TSV file instead of text.
The TSV columns vary by mode and include verbose annotation fields when ``-v``
is also used. Without ``-o``, all output is formatted text to the terminal.

Note:
    This script is part of ongoing research and is associated with an upcoming
    publication. Correct attribution is requested when used in derivative works.
    Released under the MIT License. See the LICENSE file in the repository root.

Example Usage:
    Display the sequences contained in the file::

        python3 find_gbk_features.py -i genome.gbff --list-sequences
        python3 find_gbk_features.py -i genome.gbff --list-sequences -o seqs.tsv

    Discover functional categories in an unfamiliar genome::

        python3 find_gbk_features.py -i genome.gbff --list-products
        python3 find_gbk_features.py -i genome.gbff --list-products --min-count 3

    Search by one or more product keywords::

        python3 find_gbk_features.py -i genome.gbff -q "bacteriocin" --seq NZ_CP134351.1
        python3 find_gbk_features.py -i genome.gbff -q "bacteriocin" "transporter" "regulator"
        python3 find_gbk_features.py -i genome.gbff -q "bacteriocin" "transporter" -o hits.tsv

    Search by locus tag::

        python3 find_gbk_features.py -i genome.gbff -q "RHP56_RS003" --search-type locus
        python3 find_gbk_features.py -i genome.gbff -q "RHP56_RS00345" --search-type locus-exact

    List all features on a specific contig::

        python3 find_gbk_features.py -i genome.gbff --seq NZ_CP134351.1
        python3 find_gbk_features.py -i genome.gbff --seq NZ_CP134351.1 -o all_genes.tsv

    Explore the neighbourhood of a known anchor gene::

        python3 find_gbk_features.py -i genome.gbff \\
            --context RHP56_RS00345 --window 15000 --seq NZ_CP134351.1

    List features in a coordinate range (1-based)::

        python3 find_gbk_features.py -i genome.gbff --c1 50000 --c2 80000 --seq NZ_CP134351.1
        python3 find_gbk_features.py -i genome.gbff --c1 50000 --c2 80000 -o region.tsv
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "1.5.0"

import sys
import argparse
from collections import Counter
from pathlib import Path
from typing import Any

try:
    from Bio import SeqIO
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )


# ── CLI ───────────────────────────────────────────────────────────────────────


def get_args() -> argparse.Namespace:
    """Configures the CLI parser and returns parsed arguments.

    Returns:
        An ``argparse.Namespace`` object containing all parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description=(
            "A human-readable viewer for GenBank files. "
            "Use --list-sequences first, then --list-products to discover "
            "keywords, then -q / --context to explore specific genes. "
            "Add -o FILE.tsv to save results as a tab-delimited file."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the GenBank file (.gbk or .gbff).",
    )

    # ── Discovery modes ───────────────────────────────────────────────────────
    disc = parser.add_argument_group(
        "Discovery Modes",
        "Use these first on an unfamiliar genome to understand its content.",
    )
    disc.add_argument(
        "--list-sequences",
        action="store_true",
        help=(
            "List all sequences in the file with CDS statistics: "
            "length, gene count, locus tag range, and gene coordinate span. "
            "Add -o FILE.tsv to export as a table."
        ),
    )
    disc.add_argument(
        "--list-products",
        action="store_true",
        help=(
            "List all unique product annotations sorted by frequency. "
            "Hypothetical proteins are excluded by default; "
            'use -q "hypothetical" to find them. '
            "Add -o FILE.tsv to export as a table."
        ),
    )
    disc.add_argument(
        "--min-count",
        type=int,
        default=1,
        metavar="N",
        help=(
            "With --list-products: only show products appearing N or more times. "
            "Default: 1 (show all)."
        ),
    )

    # ── Search modes ──────────────────────────────────────────────────────────
    search = parser.add_argument_group(
        "Search Modes",
        "Find specific features. Use -q alone or combine with other flags.",
    )
    search.add_argument(
        "-q",
        "--query",
        type=str,
        nargs="+",
        default=None,
        metavar="KEYWORD",
        help=(
            "One or more case-insensitive search terms. "
            "Matched against /product=, /gene_kind=, /gene_functions=, "
            "/sec_met_domain=, and /note=. "
            'Enclose multi-word terms in quotes: -q "ABC transporter" "response regulator". '
            "For locus tag searches use --search-type locus or locus-exact."
        ),
    )
    search.add_argument(
        "--search-type",
        type=str,
        choices=["product", "locus", "locus-exact"],
        default="product",
        help=(
            "'product' searches annotation text qualifiers (default). "
            "'locus' does a substring match on /locus_tag=. "
            "'locus-exact' matches a single /locus_tag= exactly."
        ),
    )
    search.add_argument(
        "--context",
        type=str,
        default=None,
        metavar="LOCUS_TAG",
        help=(
            "Show all genes within --window bp of the given locus tag. "
            "The no-antiSMASH way to define cluster boundaries. "
            "Combine with --seq on multi-record files."
        ),
    )
    search.add_argument(
        "--window",
        type=int,
        default=10000,
        metavar="BP",
        help="Neighbourhood window in bp for --context mode. Default: 10000.",
    )
    search.add_argument(
        "--c1",
        type=int,
        metavar="START_BP",
        default=None,
        help="Coordinate-range start position (1-based, inclusive).",
    )
    search.add_argument(
        "--c2",
        type=int,
        metavar="END_BP",
        default=None,
        help="Coordinate-range end position (1-based, inclusive).",
    )

    # ── Filters and output ────────────────────────────────────────────────────
    parser.add_argument(
        "--seq",
        type=str,
        default=None,
        help=(
            "Sequence ID to operate on (e.g. 'NZ_CP134351.1'). "
            "Used alone: lists every feature on that contig. "
            "Used with any search mode: restricts results to that contig. "
            "Run --list-sequences to find valid IDs."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help=(
            "Save results as a TSV file instead of printing text to the terminal. "
            "The columns depend on the active mode; "
            "verbose annotation fields are included when -v is also used."
        ),
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help=(
            "Show extended annotation details per feature: "
            "protein accession, gene name, protein length, GO terms, "
            "annotation notes, and antiSMASH biosynthetic role where present. "
            "Adds corresponding columns to TSV output when -o is used."
        ),
    )

    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    """Validates argument combinations and exits with a helpful message on error.

    Args:
        args: Parsed argument namespace from ``get_args()``.

    Raises:
        SystemExit: If a required argument combination is missing or invalid.
    """
    if args.list_sequences or args.list_products:
        return

    has_query = args.query is not None
    has_context = args.context is not None
    has_coord = args.c1 is not None or args.c2 is not None
    has_seq = args.seq is not None

    if not any([has_query, has_context, has_coord, has_seq]):
        sys.exit(
            "\n[!] Error: specify a mode. Start here if the genome is unfamiliar:\n"
            "    --list-sequences             what sequences are in this file\n"
            "    --list-products              what keywords exist\n"
            "    --seq SEQ_ID                 list all features on a contig\n\n"
            "    Then use one of these to find specific features:\n"
            "    -q 'keyword' [keyword2 ...]  keyword search (one or more terms)\n"
            "    --context LOCUS_TAG          neighbourhood explorer\n"
            "    --c1 START --c2 END          coordinate range\n"
        )

    if has_coord:
        if args.c1 is None or args.c2 is None:
            sys.exit("\n[!] Error: --c1 and --c2 must always be used together.\n")
        if args.c1 >= args.c2:
            sys.exit(
                f"\n[!] Error: --c1 ({args.c1:,}) must be less than "
                f"--c2 ({args.c2:,}).\n"
            )

    if has_context and args.window < 1:
        sys.exit("\n[!] Error: --window must be a positive integer.\n")


# ── Feature helpers ───────────────────────────────────────────────────────────


def get_annotation_text(feature) -> str:
    """Builds a single searchable string from all annotation qualifiers.

    Checks /product=, /gene_kind=, /gene_functions=, /sec_met_domain=, and
    /note=, covering both NCBI GBFF (/product=) and antiSMASH region GBK
    (/gene_kind=, /gene_functions=) in a single call.

    Args:
        feature: A BioPython ``SeqFeature`` object.

    Returns:
        A single lowercase string combining all annotation qualifier values.
    """
    parts = [
        feature.qualifiers.get("product", [""])[0],
        feature.qualifiers.get("gene_kind", [""])[0],
        " ".join(feature.qualifiers.get("gene_functions", [])),
        " ".join(feature.qualifiers.get("sec_met_domain", [])),
        feature.qualifiers.get("note", [""])[0],
    ]
    return " ".join(p for p in parts if p).lower()


def get_display_product(feature) -> str:
    """Returns the best available product description for display.

    Prefers /product=. Falls back to /gene_kind= + /gene_functions= for
    antiSMASH region files where /product= is absent on CDS features.

    Args:
        feature: A BioPython ``SeqFeature`` object.

    Returns:
        A human-readable product string, or ``'(no annotation)'``.
    """
    product = feature.qualifiers.get("product", [""])[0]
    if product:
        return product

    kind = feature.qualifiers.get("gene_kind", [""])[0]
    funcs = feature.qualifiers.get("gene_functions", [])
    func = funcs[0].replace("\n", " ").replace("  ", " ") if funcs else ""

    if kind and func:
        return f"{kind} [{func}]"
    if kind:
        return f"{kind} (antiSMASH)"
    if func:
        return func
    return "(no annotation)"


def feature_to_dict(feature, record_id: str) -> dict:
    """Converts a BioPython CDS feature into a standardized display dictionary.

    Collects all fields needed for both standard and verbose output in a
    single pass. BioPython's ``feature.location.start`` is 0-indexed; this
    function adds 1 to convert to the 1-based coordinate system used in
    GBFF files and genome viewers.

    Args:
        feature:   A BioPython ``SeqFeature`` object of type ``CDS``.
        record_id: The ID of the parent sequence record.

    Returns:
        A dictionary with keys: ``locus``, ``product``, ``start``, ``end``,
        ``strand``, ``record_id``, ``protein_id``, ``gene_name``,
        ``aa_length``, ``go_function``, ``go_process``, ``note``,
        ``gene_kind``, ``gene_funcs``.
    """
    translation = feature.qualifiers.get("translation", [""])[0]
    return {
        "locus": feature.qualifiers.get("locus_tag", [""])[0],
        "product": get_display_product(feature),
        "start": int(feature.location.start) + 1,
        "end": int(feature.location.end),
        "strand": "+" if feature.location.strand == 1 else "-",
        "record_id": record_id,
        "protein_id": feature.qualifiers.get("protein_id", [""])[0],
        "gene_name": feature.qualifiers.get("gene", [""])[0],
        "aa_length": len(translation) if translation else 0,
        "go_function": feature.qualifiers.get("GO_function", []),
        "go_process": feature.qualifiers.get("GO_process", []),
        "note": feature.qualifiers.get("note", [""])[0],
        "gene_kind": feature.qualifiers.get("gene_kind", [""])[0],
        "gene_funcs": feature.qualifiers.get("gene_functions", []),
    }


def collect_cds(record) -> list[dict]:
    """Returns all CDS features from a record as standardized feature dicts.

    Args:
        record: A BioPython ``SeqRecord`` object.

    Returns:
        A list of feature dicts produced by ``feature_to_dict()``,
        one per CDS feature in the record.
    """
    return [
        feature_to_dict(feat, record.id)
        for feat in record.features
        if feat.type == "CDS"
    ]


# ── Proximity grouping ────────────────────────────────────────────────────────


def group_by_proximity(
    results: list[dict],
    gap_threshold: int = 50_000,
) -> list[list[dict]]:
    """Groups features into spatially adjacent clusters.

    Features separated by more than ``gap_threshold`` bp are placed in
    separate groups, preventing extraction suggestions from spanning the
    whole chromosome when hits are far apart.

    Args:
        results:       Feature dicts sortable by ``'start'``.
        gap_threshold: Maximum gap in bp between consecutive hits in one group.
                       Default: 50000.

    Returns:
        A list of groups, each containing spatially adjacent feature dicts.
    """
    if not results:
        return []

    sorted_hits = sorted(results, key=lambda r: r["start"])
    groups = [[sorted_hits[0]]]

    for feat in sorted_hits[1:]:
        if feat["start"] - groups[-1][-1]["end"] <= gap_threshold:
            groups[-1].append(feat)
        else:
            groups.append([feat])

    return groups


# ── TSV helpers ───────────────────────────────────────────────────────────────


def _build_feature_tsv_row(
    feat: dict,
    keyword: str | None = None,
    is_anchor: bool | None = None,
    verbose: bool = False,
) -> dict[str, Any]:
    """Builds a flat row dict for TSV export from a feature dict.

    Args:
        feat:      Feature dict from ``feature_to_dict()``.
        keyword:   If set, prepends a ``Keyword`` column.
        is_anchor: If set, prepends an ``Is_Anchor`` column (``'yes'``/``'no'``).
        verbose:   If ``True``, appends GO terms, annotation note, and
                   antiSMASH-specific fields.

    Returns:
        An ordered dict with all requested columns. Suitable for writing
        as a TSV row.
    """
    row: dict[str, Any] = {}

    if keyword is not None:
        row["Keyword"] = keyword
    if is_anchor is not None:
        row["Is_Anchor"] = "yes" if is_anchor else "no"

    row.update(
        {
            "Sequence_ID": feat["record_id"],
            "Locus_Tag": feat["locus"],
            "Gene": feat["gene_name"],
            "Start": feat["start"],
            "End": feat["end"],
            "Strand": feat["strand"],
            "Product": feat["product"],
            "Protein_ID": feat["protein_id"],
            "Length_aa": feat["aa_length"],
        }
    )

    if verbose:
        row["GO_Function"] = "; ".join(feat["go_function"])
        row["GO_Process"] = "; ".join(feat["go_process"])
        row["Note"] = feat["note"]
        row["Gene_Kind"] = feat["gene_kind"]
        row["Gene_Functions"] = "; ".join(feat["gene_funcs"])

    return row


def _write_tsv(rows: list[dict], output_path: Path) -> None:
    """Writes a list of row dicts to a tab-delimited file.

    Args:
        rows:        List of dicts. All dicts must have identical keys.
        output_path: Destination file path.

    Note:
        Uses UTF-8 with BOM (``utf-8-sig``) for compatibility with Excel.
        Prints a confirmation or error message to ``stderr``.
    """
    if not rows:
        print("[!] No results to write.", file=sys.stderr)
        return

    headers = list(rows[0].keys())
    try:
        with open(output_path, "w", encoding="utf-8-sig", newline="") as f:
            f.write("\t".join(headers) + "\n")
            for row in rows:
                f.write("\t".join(str(row.get(h, "")) for h in headers) + "\n")
        print(f"[*] TSV written \u2192 {output_path.resolve()}", file=sys.stderr)
    except OSError as exc:
        print(f"[!] Could not write TSV to '{output_path}': {exc}", file=sys.stderr)


# ── Text output formatters ────────────────────────────────────────────────────


def write_feature_table(features: list[dict], args, out) -> None:
    """Writes a formatted feature list to ``out``.

    Standard mode shows locus tag, product, and position. Verbose mode
    (``-v``) adds protein accession, gene name, protein length, GO terms,
    annotation notes, and antiSMASH biosynthetic role where present.

    Args:
        features: Feature dicts from ``feature_to_dict()``.
        args:     Parsed argument namespace (read for ``args.verbose``).
        out:      An open, writable file-like object.
    """
    for i, r in enumerate(features, 1):
        out.write(f"Feature {i}:\n")
        out.write(f"  Locus tag  : {r['locus']}\n")
        out.write(f"  Product    : {r['product']}\n")
        out.write(
            f"  Position   : {r['start']:,}\u2013{r['end']:,} bp "
            f"({r['strand']} strand)\n"
        )

        if args.verbose:
            if r["protein_id"]:
                out.write(f"  Protein ID : {r['protein_id']}\n")
            if r["gene_name"]:
                out.write(f"  Gene name  : {r['gene_name']}\n")
            if r["aa_length"]:
                out.write(f"  Length     : {r['aa_length']} aa\n")
            for go in r["go_function"][:2]:
                out.write(f"  GO func    : {go.strip()[:72]}\n")
            for go in r["go_process"][:2]:
                out.write(f"  GO proc    : {go.strip()[:72]}\n")
            if r["note"]:
                note_display = r["note"][:80] + ("..." if len(r["note"]) > 80 else "")
                out.write(f"  Note       : {note_display}\n")
            if r["gene_kind"]:
                out.write(f"  Gene kind  : {r['gene_kind']}\n")
            for gf in r["gene_funcs"][:2]:
                out.write(f"  Gene func  : {gf.strip()[:72]}\n")
            if not args.seq:
                out.write(f"  Sequence   : {r['record_id']}\n")

        out.write("\n")


def write_extraction_suggestions(
    results: list[dict],
    input_path: Path,
    seq_filter: str | None,
    out,
) -> None:
    """Writes proximity-grouped extraction command suggestions.

    Groups hits by a 50 kb gap threshold so spatially distant hits produce
    separate commands rather than one chromosome-spanning range.

    Args:
        results:     All matching feature dicts.
        input_path:  Path to the source GenBank file (used in the command).
        seq_filter:  The ``--seq`` value if active, otherwise ``None``.
        out:         An open, writable file-like object.
    """
    if len(results) < 2:
        return

    groups = group_by_proximity(results)
    seq_flag = f" \\\n    --seq {seq_filter}" if seq_filter else ""

    if len(groups) == 1:
        c1 = min(r["start"] for r in groups[0])
        c2 = max(r["end"] for r in groups[0])
        out.write("=" * 70 + "\n")
        out.write("SUGGESTED EXTRACTION COORDINATES (1-based):\n")
        out.write("=" * 70 + "\n\n")
        out.write(
            f"  python3 extract_genome_region.py \\\n"
            f"    -i {input_path.name} \\\n"
            f"    --c1 {c1} --c2 {c2}{seq_flag} \\\n"
            f"    --faa proteins.faa --fna region.fna\n"
        )
    else:
        out.write("=" * 70 + "\n")
        out.write(
            f"SUGGESTED EXTRACTION COORDINATES ({len(groups)} SEPARATE CLUSTERS):\n"
            f"  Hits span {len(groups)} distinct genomic locations.\n"
        )
        out.write("=" * 70 + "\n\n")
        for idx, group in enumerate(groups, 1):
            c1 = min(r["start"] for r in group)
            c2 = max(r["end"] for r in group)
            span_kb = (c2 - c1) / 1000
            out.write(
                f"  Cluster {idx}  ({len(group)} hit(s), ~{span_kb:.1f} kb "
                f"at {c1:,}\u2013{c2:,} bp):\n"
                f"    python3 extract_genome_region.py \\\n"
                f"      -i {input_path.name} \\\n"
                f"      --c1 {c1} --c2 {c2}{seq_flag} \\\n"
                f"      --faa cluster{idx}_proteins.faa "
                f"--fna cluster{idx}_region.fna\n\n"
            )


# ── Mode runners ──────────────────────────────────────────────────────────────


def run_list_sequences(args) -> None:
    """Lists every sequence in the file with CDS statistics.

    Shows sequence ID, length, organism, CDS count, first and last locus tag,
    and the coordinate span of the coding region.

    When ``-o`` is specified, writes a TSV with columns: ``Sequence_ID``,
    ``Length_bp``, ``Organism``, ``CDS_Count``, ``First_Tag``, ``Last_Tag``,
    ``Gene_Start``, ``Gene_End``.

    Args:
        args: Parsed argument namespace.
    """
    rows: list[dict] = []

    print(f"\n[*] Sequences in '{args.input.name}':\n")

    for record in SeqIO.parse(args.input, "genbank"):
        org = record.annotations.get("organism", "")
        cds_features = [f for f in record.features if f.type == "CDS"]
        n_cds = len(cds_features)

        if cds_features:
            locus_tags = [
                f.qualifiers.get("locus_tag", [""])[0]
                for f in cds_features
                if f.qualifiers.get("locus_tag", [""])[0]
            ]
            first_tag = locus_tags[0] if locus_tags else "(none)"
            last_tag = locus_tags[-1] if locus_tags else "(none)"
            gene_start = min(int(f.location.start) + 1 for f in cds_features)
            gene_end = max(int(f.location.end) for f in cds_features)
        else:
            first_tag = last_tag = "(no CDS)"
            gene_start = gene_end = 0

        # Terminal output
        print(
            f"  {record.id:<35}  {len(record.seq):>12,} bp"
            + (f"  [{org}]" if org else "")
        )
        print(
            f"    {n_cds:,} CDS  |  "
            f"tags: {first_tag} \u2192 {last_tag}  |  "
            f"genes: {gene_start:,}\u2013{gene_end:,} bp"
        )
        print()

        # Collect row for TSV
        rows.append(
            {
                "Sequence_ID": record.id,
                "Length_bp": len(record.seq),
                "Organism": org,
                "CDS_Count": n_cds,
                "First_Tag": first_tag,
                "Last_Tag": last_tag,
                "Gene_Start": gene_start,
                "Gene_End": gene_end,
            }
        )

    print(f"[*] {len(rows)} sequence(s) found.\n")
    print(
        "  Next steps:\n"
        "    --list-products --seq SEQ_ID   see what keywords exist on a contig\n"
        "    --seq SEQ_ID                   list all features on a contig\n"
        "    -q 'keyword' --seq SEQ_ID      find specific features\n"
    )

    if args.output:
        _write_tsv(rows, args.output)


def run_list_products(args) -> None:
    """Lists all unique product annotations sorted by frequency.

    Hypothetical proteins are always excluded; use ``-q "hypothetical"``
    to find them explicitly.

    When ``-o`` is specified, writes a TSV with columns: ``Product``, ``Count``.

    Args:
        args: Parsed argument namespace.

    Raises:
        SystemExit: If ``--seq`` is specified but the sequence ID is not found.
    """
    counter: Counter = Counter()
    total_cds = 0
    hypo_count = 0
    seq_scanned = 0

    for record in SeqIO.parse(args.input, "genbank"):
        if args.seq and record.id != args.seq:
            continue
        seq_scanned += 1

        for feature in record.features:
            if feature.type != "CDS":
                continue
            total_cds += 1
            product = get_display_product(feature)

            if "hypothetical" in product.lower():
                hypo_count += 1
                continue

            counter[product] += 1

    if args.seq and seq_scanned == 0:
        sys.exit(
            f"\n[!] Sequence '{args.seq}' not found in '{args.input.name}'.\n"
            f"    Run --list-sequences to see available IDs.\n"
        )

    filtered = {p: c for p, c in counter.items() if c >= args.min_count}
    sorted_products = sorted(filtered.items(), key=lambda x: (-x[1], x[0]))

    scope = f"sequence '{args.seq}'" if args.seq else f"'{args.input.name}'"
    print(f"\n[*] Product annotations in {scope}")
    print(f"[*] Total CDS scanned: {total_cds:,}")
    print(
        f"[*] Excluded {hypo_count:,} 'hypothetical protein' entries "
        f'(use -q "hypothetical" to search for them)'
    )
    if args.min_count > 1:
        hidden = len(counter) - len(filtered)
        print(
            f"[*] Showing products with count >= {args.min_count} "
            f"({hidden:,} rare/singleton entries hidden)"
        )
    print(f"[*] {len(sorted_products):,} unique products:\n")
    print(f"  {'Count':>5}  Product")
    print(f"  {'-----':>5}  {'-' * 65}")

    for product, count in sorted_products:
        display = (product[:68] + "...") if len(product) > 68 else product
        print(f"  {count:>5,}  {display}")

    print(
        f"\n[*] Use any keyword above with -q, e.g.:\n"
        f'    python3 find_gbk_features.py -i {args.input.name} -q "bacteriocin"\n'
    )

    if args.output:
        tsv_rows = [{"Product": p, "Count": c} for p, c in sorted_products]
        _write_tsv(tsv_rows, args.output)


def run_sequence_dump(args) -> None:
    """Lists every CDS feature on a contig when ``--seq`` is used alone.

    The human-readable alternative to scrolling through a raw GBFF file.
    Without ``-o``, prints a formatted table to the terminal. With ``-o``,
    writes a TSV with all base feature columns (plus verbose columns if
    ``-v`` is also used).

    Args:
        args: Parsed argument namespace. ``args.seq`` must be set.

    Raises:
        SystemExit: If the specified sequence ID is not found.
    """
    for record in SeqIO.parse(args.input, "genbank"):
        if record.id != args.seq:
            continue

        all_cds = collect_cds(record)
        org = record.annotations.get("organism", "")

        if args.output:
            # TSV mode: write file, print minimal status to terminal
            print(
                f"\n[*] {record.id}  ({len(record.seq):,} bp"
                + (f", {org})" if org else ")")
                + f"\n[*] {len(all_cds):,} CDS features.",
                file=sys.stderr,
            )
            tsv_rows = [
                _build_feature_tsv_row(f, verbose=args.verbose) for f in all_cds
            ]
            _write_tsv(tsv_rows, args.output)
        else:
            # Text mode: write full table to terminal
            print(
                f"\n[*] Sequence: {record.id}  ({len(record.seq):,} bp)"
                + (f"  [{org}]" if org else ""),
                file=sys.stderr,
            )
            print(f"[*] {len(all_cds):,} CDS features.\n", file=sys.stderr)
            for f in all_cds:
                sys.stdout.write(
                    f"    {f['locus']:<25}  "
                    f"{f['start']:>10,}\u2013{f['end']:<10,}  "
                    f"({f['strand']})  {f['product'][:60]}\n"
                )
            sys.stdout.write("\n")
            sys.stdout.write("=" * 70 + "\n")
            sys.stdout.write("EXTRACTION COMMANDS:\n")
            sys.stdout.write("=" * 70 + "\n\n")
            sf = f" \\\n      --seq {record.id}"
            sys.stdout.write(
                f"  Full contig:\n"
                f"    python3 extract_genome_region.py \\\n"
                f"      -i {args.input.name} \\\n"
                f"      --c1 1 --c2 {len(record.seq)}{sf} \\\n"
                f"      --faa full_contig.faa --fna full_contig.fna\n\n"
                f"  Specific region (fill in coordinates):\n"
                f"    python3 extract_genome_region.py \\\n"
                f"      -i {args.input.name} \\\n"
                f"      --c1 XXXXX --c2 XXXXX{sf} \\\n"
                f"      --faa region.faa --fna region.fna\n\n"
                f"  Tip: use --context LOCUS_TAG --window N to narrow down "
                f"coordinates.\n"
            )
        return

    sys.exit(
        f"\n[!] Sequence '{args.seq}' not found in '{args.input.name}'.\n"
        f"    Run --list-sequences to see available IDs.\n"
    )


def run_context_search(args) -> None:
    """Shows all genes within ``--window`` bp of an anchor locus tag.

    Without ``-o``, prints a neighbourhood table to the terminal. With ``-o``,
    writes a TSV with an ``Is_Anchor`` column plus standard feature columns.

    Args:
        args: Parsed argument namespace. ``args.context`` must be set.

    Raises:
        SystemExit: If the anchor locus tag is not found.
    """
    for record in SeqIO.parse(args.input, "genbank"):
        if args.seq and record.id != args.seq:
            continue

        all_cds = collect_cds(record)
        anchor = next((f for f in all_cds if f["locus"] == args.context), None)
        if anchor is None:
            continue

        anchor_mid = (anchor["start"] + anchor["end"]) // 2
        win_c1 = max(1, anchor_mid - args.window)
        win_c2 = anchor_mid + args.window
        in_window = [f for f in all_cds if win_c1 <= f["start"] <= win_c2]

        print(
            f"\n[*] Context: '{args.context}' "
            f"({anchor['start']:,}\u2013{anchor['end']:,} bp)  "
            f"window \u00b1{args.window:,} bp  "
            f"\u2192 {win_c1:,}\u2013{win_c2:,}\n"
            f"[*] {len(in_window)} feature(s) in window:\n",
            file=sys.stderr,
        )

        if args.output:
            tsv_rows = [
                _build_feature_tsv_row(
                    f,
                    is_anchor=(f["locus"] == args.context),
                    verbose=args.verbose,
                )
                for f in in_window
            ]
            _write_tsv(tsv_rows, args.output)
        else:
            sf = f" \\\n    --seq {record.id}" if args.seq else ""
            sys.stdout.write(
                f"Neighbourhood of '{args.context}' "
                f"(\u00b1{args.window:,} bp, "
                f"{win_c1:,}\u2013{win_c2:,} bp)\n\n"
            )
            for f in in_window:
                mark = ">>> " if f["locus"] == args.context else "    "
                sys.stdout.write(
                    f"{mark}  {f['locus']:<25}  "
                    f"{f['start']:>10,}\u2013{f['end']:<10,}  "
                    f"({f['strand']})  {f['product'][:60]}\n"
                )
            sys.stdout.write(
                f"\n{'=' * 70}\nTO EXTRACT THIS NEIGHBOURHOOD:\n{'=' * 70}\n\n"
                f"  python3 extract_genome_region.py \\\n"
                f"    -i {args.input.name} \\\n"
                f"    --c1 {win_c1} --c2 {win_c2}{sf} \\\n"
                f"    --faa neighbourhood.faa --fna neighbourhood.fna\n\n"
                f"  Tip: adjust --c1/--c2 to the actual cluster boundaries.\n"
                f"  Increase --window if you need more context.\n"
            )
        return

    sys.exit(
        f"\n[!] Locus tag '{args.context}' not found"
        + (f" in sequence '{args.seq}'" if args.seq else "")
        + ".\n    Use --list-sequences to verify IDs, or\n"
        + "    -q with --search-type locus-exact to check the tag exists.\n"
    )


def run_keyword_search(args) -> None:
    """Searches CDS features for one or more query keywords.

    Supports any number of keywords supplied to ``-q``. For a single keyword,
    the output is identical to the original single-keyword format. For
    multiple keywords, the terminal shows results grouped per keyword with
    section headers, and the TSV (when ``-o`` is specified) has a leading
    ``Keyword`` column with rows sorted by keyword order.

    Annotation text search (``--search-type product``) checks /product=,
    /gene_kind=, /gene_functions=, /sec_met_domain=, and /note=.
    Locus searches check only /locus_tag=.

    Args:
        args: Parsed argument namespace. ``args.query`` must be set.
    """
    queries = args.query  # list[str] thanks to nargs='+'

    # Collect results per keyword
    results_by_query: dict[str, list[dict]] = {q: [] for q in queries}

    for record in SeqIO.parse(args.input, "genbank"):
        if args.seq and record.id != args.seq:
            continue

        for feature in record.features:
            if feature.type != "CDS":
                continue

            locus = feature.qualifiers.get("locus_tag", [""])[0]

            for query in queries:
                if args.search_type == "product":
                    matched = query.lower() in get_annotation_text(feature)
                elif args.search_type == "locus":
                    matched = query.lower() in locus.lower()
                else:  # locus-exact
                    matched = query.lower() == locus.lower()

                if matched:
                    results_by_query[query].append(feature_to_dict(feature, record.id))

    total_hits = sum(len(v) for v in results_by_query.values())

    if total_hits == 0:
        msg = (
            f"\n[!] No features found for " f"{', '.join(repr(q) for q in queries)}.\n"
        )
        if args.search_type == "product":
            msg += (
                "\n    HINT: Run --list-products to see available keywords.\n"
                "\n    If searching an antiSMASH region GBK, those files store\n"
                "    annotations in /gene_kind= and /gene_functions= rather\n"
                "    than /product=. This script checks both automatically —\n"
                "    verify your query matches antiSMASH phrasing such as\n"
                "    'biosynthetic' or 'RiPP-like'.\n"
            )
        print(msg, file=sys.stderr)
        return

    if args.output:
        # ── TSV mode ──────────────────────────────────────────────────────────
        print(
            f"\n[*] {total_hits} total hit(s) for "
            f"{', '.join(repr(q) for q in queries)}\n",
            file=sys.stderr,
        )
        tsv_rows = [
            _build_feature_tsv_row(fd, keyword=q, verbose=args.verbose)
            for q in queries
            for fd in results_by_query[q]
        ]
        _write_tsv(tsv_rows, args.output)

    else:
        # ── Terminal text mode ────────────────────────────────────────────────
        if len(queries) == 1:
            # Single keyword: existing compact format
            results = results_by_query[queries[0]]
            print(
                f"\n[*] {len(results)} feature(s) matching {repr(queries[0])}:\n",
                file=sys.stderr,
            )
            write_feature_table(results, args, sys.stdout)
            write_extraction_suggestions(results, args.input, args.seq, sys.stdout)

        else:
            # Multiple keywords: grouped display
            bar = "\u2550" * 70
            print(
                f"\n[*] {len(queries)} keyword(s): "
                f"{', '.join(repr(q) for q in queries)}  "
                f"\u2502  {total_hits} total hit(s)\n",
                file=sys.stderr,
            )
            for i, query in enumerate(queries, 1):
                results = results_by_query[query]
                sys.stdout.write(f"\n{bar}\n")
                sys.stdout.write(
                    f"  KEYWORD {i}/{len(queries)}: {repr(query)}"
                    f"  \u2502  {len(results)} hit(s)\n"
                )
                sys.stdout.write(f"{bar}\n\n")
                if results:
                    write_feature_table(results, args, sys.stdout)
                else:
                    sys.stdout.write(
                        f"  [\u2013] No features found for this keyword.\n\n"
                    )

            # Suggestions span all hits combined
            all_results = [fd for q in queries for fd in results_by_query[q]]
            write_extraction_suggestions(all_results, args.input, args.seq, sys.stdout)


def run_coordinate_search(args) -> None:
    """Shows all CDS features whose full coordinates fall within [c1, c2].

    Without ``-o``, prints a formatted feature table to the terminal. With
    ``-o``, writes a TSV with base feature columns (plus verbose columns if
    ``-v`` is used).

    Args:
        args: Parsed argument namespace. ``args.c1`` and ``args.c2`` must be set.
    """
    results = []

    for record in SeqIO.parse(args.input, "genbank"):
        if args.seq and record.id != args.seq:
            continue

        for feature in record.features:
            if feature.type != "CDS":
                continue

            start = int(feature.location.start) + 1
            end = int(feature.location.end)

            if start >= args.c1 and end <= args.c2:
                results.append(feature_to_dict(feature, record.id))

    if not results:
        print(
            f"\n[!] No features found in range {args.c1:,}\u2013{args.c2:,}.\n"
            "\n    HINT: --c1 and --c2 expect 1-based genomic coordinates.\n"
            "    Use --list-sequences to confirm the correct --seq value.\n",
            file=sys.stderr,
        )
        return

    print(
        f"\n[*] {len(results)} feature(s) in range " f"{args.c1:,}\u2013{args.c2:,}:\n",
        file=sys.stderr,
    )

    if args.output:
        tsv_rows = [_build_feature_tsv_row(r, verbose=args.verbose) for r in results]
        _write_tsv(tsv_rows, args.output)
    else:
        write_feature_table(results, args, sys.stdout)
        write_extraction_suggestions(results, args.input, args.seq, sys.stdout)


# ── Entry point ───────────────────────────────────────────────────────────────


def main() -> None:
    """Parses arguments and dispatches to the appropriate mode runner.

    Mode priority order: discovery modes first (list-sequences,
    list-products), then search modes (context, query, coordinate, seq-dump).
    """
    args = get_args()
    validate_args(args)

    if not args.input.exists():
        sys.exit(f"\n[!] File not found: {args.input}\n")

    if args.list_sequences:
        run_list_sequences(args)
        return

    if args.list_products:
        run_list_products(args)
        return

    if args.context:
        run_context_search(args)
        return

    if args.query:
        run_keyword_search(args)
        return

    if args.c1 is not None:
        run_coordinate_search(args)
        return

    if args.seq:
        run_sequence_dump(args)


if __name__ == "__main__":
    main()
