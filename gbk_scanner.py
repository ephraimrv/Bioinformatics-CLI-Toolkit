# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Jan Ephraim R. Vallente

"""Gene feature scanner and BGC explorer for GenBank files.

Scans one or more .gbk or .gbff files for CDS features matching one or more
keyword queries. Accepts a single file or an entire directory as input.
Designed to work on any GenBank-format genome without restriction on
organism type or annotation source.

Supported input formats:
    - NCBI GBFF files (RefSeq or GenBank division, prokaryote or eukaryote)
    - antiSMASH v8 region .gbk files
    - Prokka-annotated .gbk files
    - Any standard GenBank-format file

Multiple keywords can be supplied in one run. The TSV output groups results
by keyword in the order the keywords were given on the command line, so the
leftmost ``Keyword`` column runs key1 rows first, key2 rows next, and so on.
A CDS feature that matches more than one keyword appears once per matching
keyword so that no hit is missed.

Note:
    This script is part of ongoing research and is associated with an upcoming
    publication. Correct attribution is requested when used in derivative works.
    Released under the MIT License. See LICENSE in the repository root.

Example:
    Search a single NCBI GBFF file for bacteriocin genes::

        python3 bgc_explorer.py -i genome.gbff -q "bacteriocin"

    Search a directory with three keywords and export results::

        python3 bgc_explorer.py -i ./genomes/ \\
            -q "bacteriocin" "ABC transporter" "response regulator" \\
            -o results.tsv

    Search an antiSMASH region file for biosynthetic roles::

        python3 bgc_explorer.py \\
            -i C5_gnlProkkaNHJNNNGJ_1.region001.gbk \\
            -q "biosynthetic" "RiPP"
"""

__author__ = "Jan Ephraim R. Vallente"
__email__ = "ephrvallente@gmail.com"
__version__ = "2.0.0"

import re
import sys
import argparse
from collections import defaultdict
from collections.abc import Iterator
from pathlib import Path
from typing import Any

try:
    from Bio import SeqIO
except ImportError:
    sys.exit(
        "ERROR: Biopython is required but not installed.\n"
        "       Install it with: pip install biopython"
    )

# ── Constants ─────────────────────────────────────────────────────────────────

# Captures the meaningful suffix of antiSMASH /gene_functions= qualifiers.
# Input : "biosynthetic (rule-based-clusters) RiPP-like: Bacteriocin_IIc"
# Output: "RiPP-like: Bacteriocin_IIc"
_FUNC_PATTERN = re.compile(r"\)\s*([^)]+)$")

_GBK_EXTENSIONS: frozenset[str] = frozenset({".gbk", ".gbff", ".gb"})


# ── File iteration ─────────────────────────────────────────────────────────────


def iter_genbank_files(input_path: Path) -> Iterator[Path]:
    """Yields every GenBank file reachable from ``input_path``.

    Single-file inputs are yielded directly after an extension check.
    Directory inputs yield every .gbk/.gbff file in the directory (not
    recursive), sorted alphabetically for reproducible run order.

    Args:
        input_path: A ``Path`` to a .gbk/.gbff file or a directory.

    Yields:
        ``Path`` objects for each valid GenBank file.

    Raises:
        SystemExit: If the path does not exist, is a file with an
            unsupported extension, or is a directory with no .gbk/.gbff files.
    """
    if not input_path.exists():
        sys.exit(f"\n[!] Input path not found: {input_path}\n")

    if input_path.is_file():
        if input_path.suffix.lower() not in _GBK_EXTENSIONS:
            sys.exit(
                f"\n[!] '{input_path.name}' is not a recognised GenBank file.\n"
                f"    Supported extensions: {', '.join(sorted(_GBK_EXTENSIONS))}\n"
            )
        yield input_path
        return

    if input_path.is_dir():
        found = sorted(
            p
            for ext in _GBK_EXTENSIONS
            for p in input_path.glob(f"*{ext}")
            if p.is_file()
        )
        if not found:
            sys.exit(
                f"\n[!] No GenBank files found in '{input_path}'.\n"
                f"    Supported extensions: {', '.join(sorted(_GBK_EXTENSIONS))}\n"
            )
        yield from found


# ── Annotation helpers ────────────────────────────────────────────────────────


def _annotation_text(feature) -> str:
    """Builds a single lowercase searchable string from all annotation qualifiers.

    Checks /product=, /gene_kind=, /gene_functions=, /sec_met_domain=, and
    /note=, covering both NCBI GBFF (/product=) and antiSMASH region GBK
    (/gene_kind=, /gene_functions=) in one call.

    Args:
        feature: A BioPython ``SeqFeature`` object.

    Returns:
        A lowercase string joining all annotation qualifier values.
    """
    parts = [
        feature.qualifiers.get("product", [""])[0],
        feature.qualifiers.get("gene_kind", [""])[0],
        " ".join(feature.qualifiers.get("gene_functions", [])),
        " ".join(feature.qualifiers.get("sec_met_domain", [])),
        feature.qualifiers.get("note", [""])[0],
    ]
    return " ".join(p for p in parts if p).lower()


def _display_product(feature) -> str:
    """Returns the best available product description for display.

    Prefers /product=. Falls back to /gene_kind= combined with the first
    /gene_functions= entry for antiSMASH region GBK files where /product=
    is absent on CDS features.

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


def _role_str(feature) -> str:
    """Builds a biosynthetic role string from antiSMASH-specific qualifiers.

    Combines /gene_kind= with the meaningful suffix of /gene_functions= to
    produce strings such as ``"biosynthetic [RiPP-like: Bacteriocin_IIc]"``.
    Returns ``"unassigned"`` for NCBI GBFF features that lack antiSMASH
    qualifiers.

    Args:
        feature: A BioPython ``SeqFeature`` object.

    Returns:
        A role string for the terminal table and TSV ``Role_in_BGC`` column.
    """
    gene_kind = feature.qualifiers.get("gene_kind", ["unassigned"])[0]
    gene_funcs = feature.qualifiers.get("gene_functions", [])

    if not gene_funcs:
        return gene_kind

    match = _FUNC_PATTERN.search(gene_funcs[0])
    detail = match.group(1).strip() if match else gene_funcs[0].strip()
    return f"{gene_kind} [{detail}]" if detail else gene_kind


# ── Core scanner ──────────────────────────────────────────────────────────────


def scan_file(
    gbk_path: Path,
    queries: list[str],
) -> list[dict[str, Any]]:
    """Scans one GenBank file for CDS features matching any query keyword.

    Both prokaryotic and eukaryotic GenBank files are supported. For features
    with compound (multi-exon) locations, the reported start and end span the
    full extent of the coding region.

    A feature that matches multiple queries produces one hit dict per query,
    ensuring complete capture under every matching keyword.

    Args:
        gbk_path: Path to the file to scan.
        queries:  Query strings to search for. Each should already be
                  lowercased; matching is case-insensitive by construction.

    Returns:
        A list of hit dicts. Each dict contains the following keys:
        ``keyword``, ``source_file``, ``sequence_id``, ``organism``,
        ``cluster_type``, ``locus_tag``, ``gene_name``, ``start``, ``end``,
        ``strand``, ``product``, ``protein_id``, ``aa_length``, ``gene_kind``,
        ``role``.

    Raises:
        ValueError: If BioPython fails to parse the file.
    """
    hits: list[dict[str, Any]] = []

    try:
        for record in SeqIO.parse(gbk_path, "genbank"):
            organism = record.annotations.get("organism", "")
            cluster_type = ""

            # Detect antiSMASH cluster type from the 'region' feature.
            # This feature is absent in NCBI GBFF files; cluster_type
            # remains an empty string in that case.
            for feat in record.features:
                if feat.type == "region":
                    cluster_type = feat.qualifiers.get("product", [""])[0]
                    break

            for feature in record.features:
                if feature.type != "CDS":
                    continue

                text = _annotation_text(feature)
                translation = feature.qualifiers.get("translation", [""])[0]

                for query in queries:
                    if query not in text:
                        continue

                    hits.append(
                        {
                            "keyword": query,
                            "source_file": gbk_path.name,
                            "sequence_id": record.id,
                            "organism": organism,
                            "cluster_type": cluster_type,
                            "locus_tag": feature.qualifiers.get("locus_tag", [""])[0],
                            "gene_name": feature.qualifiers.get("gene", [""])[0],
                            "start": int(feature.location.start) + 1,
                            "end": int(feature.location.end),
                            "strand": "+" if feature.location.strand == 1 else "-",
                            "product": _display_product(feature),
                            "protein_id": feature.qualifiers.get("protein_id", [""])[0],
                            "aa_length": len(translation) if translation else 0,
                            "gene_kind": feature.qualifiers.get("gene_kind", [""])[0],
                            "role": _role_str(feature),
                            "protein_sequence": translation,
                        }
                    )

    except (FileNotFoundError, OSError, UnicodeDecodeError, ValueError) as exc:
        raise ValueError(f"Failed to parse '{gbk_path.name}': {exc}") from exc

    return hits


def scan_inputs(
    input_path: Path,
    queries: list[str],
) -> tuple[list[dict[str, Any]], int]:
    """Scans all GenBank files at ``input_path`` for all query keywords.

    Iterates over files from ``iter_genbank_files``, calls ``scan_file``
    on each, and prints a per-file progress line to stderr. Files that
    fail to parse are skipped with a warning rather than aborting the run.

    Args:
        input_path: Path to a single .gbk/.gbff file or a directory.
        queries:    Lowercased query strings.

    Returns:
        A tuple of ``(all_hits, files_scanned)`` where ``all_hits`` is a
        flat list of every hit dict from every file and ``files_scanned``
        is the count of files successfully processed.
    """
    all_hits: list[dict[str, Any]] = []
    files_scanned = 0

    for gbk_path in iter_genbank_files(input_path):
        print(f"  [\u2192] {gbk_path.name}", file=sys.stderr)
        try:
            hits = scan_file(gbk_path, queries)
            all_hits.extend(hits)
        except ValueError as exc:
            print(f"  [!] Skipped \u2014 {exc}", file=sys.stderr)
        files_scanned += 1

    return all_hits, files_scanned


# ── Terminal display ──────────────────────────────────────────────────────────


def _row_marker(hit: dict[str, Any]) -> str:
    """Returns the 4-character prefix used to visually mark a terminal row.

    ``>>>`` marks antiSMASH biosynthetic core genes.
    `` -> `` marks transport genes.
    ``    `` is used for all other features, including all NCBI GBFF hits
    (which carry no gene_kind qualifier).

    Args:
        hit: A hit dict from ``scan_file``.

    Returns:
        A 4-character string.
    """
    role = hit["role"].lower()
    if "biosynthetic" in role:
        return ">>> "
    if "transport" in role:
        return " -> "
    return "    "


def print_terminal_report(
    all_hits: list[dict[str, Any]],
    queries: list[str],
    files_scanned: int,
) -> None:
    """Prints a query-grouped terminal report of all matching features.

    For each query keyword, displays a sub-table of hits grouped by source
    file, then shows a global summary across all queries.

    Args:
        all_hits:      All hit dicts from ``scan_inputs``.
        queries:       Original query strings in the order they were given.
        files_scanned: Total number of files that were processed.
    """
    WIDE = "\u2550" * 110
    THIN = "\u2500" * 110

    for q_idx, query in enumerate(queries, 1):
        q_hits = [h for h in all_hits if h["keyword"] == query]

        print(f"\n{WIDE}")
        print(
            f'  QUERY {q_idx}/{len(queries)}: "{query}"'
            f"  \u2502  {len(q_hits)} hit(s)  \u2502  "
            f"{files_scanned} file(s) scanned"
        )
        print(WIDE)

        if not q_hits:
            print("  [\u2013] No features found.\n")
            continue

        # Group hits by source file, preserving scan order
        by_file: dict[str, list[dict]] = defaultdict(list)
        for h in q_hits:
            by_file[h["source_file"]].append(h)

        for filename, file_hits in by_file.items():
            s = file_hits[0]
            org_str = f"  [{s['organism']}]" if s["organism"] else ""
            bgc_str = f"  Cluster: {s['cluster_type']}" if s["cluster_type"] else ""
            print(
                f"\n  File : {filename}  \u2502  Seq: {s['sequence_id']}{org_str}{bgc_str}"
            )
            print(f"  {THIN}")
            print(
                f"  {'':4}{'Locus Tag':<22}{'Coordinates':<22}"
                f"{'Str':5}{'Product':<42}Protein ID"
            )
            print(f"  {THIN}")

            for h in sorted(file_hits, key=lambda x: x["start"]):
                mark = _row_marker(h)
                coord = f"{h['start']:,}\u2013{h['end']:,}"
                prod = (
                    (h["product"][:38] + "...")
                    if len(h["product"]) > 41
                    else h["product"]
                )
                print(
                    f"  {mark}{h['locus_tag']:<22}{coord:<22}"
                    f"{h['strand']:<5}{prod:<42}{h['protein_id']}"
                )

            print(f"\n  {len(file_hits)} hit(s) in this file.")

    # ── Global summary ────────────────────────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("  SUMMARY")
    print(f"{'=' * 60}")
    print(f"  Files scanned : {files_scanned}")
    print(f"  Total hits    : {len(all_hits)}\n")

    col_w = max(len(q) for q in queries) + 2
    for query in queries:
        q_hits = [h for h in all_hits if h["keyword"] == query]
        files_hit = len({h["source_file"] for h in q_hits})
        print(
            f"  {repr(query):<{col_w}}  "
            f"{len(q_hits):>4} hit(s)  in  "
            f"{files_hit}/{files_scanned} file(s)"
        )
    print()
    print(
        "  [i] Protein sequences are included in the TSV output\n"
        "      (Protein_Sequence column, rightmost). Use -o FILE.tsv to save.\n"
    )


# ── TSV export ────────────────────────────────────────────────────────────────

# Column names in left-to-right display order.
# Keyword is first so the TSV sorts visually by query group.
_TSV_COLUMNS: list[str] = [
    "Keyword",
    "Source_File",
    "Sequence_ID",
    "Organism",
    "Cluster_Type",
    "Locus_Tag",
    "Gene",
    "Start",
    "End",
    "Strand",
    "Product",
    "Protein_ID",
    "Length_aa",
    "Gene_Kind",
    "Role_in_BGC",
    "Protein_Sequence",
]

# Maps TSV column name → hit dict key
_COL_TO_KEY: dict[str, str] = {
    "Keyword": "keyword",
    "Source_File": "source_file",
    "Sequence_ID": "sequence_id",
    "Organism": "organism",
    "Cluster_Type": "cluster_type",
    "Locus_Tag": "locus_tag",
    "Gene": "gene_name",
    "Start": "start",
    "End": "end",
    "Strand": "strand",
    "Product": "product",
    "Protein_ID": "protein_id",
    "Length_aa": "aa_length",
    "Gene_Kind": "gene_kind",
    "Role_in_BGC": "role",
    "Protein_Sequence": "protein_sequence",
}


def write_tsv(
    all_hits: list[dict[str, Any]],
    queries: list[str],
    output_path: Path,
) -> None:
    """Writes all hits to a TSV file grouped by query keyword order.

    Rows are sorted first by the position of each keyword in ``queries``
    (so all hits for the first keyword appear before hits for the second,
    etc.), then by source file name, then by start coordinate within each
    file. This ensures the leftmost ``Keyword`` column shows query keywords
    in the order they were specified on the command line.

    The ``Cluster_Type`` and ``Gene_Kind`` columns are populated only for
    antiSMASH region GBK hits. They are empty strings for NCBI GBFF hits,
    which carry no antiSMASH-specific qualifiers.

    Args:
        all_hits:    All hit dicts from ``scan_inputs``.
        queries:     Original query strings in the order they were given.
        output_path: Destination file path for the TSV.

    Raises:
        SystemExit: If the output file cannot be written.
    """
    query_rank = {q: i for i, q in enumerate(queries)}

    sorted_hits = sorted(
        all_hits,
        key=lambda h: (query_rank[h["keyword"]], h["source_file"], h["start"]),
    )

    try:
        with open(output_path, "w", encoding="utf-8-sig", newline="") as out:
            out.write("\t".join(_TSV_COLUMNS) + "\n")
            for hit in sorted_hits:
                out.write(
                    "\t".join(str(hit[_COL_TO_KEY[col]]) for col in _TSV_COLUMNS) + "\n"
                )
        print(f"[*] TSV written \u2192 {output_path.resolve()}")

    except OSError as exc:
        sys.exit(f"\n[!] Could not write '{output_path}': {exc}\n")


# ── CLI ───────────────────────────────────────────────────────────────────────


def get_args() -> argparse.Namespace:
    """Configures and returns the CLI argument parser.

    Returns:
        Parsed argument namespace with ``input``, ``query``, and ``output``.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Scan GenBank files for genes matching one or more keywords. "
            "Accepts a single file or a directory. Works on NCBI GBFF, "
            "antiSMASH GBK, Prokka GBK, and eukaryotic assemblies."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help=(
            "Path to a single .gbk or .gbff file, or a directory of such files. "
            "All .gbk/.gbff files found in a directory are scanned."
        ),
    )
    parser.add_argument(
        "-q",
        "--query",
        type=str,
        nargs="+",
        required=True,
        metavar="KEYWORD",
        help=(
            "One or more search keywords. Matched case-insensitively against "
            "/product=, /gene_kind=, /gene_functions=, /sec_met_domain=, and "
            "/note= qualifiers. Enclose multi-word keywords in quotes: "
            '-q "ABC transporter" "response regulator". '
            "TSV output groups results in the keyword order given here."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Path to save results as a TSV file (optional).",
    )

    return parser.parse_args()


def main() -> None:
    """CLI entry point: parses arguments, runs the scan, and reports results."""
    args = get_args()

    # Deduplicate queries while preserving the order given by the user.
    # Matching is always case-insensitive; original case is preserved for display.
    queries: list[str] = list(dict.fromkeys(args.query))
    queries_lower = [q.lower() for q in queries]

    print(
        f"\n[*] BGC Explorer v{__version__}\n"
        f"[*] Input  : {args.input}\n"
        f"[*] Queries: {', '.join(repr(q) for q in queries)}\n"
        f"[*] Scanning...\n",
        file=sys.stderr,
    )

    try:
        all_hits, files_scanned = scan_inputs(args.input, queries_lower)
    except KeyboardInterrupt:
        sys.exit("\n[!] Interrupted by user.\n")

    if not all_hits:
        print(
            f"\n[!] No matches found for "
            f"{', '.join(repr(q) for q in queries)} "
            f"across {files_scanned} file(s).\n"
        )
        return

    # Restore original-case query labels in hit dicts before display/export.
    # scan_file stored lowercase keys; replace them with the display versions.
    lower_to_display = dict(zip(queries_lower, queries))
    for hit in all_hits:
        hit["keyword"] = lower_to_display[hit["keyword"]]

    print_terminal_report(all_hits, queries, files_scanned)

    if args.output:
        write_tsv(all_hits, queries, args.output)


if __name__ == "__main__":
    main()
