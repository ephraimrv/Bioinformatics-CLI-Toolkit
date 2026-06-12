"""
FastANI Result Parser

Parses FastANI output files to rank genomic similarity metrics.

This script reads raw FastANI output, sorts hits by nucleotide identity,
and provides both human-readable terminal output and machine-readable
TSV exports. It includes automated classification for species-level
matches based on user-defined ANI and alignment fraction thresholds.

Author: Jan Ephraim R. Vallente (ephrvallente@gmail.com)
Date: 2026-06-07
License: MIT
Reproducibility: Associated with upcoming research (manuscript in preparation).

Example Usage:
    $ python3 fastani_parser.py -i fastani_out.txt -o results.tsv
"""

import sys
import argparse
from pathlib import Path
from typing import Iterator
from utils import base_parser

# Standard species-level match thresholds
ANI_THRESHOLD = 95.0
AF_THRESHOLD = 60.0


def stream_fastani_results(
    file_path: Path,
) -> Iterator[tuple[str, str, float, int, int, float]]:
    """Parses a FastANI output file.

    Yields:
        A tuple of (query, reference, ani, matched, total, alignment_fraction).

    Raises:
        ValueError: If file content is malformed.
    """
    try:
        with open(file_path, "r", encoding="utf-8") as f:
            for line_number, line in enumerate(f, start=1):
                if not line.strip():
                    continue

                parts = line.strip().split()
                if len(parts) < 5:
                    raise ValueError(f"Malformed data on line {line_number}")

                query = Path(parts[0]).name
                reference = Path(parts[1]).name

                try:
                    ani = float(parts[2])
                    matched = int(parts[3])
                    total = int(parts[4])
                except ValueError as e:
                    raise ValueError(
                        f"Invalid numeric data on line {line_number}"
                    ) from e

                af = (matched / total) * 100.0 if total > 0 else 0.0
                yield query, reference, ani, matched, total, af

    except (FileNotFoundError, OSError) as e:
        raise ValueError(f"Could not read FastANI file: {file_path}") from e


def get_args() -> argparse.Namespace:
    """Configures CLI arguments."""
    parser = base_parser("Parse and sort FastANI results.")
    return parser.parse_args()


def main() -> None:
    args = get_args()
    output_path = Path(args.output) if args.output else None

    try:
        results = sorted(
            list(stream_fastani_results(args.input)), key=lambda x: x[2], reverse=True
        )

        if not results:
            sys.exit(f"[!] The file {args.input.name} is empty.")

        # Terminal Output
        print(f"\n--- FastANI Results: {results[0][0]} ---")
        print(
            f"{'Reference Genome':<40} | {'ANI (%)':<8} | {'AF (%)':<8} | {'Fragments'}"
        )
        print("-" * 80)

        for _, ref, ani, match, tot, af in results:
            match_tag = (
                "<-- SPECIES MATCH"
                if (ani >= ANI_THRESHOLD and af >= AF_THRESHOLD)
                else ""
            )
            print(f"{ref:<40} | {ani:>7.2f}  | {af:>7.2f}  | {match}/{tot} {match_tag}")
        print("-" * 80)

        # File Output
        if output_path:
            with open(output_path, "w", encoding="utf-8") as f:
                f.write("Query\tReference\tANI\tAF\tMatched\tTotal\tIs_Match\n")
                for q, r, ani, m, t, af in results:
                    is_match = ani >= ANI_THRESHOLD and af >= AF_THRESHOLD
                    f.write(f"{q}\t{r}\t{ani:.2f}\t{af:.2f}\t{m}\t{t}\t{is_match}\n")
            print(f"[*] Results written to: {output_path.name}")

    except (ValueError, OSError) as e:
        if output_path and output_path.exists():
            output_path.unlink()
        sys.exit(f"\n[!] Pipeline Error: {e}")
    except KeyboardInterrupt:
        if output_path and output_path.exists():
            output_path.unlink()
        sys.exit("\n[!] Pipeline interrupted.")


if __name__ == "__main__":
    main()
