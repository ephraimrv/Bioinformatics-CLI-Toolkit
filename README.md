# Bioinformatics CLI Toolkit

A set of command-line tools for processing and analyzing genomic sequence data. 

This toolkit is written entirely in standard Python. It is designed to process multi-FASTA files using lazy evaluation, allowing it to handle large sequence assemblies without requiring third-party libraries like `Biopython` or `Pandas`.

## Features
* **Standard Library Only:** Runs on standard Python 3.9+ without additional installations.
* **Controlled Memory Usage:** Uses line-by-line generators to parse FASTA files, preventing large sequences from overwhelming system RAM.
* **Consistent CLI:** All tools share a standard `argparse` interface for predictable input/output routing.
* **Data Validation:** Includes checks for ambiguous bases (e.g., `N`) and malformed FASTA headers, with clean error handling.

## The Toolkit

* **`utils.py`**: The core utility module containing the lazy FASTA parser, standardized CLI arguments, and biological reference dictionaries.
* **`stream_nuc_count.py`**: Calculates A, C, G, T, and anomalous (`N`) base counts for each sequence.
* **`gc_genome.py`**: Calculates per-contig GC percentage and the overall whole-genome GC content. Outputs a TSV matrix and a summary text file.
* **`dna_rev_genome.py`**: Reads a DNA genome and outputs the 3'-to-5' reverse complement.
* **`DNA_transcripton_genome.py`**: Transcribes the coding strand of a DNA assembly into mRNA.

## Usage Example

Each tool requires an input file path (`-i`) and an output file path (`-o`). 

**Running the GC Content Analyzer:**
```bash
python3 gc_genome.py -i /path/to/assembly.fasta -o /path/to/results.tsv
