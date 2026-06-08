# Bioinformatics CLI Toolkit

A command-line toolkit for genomic sequence processing, comparative genomics, and regulatory footprinting. 

This repository contains general-purpose sequence manipulation utilities alongside a specialized comparative pipeline designed to analyze biosynthetic gene clusters (BGCs) and their transcriptional regulatory architectures across reference genomes.

## Repository Structure

The toolkit is divided into two primary modules:

* **`c5_research_pipeline/`**: Advanced comparative genomics scripts developed for the structural and regulatory analysis of bacteriocin loci in *Leuconostoc mesenteroides* C5.
* **`core_utilities/`**: General-purpose sequence parsing tools for basic FASTA manipulation, open reading frame (ORF) extraction, and sequence statistics.

## Requirements and Installation

The toolkit requires Python 3.9+ and relies on the `Biopython` library for GenBank file parsing. 

Clone the repository and install the required dependencies:

```bash
git clone [https://github.com/yourusername/Bioinformatics-CLI-Toolkit.git](https://github.com/yourusername/Bioinformatics-CLI-Toolkit.git)
cd Bioinformatics-CLI-Toolkit
pip install -r requirements.txt
```

## Key Tools Overview
### Comparative Research Pipeline (c5_research_pipeline/)
- conserved_annotation_scanner.py: A core proteome profiler that aggregates CDS annotations across multiple GenBank files to identify conserved gene groups based on an exact or minimum genome frequency threshold.

- cross_genome_keyword_scanner.py: Scans multiple genomes for targeted functional keywords and extracts the corresponding locus tags and protein sequences.

- ortholog_extractor.py: Calculates mature, membrane-inserting core peptides from pre-peptide sequences and uses them to identify and extract full-length orthologs from target assemblies.

- upstream_sequence_extractor.py: Extracts exact upstream DNA coordinate blocks (promoter regions) using locus tags to evaluate regulatory structural conservation.

- (Additional tools include BGC exploration, k-mer counting, regulon scanning, and consensus profiling).

### Core Sequence Utilities (core_utilities/)
- gc_genome.py: Calculates per-contig and whole-genome GC content, outputting both summary statistics and a TSV matrix.

- stream_nuc_count.py: Calculates precise nucleotide frequencies (A, C, G, T) and anomalous base counts (N) using memory-efficient lazy evaluation.

- ORF_pipeline.py / ORF_genome_exhaustive.py: Identifies and extracts open reading frames from nucleotide sequences.

- dna_rev_genome.py: Generates the 3'-to-5' reverse complement of an input FASTA.

- DNA_transcripton_genome.py: Transcribes the coding strand of a DNA assembly into mRNA.

### Usage Examples
Most tools share a consistent argparse interface requiring an input directory or file (-i) and an output file (-o).

Example 1: Identifying conserved genes across multiple reference genomes
```
python3 c5_research_pipeline/conserved_annotation_scanner.py -i references/ --min_genomes 3 -o core_proteome.tsv -f
```

Example 2: Extracting a 150bp upstream regulatory region
```
python3 c5_research_pipeline/upstream_sequence_extractor.py -i genome.gbff -t LEUM_RS10400 -u 150 -o upstream.fasta
```
Example 3: Calculating whole-genome GC content
```
python3 core_utilities/gc_genome.py -i assembly.fasta -o gc_results.tsv
```
### Citation
If you use this toolkit in your research, please cite the repository using the provided CITATION.cff metadata or the associated Zenodo DOI.

Vallente, J. E. R. (2026). Bioinformatics-CLI-Toolkit (v1.0.0). Zenodo. 10.5281/zenodo.20586542

### License
This project is licensed under the MIT License.
