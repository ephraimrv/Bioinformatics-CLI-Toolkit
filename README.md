# Bioinformatics CLI Toolkit

A command-line toolkit for genome inspection, comparative genomics, and regulatory
motif discovery, built around GenBank-format assemblies. Originally developed
during comparative genomic analysis of a bacteriocin locus in *Leuconostoc
mesenteroides* C5, the tools are general-purpose and work on any GenBank file
(NCBI GBFF, antiSMASH region files, Prokka/BAKTA-annotated assemblies).

This is solo, ongoing research software associated with unpublished work in
progress. Interfaces and outputs may change between commits.

## Conventions

Every script in this toolkit follows the same house style: `#!/usr/bin/env
python3` shebang, an SPDX license header, a PEP 257 / Google-style module
docstring (`Note:`, `Example(s):`, etc.), and `__author__`/`__email__`/
`__version__` metadata only. All CLI tools use `python3` (never `python`)
and a consistent `-i`/`-o` argparse interface where applicable.

## Repository Structure

```
Bioinformatics-CLI-Toolkit/
├── genomics_toolkit/      # Comparative genomics & regulatory CLI tools
│   ├── utils.py                # Shared parsing/validation helpers (not standalone)
│   └── visuals/                 # Manuscript figure-generation scripts
└── core_utilities/             # Legacy basic-sequence utilities (pending review)
```

## Requirements and Installation

The toolkit requires **Python 3.10+** and Biopython. Two installation paths:

### Option A: Minimal install (pip, recommended for using the toolkit)

```bash
git clone https://github.com/ephraimrv/Bioinformatics-CLI-Toolkit.git
cd Bioinformatics-CLI-Toolkit
pip install -r requirements.txt
```

### Option B: Full reproducible environment (conda, for developing/extending)

If you want the exact research environment (including BLAST, MAFFT, HMMER, and
other system bioinformatics tools), use conda/mamba:

```bash
git clone https://github.com/ephraimrv/Bioinformatics-CLI-Toolkit.git
cd Bioinformatics-CLI-Toolkit
mamba env create -f environment-full-rosalind.yml
mamba activate rosalind
```

## Tool Overview (genomics_toolkit/)

### Discovery & Inspection
- **find_gbk_features.py** — Search/browse GenBank features by keyword, locus
  tag, coordinate range, or context window around an anchor gene. Suggests
  ready-to-run extraction commands for hits, including circular-origin
  wraparound clusters.
- **gbk_scanner.py** — Keyword-based CDS scanner across a single file or a
  directory of GenBank files (NCBI GBFF, antiSMASH, Prokka/BAKTA).
- **contig_gene_profiler.py** — List/extract all genes on a specific contig
  from a BAKTA-annotated assembly, with optional protein sequences.

### Extraction
- **extract_genome_region.py** — Extract a genomic region by coordinate,
  locus-tag range, or whole contig, to GBK/FAA/FNA. Supports circular
  origin-spanning extraction and pipeline integration via scanner TSV output.
- **universal_promoter_extractor.py** — Extract upstream regulatory regions
  for motif discovery. Auto-detects prokaryote (CDS-anchored) vs. eukaryote
  (TSS-anchored via mRNA features, isoform-aware) genomes.
- **gbk_promoter_finder.py** — Single-gene upstream promoter extraction with
  strand handling and an optional quick motif scan.

### Comparative & Ortholog Analysis
- **gbk_ortholog_finder.py** — Pairwise protein homolog detection via
  Smith-Waterman/BLOSUM62 alignment, with eukaryote-aware isoform
  deduplication, bacteriocin-aware coverage modes, and signal-peptide
  trimming. (See script docstring for the homolog-vs-ortholog
  distinction — confirming true orthology requires reciprocal-best-hit
  or phylogenetic analysis.)
  - **exact_match_ortholog_finder.py** — Exact-substring peptide search
  across reference genomes (100% identity required). Faster and
  stricter than alignment-based search — confirms presence of a known,
  exact peptide core. For divergent homolog detection, use
  gbk_ortholog_finder.py instead.
- **conserved_annotation_scanner.py** — Core-proteome profiler: aggregates
  `/product` annotations across genomes and reports conserved genes above a
  genome-frequency threshold.
- **cross_genome_keyword_scanner.py** — Targeted keyword conservation search
  across reference genomes, with TSV + matching FASTA export.
- **protein_presence_scanner.py** — Interactive exact-substring presence/
  absence scanner for a pasted peptide across a directory of reference
  genomes. Bacterial signal-peptide trimming only — use `--raw` for
  non-bacterial query peptides.
- **target_promoter_pipeline.py** — Bridges ortholog detection and promoter
  extraction into one workflow (alignment-based, not exact-match).
- **remote_blast.py** — Automated NCBI remote BLAST runner (blastp/blastn/
  blastx) with TSV output and rate-limit handling.

### Regulatory & Motif Analysis
- **regulon_scanner.py** — Regex/IUPAC operator motif scanner across upstream
  regions, with Benjamini-Hochberg FDR-corrected significance. Prokaryote-only
  (anchors on CDS start, not the transcription start site).
- **motif_discovery.py** — Simplified MEME-style Expectation-Maximization
  motif discovery (bidirectional strand scanning, diversity-aware seeding).
- **comparative_kmer_analyzer.py** — Canonical-k-mer, strand-aware enrichment
  comparison (Log2 fold change) between two genes' upstream regions.
- **alignment_conservation_profiler.py** — Builds a Position Probability
  Matrix and Information Content profile from a pre-built multiple sequence
  alignment.

### Quick Utilities
- **parse_ani.py** — Reformats raw FastANI TSV output into a clean,
  headered TSV for easier review.
- **quick_upstream.py** — Quick single-locus upstream-sequence harvester
  for one GBFF file and one target locus. Locus tag, input file, and
  upstream length are currently hardcoded for fast one-off use. For
  flexible, multi-genome, or comparative upstream extraction, use
  universal_promoter_extractor.py or gbk_promoter_finder.py instead.
  
### Visualization (visuals/)
Manuscript figure-generation scripts (e.g. functional-enrichment plots).
Each is single-purpose and tied to a specific figure, not a general CLI tool.

## Core Utilities (core_utilities/)

Legacy general-purpose sequence scripts (GC content, ORF extraction, reverse
complement, transcription). Pending review — not yet brought up to the
toolkit's current conventions.

## Usage Examples

**Search a genome for features by keyword:**
```bash
python3 genomics_toolkit/gbk_scanner.py -i genome.gbff -q "bacteriocin"
```

**Extract a 150bp upstream regulatory region:**
```bash
python3 genomics_toolkit/gbk_promoter_finder.py -i genome.gbk -l ctg1_50 -u 150 -o promoter.fasta
```

**Find conserved genes across reference genomes (≥3 genomes, TSV + FASTA):**
```bash
python3 genomics_toolkit/conserved_annotation_scanner.py -i references/ --min_genomes 3 -o core_proteome.tsv -f
```

### Citation
Citation metadata is provided via `CITATION.cff` (GitHub renders a "Cite this
repository" button from it automatically). A versioned Zenodo DOI will be
added once the first official release is published — see the project's
release notes for the current citable version.

### License
This project is licensed under the MIT License — see `LICENSE`.
