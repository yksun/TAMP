<img align="right" src="/docs/taco-icon.png">

# TACO

**Telomere-Aware Contig Optimization**

TACO is a telomere-aware genome assembly workflow for benchmarking, comparing, and refining assemblies into improved chromosome-scale candidates. The pipeline was developed for small eukaryotic genomes, with a focus on fungal genomes and PacBio HiFi reads. It runs multiple assemblers, standardizes their outputs, evaluates assembly quality, detects telomere-supported contigs, and refines a selected backbone assembly by preserving protected telomere-supported contigs while reducing redundant non-telomeric sequence.

TACO was developed at the **Grainger Bioinformatics Center, Field Museum of Natural History**.

![Latest Version](https://img.shields.io/github/v/tag/yksun/TACO?label=Latest%20Version)
![Last Commit](https://img.shields.io/github/last-commit/yksun/TACO)
![Issues](https://img.shields.io/github/issues/yksun/TACO)
![BioConda](https://img.shields.io/badge/BioConda-coming_soon-lightgrey)
![Docker](https://img.shields.io/badge/Docker-coming_soon-lightgrey)
![Singularity](https://img.shields.io/badge/Singularity-coming_soon-lightgrey)

---

<a id="table-of-contents"></a>
## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Workflow](#workflow)
- [Installation](#installation)
  - [Conda Environment](#conda-environment)
  - [Manual Installation (Recommended for Canu and Redundans)](#manual-installation)
- [Usage](#usage)
  - [Basic Command](#basic-command)
  - [Parameters](#parameters)
  - [Example Run](#example-run)
- [Pipeline Steps](#pipeline-steps)
- [Input Data Recommendations](#input-data-recommendations)
- [Dependencies](#dependencies)
- [Output Structure](#output-structure)
- [Benchmark Logging](#benchmark-logging)
- [Troubleshooting](#troubleshooting)
- [Citation and Archiving](#citation-and-archiving)
- [Acknowledgements](#acknowledgements)

---

<a id="overview"></a>
## Overview

Genome assemblers often produce different results from the same long-read dataset. One assembler may recover longer contigs, another may preserve more complete chromosome ends, and another may provide a better balance of completeness, contiguity, and redundancy. TACO was built to make these comparisons systematic, interpretable, and reproducible.

Rather than focusing only on recovering strict telomere-to-telomere (T2T) contigs, TACO is designed as a multi-assembler decision and refinement workflow. It runs several assemblers, normalizes their outputs, summarizes assembly quality metrics, detects telomere-supported contigs, and builds a unified comparison table so users can decide which assembler performs best for their dataset.

TACO then refines the selected backbone assembly by prioritizing protected telomere-supported contigs, including strict T2T contigs when present, while reducing redundant non-telomeric sequence. This makes the pipeline useful both for assembler comparison and for improving chromosome-end representation in the final assembly.

<a id="features"></a>
## Features

- Run multiple long-read assemblers from one workflow.
- Standardize assembly outputs for direct comparison across assemblers.
- Benchmark assemblies with QUAST, BUSCO, and telomere-support summaries.
- Build a unified `assembly_info.csv` table to support assembler selection.
- Separate strict T2T contigs from single-end telomeric contigs.
- Preserve protected telomere-supported contigs during final backbone refinement.
- Reduce redundant non-telomeric contigs while retaining chromosome-end support.
- Generate machine-readable benchmark logs and software version records for reproducibility.

<a id="workflow"></a>
## Workflow

TACO follows this high-level order:

1. Run one or more assemblers on the same long-read dataset.
2. Normalize assembly outputs and contig names for comparison.
3. Detect telomere-supported contigs for each assembly.
4. Evaluate assemblies with QUAST and BUSCO.
5. Build a combined summary table across assemblies to support assembler choice.
6. Select a preferred backbone assembly based on the comparison results.
7. Build a telomere-supported contig pool from all assemblies.
8. Refine the selected backbone by replacing redundant backbone contigs with protected telomere-supported contigs.
9. Run final telomere checks, final BUSCO/QUAST, and summary reporting.

A schematic workflow figure should emphasize three concepts: multi-assembler benchmarking, telomere-supported contig classification, and backbone refinement by protected contig replacement. Place the updated TACO figure in `docs/taco-workflow.png` or `docs/taco-icon.png` and link it here when ready.

<a id="installation"></a>
## Installation

TACO uses a primary Conda or Micromamba environment for most tools. At present, **Redundans is installed manually** because of dependency conflicts with the modern Conda stack used by other tools. **Canu may also be installed manually** if your local package resolution selects an unstable development build or if your environment has persistent solver conflicts.

<a id="conda-environment"></a>
### Conda Environment

Create the main environment from `tamp-env.yml`:

```bash
micromamba create -n tamp -f tamp-env.yml
micromamba activate tamp
```

Or with Conda:

```bash
conda env create -n tamp -f tamp-env.yml
conda activate tamp
```

If you use the provided installer:

```bash
bash install_tamp.sh
```

The installer is intended to place TACO under `~/opt/TACO` by default and create launchers in `~/opt/bin`.

<a id="manual-installation"></a>
### Manual Installation (Recommended for Canu and Redundans)

#### Redundans

Clone and install Redundans manually inside the TACO installation directory:

```bash
cd ~/opt/TACO
git clone --recursive https://github.com/Gabaldonlab/redundans.git
cd redundans
bash INSTALL.sh
```

If `INSTALL.sh` is not suitable on your system, use the upstream compile script:

```bash
bash bin/.compile.sh
```

Then make sure `redundans.py` is available on your `PATH`, or create a launcher script in `~/opt/bin`.

#### Canu

If `canu -version` reports a development build such as `master +XX changes`, do not use it for production runs. Install a stable Canu release instead. If the Conda solver cannot provide a stable build in your environment, compile Canu manually:

```bash
cd ~/opt
curl -LRO https://github.com/marbl/canu/releases/download/v2.3/canu-2.3.Linux-amd64.tar.xz
tar -xJf canu-2.3.*.tar.xz
```

Then expose the stable binary:

```bash
export PATH="$HOME/opt/canu-2.3/build/bin:$PATH"
```

<a id="usage"></a>
## Usage

<a id="basic-command"></a>
### Basic Command

Run TACO in a dedicated working directory and provide the input FASTQ file by absolute path:

```bash
TACO.sh -g 12m -t 16 --fastq /absolute/path/to/reads.fastq -m TGTG
```

<a id="parameters"></a>
### Parameters

| Parameter | Meaning |
|---|---|
| `-g`, `--genomesize` | Estimated haploid genome size, such as `12m`, `40m`, or `2g`. |
| `-t`, `--threads` | Number of CPU threads to use. |
| `--fastq` | Input FASTQ file. Use an absolute path if the working directory is separate from the data directory. |
| `-m`, `--motif` | Telomere motif or telomere seed motif used by the pipeline. |
| `-s`, `--steps` | Run only selected steps, for example `1,3-5`. |
| `--fasta` | External pre-assembled FASTA to include in the comparison. |
| `--busco` | Run BUSCO and optionally set the lineage dataset. |
| `--choose` | Prompt for the assembly to use in the final merge. |

<a id="example-run"></a>
### Example Run

For *Saccharomyces cerevisiae* PacBio HiFi data:

```bash
mkdir -p ~/Storage/projects/TACO/SRR13577847
cd ~/Storage/projects/TACO/SRR13577847

TACO.sh \
  -g 12m \
  -t 16 \
  --fastq ~/Storage/projects/TACO/data/SRR13577847.fastq \
  -m TGTG
```

For fungi with vertebrate-like telomere repeats, such as *Neurospora crassa*, use a different motif as appropriate, for example `TTAGGG`.

<a id="pipeline-steps"></a>
## Pipeline Steps

TACO currently implements the following major steps:

1. HiCanu assembly
2. NextDenovo assembly
3. Peregrine assembly
4. IPA assembly
5. Flye assembly
6. Hifiasm assembly
7. Copy and normalize all assemblies
8. BUSCO on all assemblies
9. Telomere contig detection and telomere metrics
10. Build a telomere-supported contig pool across assemblies
11. QUAST for assembler results
12. Final assembly refinement with telomere-supported contig replacement
13. BUSCO analysis of final assembly
14. Telomere analysis of final assembly
15. QUAST analysis of final assembly
16. Final comparison report
17. Cleanup into structured output folders

<a id="input-data-recommendations"></a>
## Input Data Recommendations

TACO is best suited to long-read data, especially PacBio HiFi reads, for small eukaryotic genomes. For publication-quality benchmarking:

- Use public datasets with stable accession numbers.
- Report the source of all third-party datasets.
- Use machine-readable input and output files.
- Archive the exact code release used for the analysis.

For fungal benchmarking, it is helpful to include a mix of:

- a small yeast genome,
- a medium filamentous fungal genome,
- and, if possible, a more repeat-rich genome.

<a id="dependencies"></a>
## Dependencies

### Core pipeline environment

The main environment is expected to provide:

- python
- canu
- nextdenovo
- peregrine-2021
- pbipa
- flye
- hifiasm
- minimap2
- bwa
- samtools
- seqtk
- miniasm
- busco
- quast
- funannotate
- quickmerge
- pigz
- curl
- wget
- git
- make
- gxx_linux-64
- parallel
- numpy
- pandas

### Additional tools used outside the main environment

#### Redundans requirements

Redundans can use several external tools. Not all of them are required for every TACO run, but the upstream software lists the following resources:

- Python
- Platanus
- Miniasm
- Minimap2
- LAST
- BWA
- SNAP aligner
- SSPACE3
- GapCloser
- GFAstats
- Meryl
- Merqury
- k8
- R
- ggplot2
- scales
- argparse for R

In practice, the TACO workflow mainly depends on the subset needed for the specific Redundans operations invoked during the final merge.

<a id="output-structure"></a>
## Output Structure

Typical output files and folders include:

- `assemblies/` — normalized assembly FASTA files and summary tables
- `logs/` — per-step log files
- `benchmark_logs/` — machine-readable benchmark and timing outputs
- `version.txt` — software versions recorded for the run
- final merged assembly files
- telomere summary files
- QUAST and BUSCO result directories

<a id="benchmark-logging"></a>
## Benchmark Logging

TACO writes benchmark-oriented logs to support reproducible testing and publication reporting. These can include:

- `benchmark_logs/run_metadata.tsv`
- `benchmark_logs/step_benchmark.tsv`
- `benchmark_logs/run_summary.txt`

These files are intended to make it easier to report runtime, software versions, and step status in benchmark tables and supplementary material.

<a id="troubleshooting"></a>
## Troubleshooting

### Canu reports `master +XX changes`

You are using a development build, not a stable release. Replace it with a stable Canu installation before benchmarking.

### Micromamba activates the wrong environment root

Set the correct root prefix in your shell configuration, then reopen the shell:

```bash
export MAMBA_ROOT_PREFIX="$HOME/.local/share/mamba"
```

### Telomere motif appears incorrect

Do not assume the same motif for all fungi. Choose a motif or seed motif appropriate for the target organism.

### `TACO.sh: command not found`

Add the install directory or launcher directory to your `PATH`, or run TACO using the full path.

<a id="citation-and-archiving"></a>
## Citation and Archiving

If you use TACO in a publication, please cite the software and archive the exact release used for the analysis in a persistent public repository such as Zenodo. Citing the archived release helps ensure that the software version used in the study remains accessible and reproducible.

TACO was developed at the Grainger Bioinformatics Center, Field Museum of Natural History, Chicago, Illinois, USA.

<a id="acknowledgements"></a>
## Acknowledgements

TACO was developed in the context of genome assembly benchmarking and telomere-aware assembly refinement workflows at the Grainger Bioinformatics Center, Field Museum of Natural History.
