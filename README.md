<img align="right" src="/docs/tamp-icon.png">

# TAMP

**Telomere Assembly Merge Pipeline**

TAMP is a telomere-aware genome assembly workflow for benchmarking, comparing, and merging assemblies into improved chromosome-scale candidates. The pipeline was developed for small eukaryotic genomes, with a focus on fungal genomes and PacBio HiFi reads. It runs multiple assemblers, standardizes their outputs, evaluates assembly quality, detects telomeric contigs, and builds a final merged assembly while protecting telomere-to-telomere (T2T) candidates.

TAMP was developed at the **Grainger Bioinformatics Center, Field Museum of Natural History**.

![Latest Version](https://img.shields.io/github/v/tag/yksun/TAMP?label=Latest%20Version)
![Last Commit](https://img.shields.io/github/last-commit/yksun/TAMP)
![Issues](https://img.shields.io/github/issues/yksun/TAMP)
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

Genome assemblers often produce different results from the same long-read dataset. One assembler may recover longer contigs, another may preserve more complete chromosome ends, and another may recover better consensus in repeat-rich regions. TAMP was built to make these comparisons systematic and reproducible.

The pipeline runs several assemblers, normalizes their outputs, measures assembly quality, detects telomeric contigs, and combines the results into a unified comparison table. It then supports a final merge step that preserves telomere-to-telomere contigs and reduces redundancy among the remaining contigs.

TAMP is especially useful when the goal is not only to maximize contiguity, but also to recover biologically meaningful chromosome-end structure.

<a id="features"></a>
## Features

- Run multiple assemblers from one workflow.
- Standardize output FASTA files for direct comparison.
- Detect telomeric contigs and summarize telomere support.
- Evaluate assemblies with QUAST and BUSCO.
- Build a unified `assembly_info.csv` table before final merging.
- Preserve T2T contigs during final merge.
- Generate machine-readable benchmark logs for reproducibility.
- Record software versions and step-by-step logs for each run.

<a id="workflow"></a>
## Workflow

TAMP follows this high-level order:

1. Run one or more assemblers.
2. Normalize assembly outputs and contig names.
3. Detect telomeric contigs for each assembly.
4. Evaluate assemblies with QUAST and BUSCO.
5. Build a combined summary table across assemblies.
6. Select a preferred assembly for the final merge.
7. Preserve T2T contigs and merge non-T2T contigs.
8. Run final telomere checks, final BUSCO/QUAST, and summary reporting.

A schematic workflow figure can be added in `docs/` and linked here when ready.

<a id="installation"></a>
## Installation

TAMP uses a primary Conda or Micromamba environment for most tools. At present, **Redundans is installed manually** because of dependency conflicts with the modern Conda stack used by other tools. **Canu may also be installed manually** if your local package resolution selects an unstable development build or if your environment has persistent solver conflicts.

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

The installer is intended to place TAMP under `~/opt/TAMP` by default and create launchers in `~/opt/bin`.

<a id="manual-installation"></a>
### Manual Installation (Recommended for Canu and Redundans)

#### Redundans

Clone and install Redundans manually inside the TAMP installation directory:

```bash
cd ~/opt/TAMP
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
git clone https://github.com/marbl/canu.git
cd canu/src
make -j 16
```

Then expose the stable binary:

```bash
export PATH="$HOME/opt/canu/build/bin:$PATH"
```

<a id="usage"></a>
## Usage

<a id="basic-command"></a>
### Basic Command

Run TAMP in a dedicated working directory and provide the input FASTQ file by absolute path:

```bash
TAMP.sh -g 12m -t 16 --fastq /absolute/path/to/reads.fastq -m TGTG
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
mkdir -p ~/Storage/projects/TAMP/SRR13577847
cd ~/Storage/projects/TAMP/SRR13577847

TAMP.sh \
  -g 12m \
  -t 16 \
  --fastq ~/Storage/projects/TAMP/data/SRR13577847.fastq \
  -m TGTG
```

For fungi with vertebrate-like telomere repeats, such as *Neurospora crassa*, use a different motif as appropriate, for example `TTAGGG`.

<a id="pipeline-steps"></a>
## Pipeline Steps

TAMP currently implements the following major steps:

1. HiCanu assembly
2. NextDenovo assembly
3. Peregrine assembly
4. IPA assembly
5. Flye assembly
6. Hifiasm assembly
7. Copy and normalize all assemblies
8. BUSCO on all assemblies
9. Telomere contig detection and telomere metrics
10. Merge all assemblies
11. QUAST for assembler results
12. Final merge using the selected assembly
13. BUSCO analysis of final assembly
14. Telomere analysis of final assembly
15. QUAST analysis of final assembly
16. Final comparison report
17. Cleanup into structured output folders

<a id="input-data-recommendations"></a>
## Input Data Recommendations

TAMP is best suited to long-read data, especially PacBio HiFi reads, for small eukaryotic genomes. For publication-quality benchmarking:

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

Redundans can use several external tools. Not all of them are required for every TAMP run, but the upstream software lists the following resources:

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

In practice, the TAMP workflow mainly depends on the subset needed for the specific Redundans operations invoked during the final merge.

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

TAMP writes benchmark-oriented logs to support reproducible testing and publication reporting. These can include:

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

### `TAMP.sh: command not found`

Add the install directory or launcher directory to your `PATH`, or run TAMP using the full path.

<a id="citation-and-archiving"></a>
## Citation and Archiving

If you use TAMP in a publication, please cite the software and archive the exact release used for the analysis in a persistent public repository such as Zenodo. Citing the archived release helps ensure that the software version used in the study remains accessible and reproducible.

TAMP was developed at the Grainger Bioinformatics Center, Field Museum of Natural History, Chicago, Illinois, USA.

<a id="acknowledgements"></a>
## Acknowledgements

TAMP was developed in the context of genome assembly benchmarking and telomere-aware merging workflows at the Grainger Bioinformatics Center, Field Museum of Natural History.
