![Latest Version](https://img.shields.io/github/v/tag/yksun/TAMP?label=Latest%20Version)
[![BioConda Install](https://img.shields.io/conda/dn/bioconda/tamp.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/tamp)
[![GitHub Clones](https://img.shields.io/badge/dynamic/json?color=success&label=Clones&query=count&url=https://gist.githubusercontent.com/Dfupa/0fc9a42bb90e0b6c38767174bce725db/raw/clone.json&logo=github)](https://github.com/MShawon/github-clone-count-badge)
![Docker Pulls](https://img.shields.io/docker/pulls/yksun/tamp)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)]()
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)]()

<img align="center" src="/docs/tamp-icon.png">

## Table of Contents

- [TAMP](#tamp)
  - [Overview](#overview)
  - [Features](#features)
  - [Workflow](#workflow)
  - [Installation](#installation)
    - [Conda Environment](#conda-environment)
    - [Manual Installation (Recommended for Canu & Redundans)](#manual-installation-recommended-for-canu--redundans)
    - [Docker (coming soon)](#docker-coming-soon)
    - [Singularity (coming soon)](#singularity-coming-soon)
  - [Usage](#usage)
    - [Basic Command](#basic-command)
    - [Parameters](#parameters)
    - [Example Run](#example-run)
  - [Pipeline Steps](#pipeline-steps)
  - [Dependencies](#dependencies)
  - [Output Structure](#output-structure)
  - [Troubleshooting](#troubleshooting)
  - [Citation](#citation)
  - [Acknowledgements](#acknowledgements)


# TAMP

**Telomere Assembly Merge Pipeline**

TAMP is a modular workflow for **genome assembly benchmarking, telomere-aware evaluation, and final assembly merging**. It is designed for **PacBio HiFi reads** and is especially useful for **small eukaryotic genomes**, with a strong initial focus on **fungal genomes**. The pipeline runs multiple assemblers, standardizes their outputs, compares their quality, identifies telomeric contigs, and builds a final merged assembly while protecting contigs that appear to be telomere-to-telomere (T2T).

TAMP was built to support both day-to-day assembly work and formal benchmarking for publication. It writes step-level logs, software version records, and machine-readable benchmark summaries that can be used directly in supplementary materials, methods documentation, and reproducibility reports.

This repository uses generic file names for the main entry points:

- `TAMP.sh`
- `install_tamp.sh`
- `tamp-env.yml`

These generic names make the workflow easier to install, document, and cite across releases.

---

## Why use TAMP?

Genome assembly often requires testing more than one assembler because different tools can produce different results from the same read set. One assembler may produce fewer contigs, another may recover more complete genes, and another may better preserve chromosome ends. TAMP helps organize this process into one workflow.

TAMP is especially useful when you want to:

- compare several long-read assemblers on the same dataset
- standardize and rank assemblies before choosing a final result
- detect telomeric contigs and estimate chromosome completeness
- merge assemblies while protecting high-confidence T2T contigs
- generate benchmark logs for internal validation or publication

---

## Main features

TAMP can:

- run multiple assemblers: **HiCanu, NextDenovo, Peregrine, IPA, Flye, and hifiasm**
- copy and normalize assemblies into a common structure
- rename contigs consistently and sort them by length
- detect telomeric contigs and summarize telomere signals
- run **BUSCO** and **QUAST** for assembly evaluation
- build a unified comparison table across assemblers
- preserve T2T contigs during final merging
- produce machine-readable benchmark logs and version tracking files

---

## System requirements

TAMP is intended for **Linux**. A recent Linux distribution with standard shell utilities is recommended. The workflow expects the required software to be available on `PATH` through the installed environment and helper launchers.

Recommended resources depend on genome size, read depth, and which assemblers are enabled. For small fungal genomes, a multi-core workstation or server is usually sufficient. Larger genomes or very deep HiFi datasets will require more RAM, more CPU time, and more temporary storage.

---

## Software dependencies

TAMP depends on a conda environment plus one manually installed component.

### Core conda environment

The main environment is defined in `tamp-env.yml`. It should include the following software families.

#### Assemblers

- `canu`
- `nextdenovo`
- `peregrine-2021`
- `pbipa`
- `flye`
- `hifiasm`

#### Alignment and assembly utilities

- `minimap2`
- `bwa`
- `samtools`
- `seqtk`
- `miniasm`
- `quickmerge`

#### Quality assessment and annotation support

- `busco`
- `quast`
- `funannotate`

#### General tools and build support

- `python=3.10`
- `numpy`
- `pandas`
- `parallel`
- `pigz`
- `curl`
- `wget`
- `git`
- `make`
- `gxx_linux-64`

### Manual dependency: Redundans

`redundans` is installed **outside the conda environment** by the installer because current conda builds of `redundans` can conflict with modern builds of tools such as `minimap2`, `seqtk`, and `hifiasm`. TAMP therefore installs the `redundans` repository directly into the TAMP installation folder and creates a launcher for it.

This manual install approach preserves a clean modern conda environment while keeping `redundans` available for the final merge workflow.

### Redundans-related tools

The upstream `redundans` project documents a broad set of tools. Not all of them are required for every TAMP run. TAMP mainly relies on the core tools needed for redundancy reduction and related alignment tasks. The most relevant runtime requirements are:

- `python` (3.8 to <3.11)
- `minimap2`
- `bwa`
- `miniasm`
- standard Unix shell tools

The upstream `redundans` repository also documents additional tools that may be used in some configurations or advanced modes:

- Platanus
- LAST
- SNAP aligner
- SSPACE3
- GapCloser
- GFAstats
- Meryl
- Merqury
- k8
- R
- `r-ggplot2`
- `r-scales`
- `r-argparse`

These extra tools are **not required for the main TAMP workflow unless you extend the redundans usage beyond the default pipeline behavior**.

---

## Installation

### Recommended installation

Run the installer from the repository directory:

```bash
bash install_tamp.sh
```

The installer will:

1. detect the latest `TAMP-*.sh` script in the folder and install it as `TAMP.sh`
2. copy `TAMP.sh`, `install_tamp.sh`, and `tamp-env.yml` into the installation directory
3. create or update the `tamp` conda environment from `tamp-env.yml`
4. clone and install `redundans` manually into the TAMP installation directory
5. create command launchers so `TAMP.sh` can be run without manually activating the environment each time
6. add the launcher path to `~/.bashrc` if needed

By default, the installer uses `~/opt` as the installation base and installs TAMP under:

```text
~/opt/TAMP
```

### Installation layout

A typical installation looks like this:

```text
~/opt/
├── TAMP/
│   ├── TAMP.sh
│   ├── TAMP-0.30.sh
│   ├── install_tamp.sh
│   ├── tamp-env.yml
│   ├── redundans/
│   ├── logs/
│   └── benchmark_logs/
└── bin/
    ├── TAMP.sh
    └── redundans.py
```

### Environment notes

TAMP can be used with `conda`, `mamba`, or `micromamba`, but it is best to keep your shell initialization simple and consistent. If you already use conda or mamba through Miniforge, that is usually the easiest option.

---

## Input requirements

TAMP uses long-read sequence data as input.

### Required inputs

- PacBio HiFi FASTQ file
- genome size estimate (for example `12m`, `40m`, or `2g`)
- telomere motif or telomere seed motif appropriate for the target species

### Optional inputs

- external pre-assembled FASTA file
- BUSCO lineage dataset
- selected steps only
- explicit assembler choice for the final merge

---

## Choosing genome size and telomere motif

The genome size and telomere motif should match the organism you are assembling.

### Examples

- *Saccharomyces cerevisiae*: genome size about **12 Mb**
- *Neurospora crassa*: genome size about **40 Mb**

Telomere motifs also vary by species.

- `TTAGGG` is appropriate for organisms with vertebrate-like telomeres, including some fungi such as *Neurospora crassa*
- *Saccharomyces cerevisiae* does **not** use a fixed `TTAGGG` repeat. Its telomeres are heterogeneous **TG1-3** repeats, so a short TG-rich seed motif such as `TGTG` is a better practical choice than `TTAGGG`

Do **not** assume one motif works for all fungi.

---

## Basic usage

### Example 1: full run

```bash
TAMP.sh \
  -g 40m \
  -t 16 \
  --fastq /path/to/reads.fastq \
  -m TTAGGG
```

### Example 2: selected steps only

```bash
TAMP.sh \
  -g 40m \
  -t 16 \
  --fastq /path/to/reads.fastq \
  -m TTAGGG \
  -s 1,3-5
```

### Example 3: BUSCO on assemblies and interactive final assembler choice

```bash
TAMP.sh \
  -g 40m \
  -t 16 \
  --fastq /path/to/reads.fastq \
  -m TTAGGG \
  --busco fungi_odb10 \
  --choose
```

### Example 4: *Saccharomyces cerevisiae*

```bash
TAMP.sh \
  -g 12m \
  -t 16 \
  --fastq /path/to/SRR13577847.fastq \
  -m TGTG
```

---

## Command-line options

```text
-g, --genomesize   Genome size (required), for example 40m or 2g
-t, --threads      Number of CPU threads (required)
--fastq            Input FASTQ file (required)
-m, --motif        Telomere motif or seed motif (required)
-s, --steps        Steps to run, for example 1,3-5
--fasta            External pre-assembled FASTA file or URL
--busco            Run BUSCO on each assembly; optional lineage can be given
--choose           Choose the final assembler interactively
```

---

## Pipeline steps

1. **HiCanu assembly**  
2. **NextDenovo assembly**  
3. **Peregrine assembly**  
4. **IPA assembly**  
5. **Flye assembly**  
6. **Hifiasm assembly**  
7. **Copy and standardize assemblies**  
8. **BUSCO on assemblies**  
9. **Telomere contig detection and telomere metrics**  
10. **Merge all assemblies**  
11. **QUAST on assembler results**  
12. **Final merge using the selected assembler**  
13. **BUSCO on final assembly**  
14. **Telomere analysis on final assembly**  
15. **QUAST on final assembly**  
16. **Final comparison report**  
17. **Cleanup and file organization**

---

## Output files

TAMP writes its results into the current working directory.

Important output locations include:

- `assemblies/` — assembly FASTA files, telomere lists, and summary tables
- `logs/` — per-step logs
- `benchmark_logs/` — machine-readable benchmark and runtime summaries
- `version.txt` — software version log
- final merged assembly FASTA files and telomere summary files

### Benchmark logging outputs

TAMP writes benchmark-friendly output files such as:

- `benchmark_logs/run_metadata.tsv`
- `benchmark_logs/step_benchmark.tsv`
- `benchmark_logs/run_summary.txt`

These files are intended to support:

- internal pipeline testing
- runtime benchmarking
- manuscript preparation
- supplementary materials for reproducibility

---

## Recommended project layout

A clean working structure is:

```text
TAMP/
├── data/
│   ├── SRR13577847.fastq
│   └── SRR18210286.fastq
├── SRR13577847/
├── SRR18210286/
├── TAMP.sh
├── install_tamp.sh
├── tamp-env.yml
└── README.md
```

You can keep the raw data in `data/` and create one working directory per dataset. This makes reruns easier and keeps raw input separate from generated output.

---

## Notes on reproducibility

TAMP was designed with reproducibility in mind.

The pipeline records:

- software versions
- per-step logs
- benchmark-ready timing tables
- standardized output tables for assembly comparison

For publication, it is recommended to archive the exact code release and environment files used for the study.

---

## Citation and archiving

If you use TAMP in a publication, please cite the software and archive the exact release used for the analysis in a persistent repository such as Zenodo. Citing the archived release in the manuscript helps ensure transparency, reproducibility, and long-term access to the specific version used in the study.

TAMP was developed at the Grainger Bioinformatics Center, Field Museum of Natural History, Chicago, Illinois, USA.

---

## Troubleshooting

### `TAMP.sh: command not found`

Run the installed launcher directly:

```bash
~/opt/TAMP/TAMP.sh -h
```

If needed, add the install path to `PATH`.

### `redundans.py` not found

Re-run the installer and check that the `redundans` repository was cloned into the installation folder and that the launcher exists in `~/opt/bin`.

### Wrong telomere motif

If telomere detection looks wrong, check the species-specific telomere repeat. Do not reuse one motif for all fungi.

---

## License


