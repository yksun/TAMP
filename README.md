<img align="center" src="/docs/taco-icon.png">

# TACO

**Telomere-Aware Contig Optimization**

TACO is a telomere-aware **all-in-one multi-assembler comparison and refinement pipeline** for genome assembly benchmarking, decision-making, and chromosome-end improvement. The pipeline was developed for small eukaryotic genomes, with a focus on fungal genomes and PacBio HiFi reads. It runs multiple assemblers, standardizes their outputs, evaluates assembly quality, detects telomere-supported contigs, and can either: **(1)** stop after generating a unified assembler comparison table for convenient benchmarking, or **(2)** continue into telomere-aware backbone refinement for an improved chromosome-scale candidate assembly.

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
  - [Assembly-only Mode](#assembly-only-mode)
  - [Parameters](#parameters)
  - [Parameter Details](#parameter-details)
  - [Example Run](#example-run)
- [Assembly Selection Strategy](#assembly-selection-strategy)
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

Rather than focusing only on recovering strict telomere-to-telomere (T2T) contigs, TACO is designed as a multi-assembler decision and refinement workflow. It runs several assemblers, normalizes their outputs, summarizes assembly quality metrics, detects telomere-supported contigs, and builds a unified comparison table so users can decide which assembler performs best for their dataset. This makes TACO useful not only as a merge/refinement workflow, but also as a convenient **all-in-one assembler comparison pipeline**.

TACO can then refine the selected backbone assembly by prioritizing protected telomere-supported contigs, including strict T2T contigs when present, while reducing redundant non-telomeric sequence. In other words, users can treat TACO as either a comparison-first workflow or a comparison-plus-refinement workflow, depending on whether `--assembly-only` is used.

<a id="features"></a>
## Features

- Run multiple long-read assemblers from one workflow.
- Standardize assembly outputs for direct comparison across assemblers.
- Benchmark assemblies with QUAST, BUSCO, telomere-support summaries, and optional Merqury.
- Build a unified `assemblies/assembly_info.csv` table for convenient assembler comparison.
- Support an **assembly-only mode** for users who want benchmarking and comparison without refinement.
- Separate strict T2T contigs from single-end telomeric contigs.
- Build an optimized telomere pool before final refinement.
- Preserve protected telomere-supported contigs during final backbone refinement.
- Reduce redundant non-telomeric contigs while retaining chromosome-end support.
- Generate machine-readable benchmark logs and software version records for reproducibility.

<a id="workflow"></a>
## Workflow

TACO follows this high-level order:

1. Run one or more assemblers on the same long-read dataset.
2. Normalize assembly outputs and contig names for comparison.
3. Detect telomere-supported contigs for each assembly.
4. Evaluate assemblies with QUAST, BUSCO, and optional Merqury.
5. Build a combined summary table across assemblies to support assembler choice.

At this point, the workflow can branch:

- **Assembly-only mode (`--assembly-only`)**: stop after benchmarking and write `assemblies/assembly_info.csv` as the main comparison table.
- **Full refinement mode**: continue into telomere-pool optimization and final backbone refinement.

Full refinement mode then continues with:

6. Select a preferred backbone assembly based on the comparison results.
7. Build an optimized telomere pool from all assemblies.
8. Refine the selected backbone by replacing redundant backbone contigs with protected telomere-supported contigs and rescuing useful telomeric ends.
9. Run final telomere checks, final BUSCO/QUAST, and summary reporting.

A schematic workflow figure should emphasize four concepts: multi-assembler benchmarking, assembly-only comparison mode, telomere-supported contig classification, and backbone refinement by protected contig replacement. Place the updated TACO figure in `docs/taco-workflow.png` or `docs/taco-icon.png` and link it here when ready.

<a id="installation"></a>
## Installation

TACO uses a primary Conda or Micromamba environment for most tools. At present, **Redundans is installed manually** because of dependency conflicts with the modern Conda stack used by other tools. **Canu may also be installed manually** if your local package resolution selects an unstable development build or if your environment has persistent solver conflicts.

<a id="conda-environment"></a>
### Conda Environment

Create the main environment from `taco-env.yml`:

```bash
micromamba create -n taco -f taco-env.yml
micromamba activate taco
```

Or with Conda:

```bash
conda env create -n taco -f taco-env.yml
conda activate taco
```

If you use the provided installer:

```bash
bash install_taco.sh
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

Run TACO in a dedicated working directory and provide the input FASTQ file by absolute path.

**Most users will use TACO in one of two ways:**

1. **Assembler comparison only** with `--assembly-only`, when the main goal is to compare assemblers and produce `assemblies/assembly_info.csv`.
2. **Full telomere-aware refinement**, when the goal is to compare assemblers and then build a final refined assembly.

#### Full refinement run

```bash
TACO.sh -g 12m -t 16 --fastq /absolute/path/to/reads.fastq -m TGTG
```

#### Assembly-only comparison run

```bash
TACO.sh -g 12m -t 16 --fastq /absolute/path/to/reads.fastq -m TGTG --assembly-only
```

#### Assembly-only comparison with Merqury

```bash
TACO.sh -g 12m -t 16 --fastq /absolute/path/to/reads.fastq -m TGTG --assembly-only --merqury-db reads.meryl
```

<a id="assembly-only-mode"></a>
### Assembly-only Mode

Use `--assembly-only` when you want TACO primarily as an **all-in-one assembler comparison pipeline**.

In this mode, TACO runs the assembler generation and comparison steps, including:
- assembler execution
- assembly standardization
- BUSCO
- telomere summaries
- QUAST
- optional Merqury for all assemblies

It then writes the combined comparison table to:

- `assemblies/assembly_info.csv`

and a copied summary table to:

- `final_results/assembly_only_result.csv`

This mode is especially useful when you want to:
- benchmark multiple assemblers on the same dataset,
- decide which assembler is best before doing any refinement,
- or use TACO mainly as a convenient one-command comparison workflow.


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
| `--choose` | Manually choose the assembler to use as the final refinement backbone. |
| `--assembly-only` | Run assembler comparison only and stop after generating `assemblies/assembly_info.csv` and the assembly-only summary outputs. |
| `--auto-mode smart` | Use the biologically informed automatic backbone scoring system that combines BUSCO, telomere support, optional Merqury metrics, contig count, and N50. |
| `--auto-mode n50` | Use legacy automatic selection based only on the highest N50. |
| `--merqury` | Enable optional Merqury scoring and reporting using an auto-detected `.meryl` database if available. |
| `--merqury-db <reads.meryl>` | Enable Merqury and explicitly provide the path to a read-derived `.meryl` database. |
| `--no-merqury` | Explicitly disable Merqury, even if `merqury.sh` and a `.meryl` database are present. |

<a id="parameter-details"></a>
### Parameter Details

#### `--auto-mode smart`
This is the default automatic selection mode. TACO scores candidate backbone assemblies using a composite metric that prioritizes biological completeness and chromosome-end support over contiguity alone:

```text
score =
  BUSCO_C (%) x 1000
+ Telomere_double_end_contigs x 500
+ Telomere_single_end_contigs x 100
+ Merqury_completeness (%) x 200
+ Merqury_QV x 20
- Number_of_contigs x 10
+ log10(N50) x 100
```

Use this mode when you want TACO to choose the backbone assembly automatically in a biologically informed way.

#### `--auto-mode n50`
This mode uses the highest N50 only. It is useful for reproducing older behavior or for quick continuity-only ranking, but it can favor assemblies with long contigs despite lower completeness.

#### `--assembly-only`
Use this mode when your main goal is assembler benchmarking and comparison rather than telomere-aware refinement. TACO will run the comparison-oriented steps, build `assemblies/assembly_info.csv`, optionally add Merqury metrics if enabled, and stop before final backbone refinement.

#### `--merqury`
This enables optional Merqury integration. TACO will try to detect a usable read-derived `.meryl` database automatically from common locations such as `reads.meryl`, `meryl/reads.meryl`, or `merqury/reads.meryl`.

When enabled, TACO:
- runs Merqury on available assembler outputs before backbone selection,
- adds Merqury QV and completeness to `assemblies/assembly_info.csv`,
- uses those metrics in the smart backbone scoring system,
- runs Merqury again on the final refined assembly,
- and reports final Merqury metrics in `final_results/final_result.csv`.

#### `--merqury-db <reads.meryl>`
Use this when you want Merqury enabled and already know the exact path to the `.meryl` database that should be used. This is the most reproducible way to run Merqury within TACO.

#### `--no-merqury`
Use this to force TACO to skip Merqury, even if `merqury.sh` is installed and a `.meryl` database is present. This is useful for lightweight runs or when no trustworthy read-derived k-mer database is available.


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
10. Build an optimized telomere pool across assemblies
11. QUAST for assembler results
12. Final assembly refinement with optimized telomere-end replacement
13. BUSCO analysis of final assembly
14. Telomere analysis of final assembly
15. QUAST analysis of final assembly
16. Final comparison report
17. Cleanup into structured output folders
18. Assembly-only comparison summary

If `--assembly-only` is used, TACO follows the comparison path and stops after generating the assembly comparison outputs instead of continuing into final refinement.

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

- `assemblies/` — normalized assembly FASTA files and summary tables, including `assembly_info.csv`
- `final_results/` — copied comparison tables and final refined outputs
- `logs/` — per-step log files
- `benchmark_logs/` — machine-readable benchmark and timing outputs
- `version.txt` — software versions recorded for the run
- final merged assembly files
- telomere summary files
- QUAST, BUSCO, and optional Merqury result directories

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


## Assembly Selection Strategy

When `--choose` is not provided, TACO can automatically select the backbone assembly for final refinement.

### Auto-selection modes

- `smart` (default): biologically informed multi-factor scoring
- `n50`: legacy behavior using the highest N50 only

### Smart scoring system

In full refinement mode, TACO ranks assemblies using a composite score that prioritizes biological completeness and chromosome-end support over contiguity alone:

```text
score =
  BUSCO_C (%) x 1000
+ Telomere_double_end_contigs x 600
+ Telomere_single_end_contigs x 250
+ Merqury_completeness (%) x 200
+ Merqury_QV x 20
- Number_of_contigs x 10
+ log10(N50) x 100
```

### Rationale

- **BUSCO completeness** is the primary driver and strongly favors biologically complete assemblies.
- **Telomere double-end contigs** provide the strongest evidence for chromosome-level continuity.
- **Telomere single-end contigs** contribute partial chromosome-end support.
- **Contig count** penalizes fragmented assemblies.
- **N50** contributes as a secondary contiguity term rather than dominating selection by itself.

### Example behavior

This scoring strategy prevents assemblies with high N50 but poor completeness from being selected over more biologically complete assemblies with fewer contigs and stronger telomere support.

### Command-line control

```bash
# Default biologically informed mode
./TACO.sh --auto-mode smart ...

# Legacy N50-only behavior
./TACO.sh --auto-mode n50 ...
```

### Debug output

During automatic selection, TACO reports the values used for scoring for each assembler, for example:

```text
[DEBUG] canu: BUSCO=97.5 T2T=2 single=1 MerquryQV=41.2 MerquryComp=99.1 contigs=56 N50=785256 score=...
```
