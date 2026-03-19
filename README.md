# TAMP

**Telomere Assembly Merge Pipeline**

TAMP is a modular workflow for benchmarking, selecting, and merging genome assemblies. It is built for fungal haploid genomes and PacBio HiFi reads, and it can also be used for other small eukaryotic genomes.

This README uses **generic file names** for the main pipeline and install files:

- `TAMP.sh`
- `install_tamp.sh`
- `tamp-env.yml`

It does **not** tie the commands in this README to a specific versioned script name.

## Main features

TAMP can:

- run multiple assemblers: HiCanu, NextDenovo, Peregrine, IPA, Flye, and RAFT-hifiasm
- standardize assembly FASTA files
- detect telomeric contigs
- evaluate assemblies with QUAST and BUSCO
- build a combined metrics table for assembly comparison
- protect telomere-to-telomere contigs during final merging
- write benchmark logs for reproducibility and paper-ready reporting

## Requirements

TAMP is intended for Linux systems. The workflow expects the tools in the conda environment to be available on `PATH`.

## Installation

Create the environment and install TAMP with the generic installer name:

```bash
bash install_tamp.sh tamp tamp-env.yml
conda activate tamp
```

If you use `micromamba` or `mamba`, keep the same file names and workflow, but use the solver available on your system.

## Input data

TAMP uses long-read sequencing data as input.

Required input:

- PacBio HiFi FASTQ file
- genome size estimate
- telomere motif for the target species

Optional input:

- pre-assembled FASTA file
- BUSCO lineage
- selected pipeline steps
- explicit assembler choice for final merge

## Basic usage

Run the pipeline with:

```bash
bash TAMP.sh -g 40m -t 16 --fastq reads.fastq.gz -m TTAGGG
```

Run only selected steps:

```bash
bash TAMP.sh -g 40m -t 16 --fastq reads.fastq.gz -m TTAGGG -s 1,3-5
```

Run with BUSCO and interactive final assembler choice:

```bash
bash TAMP.sh -g 40m -t 16 --fastq reads.fastq.gz -m TTAGGG --busco fungi_odb10 --choose
```

## Command-line options

```text
-g, --genomesize   Genome size, for example 40m or 2g
-t, --threads      Number of CPU threads
--fastq            Input FASTQ file
-m, --motif        Telomere motif for the target species
-s, --steps        Steps to run, for example 1,3-5
--fasta            External pre-assembled FASTA file or URL
--busco            Run BUSCO; optional lineage can be given
--choose           Choose the final assembler interactively
```

## Pipeline steps

1. HiCanu assembly  
2. NextDenovo assembly  
3. Peregrine assembly  
4. IPA assembly  
5. Flye assembly  
6. Hifiasm assembly  
7. Copy and standardize assemblies  
8. BUSCO on assemblies  
9. Telomere contig detection and telomere metrics  
10. Merge all assemblies  
11. QUAST on assembler results  
12. Final merge using the selected assembler  
13. BUSCO on final assembly  
14. Telomere analysis on final assembly  
15. QUAST on final assembly  
16. Final comparison report  
17. Cleanup and file organization  

## Output files

Main outputs include:

- `assemblies/` for assembly FASTA files and summary tables
- `logs/` for per-step logs
- `benchmark_logs/` for benchmarking and reproducibility logs
- `version.txt` for software versions
- final assembly FASTA files and telomere summary files

## Benchmark logging

TAMP writes machine-readable benchmark files that can be used for internal testing and manuscript preparation.

Expected benchmark outputs include:

- `benchmark_logs/run_metadata.tsv`
- `benchmark_logs/step_benchmark.tsv`
- `benchmark_logs/run_summary.txt`

These files record runtime information, step status, and run metadata.

## Telomere motif note

The telomere motif must match the species being assembled. Do not assume one motif works for all fungi.

Example:

- `TTAGGG` is suitable for *Neurospora crassa*
- other fungi may require a different repeat

## Reproducibility and manuscript use

TAMP is structured to support reproducible benchmarking.

For a PeerJ paper, it is helpful to provide:

- the exact `TAMP.sh` script used for the analysis
- the exact `install_tamp.sh` and `tamp-env.yml` files used for installation
- the benchmark output files from `benchmark_logs/`
- the raw or public input data accession numbers
- an archived repository release with a DOI

PeerJ states that code and raw data are almost always required, that files should be machine-readable, and that repository snapshots with a DOI are preferred for reproducibility. citeturn746666search0turn746666search3

For a methods or software paper, you should also report the source of third-party datasets, the preprocessing steps, the computing environment, and the evaluation metrics used in benchmarking. citeturn746666search4turn746666search0

## Recommended repository layout

```text
TAMP/
├── TAMP.sh
├── install_tamp.sh
├── tamp-env.yml
├── README.md
├── logs/
├── benchmark_logs/
└── assemblies/
```

## Citation and archiving

If you publish results generated with TAMP, archive the exact release used for the analysis in Zenodo or another DOI-minting repository. Then cite that archived release in the manuscript.

## Notes

- Keep helper scripts integrated into `TAMP.sh` when possible.
- Keep tool names on `PATH` through one activated conda environment.
- Keep benchmark files with the final results so the run can be checked later.
