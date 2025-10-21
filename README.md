# 🧬 TAMP: Telomere Assembly Merge Pipeline (`TAMP_pre-0.1.sh`)

**TAMP (Telomere Assembly Merge Pipeline)** is a modular, automated pipeline for benchmarking and merging genome assemblies.  
It was tested and optimized for **fungal haploid genomes**, but can be easily adapted for any **PacBio HiFi dataset**.  
TAMP serves as both a **benchmarking framework**—automatically evaluating multiple assemblers—and a **decision-support tool** that simplifies the process of identifying the best-performing assembly for your genome.

---

## ✨ Key Features

- **Purpose-built for fungal genomes** (haploid assemblies).  
  - Tuned parameters for high accuracy and telomere-complete assemblies.
- **Flexible input:**  
  - Accepts either PacBio HiFi reads (`--fastq`) or pre-assembled genomes (`--fasta`).
- **Multi-assembler integration:**  
  - Runs or imports assemblies from **HiCanu**, **NextDenovo**, **Peregrine**, **IPA**, **Flye**, **RAFT-hifiasm**, and user-provided assemblies.
- **Automated evaluation & ranking:**
  - Runs **BUSCO** for completeness (`--busco <lineage>`).  
  - Runs **QUAST** for contiguity, accuracy, and GC metrics.  
  - Detects **telomere-containing contigs** and counts single/double-end contigs.
- **Unified reporting:**
  - Aggregates all statistics in `assembly_info.csv` and a Markdown summary `assembly_info.md`.
  - Generates human-readable pre-merge summaries (pretty-printed tables + markdown).
- **Version tracking:**  
  - Captures assembler versions and records them automatically.

---

## ⚙️ Installation & Requirements

### Dependencies

You’ll need the following tools available in `$PATH` or via conda:

| Tool | Purpose |
|------|----------|
| `bash >= 4.2` | Core shell interpreter |
| `conda` | Environment management |
| `seqtk` | Telomere motif scanning |
| `busco >= 5.x` | Completeness evaluation |
| `quast >= 5.x` | Assembly statistics |
| `funannotate` | Sorting & cleaning assemblies |
| `awk`, `sed`, `grep`, `column` | Core utilities |

Each assembler should also be installed and callable:

- `canu`, `nextDenovo`, `pg_asm` (Peregrine), `ipa`, `flye`, `hifiasm`, and optionally `RAFT`.

---

## 🚀 Usage

```bash
bash TAMP.sh   -g 2g   -t 32   --fastq reads.fastq   --fasta preassembly.fa   -m AACCCT   --busco ascomycota_odb10   --choose
```

### Arguments

| Flag | Description |
|------|--------------|
| `-g` | Genome size (required). e.g. `2g`, `500m` |
| `-t` | Threads (required). e.g. `32` |
| `--fastq` | Path to PacBio HiFi reads (optional if `--fasta` provided) |
| `--fasta` | Path to pre-assembled genome FASTA (optional or combined) |
| `-m` | Telomere motif (required). e.g. `AACCCT` |
| `-s` | Steps to run (comma/range). e.g. `1,3-5` (default: all) |
| `--busco` | Lineage for BUSCO evaluation (default: `ascomycota_odb10`) |
| `--choose` | Prompt to choose assembler for final merge interactively |

---

## 🔢 Pipeline Steps Overview

| Step | Description |
|------|--------------|
| **1–6** | Run supported assemblers (HiCanu, NextDenovo, Peregrine, IPA, Flye, RAFT-hifiasm) |
| **7** | Copy all assemblies → logs to `assembly_info.csv`, runs BUSCO if `--busco` |
| **8** | Detect telomere-containing contigs → counts single/double-end telomeres |
| **9** | Merge assemblies → generates `allmerged_telo.fasta` |
| **10** | Run QUAST → extracts all available metrics (`report.tsv` or fallback `report.txt`) |
| **11** | Previews QUAST report + full Markdown summary before final merge |
| **12–13** | Post-merge BUSCO & telomere cleanup |

---

## 📊 Unified Reporting System

After execution, TAMP outputs the following summaries:

### **assembly_info.csv**
Comprehensive summary of all assemblies, including:

```
Name, Version, Lineage, Timestamp,
BUSCO metrics (C, S, D, F, M),
Telomere counts (DoubleEndTelo, SingleEndTelo),
QUAST metrics (Total length, #contigs, N50, L50, GC%, etc.)
```

### **assembly_info.md**
Markdown-formatted summary generated automatically — easy to view in GitHub or share with collaborators.

**Pre-merge preview:**  
Before final merging, the script prints a formatted table and exports `assembly_info.md`, allowing you to choose the best assembly visually.

---

## 🧩 Outputs

| File | Description |
|------|--------------|
| `assembly_info.csv` | Incremental statistics log (updated after each step) |
| `assembly_info.md` | Markdown version for GitHub and sharing |
| `quast_results/` | QUAST output files |
| `busco_*` | BUSCO result directories |
| `t2t_clean.fasta` | Cleaned telomere-to-telomere assembly |
| `merged_<assembler>_sort.fa` | Final merged assembly output |

---

## 🧠 Example Workflows

### Evaluate pre-assembled FASTA only
```bash
./TAMP.sh -g 2g -t 32 --fasta genome.fa -m AACCCT --busco ascomycota_odb10 -s 7-13
```

### Run full benchmarking and select best assembly
```bash
./TAMP.sh -g 2g -t 32 --fastq reads.fastq -m AACCCT --busco fungi_odb10 --choose
```

### Resume evaluation from step 7
```bash
./TAMP.sh -g 2g -t 32 -m AACCCT --fasta existing.fa -s 7-10 --busco ascomycota_odb10
```

---

## 🧾 Citation

If you use **TAMP** in your research, please cite this repository and the underlying tools:
- *HiCanu, NextDenovo, Peregrine, IPA, Flye, Hifiasm, RAFT, BUSCO, QUAST, Funannotate, Seqtk.*

---

## 📁 Repository Layout

```
├── TAMP.sh                       # Main pipeline script
├── assembly_info.csv             # Auto-generated run summary
├── assembly_info.md              # Markdown summary for GitHub
├── quast_results/                # QUAST outputs
└── README.md                     # This documentation
```

---

## 🧬 About

TAMP was developed for **fungal genome assembly benchmarking and telomere-to-telomere merging** workflows.  
It is ideal for researchers seeking a reproducible and automated evaluation framework for **PacBio HiFi** datasets.

