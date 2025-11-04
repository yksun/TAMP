# 🧬 TAMP: Telomere Assembly Merge Pipeline (v0.1)

**TAMP (Telomere Assembly Merge Pipeline)** is a modular, automated pipeline for benchmarking, selecting, and merging genome assemblies—optimized for **fungal haploid genomes** and **PacBio HiFi** reads, but broadly applicable to other small eukaryotic genomes.

---

## ✨ What’s new in v0.1
- Script renamed to **`TAMP-0.1.sh`**.
- **Per‑step logging**: each step writes timestamped stdout/stderr to `logs/step_<N>.log`. 
- **Version snapshot**: a consolidated `version.log` is generated at run start (records versions of key tools found on `$PATH`).

---

## 📁 Repository Layout

> All conda YAML environments live under `TAMP/dependency/`.  
> All helper scripts (e.g., `extract_contig_T_V3.sh`) live under `TAMP/dependency/scripts/`.

```
TAMP/
├── TAMP-0.1.sh                 # Main pipeline
├── dependency/
│   ├── envs and tool YAMLs (*.yml)
│   └── scripts/
│       └── extract_contig_T_V3.sh
├── README.md                   # This file
└── (generated at runtime)
    ├── logs/                   # Per-step logs: step_1.log, step_2.log, ...
    ├── version.log             # Software versions snapshot
    ├── quast_results/          # QUAST outputs (if enabled)
    ├── busco                   # BUSCO outputs (if enabled)
    ├── assembly_info.csv       # Aggregated metrics
    └── assembly_info.md        # Markdown summary
```

---

## ⚙️ Installation & Environments

All conda environment recipes are provided in **`TAMP/dependency/`** as YAML files. You can create one consolidated environment or separate ones for assembly/evaluation depending on your setup.


### Option: Multiple environments
```bash
conda env create -f dependency/busco.yml
conda env create -f dependency/funannotate.yml
conda env create -f dependency/pacbiohifi.yml
conda env create -f dependency/redundans.yml
# activate the one relevant to your step(s)
conda activate pacbiohifi # pipeline will automatically activate the env
```

> **Tip:** List available YAMLs with `ls dependency/*.yml` and choose the one(s) that match your site/tooling.

### External data requirements
- **BUSCO** requires a lineage dataset (e.g., `ascomycota_odb10`). See BUSCO docs on downloading lineages.
- **QUAST** may need Python and R dependencies depending on modules used.

---

## 🚀 Quick Start

```bash
# Make the pipeline executable
chmod +x TAMP-0.1.sh

# Run on HiFi reads (with BUSCO):
./TAMP-0.1.sh   -g 90m -t 20   --fastq Pste_1.0.dup.fastq.gz   -m AACCCT   --busco ascomycota_odb10

# Evaluate an existing FASTA (steps 7–13 only):
./TAMP-0.1.sh -g 90m -t 20 --fasta Physcia_stellaris.scaffolds.fa -m AACCCT -s 7-13 --busco ascomycota_odb10
```

### Common flags
| Flag | Description |
|------|-------------|
| `-g` | Genome size, e.g., `2g`, `500m`, `90m` (required) |
| `-t` | Threads, e.g., `20` (required) |
| `--fastq` | PacBio HiFi reads (optional if `--fasta` provided) |
| `--fasta` | Pre-assembled genome (optional; may be combined with `--fastq`) |
| `-m` | Telomere motif, e.g., `AACCCT` (required) |
| `-s` | Steps to run (comma/range), e.g., `1,3-5` (default: all) |
| `--busco` | BUSCO lineage for completeness (e.g., `ascomycota_odb10`) |
| `--choose` | Interactively choose the best assembly before merge |

---

## 🔢 Pipeline Overview (Typical Steps)

| Step | Description |
|------|-------------|
| **1–6** | Run supported assemblers (HiCanu, NextDenovo, Peregrine, IPA, Flye, RAFT‑hifiasm) |
| **7** | Aggregate assemblies and baseline stats; run BUSCO if `--busco` provided |
| **8** | Detect telomere‑containing contigs; count single/double‑end telomeres |
| **9** | Merge assemblies into consensus (`allmerged_telo.fasta`) |
| **10** | QUAST: contiguity/accuracy/GC metrics |
| **11** | Pre‑merge preview & Markdown summary (`assembly_info.csv`) |
| **12–13** | Post‑merge BUSCO & telomere cleanup |

> The helper script **`dependency/scripts/extract_contig_T_V3.sh`** is invoked by the pipeline for telomere‑related extraction tasks where required.

---

## 📝 Logging & Reproducibility

- **Per‑step logs:** created under `logs/` as `step_<N>.log`, each line timestamped.
- **Versions snapshot:** `version.log` is written at the start of the run, capturing versions of key tools in your `$PATH` (e.g., canu, nextDenovo, pg_asm, ipa, flye, hifiasm, seqtk, BUSCO, QUAST, minimap2, bwa, samtools, python3, etc.).
- **Run summaries:** consolidated in `assembly_info.csv` and a human‑readable `assembly_info.md` for GitHub.

---

## 📤 Outputs

| File/Folder | Description |
|-------------|-------------|
| `assembly_info.csv` | Incremental statistics log updated after each step |
| `assembly_info.md` | Markdown summary that’s easy to review on GitHub |
| `quast_results/` | QUAST reports and tables |
| `busco_*` | BUSCO result directories |
| `t2t_clean.fasta` | Cleaned telomere‑to‑telomere assembly |
| `merged_<assembler>_sort.fa` | Final merged assembly output |

---

## 🧪 Examples

**Full benchmark & choose best assembly interactively**
```bash
./TAMP-0.1.sh -g 2g -t 32 --fastq reads.fastq -m AACCCT --busco fungi_odb10 --choose
```

**Resume evaluation from step 7**
```bash
./TAMP-0.1.sh -g 2g -t 32 -m AACCCT --fasta existing.fa -s 7-10 --busco ascomycota_odb10
```

---

## 🧯 Troubleshooting

- If timestamping fails, ensure you have `gawk` available. A POSIX‑`awk` fallback can be enabled if needed.
- For BUSCO, verify lineage data are downloaded and accessible.
- QUAST output paths vary slightly with versions; the pipeline handles `report.tsv`/`report.txt` gracefully.

---

## 🧾 Citation

If you use **TAMP** in your research, please cite this repository and the underlying tools:
HiCanu, NextDenovo, Peregrine, IPA, Flye, Hifiasm, RAFT, BUSCO, QUAST, Funannotate, Seqtk.
