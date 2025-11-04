# 🧬 TAMP: Telomere Assembly Merge Pipeline (v0.1)

**TAMP** is a modular pipeline for benchmarking multiple assemblies, selecting a best candidate, and producing a telomere-aware final merge with unified reports.

> This release renames the main script to **`TAMP-0.1.sh`**, adds per-step logs, tool/version capture, and a canonical final assembly path: `assemblies/final.merged.fasta`.

---

## ✨ Highlights

- Evaluate assemblies from multiple tools (canu, nextDenovo, peregrine, ipa, flye, RAFT-hifiasm, + optional external).
- Compute **BUSCO**, **QUAST**, and **Telomere** metrics.
- Final merge + redundancy reduction with **canonical** output: `assemblies/final.merged.fasta`.
- Reproducible logging/provenance: per-step logs, tool versions, and conda env tracking.

---

## 📦 Repository layout (key)

```
.
├── TAMP-0.1.sh                 # main pipeline (renamed & instrumented)
├── README.md                   # this file
├── tamp-env.yml                # monolithic environment (optional)
├── assemblies/                 # per-assembler inputs & outputs
├── dependency/                 # <<< ENV YMLs and helper scripts live here
│   ├── busco.yml
│   ├── funannotate.yml
│   ├── pacbiohifi.yml
│   ├── quast.yml
│   ├── redundans.yml
│   ├── seqtk.yml
│   ├── scripts/
│   │   └── extract_contig_T_V3.sh
│   └── ... (any other .yml / .sh you added)
└── log/                        # logs written at runtime
```

> Put any additional environments (`*.yml`) and helper scripts (`*.sh`) under **`dependency/`** and they will be documented by this README’s guidance below.

---

## ⚙️ Installation

### Option A — single environment
```bash
mamba env create -f tamp-env.yml
mamba activate tamp
```

### Option B — per-tool environments (recommended for long-term reproducibility)
Create only what you need; activate before the corresponding steps.

| Env file | Purpose | Main tools | Create | Activate | Quick check |
|---|---|---|---|---|---|
| `dependency/busco.yml` | BUSCO completeness | `busco` | `mamba env create -f dependency/busco.yml -n busco` | `mamba activate busco` | `busco --version` |
| `dependency/quast.yml` | QUAST assembly stats | `quast.py`/`quast` | `mamba env create -f dependency/quast.yml -n quast` | `mamba activate quast` | `quast.py --version || quast -v` |
| `dependency/funannotate.yml` | Sorting final contigs; `seqtk` lives here in our setup | `funannotate`, `seqtk` | `mamba env create -f dependency/funannotate.yml -n funannotate` | `mamba activate funannotate` | `funannotate --version && seqtk 2>&1 | head -n1` |
| `dependency/pacbiohifi.yml` | Runtime for QUAST and other HiFi tasks (if you prefer this env) | `quast`, `samtools` | `mamba env create -f dependency/pacbiohifi.yml -n pacbiohifi` | `mamba activate pacbiohifi` | `quast.py --version || quast -v` |
| `dependency/redundans.yml` | Redundancy reduction | `redundans.py` | `mamba env create -f dependency/redundans.yml -n redundans` | `mamba activate redundans` | `redundans.py --help | head -n1` |
| `dependency/seqtk.yml` | Minimal env for seqtk-only use | `seqtk` | `mamba env create -f dependency/seqtk.yml -n seqtk` | `mamba activate seqtk` | `seqtk 2>&1 | head -n1` |

> If your `.yml` files already contain a `name:` field, the `-n <name>` flag is optional.

---

## 🚀 Quick start

```bash
# Example: run only Step 12 (final merge & publish)
bash TAMP-0.1.sh -g 50m -t 20 --fasta assemblies/canu.result.fasta -m AACCCT --busco ascomycota_odb10 -s 12

# Full compare + choose final (interactive prompt unless you pass --choose)
bash TAMP-0.1.sh -g 50m -t 20 --fastq reads.fastq.gz -m AACCCT --busco ascomycota_odb10
```

### Key arguments
- `-g <SIZE>`: genome size, e.g. `50m`, `2g`
- `-t <INT>`: threads
- `--fastq <reads.fastq.gz>`: PacBio HiFi reads (optional)
- `--fasta <assembly.fa>`: external/prebuilt assembly to include
- `-m <MOTIF>`: telomere motif (default `TTAGGG`; reverse is `AACCCT`)
- `--busco <LINEAGE>`: BUSCO lineage (e.g. `ascomycota_odb10`)
- `--choose [NAME]`: pick final (`canu|external|flye|ipa|nextDenovo|peregrine|RAFT-hifiasm`)
- `-s <N>`: run **only** step N (e.g., `12..17`)

---

## 🔢 Steps (final stages)

- **Step 12** — Final merge + redundans; writes `assemblies/final.merged.fasta`. If `--choose` is not provided, prints `assemblies/assembly_info.csv` (if present) and prompts for selection.
- **Step 13** — BUSCO on final; writes `assemblies/final.busco.csv` (metrics-as-rows; single column `final`).
- **Step 14** — Telomere analysis on final; writes `assemblies/final.telo.csv` plus `assemblies/final.telo.fasta`.
- **Step 15** — QUAST on final; writes `assemblies/final.quast.csv` and `assemblies/final-quast.tsv` (from `quast_final/`).
- **Step 16** — Merge final metrics into **`final_result.csv`** (in CWD), optionally appending to `assemblies/assembly_info.csv` columns.
- **Step 17** — Cleanup; organizes temp artifacts into `temp/*` and logs into `log/`.

Per-step logs: `log/step_<N>_<runid>.log`. Main aggregated log: `log/TAMP_<runid>.log`.  
Provenance: `version.log` (tool versions), `version_info.txt` (conda envs used).

---

## 🧰 Helper scripts (in `dependency/scripts/`)

### `extract_contig_T_V3.sh`
Extracts contigs from a FASTA by ID list (e.g., telomere-positive IDs derived from `seqtk telo`).

**Usage (examples):**
```bash
# Using telomere ID list generated in pipeline
bash dependency/scripts/extract_contig_T_V3.sh   -i assemblies/canu.result.fasta   -l assemblies/canu.result.telo.list.ids   -o canu.telo.fasta

# Generic extraction by IDs
bash dependency/scripts/extract_contig_T_V3.sh   -i assemblies/some_assembly.fasta   -l ids.txt   -o subset.fasta
```

**Common flags (may vary; run with -h to confirm):**
- `-i <FASTA>`: input assembly fasta
- `-l <LIST>`: file with one contig ID per line
- `-o <FASTA>`: output fasta

> Add additional helper scripts here (one subsection per script). Provide 1–2 sentence intros + a minimal example invocation.

---

## 🧪 Reproducibility & Logging

- Each step writes a dedicated log with timestamps and `set -x` tracing.
- The pipeline writes:
  - `version.log` — bash/python/conda + tool versions (BUSCO, QUAST, seqtk, funannotate, redundans, merge_wrapper where available)
  - `version_info.txt` — unique conda environments **actually activated** during the run + `conda info --envs`

---

## 🧯 Troubleshooting

- **QUAST not found:** activate the `quast` or `pacbiohifi` environment.
- **seqtk not found:** activate the `funannotate` (or `seqtk`) environment.
- **BUSCO rerun refused:** remove stale `busco/final` or run with force inside the pipeline; the pipeline already tries to reuse existing metrics when present.
- **No `assemblies/` folder:** create it and place your `*.result.fasta` files.

---

## 📜 License & Citation

Please cite this repository and the tools utilized (BUSCO, QUAST, Seqtk, Funannotate, MUMmer/Redundans, and the assemblers).

---

## 🙌 Acknowledgements

Developed to streamline fungal genome assembly evaluation and telomere-to-telomere merging. Contributions & PRs welcome!
