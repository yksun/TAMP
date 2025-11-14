# 🧬 TAMP: Telomere Assembly Merge Pipeline (v0.2.7)

**TAMP** is a modular, automated pipeline for benchmarking, selecting, and merging genome assemblies. It is optimized for **fungal haploid genomes** and **PacBio HiFi** reads, but is generally applicable to small eukaryotic genomes.

---

## What TAMP Does (at a glance)

- **Runs multiple assemblers** (1–6): HiCanu, NextDenovo, Peregrine, IPA, Flye, RAFT-hifiasm.
- **Normalizes and copies assemblies** (7): consistent contig headers and length-sorted FASTA per assembler.
- **Telomere discovery & per-assembler telomeric subsets** (8–9): finds telomeric contigs (via `seqtk telo`) and prepares per-assembler telomeric FASTAs.
- **Quality assessment** (10, 13): QUAST for structure metrics; BUSCO for completeness.
- **Unified assembly metrics** (12): *before merging*, TAMP **rebuilds** `assemblies/assembly_info.csv` by combining BUSCO/QUAST/TELO CSVs so you can compare candidates in one table.
- **Final merge with T2T protection** (12):
  - You choose an assembly interactively or with `--choose`.
  - T2T contigs from `t2t_clean.fasta` are **always preserved**.
  - Redundancy reduction (`redundans`) runs **only** on non‑T2T contigs, with an optional `minimap2` filter to drop non‑T2T contigs redundant to T2T.
  - Contigs are recombined and length‑sorted for the final FASTA.
- **Telomere QC & final summary** (14–17): telomere counts, telomeric FASTA, and consolidated reports.
- **Reproducibility**: each run writes tool versions to `version.txt` and step-by-step logs to `logs/`.

---

## Installation & Environments

```bash
conda env create -f dependency/pacbiohifi.yml
conda env create -f dependency/busco.yml
conda env create -f dependency/quast.yml
```
TAMP activates the appropriate environment automatically during steps.

---

## Quick Start

**Full run (HiFi reads):**
```bash
bash TAMP-0.2.7.sh --fastq reads.fastq.gz -g 90m -t 20 -m AACCCT --busco ascomycota_odb10
```

**Resume from Step 7:**
```bash
bash TAMP-0.2.7.sh --fasta genome.fa -g 90m -t 20 -m AACCCT -s 7-17
```

**Run merge only (Step 12+) with a prebuilt T2T:**
```bash
bash TAMP-0.2.7.sh -g 90m -t 32 --fasta myassembly.fa -s 12-17
# Non-interactive choose:
bash TAMP-0.2.7.sh -g 90m -t 32 --fasta myassembly.fa -s 12 --choose flye
```

---

## Pipeline Steps (summary)

| Step | Purpose |
|------|--------|
| 1–6  | Assemble with HiCanu, NextDenovo, Peregrine, IPA, Flye, RAFT‑Hifiasm |
| 7    | Copy to `assemblies/` and normalize headers (stable, length‑sorted IDs) |
| 8–9  | Detect telomeres per assembler (`seqtk telo`) and create telomeric subsets |
| 10   | QUAST metrics (`quast_final/`) |
| 11   | (If present) Create assembly summaries |
| 12   | **Rebuild `assemblies/assembly_info.csv` and perform final merge with T2T protection** |
| 13   | BUSCO completeness |
| 14–17| Telomere QC, counts, telomeric FASTA, and summary tables |

---

## Step 12 Details — Rebuild, Compare, Merge

1. **Rebuild summary table first**  
   At the start of Step 12, TAMP **rebuilds** `assemblies/assembly_info.csv` by merging any of these that exist:
   - `assemblies/assembly.busco.csv` (or common alternates)
   - `assemblies/assembly.quast.csv`
   - `assemblies/assembly.telo.csv`  
   The builder is Python‑based in v0.2.7+, so it’s robust to CRLF line endings and awk variants.

2. **Choose an assembly**  
   - By default (v0.2.7), Step 12 **auto-selects** the assembler with the highest N50 in `assemblies/assembly_info.csv`.
   - You can override this with `--choose <assembler>`.
   - If auto-selection cannot decide (e.g., missing N50 row), an interactive prompt lists assemblies present in `assemblies/*.result.fasta`.

3. **T2T protection and merge**  
   - **All** contigs from `t2t_clean.fasta` are protected and **kept as‑is**.
   - The chosen assembly’s non‑T2T contigs are reduced with `redundans` (T2T is never passed into `redundans`).  
   - (Recommended) `minimap2 -x asm20` can drop non‑T2T contigs that are redundant vs. T2T (e.g., identity ≥ 0.95; covered fraction ≥ 0.95).  
   - Protected T2T + reduced “others” → recombined and length‑sorted → `assemblies/final.merged.fasta`.

---

## Telomere Metrics: Definitions (used in Steps 9 & 14)

- **Double‑end contig**: a contig with telomeric signal at **both** ends (has at least one hit where `start == 0` and at least one hit where `end == length`).  
- **Single‑end contig**: telomeric signal at **exactly one** end.

These definitions are enforced in v0.2.2+ to avoid under‑counting double‑ended contigs.

---

## Outputs

- `assemblies/*.result.fasta` — normalized per‑assembler FASTAs  
- `assemblies/final.merged.fasta` — final assembly  
- `assemblies/assembly_info.csv` — unified matrix of BUSCO/QUAST/TELO metrics  
- `assemblies/merged.telo.csv` — final telomere metrics table  
- `assemblies/merged.busco.csv` — final BUSCO metrics table  
- `assemblies/merged.quast.csv` — final QUAST metrics table  
- `merged_result.csv` — consolidated comparison of assemblies (from Step 16)  
- `quast_final/` — QUAST outputs  
- `logs/step_<N>.log` — logs per step  
- `version.txt` — tool versions captured during the run

---

## Troubleshooting

- **No assembler found / empty files**: check paths and that each assembler completed.  
- **Step 12 table looks empty**: confirm at least one of BUSCO / QUAST / TELO CSVs exists in `assemblies/`.  
- **Seqtk not found**: ensure the `funannotate` (or appropriate) env is active; TAMP will try to auto‑activate.  
- **Line‑ending issues**: v0.2.7+ Python builder strips CRs; previous awk errors (e.g., `gsub(/` at Step 12) are resolved.  
- **T2T contig names changed by a merger**: use the alignment‑based T2T split path to protect T2T contigs by alignment rather than name.

---

## Citation

If you use **TAMP v0.2.7**, please cite:

> Sun, Y. (2025). *TAMP: Telomere Assembly Merge Pipeline v0.2.7.*  
> Grainger Bioinformatics Center, Field Museum of Natural History.

---

## Changelog

- **v0.2.7** — *Final merge & summaries*: Auto-selects the assembler with the highest N50 in Step 12 (unless `--choose` is given); runs `redundans` with minimap2-based reduction on non‑T2T contigs; writes telomere/BUSCO/QUAST summary tables as `assemblies/merged.*.csv` and a consolidated `merged_result.csv`.
- **v0.2.6** — *Step 12 builder hardened*: Replaced awk table builder with **Python‑based** `build_assembly_info_v2` (CRLF‑safe, header‑normalized) and invoked it at the start of Step 12.  
- **v0.2.5** — *Step 12 robustness*: Made CR removal explicit to avoid awk `/.../` parse errors on some platforms.  
- **v0.2.4** — *Step 12 flow*: Always rebuild `assemblies/assembly_info.csv` from BUSCO/QUAST/TELO prior to prompting/`--choose`; print the matrix.  
- **v0.2.3** — *Step 7 bugfix*: Removed stray Bash call embedded in the Python heredoc of `rename_and_sort_fasta` (no more `SyntaxError` in `/tmp/rename_fa*.py`).  
- **v0.2.2** — *Telomere counts*: Fixed double‑end logic in Steps 9 & 14 to require signal at **both** ends; single‑end unchanged.  
- **v0.2.1** — *T2T protection*: Step 12 ensures **T2T contigs are preserved**; run `redundans` only on non‑T2T contigs; keep `--choose` and interactive prompt.
