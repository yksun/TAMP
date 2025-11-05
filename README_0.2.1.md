# 🧬 TAMP: Telomere Assembly Merge Pipeline (v0.2.1)

**TAMP v0.2.1** is a maintenance update to the v0.2 release that clarifies and hardens **Step 12** so that *all contigs from `t2t_clean.fasta` are preserved*, while redundancy reduction is applied only to non‑T2T contigs. The interactive assembler selection and `--choose` flag remain supported.

---

## 🚀 What’s new in v0.2.1

- **T2T protection:** contigs from `t2t_clean.fasta` are *never* passed into `redundans`. Redundancy reduction runs only on the “others” set; T2T is re‑added unchanged.
- **Optional alignment filter:** non‑T2T contigs that are redundant vs. T2T (e.g., `minimap2 -x asm20`, identity ≥ 0.95 and covered fraction ≥ 0.95) are dropped before recombining.
- **Choice preserved:** interactive prompt or `--choose <assembler>` still selects the non‑T2T assembly to merge with T2T.
- **Docs/usage refresh:** examples updated to use `TAMP-0.2.1.sh`. Changelog header embedded in the script.

> No other pipeline behavior changes: outputs, dependency management, and step numbering remain consistent with v0.2.

---

## 📁 Output Layout

```
final_results/
├── final.merged.fasta
├── assembly_info.csv
├── assembly_info.md
└── t2t_clean.fasta
```
Other folders:
- `logs/step_*.log` — detailed per-step logs  
- `version.txt` — software version tracking (auto-generated)  
- `assemblies/` — per-assembler outputs

---

## ⚙️ Dependencies

Set up environments:
```bash
conda env create -f dependency/pacbiohifi.yml
conda env create -f dependency/busco.yml
conda env create -f dependency/quast.yml
```
The script activates environments automatically as needed.

---

## 🧩 Workflow Summary

| Step | Function |
|------|----------|
| 1–6  | Run assemblers: HiCanu, NextDenovo, Peregrine, IPA, Flye, RAFT‑Hifiasm |
| 7    | Copy assemblies, normalize contigs |
| 8–9  | Telomere detection & preliminary merge |
| 10   | QUAST quality assessment |
| 11   | Generate `assembly_info.csv` & `assembly_info.md` |
| 12   | **Final assembly merge (T2T contigs protected)** |
| 13   | BUSCO completeness evaluation |
| 14–17| Telomere QC & final summary |

---

## 🔐 Step 12 (v0.2.1): “Protect T2T, reduce others, recombine”

1. **Select the non‑T2T assembly**  
   - Interactive prompt lists available assemblies in `assemblies/` or pass `--choose <assembler>`.

2. **Merge chosen assembly with T2T (name‑stable)**  
   - Merge while preserving original T2T contig names (required for name‑based split).

3. **Split by contig origin**  
   - **Protected:** contigs whose names occur in `t2t_clean.fasta` → `assemblies/protected.t2t.fa`  
   - **Others:** everything else → `assemblies/others.fa`

4. **Reduce only “others” with redundans**  
   ```bash
   redundans.py --noscaffolding --nogapclosing -t <threads> \
     -f assemblies/others.fa --identity 0.50 --overlap 0.80
   # → redundans/scaffolds.reduced.fa
   ```

5. **(Recommended) Drop “others” redundant to T2T via alignment**  
   ```bash
   minimap2 -x asm20 -t <threads> protected.t2t.fa redundans/scaffolds.reduced.fa \
     | <parse-and-keep-only-nonredundant>  # keep others not covering T2T at ≥0.95/0.95
   # → assemblies/others.filtered.fa
   ```

6. **Recombine & sort**  
   ```bash
   cat assemblies/protected.t2t.fa assemblies/others.filtered.fa > assemblies/merged_protected_priority.fa
   funannotate sort -i assemblies/merged_protected_priority.fa -b contig -o final_results/final.merged.fasta --minlen 500
   ```

> If the merger renames contigs, use an alignment‑based split against `t2t_clean.fasta` instead of name matching.

---

## 🧾 Version Tracking

Example `version.txt`:
```
Pipeline: TAMP-0.2.1.sh
Run ID:   20251105_104139
Date:     2025-11-05 10:41:39
Host:     <hostname>

Software versions:
canu: v2.3
nextDenovo: v3.6
peregrine: v0.5.3
ipa: v1.3
flye: v2.9
hifiasm: v0.19
busco: v5.7.1
quast: v5.2.0
minimap2: 2.28-r1209
bwa: v0.7.17
samtools: v1.13
python3: 3.10.14
```

---

## 🧪 Example Usage

Full pipeline:
```bash
bash TAMP-0.2.1.sh --fastq reads.fastq.gz -g 90m -t 20 -m AACCCT --busco ascomycota_odb10
```

Continue from Step 7:
```bash
bash TAMP-0.2.1.sh --fasta genome.fa -g 90m -t 20 -m AACCCT -s 7-17
```

Merge & finalize (Step 12+):
```bash
bash TAMP-0.2.1.sh -g 90m -t 32 --fasta myassembly.fa --busco ascomycota_odb10 -s 12-17
# Non-interactive choice example:
bash TAMP-0.2.1.sh -g 90m -t 32 --fasta myassembly.fa -s 12 --choose flye
```

---

## 📊 Final Deliverables

| File | Description |
|------|-------------|
| `final.merged.fasta` | Final merged genome |
| `assembly_info.csv`  | Metrics per assembler |
| `assembly_info.md`   | Summary (Markdown) |
| `t2t_clean.fasta`    | Telomere-to-telomere cleaned FASTA |

---

## 🧯 Troubleshooting

- Logs in `logs/step_<N>.log`
- Check conda envs if tools show as NOT FOUND
- BUSCO lineage: verify dataset (e.g., `ascomycota_odb10`) exists
- If `version.txt` missing entries, rerun within correct environment

---

## 📜 Citation

If you use **TAMP v0.2.1**, please cite:

> Sun, Y. (2025). *TAMP: Telomere Assembly Merge Pipeline v0.2.1.*  
> Grainger Bioinformatics Center, Field Museum of Natural History.

---
Generated on 2025-11-05
