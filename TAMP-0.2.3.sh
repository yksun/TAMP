#!/bin/bash
# ------------------------------------------------------------
# v0.2.3 (2025-11-05)
# - FIX: Step 7 'rename_and_sort_fasta': removed shell call inside Python heredoc
#        and log Python version from bash before invoking Python (no more SyntaxError).
# ------------------------------------------------------------
# ------------------------------------------------------------
# v0.2.2 (2025-11-05)
# - FIX: Telomere double-end logic in Steps 9 & 14: now counts contigs with telomeres at both ends.
# ------------------------------------------------------------
# ------------------------------------------------------------
# v0.2.1 (2025-11-05)
# - Version bump from 0.2 → 0.2.1.
# - Step 12 note: t2t_clean.fasta contigs are preserved; redundans runs only on non-T2T contigs.
#   (--choose still supported; interactive prompt retained.)
# ------------------------------------------------------------

# ===== Version logging helpers =====
VERSION_FILE="${VERSION_FILE:-version.txt}"

init_version_file() {
  if [[ ! -f "$VERSION_FILE" ]]; then
    {
      echo "Pipeline: ${0##*/}"
      echo "Run ID:   $(date +%Y%m%d_%H%M%S)"
      echo "Date:     $(date '+%Y-%m-%d %H:%M:%S')"
      echo "Host:     $(hostname)"
      echo
      echo "Software versions:"
    } > "$VERSION_FILE"
  fi
}

log_version() {
  # usage: log_version <label> <cmd>
  local label="$1"; shift
  local cmd="$1"; shift || true
  init_version_file
  if ! grep -q "^${label}:" "$VERSION_FILE" 2>/dev/null; then
    local ver="NOT FOUND"
    if command -v "$cmd" >/dev/null 2>&1; then
      local out=""
      for probe in \
        "$cmd --version" \
        "$cmd -V" \
        "$cmd version" \
        "$cmd -v" \
        "$cmd 2>&1"; do
        out=$(eval $probe 2>&1 | head -n 1) || true
        if [[ -n "$out" ]]; then ver="$out"; break; fi
      done
    fi
    echo "${label}: ${ver}" >> "$VERSION_FILE"
  fi
}
# ===================================

# --- helper: rename contigs to <prefix>_1..N by descending length (wrap at 60bp) ---
rename_and_sort_fasta() {
  # usage: rename_and_sort_fasta INPUT.fa OUTPUT.fa PREFIX
  local infa="$1"; local outfa="$2"; local prefix="$3"
  if [[ ! -s "$infa" ]]; then
    echo "[warn] rename_and_sort_fasta: input '$infa' missing or empty; skipping." >&2
    return 0
  fi
  local py="$(mktemp -t rename_fa.XXXXXX.py)"
  cat > "$py" <<'PY'
import sys

def read_fasta(fp):
  name=None; seq=[]
  for line in fp:
    if line.startswith('>'):
      if name is not None:
        yield name, ''.join(seq)
      name=line[1:].strip()
      seq=[]
    else:
      seq.append(line.strip())
  if name is not None:
    yield name, ''.join(seq)

def wrap(s, w=60):
  return '\n'.join(s[i:i+w] for i in range(0, len(s), w)) if s else ''

inp, outp, prefix = sys.argv[1:4]
with open(inp, 'r') as f:
  recs = list(read_fasta(f))

# sort by length desc, stable
recs.sort(key=lambda x: len(x[1]), reverse=True)

with open(outp, 'w') as o:
  for i, (name, seq) in enumerate(recs, start=1):
    o.write(f">{prefix}_{i}\n")
    o.write(wrap(seq) + ("\n" if seq and not seq.endswith("\n") else ""))

PY
  log_version "python3" "python3"
  python3 "$py" "$infa" "$outfa" "$prefix"
  local rc=$?
  rm -f "$py"
  return $rc
}
# --- end helper ---


# ==== Interactive prompt helpers (inserted by fix) ====
# Open a dedicated read-only file descriptor to the terminal for robust prompts
if [ -z "${__TTY_FD3_OPENED__:-}" ] && [ -r /dev/tty ]; then
  exec 3</dev/tty
  __TTY_FD3_OPENED__=1
fi

# Print a prompt directly to the terminal and read from terminal (works even in pipelines/tee)
prompt_from_tty() {
  # usage: prompt_from_tty VAR "Your prompt: "
  local __outvar=$1; shift
  local __prompt="$*"
  if [ -r /dev/tty ]; then
    # show prompt immediately on real terminal
    printf "%s" "$__prompt" > /dev/tty
    local __ans=""
    if [ -e /proc/self/fd/3 ]; then
      IFS= read -r __ans <&3 || __ans=""
    else
      IFS= read -r __ans < /dev/tty || __ans=""
    fi
    # Trim leading/trailing whitespace
    __ans="$(printf "%s" "$__ans" | awk '{$1=$1};1')"
    printf -v "$__outvar" '%s' "$__ans"
    return 0
  fi
  return 1
}
# ==== end helpers ====

# Version: TAMP-0.1

# Default values
genomesize=""
threads=20
fastq=""
steps=()
motif="AACCCT"
# assembler default changed to empty to force prompt when --choose not given
assembler="" # Default assembler
external_fasta=""  # Optional external pre-assembled FASTA provided by user
run_busco=false    # Whether to run BUSCO on individual assemblies
busco_lineage="ascomycota_odb10"  # Default BUSCO lineage
# --- Logging & versioning (added by TAMP-0.1) ---
PIPELINE_NAME="TAMP-0.1.sh"
RUN_ID="$(date +'%Y%m%d_%H%M%S')"
LOG_DIR="logs"
mkdir -p "$LOG_DIR"

# awk program to prefix each line with a timestamp
LOG_AWK='{ print strftime("[%Y-%m-%d %H:%M:%S]"), $0; fflush(); }'

timestamp() { date +"%Y-%m-%d %H:%M:%S"; }

write_versions() {
  local vf="version.txt"
  {
    echo "Pipeline: ${PIPELINE_NAME}"
    echo "Run ID:   ${RUN_ID}"
    echo "Date:     $(timestamp)"
    echo "Host:     $(hostname 2>/dev/null || echo unknown)"
    echo
    echo "Software versions:"
  } > "$vf"

  getv() {
    local c="$1"
    if ! command -v "$c" >/dev/null 2>&1; then
      echo "$c: NOT FOUND"
      return
    fi
    # Try common flags
    local out=""
    for f in "--version" "-V" "-v" "version"; do
      out=$("$c" $f 2>&1 | head -n1)
      if [[ -n "$out" ]]; then
        echo "$c: $out"
        return
      fi
    done
    echo "$c: FOUND (version unknown)"
  }

  # Common tools in this pipeline; missing ones will be noted as NOT FOUND
  for c in canu nextDenovo peregrine ipa flye hifiasm seqtk busco quast.py quast minimap2 merge_wrapper.py raft bwa samtools python3; do
    getv "$c" >> "$vf"
  done

  echo "" >> "$vf"
}



# Function to display usage
usage() {
  cat <<USAGE
Usage: ${PIPELINE_NAME:-$0} -g <genomesize> -t <threads> --fastq <fastq> -m <motif> [-s <steps>] [--fasta <path>] [--busco <lineage>] [--choose]

Options:
  -g, --genomesize   Genome size (required), e.g. 2g, 500m
  -t, --threads      Number of threads (required)
  --fastq            Path to the FASTQ file (required)
  -m, --motif        Telomere motif (required), e.g. AACCCT
  -s, --steps        Steps to run (optional, default: all). Accepts comma/range list (e.g. 1,2,5-7)
  --fasta            External pre-assembled FASTA (optional, single file or URL)
  --busco            Run BUSCO on each individual assembly; optional lineage (default: ascomycota_odb10)
  --choose           Prompt to choose an assembler for the final merge (default: peregrine)

Notes:
  • Per-step logs: logs/step_<N>.log with timestamps.
  • A consolidated software version summary is written to ./version.txt at start.
  • Example: ${PIPELINE_NAME:-$0} -g 2g -t 16 --fastq reads.fastq -m AACCCT -s 1,3-5 --choose

Steps:
  1. HiCanu assembly
  2. NextDenovo assembly
  3. Peregrine assembly
  4. IPA assembly
  5. Flye assembly
  6. Hifiasm assembly
  7. Copy all assemblies
  8. BUSCO on all assemblies (including external)
  9. Telomere contigs + metrics
  10. Merge all assemblies
  11. QUAST for all assembler results
  12. Final merge using selected assembler
  13. BUSCO analysis (final)
  14. Telomere analysis (final)
  15. QUAST analysis (final)
  16. Generate final comparison report
  17. Cleanup temporary files into structured folders
USAGE
  exit 1
}

# Function to expand ranges in the step input (e.g., "1,3-5" becomes "1 3 4 5")
expand_steps() {
  local IFS=','; read -ra ranges <<< "$1"
  for range in "${ranges[@]}"; do
    if [[ $range == *"-"* ]]; then
      local start=${range%-*}
      local end=${range#*-}
      for ((i=start; i<=end; i++)); do
        steps+=("$i")
      done
    else
      steps+=("$range")
    fi
  done
}

# Parse command-line options
while getopts ":g:t:s:m:-:" opt; do
  case $opt in
    g) genomesize="$OPTARG"
    ;;
    t) threads="$OPTARG"
    ;;
    s) expand_steps "$OPTARG"
    ;;
    m) motif="$OPTARG"
    ;;
    -) case "${OPTARG}" in
         fastq) fastq="${!OPTIND}"; OPTIND=$((OPTIND + 1))
         ;;
         fasta) external_fasta="${!OPTIND}"; OPTIND=$((OPTIND + 1))
         ;;
         busco)
           run_busco=true
           next_arg="${!OPTIND}"
           if [[ -n "$next_arg" && "$next_arg" != -* ]]; then
             busco_lineage="$next_arg"
             OPTIND=$((OPTIND + 1))
           fi
         ;;
         choose) 
           echo "Please enter the assembler you want to use for the final merge (caun, nextDenovo, peregrine, ipa, flye, RAFT-hifiasm):"
           read assembler
           ;;
         *) echo "Unknown option --${OPTARG}" >&2; usage
         ;;
       esac
    ;;
    \?) echo "Invalid option: -$OPTARG" >&2
        usage
    ;;
    :) echo "Option -$OPTARG requires an argument." >&2
       usage
    ;;
  esac
done

# Ensure that required options are set
if [[ -z "$genomesize" || -z "$threads" || -z "$fastq" || -z "$motif" ]]; then
  usage
fi

# Write versions summary at the start of a run
write_versions



# Optional: resolve external FASTA (download if URL) AFTER required args check
if [[ -n "$external_fasta" ]]; then
  if [[ "$external_fasta" =~ ^https?:// || "$external_fasta" =~ ^ftp:// ]]; then
    echo "[info] --fasta looks like a URL, downloading..."
    tmp_name="external_input.fasta"
    if command -v curl >/dev/null 2>&1; then
      curl -L "$external_fasta" -o "$tmp_name" || { echo "[error] Failed to download $external_fasta"; exit 1; }
    elif command -v wget >/dev/null 2>&1; then
      wget -O "$tmp_name" "$external_fasta" || { echo "[error] Failed to download $external_fasta"; exit 1; }
    else
      echo "[error] Neither curl nor wget found for downloading $external_fasta"
      exit 1
    fi
    external_fasta="$tmp_name"
  fi
  # If gzipped, decompress to .fa
  if [[ "$external_fasta" =~ \.gz$ ]]; then
    echo "[info] Decompressing gzipped FASTA: $external_fasta"
    gz_out="external_input.fa"
    gunzip -c "$external_fasta" > "$gz_out" || { echo "[error] gunzip failed for $external_fasta"; exit 1; }
    external_fasta="$gz_out"
  fi
  # Existence check
  if [[ ! -s "$external_fasta" ]]; then
    echo "[error] --fasta provided but file not found or empty: $external_fasta"
    exit 1
  fi
fi
project=$(basename "$fastq" .fastq)

# Modify run.cfg with genomesize and threads
cat <<EOT >> ./run_${project}.cfg
[General]
job_type = local # local, slurm, sge, pbs, lsf
job_prefix = nextDenovo
task = all # all, correct, assemble
rewrite = yes # yes/no
deltmp = yes
parallel_jobs = ${threads} # number of tasks used to run in parallel
input_type = raw # raw, corrected
read_type = hifi # clr, ont, hifi
input_fofn = ./input_${project}.fofn
workdir = NextDenovo

[correct_option]
read_cutoff = 1k
genome_size = ${genomesize} # estimated genome size
sort_options = -m 20g -t ${threads}
minimap2_options_raw = -t ${threads}
pa_correction = 3 # number of corrected tasks used to run in parallel, each corrected task requires ~TOTAL_INPUT_BASES/4 bytes of memory usage.
correction_options = -p ${threads}

[assemble_option]
minimap2_options_cns = -t ${threads}
nextgraph_options = -a 1

# see https://nextdenovo.readthedocs.io/en/latest/OPTION.html for a detailed introduction about all the parameters"
EOT

echo $PWD/${fastq} > input_${project}.fofn
echo $PWD/${fastq} > reads_${project}.lst

# Function to check the success of the previous command
check_command() {
  if [ $? -ne 0 ]; then
    echo "Error: Command failed. Exiting."
    exit 1
  fi
  if [[ -n "${1:-}" && "${1#-}" = "$1" ]]; then
    if [[ "$1" != *" "* ]]; then
      log_version "$1" "$1"
    fi
  fi
}


echo "Activating assembly environment"
eval "$(conda shell.bash hook)"
conda activate pacbiohifi
check_command

# If no specific steps provided, default to running all steps (1-17)
if [ ${#steps[@]} -eq 0 ]; then
  steps=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17)
fi

# Execute the specified steps
for step in "${steps[@]}"; do
  STEP_LOG="${LOG_DIR}/step_${step}.log"
  {
    echo "===== [$(timestamp)] STEP ${step} START ====="
    case $step in
    1)
log_version "canu" "canu"
      echo "Step 1 - Assembly of the genome using HiCanu"
      canu -p canu -d hicanu genomeSize=$genomesize maxThreads=$threads -pacbio-hifi $fastq
      check_command
      ;;
    2)
log_version "nextDenovo" "nextDenovo"
      echo "Step 2 - Assembly of the genome using NextDenovo"
      nextDenovo run_${project}.cfg
      check_command
      ;;
    3)
      echo "Step 3 - Assembly of the genome using Peregrine"
      pg_asm reads_${project}.lst peregrine-2021
      check_command
      ;;
    4)
log_version "ipa" "ipa"
      echo "Step 4 - Assembly of the genome using IPA"
      ipa local --nthreads $threads --njobs 1 --run-dir ipa -i $fastq
      check_command
      ;;
    5)
log_version "flye" "flye"
      echo "Step 5 - Assembly of the genome using Flye"
      flye --pacbio-hifi $fastq --out-dir flye --threads $threads
      check_command
      ;;
    6)
      echo "Step 6 - Assembly of the genome using Hifiasm"
      mkdir hifiasm
log_version "hifiasm" "hifiasm"
      cd hifiasm
      hifiasm -o errorcorrect -t$threads --write-ec ../$fastq 2> errorcorrect.log
      check_command
      COVERAGE=$(grep "homozygous" errorcorrect.log | tail -1 | awk '{print $6}')
      hifiasm -o getOverlaps -t$threads --dbg-ovec errorcorrect.ec.fa 2> getOverlaps.log
      check_command
      cat getOverlaps.0.ovlp.paf getOverlaps.1.ovlp.paf > overlaps.paf
      ~/opt/RAFT/raft -e ${COVERAGE} -o fragmented errorcorrect.ec.fa overlaps.paf
      hifiasm -o finalasm -t$threads -r1 fragmented.reads.fasta 2> finalasm.log
      check_command
      awk '/^S/{print ">"$2;print $3}' finalasm.bp.hap1.p_ctg.gfa > RAFT-hifiasm.fasta
      cd ..
      ;;    7)
      echo "Step 7 - Copy all assemblies"
      mkdir -p assemblies
      cp ./hicanu/canu.contigs.fasta ./assemblies/canu.result.fasta
      # Normalize headers: canu_1..N by contig length
      tmp_renamed="assemblies/.canu.renamed.tmp.fasta"
      rename_and_sort_fasta "./assemblies/canu.result.fasta" "$tmp_renamed" "canu" && mv -f "$tmp_renamed" "./assemblies/canu.result.fasta"

      cp ./NextDenovo/03.ctg_graph/nd.asm.fasta ./assemblies/nextDenovo.result.fasta
      # Normalize headers: nextDenovo_1..N by contig length
      tmp_renamed="assemblies/.nextDenovo.renamed.tmp.fasta"
      rename_and_sort_fasta "./assemblies/nextDenovo.result.fasta" "$tmp_renamed" "nextDenovo" && mv -f "$tmp_renamed" "./assemblies/nextDenovo.result.fasta"

      cp ./peregrine-2021/asm_ctgs_m_p.fa ./assemblies/peregrine.result.fasta
      # Normalize headers: peregrine_1..N by contig length
      tmp_renamed="assemblies/.peregrine.renamed.tmp.fasta"
      rename_and_sort_fasta "./assemblies/peregrine.result.fasta" "$tmp_renamed" "peregrine" && mv -f "$tmp_renamed" "./assemblies/peregrine.result.fasta"

      cp ./ipa/assembly-results/final.p_ctg.fasta ./assemblies/ipa.result.fasta
      # Normalize headers: ipa_1..N by contig length
      tmp_renamed="assemblies/.ipa.renamed.tmp.fasta"
      rename_and_sort_fasta "./assemblies/ipa.result.fasta" "$tmp_renamed" "ipa" && mv -f "$tmp_renamed" "./assemblies/ipa.result.fasta"

      cp ./flye/assembly.fasta ./assemblies/flye.result.fasta
      # Normalize headers: flye_1..N by contig length
      tmp_renamed="assemblies/.flye.renamed.tmp.fasta"
      rename_and_sort_fasta "./assemblies/flye.result.fasta" "$tmp_renamed" "flye" && mv -f "$tmp_renamed" "./assemblies/flye.result.fasta"

      cp ./hifiasm/RAFT-hifiasm.fasta ./assemblies/RAFT-hifiasm.result.fasta
      # Normalize headers: RAFT-hifiasm_1..N by contig length
      tmp_renamed="assemblies/.RAFT-hifiasm.renamed.tmp.fasta"
      rename_and_sort_fasta "./assemblies/RAFT-hifiasm.result.fasta" "$tmp_renamed" "RAFT-hifiasm" && mv -f "$tmp_renamed" "./assemblies/RAFT-hifiasm.result.fasta"

      # If user provided an external FASTA, include it as an assembly, too
      if [[ -n "$external_fasta" ]]; then
        # Normalize/resolve if relative path
        ext_basename=$(basename "$external_fasta")
        cp "$external_fasta" "./assemblies/external.result.fasta"
      fi

      # === Run BUSCO on all assembled genomes (including external) ===
      # Requires helper: run_busco_on_assembly()
      ;;

8)
echo "Step 8 - Run BUSCO on all assembled genomes (including external)"
eval "$(conda shell.bash hook)"; conda deactivate 2>/dev/null || true; conda activate busco || true

shopt -s nullglob
a=(assemblies/*.result.fasta)
if [[ ${#a[@]} -eq 0 ]]; then
  echo "[error] No assemblies in ./assemblies for BUSCO."
  exit 1
fi

# Fixed output order to match other matrices
cols=("canu" "external" "flye" "ipa" "nextDenovo" "peregrine" "RAFT-hifiasm")

# Settings
lineage=${busco_lineage:-"ascomycota_odb10"}   # default lineage
threads=${threads:-8}

# Put all BUSCO results here (ok if already exists)
[[ -d busco ]] || mkdir -p busco

has_busco_metrics() {
  local base="$1"
  # Look under busco/<base> and legacy busco/run_<base> for short_summary or full_table
  if [[ -d "busco/${base}" ]]; then
    if find "busco/${base}" -maxdepth 3 -type f \( -name 'short_summary*.txt' -o -name 'full_table*.tsv' \) | grep -q .; then
      return 0
    fi
  fi
  if [[ -d "busco/run_${base}" ]]; then
    if find "busco/run_${base}" -maxdepth 3 -type f \( -name 'short_summary*.txt' -o -name 'full_table*.tsv' \) | grep -q .; then
      return 0
    fi
  fi
  return 1
}

BUSCO_BIN="$(command -v busco || true)"
if [[ -z "$BUSCO_BIN" ]]; then
  echo "[warn] 'busco' not found on PATH; will only parse existing results under ./busco" >&2
fi

for fasta in "${a[@]}"; do
  base=$(basename "$fasta"); base="${base%.result.fasta}"; base="${base%.fasta}"

  if has_busco_metrics "$base"; then
    echo "[ok] Found existing BUSCO metrics for $base → using busco/${base} (or legacy path)"
    continue
  fi

  if [[ -z "$BUSCO_BIN" ]]; then
    echo "[error] Missing BUSCO metrics for '$base' and 'busco' binary not available; cannot run." >&2
    exit 127
  fi

  # If an output dir exists but is incomplete, re-run with -f to overwrite
  force_flag=""
  if [[ -d "busco/${base}" ]]; then
    echo "[info] Incomplete BUSCO dir detected for $base → re-running with -f"
    force_flag="-f"
  fi

  echo "[run] BUSCO on $base (lineage=$lineage, threads=$threads)"
log_version "busco" "busco"
  pushd busco >/dev/null
  busco -i "../$fasta" -l "$lineage" -m genome -c "$threads" -o "$base" $force_flag 2>&1 | tee "busco_${base}.log" || true
  popd >/dev/null
done

# === Build assemblies/assembly.busco.csv with metrics as rows (use a temp Python file) ===
pyfile="$(mktemp -t busco_matrix.XXXXXX.py)"
cat > "$pyfile" <<'PY'
import os, re, glob, csv

OUT_CSV = os.path.join("assemblies", "assembly.busco.csv")
DESIRED = ["canu","external","flye","ipa","nextDenovo","peregrine","RAFT-hifiasm"]

def newest(paths):
    return max(paths, key=os.path.getmtime) if paths else None

def find_run_root(base):
    # Prefer busco/<base>, fallback to busco/run_<base>, then any close match
    candidates = []
    p1 = os.path.join("busco", base)
    p2 = os.path.join("busco", f"run_{base}")
    if os.path.isdir(p1): candidates.append(p1)
    if os.path.isdir(p2): candidates.append(p2)
    if not candidates:
        # loose match (handles RAFT-hifiasm variants)
        candidates = sorted(glob.glob(os.path.join("busco", f"{base}*"))) + \
                     sorted(glob.glob(os.path.join("busco", f"run_{base}*")))
    return newest(candidates)

def read_counts_from_anywhere(run_root):
    if not run_root: return None
    # Look for full_table first (best source), searching recursively
    fulls = glob.glob(os.path.join(run_root, "**", "full_table*.tsv"), recursive=True)
    if fulls:
        p = newest(fulls)
        with open(p, newline="") as f:
            lines = [ln.rstrip("\n") for ln in f]
        status_idx = None
        for ln in lines:
            if ln.startswith("#") and "Status" in ln:
                hdr = ln.lstrip("#").strip().split("\t")
                for i, h in enumerate(hdr):
                    if h.strip().lower() == "status":
                        status_idx = i
                break
        S = D = F = M = 0
        n = 0
        for ln in lines:
            if not ln or ln.startswith("#"): continue
            parts = ln.split("\t")
            st = parts[1].strip() if status_idx is None else (parts[status_idx].strip() if len(parts) > status_idx else "")
            n += 1
            if   st == "Complete":   S += 1
            elif st == "Duplicated": D += 1
            elif st == "Fragmented": F += 1
            elif st == "Missing":    M += 1
        if n == 0: return None
        return {"S": S, "D": D, "F": F, "M": M, "n": n}
    # Fallback: parse short_summary
    sums = glob.glob(os.path.join(run_root, "**", "short_summary*.txt"), recursive=True)
    if sums:
        txt = open(newest(sums), "r", errors="ignore").read()
        m = re.search(r"C:(\d+(?:\.\d+)?)%.*?S:(\d+(?:\.\d+)?)%.*?D:(\d+(?:\.\d+)?)%.*?F:(\d+(?:\.\d+)?)%.*?M:(\d+(?:\.\d+)?)%.*?n:(\d+)", txt, re.S)
        if not m: return None
        Cpct, Spct, Dpct, Fpct, Mpct = map(float, m.groups()[:5])
        n = int(m.group(6))
        # Approximate counts from percents
        Ccnt = round(n * Cpct / 100.0)
        Mcnt = round(n * Mpct / 100.0)
        return {"Cpct": Cpct, "Spct": Spct, "Dpct": Dpct, "Fpct": Fpct, "Mpct": Mpct, "Ccnt": Ccnt, "Mcnt": Mcnt, "n": n}
    return None

def safe_pct(x, n): return 0.0 if not n else 100.0 * float(x) / float(n)
def fmt_pct(x):
    s = f"{x:.1f}"
    return s.rstrip('0').rstrip('.') if '.' in s else s

results = {}
for asm in DESIRED:
    run_root = find_run_root(asm)
    counts = read_counts_from_anywhere(run_root)
    if counts and "S" in counts:
        S, D, F, M, n = counts["S"], counts["D"], counts["F"], counts["M"], counts["n"]
        C = S + D
        results[asm] = {
            "Cpct": fmt_pct(safe_pct(C, n)),
            "Spct": fmt_pct(safe_pct(S, n)),
            "Dpct": fmt_pct(safe_pct(D, n)),
            "Fpct": fmt_pct(safe_pct(F, n)),
            "Mpct": fmt_pct(safe_pct(M, n)),
            "Ccnt": str(C), "Mcnt": str(M), "n": str(n),
        }
    elif counts and "Cpct" in counts:
        results[asm] = {
            "Cpct": fmt_pct(counts["Cpct"]), "Spct": fmt_pct(counts["Spct"]), "Dpct": fmt_pct(counts["Dpct"]),
            "Fpct": fmt_pct(counts["Fpct"]), "Mpct": fmt_pct(counts["Mpct"]),
            "Ccnt": str(counts["Ccnt"]), "Mcnt": str(counts["Mcnt"]), "n": str(counts["n"]),
        }
    else:
        results[asm] = {"Cpct":"0","Spct":"0","Dpct":"0","Fpct":"0","Mpct":"0","Ccnt":"0","Mcnt":"0","n":"0"}

rows = []
header = ["Metric"] + DESIRED
rows.append(header)
rows.append(["BUSCO C (%)"]      + [results[a]["Cpct"] for a in DESIRED])
rows.append(["BUSCO S (%)"]      + [results[a]["Spct"] for a in DESIRED])
rows.append(["BUSCO D (%)"]      + [results[a]["Dpct"] for a in DESIRED])
rows.append(["BUSCO F (%)"]      + [results[a]["Fpct"] for a in DESIRED])
rows.append(["BUSCO M (%)"]      + [results[a]["Mpct"] for a in DESIRED])
rows.append(["BUSCO C (count)"]  + [results[a]["Ccnt"] for a in DESIRED])
rows.append(["BUSCO M (count)"]  + [results[a]["Mcnt"] for a in DESIRED])

with open(OUT_CSV, "w", newline="") as f:
    csv.writer(f).writerows(rows)

print(f"Wrote {OUT_CSV}")
PY

python3 "$pyfile"
rm -f "$pyfile"

echo "[ok] Wrote assemblies/assembly.busco.csv"
  ;;

9)
  echo "Step 9 - Extract telomere-containing contigs and compute telomere metrics"
# Ensure assemblies/ exists and is populated even if Step 7 was skipped
mkdir -p assemblies
shopt -s nullglob
existing_assemblies=(assemblies/*.result.fasta)
if [[ ${#existing_assemblies[@]} -eq 0 ]]; then
  # Opportunistically gather any assembler outputs present
  [[ -s ./hicanu/canu.contigs.fasta ]] && cp ./hicanu/canu.contigs.fasta ./assemblies/canu.result.fasta
  [[ -s ./NextDenovo/03.ctg_graph/nd.asm.fasta ]] && cp ./NextDenovo/03.ctg_graph/nd.asm.fasta ./assemblies/nextDenovo.result.fasta
  [[ -s ./peregrine-2021/asm_ctgs_m_p.fa ]] && cp ./peregrine-2021/asm_ctgs_m_p.fa ./assemblies/peregrine.result.fasta
  [[ -s ./ipa/assembly-results/final.p_ctg.fasta ]] && cp ./ipa/assembly-results/final.p_ctg.fasta ./assemblies/ipa.result.fasta
  [[ -s ./flye/assembly.fasta ]] && cp ./flye/assembly.fasta ./assemblies/flye.result.fasta
  [[ -s ./hifiasm/RAFT-hifiasm.fasta ]] && cp ./hifiasm/RAFT-hifiasm.fasta ./assemblies/RAFT-hifiasm.result.fasta
  if [[ -n "$external_fasta" && -s "$external_fasta" ]]; then
    cp "$external_fasta" ./assemblies/external.result.fasta
  fi
fi

# Re-evaluate after attempted population
existing_assemblies=(assemblies/*.result.fasta)
if [[ ${#existing_assemblies[@]} -eq 0 ]]; then
  echo "[error] No assemblies found in ./assemblies. Run step 7 first or supply --fasta." >&2
  exit 1
fi

# Fixed column order for matrix outputs
cols=("canu" "external" "flye" "ipa" "nextDenovo" "peregrine" "RAFT-hifiasm")

# Store metrics in associative arrays keyed by assembler
declare -A tdouble tsingle

# === Single loop over assemblies ===
for fasta in assemblies/*.result.fasta; do
  asm="${fasta##*/}"; asm="${asm%.result.fasta}"          # e.g., canu / RAFT-hifiasm
  list="${fasta%.result.fasta}.telo.list"
  out="${fasta%.result.fasta}.telo.fasta"

log_version "seqtk" "seqtk"
  # 1) Find telomere-containing contigs and write list
  seqtk telo -s 1 -m "$motif" "$fasta" > "$list"
  check_command

  # 2) Unique IDs -> subseq to .telo.fasta
  awk '{print $1}' "$list" | sed 's/[[:space:]].*$//' | tr -d $'\r' | sort -u > "${list}.ids"
  if [[ -s "${list}.ids" ]]; then
    seqtk subseq "$fasta" "${list}.ids" > "$out"
    [[ -s "$out" ]] && echo "[ok] Wrote $(basename "$out")" || echo "[warn] ${out} is empty"
  else
    echo "[warn] No telomere-containing contigs for $(basename "$fasta"); writing empty ${out}"
    : > "$out"
  fi

  # 3) Compute metrics for this assembler
  double=$(awk 'NF>=4 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $4~/^[0-9]+$/ {c=$1; sub(/[ \t].*$/,"",c); if($2==0 && $3==$4) print c}' "$list" | sort -u | wc -l)
  # v0.2.2 override: correct double-end contig counting (both ends flagged)
  double=$(awk '
    NF>=4 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $4~/^[0-9]+$/ {
      c=$1; sub(/[ \t].*$/,"",c);
      s[c]+=($2==0)?1:0; e[c]+=($3==$4)?1:0
    }
    END{ for(c in s){ if(s[c]>0 && e[c]>0) print c } }
  ' "$list" | sort -u | wc -l)
  single=$(awk '
    NF>=4 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $4~/^[0-9]+$/ {
      c=$1; sub(/[ \t].*$/,"",c);
      s[c]+=($2==0)?1:0; e[c]+=($3==$4)?1:0
    }
    END{ for(c in s){ if((s[c]+e[c])==1) print c } }
  ' "$list" | wc -l)

  tdouble["$asm"]="$double"
  tsingle["$asm"]="$single"
done

# 4) Write the single matrix file: assemblies/assembly.telo.csv
matrix_csv="assemblies/assembly.telo.csv"
{
  printf "Metric"
  for c in "${cols[@]}"; do printf ",%s" "$c"; done
  printf "\n"

  # exactly one row per metric
  for metric in "Telomere double-end contigs" "Telomere single-end contigs"; do
    printf "%s" "$metric"
    for c in "${cols[@]}"; do
      if [[ "$metric" == "Telomere double-end contigs" ]]; then
        val="${tdouble[$c]:-0}"
      else
        val="${tsingle[$c]:-0}"
      fi
      printf ",%s" "$val"
    done
    printf "\n"
  done
} > "$matrix_csv"
echo "[ok] Wrote $(basename "$matrix_csv")"

# 5) Write a tidy per-assembler table: assemblies/total_telo.csv
total_csv="assemblies/total_telo.csv"
{
  echo "Assembler,Telomere double-end contigs,Telomere single-end contigs"
  for c in "${cols[@]}"; do
    printf "%s,%s,%s\n" "$c" "${tdouble[$c]}" "${tsingle[$c]}"
  done
} > "$total_csv"
echo "[ok] Wrote $(basename "$total_csv")"

check_command

;;
    10)
      echo "Step 10 - Merge all telo assemblies"
      shopt -s nullglob
      fasta_files=(assemblies/*.telo.fasta)
      if [[ ${#fasta_files[@]} -eq 0 ]]; then
        echo "[info] No per-assembly *.telo.fasta files found; generating per-assembly telomere-contig FASTAs from ./assemblies/*.result.fasta using motif '${motif:-TTAGGG}'."
        shopt -s nullglob
        for f in assemblies/*.result.fasta; do
          [[ -s "$f" ]] || continue
          base=$(basename "$f"); base="${base%.result.fasta}"; base="${base%.fasta}"
          # Build per-assembly telomere FASTA
          M="${motif:-TTAGGG}"
          RC="$(echo "$M" | tr 'ACGTacgtnN' 'TGCAtgcanN' | rev)"
          awk -v W="${TELO_WINDOW:-200}" -v R="${TELO_MIN_REPEATS:-3}" -v M="$M" -v RC="$RC" '
            BEGIN{ pat="(" M "|" RC "){" R ",}" }
            /^>/ { if(seq!=""){ 
                      left=substr(seq,1,W); right=substr(seq,length(seq)-W+1,W);
                      if( (left ~ pat) || (right ~ pat) ){ print header; print seq; }
                   }
                   header=$0; seq=""; next }
            { gsub(/[ \t]/,""); seq=seq $0 }
            END{
              if(seq!=""){
                left=substr(seq,1,W); right=substr(seq,length(seq)-W+1,W);
                if( (left ~ pat) || (right ~ pat) ){ print header; print seq; }
              }
            }' "$f" > "assemblies/${base}.telo.fasta"
          if [[ ! -s "assemblies/${base}.telo.fasta" ]]; then
            rm -f "assemblies/${base}.telo.fasta"
          fi
          # Also emit per-assembly telomere coordinate list using seqtk (matches allmerged.telo.list format)
          if [[ -s "assemblies/${base}.telo.fasta" ]]; then
            eval "$(conda shell.bash hook)"; conda deactivate 2>/dev/null || true; conda activate pacbiohifi || true
            seqtk telo -s 1 -m "$motif" "assemblies/${base}.telo.fasta" > "assemblies/${base}.telo.list" || true
          fi
        done
        # Refresh list
        fasta_files=(assemblies/*.telo.fasta)
        if [[ ${#fasta_files[@]} -eq 0 ]]; then
          echo "[error] Still no *.telo.fasta files after generation. Check that Step 7 produced assemblies and that motif '$motif' is correct."
          exit 1
        fi
      fi
      for file1 in "${fasta_files[@]}"; do
        for file2 in "${fasta_files[@]}"; do
          base1=$(basename "$file1" .telo.fasta)
log_version "merge_wrapper.py" "merge_wrapper.py"
          base2=$(basename "$file2" .telo.fasta)
          merge_wrapper.py -l 1000000 "$file1" "$file2" --prefix merged_"$base1"_"$base2"
          check_command
        done
      done
      cat merged_*.fasta > allmerged_telo.fasta
      eval "$(conda shell.bash hook)"
      conda deactivate
      conda activate funannotate
      ( eval "$(conda shell.bash hook)"; conda deactivate 2>/dev/null || true; conda activate funannotate 2>/dev/null || true; if ! command -v funannotate >/dev/null 2>&1; then echo "[error] funannotate not found in env \"funannotate\"" >&2; exit 127; fi; funannotate sort -i allmerged_telo.fasta -b contig -o allmerged_telo_sort.fasta --minlen 500 )
      eval "$(conda shell.bash hook)"; conda deactivate 2>/dev/null || true; conda activate pacbiohifi || true
      seqtk telo -s 1 -m "$motif" allmerged_telo_sort.fasta > allmerged.telo.list
      bash ./t2t_list.sh -i allmerged.telo.list -o t2t.list
      ~/opt/scripts/faSomeRecords allmerged_telo_sort.fasta t2t.list t2t.fasta
      eval "$(conda shell.bash hook)"; conda deactivate 2>/dev/null || true; conda activate funannotate || true
      ( eval "$(conda shell.bash hook)"; conda deactivate 2>/dev/null || true; conda activate funannotate 2>/dev/null || true; if ! command -v funannotate >/dev/null 2>&1; then echo "[error] funannotate not found in env \"funannotate\"" >&2; exit 127; fi; funannotate clean -i  t2t.fasta -p 30 -o  t2t_clean.fasta --exhaustive )
      check_command
      ;;


11)
echo "Step 11 - QUAST metricscs for all assemblies"
eval "$(conda shell.bash hook)"; conda deactivate 2>/dev/null || true; conda activate pacbiohifi || true

shopt -s nullglob
a=(assemblies/*.result.fasta)
if [[ ${#a[@]} -eq 0 ]]; then
  echo "[error] No assemblies to run QUAST." >&2
  exit 1
fi

mkdir -p quast_out assemblies

# Try to run QUAST if available (but okay if it's not on PATH)
QUAST_BIN="$(command -v quast.py || command -v quast || true)"
if [[ -n "$QUAST_BIN" ]]; then
  "$QUAST_BIN" "${a[@]}" --threads "${threads:-8}" -o quast_out 2>&1 | tee quast.log || true
else
  echo "[warn] quast.py/quast not found; will use existing quast_out/report files if present." >&2
fi

# Build assemblies/assembly.quast.csv with metrics as rows
python3 - <<'PY'
import os, csv, sys

OUT = os.path.join("assemblies", "assembly.quast.csv")
TREPORT = os.path.join("quast_out", "transposed_report.tsv")
REPORT  = os.path.join("quast_out", "report.tsv")

DESIRED = ["canu","external","flye","ipa","nextDenovo","peregrine","RAFT-hifiasm"]

def norm(s):
    return (s or "").lower().replace("-", "").replace("_", "").replace(".result", "").replace(".fasta","").strip()

def write_csv(rows, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        csv.writer(f).writerows(rows)

def build_from_report(pth):
    """
    QUAST report.tsv (usual): metrics are ROWS, assemblies are COLUMNS.
    Header: ["Assembly", <asm1>, <asm2>, ...]
    """
    with open(pth, newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    if not rows: return None

    header = rows[0]            # ["Assembly", asm1, asm2, ...]
    asm_headers = header[1:]    # assembly labels in columns
    n_asm = len(asm_headers)

    # Map desired names -> column index in asm_headers
    mapcol = {}
    used = set()
    for d in DESIRED:
        nd = norm(d)
        idx = None
        for j, h in enumerate(asm_headers):
            if j in used: continue
            if nd in (norm(h),) or nd == norm(h) or norm(h).startswith(nd) or nd in norm(h):
                idx = j; break
        if idx is None and nd == "rafthifiasm":
            for j, h in enumerate(asm_headers):
                if j in used: continue
                nh = norm(h)
                if "raft" in nh and "hifiasm" in nh:
                    idx = j; break
        mapcol[d] = idx
        if idx is not None: used.add(idx)

    out = [["Metric"] + DESIRED]
    # Each subsequent row: [metric_name, val_asm1, val_asm2, ...]
    for r in rows[1:]:
        if not r: continue
        metric = r[0]
        # pad row to full length
        r = r + [""] * (1 + n_asm - len(r))
        vals = []
        for d in DESIRED:
            j = mapcol[d]
            vals.append(r[1 + j] if j is not None else "")
        out.append([metric] + vals)
    return out

def build_from_transposed(pth):
    """
    transposed_report.tsv: assemblies are ROWS, metrics are COLUMNS.
    Header: ["Assembly", metric1, metric2, ...]
    """
    with open(pth, newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    if not rows: return None

    header = rows[0]
    metrics = header[1:]  # metric names
    # Map desired names -> row index where first cell is assembly label
    maprow = {d: None for d in DESIRED}
    for i in range(1, len(rows)):
        asm_name = rows[i][0] if rows[i] else ""
        na = norm(asm_name)
        for d in DESIRED:
            nd = norm(d)
            if maprow[d] is not None: continue
            if nd == na or nd in na or na in nd:
                maprow[d] = i
            elif nd == "rafthifiasm" and ("raft" in na and "hifiasm" in na):
                maprow[d] = i

    # Build output with metrics as rows
    out = [["Metric"] + DESIRED]
    for j, m in enumerate(metrics, start=1):
        line = [m]
        for d in DESIRED:
            i = maprow[d]
            if i is None:
                line.append("")
            else:
                row = rows[i]
                # pad row
                row = row + [""] * (j + 1 - len(row))
                line.append(row[j] if j < len(row) else "")
        out.append(line)
    return out

rows = None
if os.path.exists(REPORT):
    rows = build_from_report(REPORT)
elif os.path.exists(TREPORT):
    rows = build_from_transposed(TREPORT)
else:
    sys.stderr.write("[error] Could not find quast_out/report.tsv or quast_out/transposed_report.tsv\n")
    sys.exit(1)

write_csv(rows, OUT)
print(f"Wrote {OUT}")
PY

echo "Done: assemblies/assembly.quast.csv"
  ;;
      12)
        # Step 12 - Final merge (choose assembler via --choose or prompt)
        echo "Step 12 - Final merge"

        # Parse optional --choose argument
        # assembler default changed to empty to force prompt when --choose not given
        assembler=""
        CHOOSE_FLAG=0
        while [[ $# -gt 0 ]]; do
          case "$1" in
            --choose)
              CHOOSE_FLAG=1
              if [ -n "${2:-}" ] && [ "${2#-}" != "$2" ] ; then
                :
              else
                if [ -n "${2:-}" ]; then assembler="$2"; shift; fi
              fi
              shift
              ;;
            --choose=*)
              CHOOSE_FLAG=1
              assembler="${1#--choose=}"
              shift
              ;;
            *)
              shift
              ;;
          esac
        done

        # If not provided, build/print assemblies/assembly_info.csv (if parts exist) and then prompt
        if [[ -z "$assembler" ]]; then
          if [[ ! -d assemblies ]]; then
            echo "[error] 'assemblies' folder not found (expected to exist)." >&2
            exit 1
          fi

          parts=()
          for cand in assemblies/assembly.telo.csv assemblies/assembly.busco.csv assemblies/assembly.quast.csv assemblies/assemblies.quast.csv; do
            [[ -s "$cand" ]] && parts+=("$cand")
          done

          if (( ${#parts[@]} > 0 )); then
            info_csv="assemblies/assembly_info.csv"
            desired_header="Metric,canu,external,flye,ipa,nextDenovo,peregrine,RAFT-hifiasm"
            echo "$desired_header" > "$info_csv"
            for f in "${parts[@]}"; do
              awk -v TARGET="$desired_header" -F',' -v OFS=',' '
                BEGIN{
                  n=split(TARGET, want, ",")
                  for(i=1;i<=n;i++){ gsub(/^ *| *$/,"", want[i]); lwant[i]=tolower(want[i]) }
                }
                NR==1{
                  for(i=1;i<=NF;i++){
                    gsub(/
/,"",$i); h=$i; gsub(/^ *| *$/,"",h)
                    l=tolower(h); hmap[l]=i
                  }
                  for(i=1;i<=n;i++){ idx[i]= (lwant[i] in hmap ? hmap[lwant[i]] : 0) }
                  next
                }
                NR>1{
                  row=""
                  for(i=1;i<=n;i++){
                    v=(idx[i]>0 && idx[i]<=NF)? $idx[i] : ""
                    gsub(/
/,"",v)
                    row = (i==1)? v : (row OFS v)
                  }
                  print row
                }' "$f" >> "$info_csv"
            done
          fi

          if [[ -s assemblies/assembly_info.csv ]]; then
            echo "==== assemblies/assembly_info.csv ===="
            if command -v column >/dev/null 2>&1; then
              column -s, -t assemblies/assembly_info.csv | sed 's/^/  /'
            else
              sed 's/^/  /' assemblies/assembly_info.csv
            fi
            echo "======================================"
          else
            echo "[warn] No summary CSVs found to build assemblies/assembly_info.csv; available assemblies:"
            shopt -s nullglob
            for f in assemblies/*.result.fasta; do
              b=$(basename "$f"); echo "  - ${b%.result.fasta}"
            done
          fi
        fi

        # --- interactive selection for assembler (Step 12) ---
        valid_assemblers="canu external flye ipa nextDenovo peregrine RAFT-hifiasm"
        if [[ -z "${assembler:-}" ]]; then
          if ! prompt_from_tty assembler "Enter the assembler to use for the final merge (e.g., ${valid_assemblers// /, }): "; then
            echo "[error] No interactive TTY available; re-run with --choose=<assembler>." >&2
            exit 2
          fi
        fi
        while :; do
          assembler="$(printf "%s" "$assembler" | awk '{$1=$1};1')"
          if printf ' %s ' "$valid_assemblers" | grep -q " $assembler "; then
            break
          fi
          echo "[warn] '$assembler' is not a valid assembler." >&2
          if ! prompt_from_tty assembler "Enter one of [${valid_assemblers// /, }]: "; then
            echo "[error] No interactive TTY available; re-run with --choose=<assembler>." >&2
            exit 2
          fi
        done
        # --- end interactive selection ---

        # Validate choice by checking the FASTA file exists
        asm_fa="assemblies/${assembler}.result.fasta"
        if [[ ! -s "$asm_fa" ]]; then
          echo "[error] Selected assembler '$assembler' does not have '$asm_fa'." >&2
          echo "       Tip: choices are based on files named assemblies/<name>.result.fasta" >&2
          exit 1
        fi

        # Merge chosen assembly with t2t_clean.fasta
        log_version "merge_wrapper.py" "merge_wrapper.py" 2>/dev/null || true
        merge_wrapper.py -l 1000000 "$asm_fa" "t2t_clean.fasta" --prefix "${assembler}"
        check_command

        # Locate merged output produced by merge_wrapper.py
        merged_fa=""
        for cand in "merged_${assembler}.fasta" "merged_${assembler}.fa" "${assembler}.fasta" "${assembler}.fa"; do
          if [[ -s "$cand" ]]; then merged_fa="$cand"; break; fi
        done
        if [[ -z "$merged_fa" ]]; then
          echo "[error] Could not find merged FASTA for prefix '${assembler}' after merge_wrapper.py" >&2
          exit 1
        fi
        echo "[ok] Using merged FASTA: $merged_fa"

        # === Protect t2t contigs during redundans ===
        awk '/^>/{gsub(/^>/,""); print $1}' t2t_clean.fasta > assemblies/.protected.ids

        python3 - <<'PY'
import sys
from pathlib import Path

merged = Path(sys.argv[1])
t2t    = Path(sys.argv[2])
prot_ids = {ln.strip().split()[0] for ln in open(t2t) if ln.startswith(">")}

def read_fa(p):
    name=None; seq=[]
    with open(p) as f:
        for ln in f:
            if ln.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name=ln[1:].strip().split()[0]
                seq=[]
            else:
                seq.append(ln.strip())
    if name is not None:
        yield name, "".join(seq)

prot_out = Path("assemblies/protected.t2t.fa")
oth_out  = Path("assemblies/others.fa")

with open(prot_out,"w") as po, open(oth_out,"w") as oo:
    for name,seq in read_fa(merged):
        tgt = po if name in prot_ids else oo
        tgt.write(f">{name}
")
        for i in range(0,len(seq),60):
            tgt.write(seq[i:i+60]+"
")
print("[split] wrote", prot_out, "and", oth_out, file=sys.stderr)
PY
        "$merged_fa" t2t_clean.fasta
        check_command

        # Run redundans ONLY on 'others'
        eval "$(conda shell.bash hook)"
        conda deactivate
        conda activate redundans
        log_version "redundans.py" "redundans.py" 2>/dev/null || true
        redundans.py --noscaffolding --nogapclosing -t "${threads:-8}" -f assemblies/others.fa --identity 0.50 --overlap 0.80 --log redundans.log
        check_command

        reduced_others="redundans/scaffolds.reduced.fa"
        if [[ ! -s "$reduced_others" ]]; then
          echo "[error] redundans reduced FASTA not found at $reduced_others" >&2
          exit 1
        fi

        # Optional filter: drop others that are redundant to t2t (id>=0.50, cov>=0.80)
        eval "$(conda shell.bash hook)"
        conda deactivate
        conda activate pacbiohifi 2>/dev/null || true
        if command -v minimap2 >/dev/null 2>&1; then
          log_version "minimap2" "minimap2" 2>/dev/null || true
          paf="assemblies/others_vs_t2t.paf"
          minimap2 -x asm20 -t "${threads:-8}" t2t_clean.fasta "$reduced_others" > "$paf"
          python3 - <<'PY'
import sys
best = {}
with open(sys.argv[1]) as f:
    for ln in f:
        if ln.startswith("#") or not ln.strip(): continue
        t = ln.rstrip("
").split("	")
        if len(t) < 12: continue
        qname, qlen, qstart, qend = t[0], int(t[1]), int(t[2]), int(t[3])
        tags = {kv.split(":")[0]: kv for kv in t[12:]}
        alnlen = abs(qend - qstart)
        if alnlen == 0: continue
        ident = None
        if 'NM' in tags:
            nm = int(tags['NM'].split(":")[-1])
            ident = max(0.0, 1.0 - (nm / max(1, alnlen)))
        elif 'de' in tags:
            try:
                de = float(tags['de'].split(":")[-1])
                ident = max(0.0, 1.0 - de)
            except: pass
        if ident is None: 
            continue
        cov = alnlen / float(qlen)
        cur = best.get(qname, (0.0, 0.0))
        if cov > cur[0] or (abs(cov-cur[0])<1e-6 and ident > cur[1]):
            best[qname] = (cov, ident)
drop = {q for q,(cov,iden) in best.items() if cov >= 0.80 and iden >= 0.50}
keep = []
with open(sys.argv[2]) as f:
    name=None; seq=[]
    for ln in f:
        if ln.startswith(">"):
            if name is not None and name not in drop:
                keep.append(name)
            name=ln[1:].strip().split()[0]
    if name is not None and name not in drop:
        keep.append(name)
print("\n".join(keep))
PY
          "$paf" "$reduced_others" > assemblies/others.keep.ids

          python3 - <<'PY'
import sys
ids = set(ln.strip() for ln in open(sys.argv[1]) if ln.strip())
inp = sys.argv[2]
outp = sys.argv[3]
with open(inp) as f, open(outp,"w") as o:
    name=None; seq=[]
    def flush():
        if name is not None and name in ids:
            o.write(f">{name}
")
            for i in range(0,len(seq),60): o.write(seq[i:i+60]+"
")
    for ln in f:
        if ln.startswith(">"):
            flush()
            name=ln[1:].strip().split()[0]
            seq=[]
        else:
            seq.append(ln.strip())
    flush()
PY
          assemblies/others.keep.ids "$reduced_others" assemblies/others.filtered.fa
        else
          echo "[warn] minimap2 not found; keeping all reduced 'others' without extra filtering." >&2
          cp -f "$reduced_others" assemblies/others.filtered.fa
        fi

        # Recombine: protected t2t + filtered others
        cat assemblies/protected.t2t.fa assemblies/others.filtered.fa > assemblies/merged_protected_priority.fa
        echo "[ok] Built assemblies/merged_protected_priority.fa (t2t protected)"

        # Sort final contigs
        eval "$(conda shell.bash hook)"
        conda deactivate
        conda activate funannotate
        check_command
        log_version "funannotate" "funannotate" 2>/dev/null || true
        ( eval "$(conda shell.bash hook)"; conda deactivate 2>/dev/null || true; conda activate funannotate 2>/dev/null || true;           if ! command -v funannotate >/dev/null 2>&1; then echo "[error] funannotate not found in env "funannotate"" >&2; exit 127; fi;           funannotate sort -i assemblies/merged_protected_priority.fa -b contig -o "merged_${assembler}_sort.fa" --minlen 500 )
        check_command

        # Publish final
        if [[ ! -d assemblies ]]; then
          echo "[error] 'assemblies' folder not found (expected to exist)." >&2
          exit 1
        fi
        cp -f "merged_${assembler}_sort.fa" "assemblies/final.merged.fasta"
        echo "[ok] Wrote assemblies/final.merged.fasta"
        ;;

13)
# Step 13 - BUSCO analysis for final merged assembly (reuse existing results; JSON/TXT/TSV parsing)
echo "Step 13 - BUSCO analysis"
eval "$(conda shell.bash hook)"; conda deactivate 2>/dev/null || true; conda activate busco || true

threads=${threads:-8}
lineage=${busco_lineage:-"ascomycota_odb10"}   # set via --busco <lineage>

# Preconditions
final_fa="assemblies/final.merged.fasta"
if [[ ! -s "$final_fa" ]]; then
  echo "[error] Final merged FASTA '$final_fa' not found." >&2
  exit 1
fi




[[ -d assemblies ]] || { echo "[error] 'assemblies' folder not found."; exit 1; }
[[ -d busco ]] || mkdir busco

# Helper: does the 'final' run already have metrics?
has_final_metrics() {
  for d in "busco/final" "busco/run_final"; do
    [[ -d "$d" ]] || continue
    if find "$d" -type f \( -name 'short_summary*.json' -o -name 'short_summary*.txt' -o -name 'full_table*.tsv' \) -print -quit | grep -q .; then
      return 0
    fi
  done
  return 1
}

# Run BUSCO only if needed
if has_final_metrics; then
  echo "[ok] Found existing BUSCO metrics for final → using busco/final (or legacy path)"
else
  BUSCO_BIN="$(command -v busco || true)"
  if [[ -z "$BUSCO_BIN" ]]; then
    echo "[error] Missing final BUSCO metrics and 'busco' is not available to run." >&2
    exit 127
  fi
  force_flag=""
  if [[ -d busco/final || -d busco/run_final ]]; then
    echo "[info] Incomplete 'final' BUSCO dir detected → re-running with -f"
    force_flag="-f"
  fi
  echo "[run] BUSCO on final (lineage=$lineage, threads=$threads)"
  pushd busco >/dev/null
  busco -o final -i "../$final_fa" -l "$lineage" -m genome -c "$threads" $force_flag 2>&1 | tee busco_final.log || true
  popd >/dev/null
fi

check_command

# Build assemblies/final.busco.csv (metrics-as-rows; single column 'final') using a temp Python script
pyfile="$(mktemp -t final_busco.XXXXXX.py)"
cat > "$pyfile" <<'PY'
import os, re, glob, json, csv

OUT_CSV = os.path.join("assemblies", "final.busco.csv")

def newest(paths):
    return max(paths, key=os.path.getmtime) if paths else None

def find_run_roots():
    roots = []
    for d in ("busco/final", "busco/run_final"):
        if os.path.isdir(d): roots.append(d)
    return roots

def parse_short_summary_json(path):
    try:
        with open(path, "r") as f:
            data = json.load(f)
        n = data.get("n") or data.get("lineage_dataset", {}).get("n")
        per = data.get("percentages") or {}
        if n and per:
            Cpct = float(per.get("C", 0)); Spct = float(per.get("S", 0))
            Dpct = float(per.get("D", 0)); Fpct = float(per.get("F", 0)); Mpct = float(per.get("M", 0))
            Ccnt = round(n * Cpct / 100.0); Mcnt = round(n * Mpct / 100.0)
            return {"Cpct":Cpct, "Spct":Spct, "Dpct":Dpct, "Fpct":Fpct, "Mpct":Mpct, "Ccnt":Ccnt, "Mcnt":Mcnt, "n":int(n)}
        S = data.get("S"); D = data.get("D"); F = data.get("F"); M = data.get("M")
        if None not in (S, D, F, M, n):
            n = int(n); S=int(S); D=int(D); F=int(F); M=int(M)
            C = S + D
            def pct(x): return 0.0 if not n else 100.0*float(x)/float(n)
            return {"Cpct":pct(C), "Spct":pct(S), "Dpct":pct(D), "Fpct":pct(F), "Mpct":pct(M), "Ccnt":C, "Mcnt":M, "n":n}
    except Exception:
        pass
    return None

def parse_short_summary_txt(path):
    try:
        txt = open(path, "r", errors="ignore").read()
        m = re.search(r"C:(\d+(?:\.\d+)?)%.*?S:(\d+(?:\.\d+)?)%.*?D:(\d+(?:\.\d+)?)%.*?F:(\d+(?:\.\d+)?)%.*?M:(\d+(?:\.\d+)?)%.*?n:(\d+)", txt, re.S)
        if not m: return None
        Cpct, Spct, Dpct, Fpct, Mpct = map(float, m.groups()[:5]); n = int(m.group(6))
        Ccnt = round(n * Cpct / 100.0); Mcnt = round(n * Mpct / 100.0)
        return {"Cpct":Cpct, "Spct":Spct, "Dpct":Dpct, "Fpct":Fpct, "Mpct":Mpct, "Ccnt":Ccnt, "Mcnt":Mcnt, "n":n}
    except Exception:
        return None

def parse_full_table_tsv(path):
    lines = [ln.rstrip("\n") for ln in open(path, newline="")]
    status_idx = None
    for ln in lines:
        if ln.startswith("#") and "Status" in ln:
            hdr = ln.lstrip("#").strip().split("\t")
            for i, h in enumerate(hdr):
                if h.strip().lower() == "status": status_idx = i
            break
    S = D = F = M = 0; n = 0
    for ln in lines:
        if not ln or ln.startswith("#"): continue
        parts = ln.split("\t")
        st = parts[1].strip() if status_idx is None else (parts[status_idx].strip() if len(parts) > status_idx else "")
        n += 1
        if   st == "Complete":   S += 1
        elif st == "Duplicated": D += 1
        elif st == "Fragmented": F += 1
        elif st == "Missing":    M += 1
    if n == 0: return None
    C = S + D
    def pct(x): return 0.0 if not n else 100.0*float(x)/float(n)
    return {"Cpct":pct(C), "Spct":pct(S), "Dpct":pct(D), "Fpct":pct(F), "Mpct":pct(M), "Ccnt":C, "Mcnt":M, "n":n}

def read_metrics():
    roots = find_run_roots()
    if not roots: return None
    # Prefer JSON → TXT → TSV
    for r in sorted(roots, key=os.path.getmtime, reverse=True):
        js = sorted(glob.glob(os.path.join(r, "**", "short_summary*.json"), recursive=True), key=os.path.getmtime, reverse=True)
        if js:
            m = parse_short_summary_json(js[0])
            if m: return m
        tx = sorted(glob.glob(os.path.join(r, "**", "short_summary*.txt"), recursive=True), key=os.path.getmtime, reverse=True)
        if tx:
            m = parse_short_summary_txt(tx[0])
            if m: return m
        ft = sorted(glob.glob(os.path.join(r, "**", "full_table*.tsv"), recursive=True), key=os.path.getmtime, reverse=True)
        if ft:
            m = parse_full_table_tsv(ft[0])
            if m: return m
    return None

def fmt_pct(x):
    s = f"{float(x):.1f}"
    return s.rstrip('0').rstrip('.') if '.' in s else s

m = read_metrics()
if not m:
    raise SystemExit("[error] Could not find BUSCO outputs in busco/final or busco/run_final")

rows = [
    ["Metric", "final"],
    ["BUSCO C (%)", fmt_pct(m["Cpct"])],
    ["BUSCO S (%)", fmt_pct(m["Spct"])],
    ["BUSCO D (%)", fmt_pct(m["Dpct"])],
    ["BUSCO F (%)", fmt_pct(m["Fpct"])],
    ["BUSCO M (%)", fmt_pct(m["Mpct"])],
    ["BUSCO C (count)", str(m["Ccnt"])],
    ["BUSCO M (count)", str(m["Mcnt"])],
]

with open(OUT_CSV, "w", newline="") as f:
    csv.writer(f).writerows(rows)

print(f"Wrote {OUT_CSV}")
PY

python3 "$pyfile"
rm -f "$pyfile"

echo "[ok] Wrote assemblies/final.busco.csv (lineage: ${lineage})"
      ;;

    14)
# Step 14 - Telomere analysis (final assembly)
echo "Step 14 - Telomere analysis (final assembly)"

# Use the env that contains seqtk
eval "$(conda shell.bash hook)"; conda deactivate 2>/dev/null || true; conda activate funannotate 2>/dev/null || true

# Preconditions
[[ -d assemblies ]] || { echo "[error] 'assemblies' folder not found."; exit 1; }
final_fa="assemblies/final.merged.fasta"
if [[ ! -s "$final_fa" ]]; then
  echo "[error] Final merged FASTA '$final_fa' not found." >&2
  exit 1
fi
if ! command -v seqtk >/dev/null 2>&1; then
  echo "[error] seqtk not found in env 'funannotate'." >&2
  exit 127
fi

# Outputs (all under assemblies/)
list="assemblies/final.telo.list"
ids="assemblies/final.telo.list.ids"
out="assemblies/final.telo.fasta"
csv="assemblies/final.telo.csv"

# 1) Find telomere-containing contigs and write list
#    Use provided $motif if set; fall back to TTAGGG
seqtk telo -s 1 -m "${motif:-TTAGGG}" "$final_fa" > "$list"
check_command

# 2) Unique IDs -> subseq to .telo.fasta
awk '{print $1}' "$list" | sed 's/[[:space:]].*$//' | tr -d $'\r' | sort -u > "$ids"
if [[ -s "$ids" ]]; then
  seqtk subseq "$final_fa" "$ids" > "$out"
  [[ -s "$out" ]] && echo "[ok] Wrote $(basename "$out")" || echo "[warn] ${out} is empty"
else
  echo "[warn] No telomere-containing contigs for final assembly; writing empty ${out}"
  : > "$out"
fi

# 3) Compute metrics (same logic as for individual assemblies)
double=$(awk 'NF>=4 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $4~/^[0-9]+$/ {c=$1; sub(/[ \t].*$/,"",c); if($2==0 && $3==$4) print c}' "$list" | sort -u | wc -l)
  # v0.2.2 override: correct double-end contig counting (both ends flagged)
  double=$(awk '
    NF>=4 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $4~/^[0-9]+$/ {
      c=$1; sub(/[ \t].*$/,"",c);
      s[c]+=($2==0)?1:0; e[c]+=($3==$4)?1:0
    }
    END{ for(c in s){ if(s[c]>0 && e[c]>0) print c } }
  ' "$list" | sort -u | wc -l)
single=$(awk '
  NF>=4 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $4~/^[0-9]+$/ {
    c=$1; sub(/[ \t].*$/,"",c);
    s[c]+=($2==0)?1:0; e[c]+=($3==$4)?1:0
  }
  END{ for(c in s){ if((s[c]+e[c])==1) print c } }
' "$list" | wc -l)

# 4) Write metrics CSV (metrics as rows; single column 'final')
{
  echo "Metric,final"
  echo "Telomere double-end contigs,${double:-0}"
  echo "Telomere single-end contigs,${single:-0}"
} > "$csv"
echo "[ok] Wrote $(basename "$csv")"
      ;;

15)
# Step 15 - QUAST metrics for FINAL assembly
echo "Step 15 - QUAST metrics for final assembly"
eval "$(conda shell.bash hook)"; conda deactivate 2>/dev/null || true; conda activate pacbiohifi || true

# Preconditions
final_fa="assemblies/final.merged.fasta"
if [[ ! -s "$final_fa" ]]; then
  echo "[error] Final merged FASTA '$final_fa' not found." >&2
  exit 1
fi
[[ -d assemblies ]] || { echo "[error] 'assemblies' folder not found."; exit 1; }

# Run QUAST → quast_final/
mkdir -p quast_final
QUAST_BIN="$(command -v quast.py || command -v quast || true)"
if [[ -n "$QUAST_BIN" ]]; then
  "$QUAST_BIN" "$final_fa" --labels final --threads "${threads:-8}" -o quast_final 2>&1 | tee quast_final.log || true
else
  echo "[warn] quast.py/quast not found; will try to use any existing quast_final/*report*.tsv" >&2
fi


# Build assemblies/final-quast.tsv and assemblies/final.quast.csv
python3 - <<'PY'
import os, csv, sys

out_dir_csv = os.path.join("assemblies")
if not os.path.isdir(out_dir_csv):
    sys.stderr.write("[error] 'assemblies' folder not found.\n")
    sys.exit(1)

report  = os.path.join("quast_final", "report.tsv")
treport = os.path.join("quast_final", "transposed_report.tsv")

final_report_tsv = os.path.join("assemblies", "final-quast.tsv")
final_csv        = os.path.join("assemblies", "final.quast.csv")

def write_tsv(rows, path):
    with open(path, "w", newline="") as f:
        csv.writer(f, delimiter="\t").writerows(rows)

def write_csv(rows, path):
    with open(path, "w", newline="") as f:
        csv.writer(f).writerows(rows)

def from_report_to_csv(report_path):
    with open(report_path, newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    if not rows or len(rows[0]) < 2:
        sys.stderr.write("[error] report.tsv appears malformed.\n"); sys.exit(1)
    out = [["Metric","final"]]
    for r in rows[1:]:
        if not r: continue
        metric = r[0]
        val = r[1] if len(r) > 1 else ""
        out.append([metric, val])
    return out

def from_treport_to_outputs(treport_path):
    with open(treport_path, newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    if not rows or len(rows[0]) < 2:
        sys.stderr.write("[error] transposed_report.tsv appears malformed.\n"); sys.exit(1)
    header  = rows[0]
    metrics = header[1:]

    row_idx = None
    for i in range(1, len(rows)):
        if rows[i] and rows[i][0].strip().lower() == "final":
            row_idx = i; break
    if row_idx is None:
        row_idx = 1 if len(rows) > 1 else None
    if row_idx is None:
        sys.stderr.write("[error] No assembly rows in transposed_report.tsv.\n"); sys.exit(1)

    vals_row = rows[row_idx]
    rep_rows = [["Assembly","final"]]
    csv_rows = [["Metric","final"]]
    for j, m in enumerate(metrics, start=1):
        v = vals_row[j] if j < len(vals_row) else ""
        rep_rows.append([m, v])
        csv_rows.append([m, v])
    return rep_rows, csv_rows

if os.path.exists(report):
    with open(report, newline="") as f:
        rep_rows = list(csv.reader(f, delimiter="\t"))
    if rep_rows:
        rep_rows[0][0] = "Assembly"
    write_tsv(rep_rows, final_report_tsv)
    csv_rows = from_report_to_csv(final_report_tsv)
    write_csv(csv_rows, final_csv)
elif os.path.exists(treport):
    rep_rows, csv_rows = from_treport_to_outputs(treport)
    write_tsv(rep_rows, final_report_tsv)
    write_csv(csv_rows, final_csv)
else:
    sys.stderr.write("[error] Could not find quast_final/report.tsv or transposed_report.tsv\n")
    sys.exit(1)

print(f"Wrote {final_report_tsv}")
print(f"Wrote {final_csv}")
PY

echo "[ok] Wrote assemblies/final-quast.tsv and assemblies/final.quast.csv"
      ;;
16)
# Step 16 - Final assembly comparison
echo "Step 16 - Final assembly comparison"

# Inputs (in assemblies/)
final_telo_csv="assemblies/final.telo.csv"
final_busco_csv="assemblies/final.busco.csv"
final_quast_csv="assemblies/final.quast.csv"
asm_info_csv="assemblies/assembly_info.csv"

out_csv="final_result.csv"   # write to current working directory

python3 - <<'PY'
import os, csv, sys

# Inputs
FINAL_FILES = [
    ("telo",  "assemblies/final.telo.csv"),
    ("busco", "assemblies/final.busco.csv"),
    ("quast", "assemblies/final.quast.csv"),
]
ASM_INFO = "assemblies/assembly_info.csv"
OUT = "final_result.csv"

def read_metric_final(path):
    """Read a (Metric,final) CSV into an ordered list of (metric, value)."""
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path, newline="") as f:
        r = csv.reader(f)
        hdr = next(r, None)
        if not hdr or len(hdr) < 2:
            return rows
        # find columns by name (robust to slight variance)
        idx_metric = None
        idx_final = None
        for i,h in enumerate(hdr):
            name = (h or "").strip().lower()
            if name == "metric":
                idx_metric = i
            if name == "final":
                idx_final = i
        # Fall back to positional if needed
        if idx_metric is None: idx_metric = 0
        if idx_final  is None: idx_final  = 1 if len(hdr) > 1 else 0
        for row in r:
            if not row: continue
            m = (row[idx_metric] if idx_metric < len(row) else "").strip()
            v = (row[idx_final]  if idx_final  < len(row) else "").strip()
            if m:
                rows.append((m, v))
    return rows

def read_assembly_info(path):
    """Read assemblies/assembly_info.csv (Metric + assemblers columns)."""
    if not os.path.exists(path):
        return None, None
    with open(path, newline="") as f:
        r = csv.reader(f)
        hdr = next(r, None)
        if not hdr: return None, None
        # normalize first header cell to 'Metric'
        hdr = list(hdr)
        hdr[0] = "Metric"
        body = [row for row in r if row]
    return hdr, body

# 1) Collect all final metrics into an ordered dict (preserve file order: telo->busco->quast)
from collections import OrderedDict
final_map = OrderedDict()
order_seen = []  # preserve metric order as we read
for tag, p in FINAL_FILES:
    for m, v in read_metric_final(p):
        if m not in final_map:
            order_seen.append(m)
        final_map[m] = v

# 2) If assembly_info exists, append a 'final' column matched by 'Metric'
hdr, body = read_assembly_info(ASM_INFO)

if hdr is not None:
    # Ensure exact target order and add 'final' at end
    if "final" not in [h.lower() for h in hdr]:
        hdr_out = list(hdr) + ["final"]
    else:
        # unlikely, but avoid duplicate
        hdr_out = list(hdr)

    # Map metric -> row index for quick lookup
    metric_col = 0
    out_rows = []
    for row in body:
        # pad row to header length if needed
        if len(row) < len(hdr_out)-1:
            row = row + [""] * (len(hdr_out)-1 - len(row))
        metric = (row[metric_col] if row else "").strip()
        final_val = final_map.get(metric, "")
        out_rows.append(row + [final_val])

    # 3) Append any final-only metrics not present in assembly_info
    present_metrics = set(r[0].strip() for r in body if r and r[0].strip())
    for m in order_seen:
        if m not in present_metrics:
            r = [""] * (len(hdr_out))
            r[0] = m
            r[-1] = final_map.get(m, "")
            out_rows.append(r)

    # Write OUT
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(hdr_out)
        w.writerows(out_rows)

else:
    # No assembly_info.csv → write a simple Metric,final sheet
    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Metric","final"])
        for m in order_seen:
            w.writerow([m, final_map.get(m, "")])

print(f"Wrote {OUT}")
PY

check_command
echo "[ok] Wrote final_result.csv"
      ;;
      
17)
# Step 17 - Cleanup temporary files
echo "Step 17 - Cleanup temporary files"
cleanup_temp() {
  echo "Cleanup: organizing temporary outputs…"
  shopt -s nullglob dotglob extglob nocaseglob
  mkdir -p temp/merge temp/merge/fasta temp/merge/param temp/busco temp/log final_results
  find . -maxdepth 1 -type f -name 'aln_summary_merged*.tsv'   -exec mv -f -- {} temp/merge/ \; 2>/dev/null || true
  find . -maxdepth 1 -type f -name 'anchor_summary_merged_*.txt' -exec mv -f -- {} temp/merge/ \; 2>/dev/null || true
  find . -maxdepth 1 -type f \( -name '.merged_*' -o -name 'merged_*.delta' -o -name 'merged_*.coords' -o -name 'merged_*.snps' -o -name 'merged_*.delta.*' -o -name 'merged_*.crunch' -o -name 'merged_*.filter' -o -name 'merged_*.qdiff' -o -name 'merged_*.rdiff' -o -name 'merged_*.mcoords' \) \
       -exec mv -f -- {} temp/merge/ \; 2>/dev/null || true
  find busco -maxdepth 2 -type f -name '*.log' -exec mv -f -- {} temp/busco/ \; 2>/dev/null || true
  find .     -maxdepth 1 -type f \( -name 'busco_*.log' -o -name '*busco*.log' \) -exec mv -f -- {} temp/busco/ \; 2>/dev/null || true
  find . -maxdepth 1 -type f \( -name 'merged_*.fasta' -o -name 'merged_*.fa' \) -exec mv -f -- {} temp/merge/fasta/ \; 2>/dev/null || true
  for f in merged_*; do
    [[ -e "$f" && -f "$f" ]] || continue
    case "$f" in merged_*.fasta|merged_*.fa) continue ;; esac
    mv -f -- "$f" temp/merge/
  done
  find . -maxdepth 1 -type f -name 'param_summary_merged_*.txt' -exec mv -f -- {} temp/merge/param/ \; 2>/dev/null || true
  for f in *.log; do
    [[ -f "$f" ]] || continue
    [[ "$f" == *busco* ]] && continue
    mv -f -- "$f" temp/log/
  done
  mv -f -- final_result.csv final_results/
  mv -f -- assemblies/final.merged.fasta final_results/
  echo "[ok] Cleanup complete."
}
cleanup_temp
      ;;
  esac
  echo "===== [$(timestamp)] STEP ${step} END ====="
  } 2>&1 | awk "$LOG_AWK" | tee -a "$STEP_LOG"

done
