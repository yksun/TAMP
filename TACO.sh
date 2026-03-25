#!/bin/bash
if [ -z "$BASH_VERSION" ]; then
  exec /usr/bin/env bash "$0" "$@"
fi

# v0.5.3
# - MAJOR CHANGE (Step 10): Reworked telomere-pool construction from seqtk telo
#   coordinate output. Contigs are now separated into biologically distinct classes:
#   strict T2T contigs, best single-end telomeric contigs, and optimized telomere-supported contigs.
#
# - MAJOR CHANGE (Step 10): Preserved the strict meaning of t2t.fasta.
#   It contains only true double-end telomere-to-telomere contigs.
#   Best single-end telomeric representatives are written to single_tel_best.fasta,
#   and the optimized combined pool is written to telomere_supported_best.fasta.
#
# - MAJOR CHANGE (Step 10): Added redundancy reduction for single-end telomeric contigs
#   using all-vs-all minimap2 clustering and longest-representative selection.
#
# - MAJOR CHANGE (Step 10): Updated protected telomere-pool priority to:
#   strict T2T > best single-end telomeric representatives > optimized telomere-supported pool.
#
# - MAJOR CHANGE (Step 12): Replaced older final merge behavior with cleaned backbone
#   refinement logic. The selected assembler output is used as the backbone assembly,
#   while redundant backbone contigs can be replaced or rescued using the optimized telomere pool.
#
# - MAJOR CHANGE (Step 12): Added support for optional automatic backbone selection modes:
#   smart scoring or legacy N50-only selection.
#
# - CHANGE (Step 12): Increased telomere contribution in smart backbone scoring to better
#   reflect telomere-end retention during final assembly refinement.
#
# - CHANGE (Step 12): Added optional Merqury integration for assembler ranking and final
#   reporting when --merqury or --merqury-db is supplied.
#
# - FIX (Step 12): Improved fungal single-end telomere rescue by relaxing terminal overhang,
#   alignment identity, and coverage thresholds.
#
# - CHANGE (Step 12): Telomere rescue now prioritizes single_tel_best_clean.fasta as the
#   preferred rescue pool before falling back to broader telomeric sets.
#
# - CHANGE (Step 14): Updated final telomere analysis to report strict T2T, single-end,
#   total telomere-supported contigs, protected telomere mode, and rescue counts.
#
# - CHANGE (Step 16): Updated final comparison reporting to include Merqury metrics,
#   telomere-pool statistics, rescue counts, selection score, selected assembler,
#   auto-selection mode, and score formula in final_results/final_result.csv.
#
# - CHANGE (Step 18): Added assembly-only mode to run assembler benchmarking/comparison
#   without final telomere-aware refinement, while still generating assemblies/assembly_info.csv
#   and optional Merqury comparison output.
#
# - FIX: Removed duplicated old parser/selection code, aligned usage text with current logic,
#   and cleaned older command fragments left from pre-0.5.x versions.
#
# - DESIGN UPDATE: TACO now follows a telomere-pool optimization and backbone-refinement
#   strategy rather than repeated structural merging, improving chromosome-end retention
#   while reducing non-telomeric redundancy.

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

# ==== Interactive prompt helpers ====
if [ -z "${__TTY_FD3_OPENED__:-}" ] && [ -r /dev/tty ]; then
  exec 3</dev/tty
  __TTY_FD3_OPENED__=1
fi

prompt_from_tty() {
  # usage: prompt_from_tty VAR "Your prompt: "
  local __outvar=$1; shift
  local __prompt="$*"
  if [ -r /dev/tty ]; then
    printf "%s" "$__prompt" > /dev/tty
    local __ans=""
    if [ -e /proc/self/fd/3 ]; then
      IFS= read -r __ans <&3 || __ans=""
    else
      IFS= read -r __ans < /dev/tty || __ans=""
    fi
    __ans="$(printf "%s" "$__ans" | awk '{$1=$1};1')"
    printf -v "$__outvar" '%s' "$__ans"
    return 0
  fi
  return 1
}

# Default values
genomesize=""
threads=20
fastq=""
steps=()
motif="AACCCT"
assembler=""
external_fasta=""
run_busco=false
busco_lineage="ascomycota_odb10"
assembly_only=false
merqury_enable=false
merqury_db=""
auto_mode="smart"
CHOOSE_FLAG=0

# --- Logging & versioning ---
PIPELINE_NAME="TACO-0.5.3.sh"
RUN_ID="$(date +'%Y%m%d_%H%M%S')"
LOG_DIR="logs"
mkdir -p "$LOG_DIR"

LOG_AWK='{ print strftime("[%Y-%m-%d %H:%M:%S]"), $0; fflush(); }'
timestamp() { date +"%Y-%m-%d %H:%M:%S"; }

BENCH_DIR="benchmark_logs"
BENCH_STEP_TSV="${BENCH_DIR}/step_benchmark.tsv"
BENCH_RUN_TSV="${BENCH_DIR}/run_metadata.tsv"
BENCH_SUMMARY_TXT="${BENCH_DIR}/run_summary.txt"
mkdir -p "$BENCH_DIR"

write_run_metadata() {
  mkdir -p "$BENCH_DIR"
  {
    echo -e "field	value"
    echo -e "pipeline	${PIPELINE_NAME}"
    echo -e "run_id	${RUN_ID}"
    echo -e "date	$(timestamp)"
    echo -e "host	$(hostname 2>/dev/null || echo unknown)"
    echo -e "working_directory	$(pwd)"
    echo -e "conda_env	${CONDA_DEFAULT_ENV:-NA}"
    echo -e "python	$(command -v python3 2>/dev/null || echo NA)"
    echo -e "threads	${threads:-NA}"
    echo -e "genome_size	${genomesize:-NA}"
    echo -e "fastq	${fastq:-NA}"
    echo -e "motif	${motif:-NA}"
    echo -e "steps	${steps[*]:-all}"
    echo -e "os	$(uname -srvmo 2>/dev/null || uname -a 2>/dev/null || echo NA)"
    echo -e "cpu_count	$(getconf _NPROCESSORS_ONLN 2>/dev/null || nproc 2>/dev/null || echo NA)"
    echo -e "memory_kb	$(awk '/MemTotal:/ {print $2; exit}' /proc/meminfo 2>/dev/null || echo NA)"
  } > "$BENCH_RUN_TSV"
}

init_benchmark_step_table() {
  mkdir -p "$BENCH_DIR"
  if [[ ! -f "$BENCH_STEP_TSV" ]]; then
    echo -e "run_id	step	step_name	start_time	end_time	runtime_sec	status	log_file" > "$BENCH_STEP_TSV"
  fi
}

step_name() {
  case "$1" in
    1) echo "HiCanu assembly" ;;
    2) echo "NextDenovo assembly" ;;
    3) echo "Peregrine assembly" ;;
    4) echo "IPA assembly" ;;
    5) echo "Flye assembly" ;;
    6) echo "Hifiasm assembly" ;;
    7) echo "Copy all assemblies" ;;
    8) echo "BUSCO for assemblies" ;;
    9) echo "Telomere contigs and metrics" ;;
    10) echo "Build optimized telomere pool" ;;
    11) echo "QUAST for assemblies" ;;
    12) echo "Final assembly refinement with optimized telomere-end replacement" ;;
    13) echo "BUSCO final" ;;
    14) echo "Telomere analysis final" ;;
    15) echo "QUAST final" ;;
    16) echo "Generate final comparison report" ;;
    17) echo "Cleanup temporary files" ;;
    18) echo "Assembly-only comparison summary" ;;
    *) echo "Unknown step" ;;
  esac
}

append_step_benchmark() {
  local step="$1"
  local start_epoch="$2"
  local end_epoch="$3"
  local status="$4"
  local log_file="$5"
  local sname
  sname="$(step_name "$step")"
  printf '%s	%s	%s	%s	%s	%s	%s	%s
'     "$RUN_ID"     "$step"     "$sname"     "$(date -d "@$start_epoch" '+%Y-%m-%d %H:%M:%S' 2>/dev/null || date -r "$start_epoch" '+%Y-%m-%d %H:%M:%S' 2>/dev/null || echo "$start_epoch")"     "$(date -d "@$end_epoch" '+%Y-%m-%d %H:%M:%S' 2>/dev/null || date -r "$end_epoch" '+%Y-%m-%d %H:%M:%S' 2>/dev/null || echo "$end_epoch")"     "$((end_epoch-start_epoch))"     "$status"     "$log_file" >> "$BENCH_STEP_TSV"
}

write_benchmark_summary() {
  mkdir -p "$BENCH_DIR"
  {
    echo "Pipeline: ${PIPELINE_NAME}"
    echo "Run ID: ${RUN_ID}"
    echo "Date: $(timestamp)"
    echo "Benchmark table: ${BENCH_STEP_TSV}"
    echo "Run metadata: ${BENCH_RUN_TSV}"
    echo
    if [[ -f "$BENCH_STEP_TSV" ]]; then
      awk -F'	' 'NR==1{next} {c[$7]++} END{for(k in c) print k ": " c[k] " step(s)"}' "$BENCH_STEP_TSV"
    fi
  } > "$BENCH_SUMMARY_TXT"
}

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
    local out=""
    for f in "--version" "-V" "-v" "version"; do
      out=$("$c" $f 2>&1 | head -n 1)
      if [[ -n "$out" ]]; then
        echo "$c: $out"
        return
      fi
    done
    echo "$c: FOUND (version unknown)"
  }

  for c in canu nextDenovo peregrine ipa flye hifiasm seqtk busco quast.py quast minimap2 merge_wrapper.py bwa samtools python3; do
    getv "$c" >> "$vf"
  done
  echo "" >> "$vf"
}

usage() {
  cat <<USAGE
Usage: ${PIPELINE_NAME:-$0} -g <genomesize> -t <threads> --fastq <fastq> -m <motif> [-s <steps>] [--fasta <path>] [--busco <lineage>] [--choose <assembler>] [--auto-mode <smart|n50>] [--assembly-only] [--merqury|--merqury-db <reads.meryl>]

Options:
  -g, --genomesize   Genome size (required), e.g. 2g, 500m
  -t, --threads      Number of threads (required)
  --fastq            Path to the FASTQ file (required)
  -m, --motif        Telomere motif (required), e.g. AACCCT
  -s, --steps        Steps to run (optional, default: all). Accepts comma/range list (e.g. 1,2,5-7)
  --fasta            External pre-assembled FASTA (optional, single file or URL)
  --busco            Run BUSCO on each individual assembly; optional lineage (default: ascomycota_odb10)
  --choose           Prompt for or provide the assembler to use as the final refinement backbone
  --auto-mode        Auto-selection mode for backbone choice: smart or n50
  --assembly-only    Run only assembler benchmarking/comparison (steps 1-9,11,18)
  --merqury          Enable optional Merqury scoring/reporting using auto-detected .meryl database
  --merqury-db       Enable Merqury and explicitly provide a read-derived .meryl directory
  --no-merqury       Disable Merqury even if present

Notes:
  • Per-step logs: logs/step_<N>.log with timestamps.
  • A consolidated software version summary is written to ./version.txt at start.
  • Example full run:
      ${PIPELINE_NAME:-$0} -g 2g -t 16 --fastq reads.fastq -m AACCCT --auto-mode smart
  • Example assembly-only run:
      ${PIPELINE_NAME:-$0} -g 2g -t 16 --fastq reads.fastq -m AACCCT --assembly-only
  • Example assembly-only run with Merqury:
      ${PIPELINE_NAME:-$0} -g 2g -t 16 --fastq reads.fastq -m AACCCT --assembly-only --merqury-db reads.meryl

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
  10. Build optimized telomere pool
  11. QUAST for all assembler results
  12. Final assembly refinement using selected assembler
  13. BUSCO analysis (final)
  14. Telomere analysis (final)
  15. QUAST analysis (final)
  16. Generate final comparison report
  17. Cleanup temporary files into structured folders
  18. Assembly-only comparison summary
USAGE
  exit 1
}

build_assembly_info_v2() {
  mkdir -p assemblies
  local info_csv="assemblies/assembly_info.csv"
  local desired_header="Metric,canu,external,flye,ipa,nextDenovo,peregrine,hifiasm"

  printf '%s\n' "$desired_header" > "$info_csv"

  cat > .merge_posix.awk <<'AWK'
BEGIN{
  FS=","; OFS=",";
  n=split(TARGET, want, ",")
  for(i=1;i<=n;i++){ gsub(/^ *| *$/,"", want[i]); lwant[i]=tolower(want[i]) }
}
NR==1{
  for(i=1;i<=NF;i++){
    h=$i
    gsub(/\r/,"",h)
    sub(/^[ \t]+/,"",h); sub(/[ \t]+$/,"",h)
    l=tolower(h)
    hmap[l]=i
  }
  for(i=1;i<=n;i++){
    tmp=hmap[lwant[i]]
    if (tmp=="") idx[i]=0
    else idx[i]=tmp
  }
  next
}
{
  row=""
  for(i=1;i<=n;i++){
    if (idx[i]>0 && idx[i]<=NF) v=$idx[i]; else v=""
    gsub(/\r/,"",v)
    if (i==1) row=v; else row=row OFS v
  }
  print row
}
AWK

  local candidates=(
    "assemblies/assembly.busco.csv" "assemblies/assembly_busco.csv" "assemblies/busco.csv"
    "assemblies/assembly.quast.csv" "assemblies/assembly_quast.csv" "assemblies/quast.csv"
    "assemblies/assembly.telo.csv"  "assemblies/assembly_telo.csv"  "assemblies/telo.csv"
  )
  for f in "${candidates[@]}"; do
    if [[ -s "$f" ]]; then
      tr -d '\r' < "$f" | awk -v TARGET="$desired_header" -f .merge_posix.awk >> "$info_csv"
    fi
  done
  rm -f .merge_posix.awk

  echo "[ok] Wrote $(basename "$info_csv")"
  if command -v column >/dev/null 2>&1; then
    echo "==== assemblies/assembly_info.csv ===="
    column -s, -t "$info_csv" | sed 's/^/  /'
    echo "======================================"
  fi
}

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
    g) genomesize="$OPTARG" ;;
    t) threads="$OPTARG" ;;
    s) expand_steps "$OPTARG" ;;
    m) motif="$OPTARG" ;;
    -)
      case "${OPTARG}" in
        fastq)
          fastq="${!OPTIND}"
          OPTIND=$((OPTIND + 1))
          ;;
        fasta)
          external_fasta="${!OPTIND}"
          OPTIND=$((OPTIND + 1))
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
          CHOOSE_FLAG=1
          next_arg="${!OPTIND}"
          if [[ -n "$next_arg" && "$next_arg" != -* ]]; then
            assembler="$next_arg"
            OPTIND=$((OPTIND + 1))
          else
            echo "Please enter the assembler you want to use for the final refinement backbone (canu, nextDenovo, peregrine, ipa, flye, hifiasm):"
            read assembler
          fi
          ;;
        auto-mode)
          auto_mode="${!OPTIND}"
          OPTIND=$((OPTIND + 1))
          ;;
        assembly-only)
          assembly_only=true
          ;;
        merqury)
          merqury_enable=true
          ;;
        no-merqury)
          merqury_enable=false
          ;;
        merqury-db)
          merqury_enable=true
          merqury_db="${!OPTIND}"
          OPTIND=$((OPTIND + 1))
          ;;
        *)
          echo "Unknown option --${OPTARG}" >&2
          usage
          ;;
      esac
      ;;
    \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
    :)  echo "Option -$OPTARG requires an argument." >&2; usage ;;
  esac
done

# Ensure required options
if [[ -z "$genomesize" || -z "$threads" || -z "$fastq" || -z "$motif" ]]; then
  usage
fi

if [[ "$assembly_only" == true ]]; then
  echo "[info] Assembly-only mode enabled"
  steps=()
  expand_steps "1-9,11,18"
fi

write_versions
check_single_env_requirements

# Resolve external FASTA if provided
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
  if [[ "$external_fasta" =~ \.gz$ ]]; then
    echo "[info] Decompressing gzipped FASTA: $external_fasta"
    gz_out="external_input.fa"
    gunzip -c "$external_fasta" > "$gz_out" || { echo "[error] gunzip failed for $external_fasta"; exit 1; }
    external_fasta="$gz_out"
  fi
  if [[ ! -s "$external_fasta" ]]; then
    echo "[error] --fasta provided but file not found or empty: $external_fasta"
    exit 1
  fi
fi

# Project name robust for .fastq.gz
project="$(basename "${fastq%.gz}" .fastq)"

# nextDenovo config
cat <<EOT > ./run_${project}.cfg
[General]
job_type = local # local, slurm, sge, pbs, lsf
job_prefix = nextDenovo
task = assemble # all, correct, assemble
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
pa_correction = 3
correction_options = -p ${threads}

[assemble_option]
minimap2_options_cns = -t ${threads}
nextgraph_options = -a 1
EOT

fastq=$(realpath "$fastq")
echo "$fastq" > input_${project}.fofn
echo "$fastq" > reads_${project}.lst

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

require_cmd() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "[error] Required command not found in the active conda environment: $cmd" >&2
    exit 127
  fi
}

t2t_list_inline() {
  local seqtk_file=""
  local output="t2t.list"
  local single_output="single_tel.list"
  local supported_output="telomere_supported.list"
  local tel_window="${TEL_WINDOW:-100}"

  OPTIND=1
  while getopts ":i:o:s:u:w:" opt; do
    case "$opt" in
      i) seqtk_file="$OPTARG" ;;
      o) output="$OPTARG" ;;
      s) single_output="$OPTARG" ;;
      u) supported_output="$OPTARG" ;;
      w) tel_window="$OPTARG" ;;
      *) echo "Usage: t2t_list_inline -i <seqtk_telo_result> [-o <t2t_list>] [-s <single_end_list>] [-u <supported_list>] [-w <window_bp>]" >&2; return 1 ;;
    esac
  done
  shift $((OPTIND-1))

  if [[ -z "$seqtk_file" ]]; then
    echo "[error] t2t_list_inline: seqtk telo result file (-i) is required." >&2
    return 1
  fi

  if [[ ! -s "$seqtk_file" ]]; then
    : > "$output"
    : > "$single_output"
    : > "$supported_output"
    echo "[warn] t2t_list_inline: input '$seqtk_file' missing or empty; wrote empty output lists" >&2
    return 0
  fi

  : > "$output"
  : > "$single_output"
  : > "$supported_output"

  awk -v w="$tel_window" -v T2T="$output" -v SINGLE="$single_output" '
    NF>=4 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $4~/^[0-9]+$/ {
      c=$1
      sub(/[ \t].*$/, "", c)

      left=$2
      right=$3
      len=$4

      has_left  = (left <= w)
      has_right = ((len - right) <= w)

      if (has_left && has_right) print c >> T2T
      else if (has_left || has_right) print c >> SINGLE
    }
  ' "$seqtk_file"

  sort -u "$output" -o "$output" 2>/dev/null || true
  sort -u "$single_output" -o "$single_output" 2>/dev/null || true
  cat "$output" "$single_output" 2>/dev/null | sort -u > "$supported_output"

  echo "[INFO] Strict T2T contigs written to $output"
  echo "[INFO] Single-end telomeric contigs written to $single_output"
  echo "[INFO] Telomere-supported contigs written to $supported_output"
}

extract_by_list() {
  local listfile="$1"
  local infa="$2"
  local outfa="$3"

  if [[ ! -s "$listfile" ]]; then
    : > "$outfa"
    return 0
  fi

  awk '
  BEGIN {
    while ((getline < LST) > 0) ids[$1]=1
  }
  /^>/ {
    h=substr($0,2)
    split(h,a,/[\t ]/)
    keep=(a[1] in ids)
  }
  keep
  ' LST="$listfile" "$infa" > "$outfa"
}

check_single_env_requirements() {
  local cmds=(python3 canu nextDenovo pg_asm ipa flye hifiasm seqtk busco minimap2 bwa samtools funannotate redundans.py merge_wrapper.py)
  local missing=()
  for cmd in "${cmds[@]}"; do
    command -v "$cmd" >/dev/null 2>&1 || missing+=("$cmd")
  done
  if ! command -v quast.py >/dev/null 2>&1 && ! command -v quast >/dev/null 2>&1; then
    missing+=("quast.py/quast")
  fi
  if [[ ${#missing[@]} -gt 0 ]]; then
    echo "[error] Missing tools in the active conda environment: ${missing[*]}" >&2
    echo "[error] Create and activate the unified TACO environment first, then rerun the pipeline." >&2
    exit 127
  fi
}

echo "[info] Using the current active conda environment (single-env mode)."
write_versions
write_run_metadata
init_benchmark_step_table

if [ ${#steps[@]} -eq 0 ]; then
  steps=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17)
fi

for step in "${steps[@]}"; do
  STEP_LOG="${LOG_DIR}/step_${step}.log"
  step_started_epoch=$(date +%s)
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
      log_version "pg_asm" "pg_asm"
      echo "Step 3 - Assembly of the genome using Peregrine"
      pg_asm reads_${project}.lst peregrine-2021
      check_command
      ;;
    4)
      log_version "ipa" "ipa"
      echo "Step 4 - Assembly of the genome using IPA"
       # Clean previous run to avoid Snakemake conflict
       if [[ -d ipa ]]; then
          echo "[INFO] Removing existing IPA run directory"
          rm -rf ipa
        fi
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
      mkdir -p hifiasm
      log_version "hifiasm" "hifiasm"
      cd hifiasm
      hifiasm -o hifiasm.asm -t "$threads" "$fastq" 2> hifiasm.log
      check_command
      if [[ -f hifiasm.asm.bp.p_ctg.gfa ]]; then
        awk '/^S/{print ">"$2;print $3}' hifiasm.asm.bp.p_ctg.gfa > hifiasm.fasta
        check_command
      else
        echo "[error] Hifiasm primary contig GFA not found: hifiasm.asm.bp.p_ctg.gfa" >&2
        exit 1
      fi
      cd ..
      ;;
    7)
      echo "Step 7 - Copy all assemblies"
      mkdir -p assemblies
      cp ./hicanu/canu.contigs.fasta ./assemblies/canu.result.fasta
      tmp_renamed="assemblies/.canu.renamed.tmp.fasta"
      rename_and_sort_fasta "./assemblies/canu.result.fasta" "$tmp_renamed" "canu" && mv -f "$tmp_renamed" "./assemblies/canu.result.fasta"

      cp ./NextDenovo/03.ctg_graph/nd.asm.fasta ./assemblies/nextDenovo.result.fasta
      tmp_renamed="assemblies/.nextDenovo.renamed.tmp.fasta"
      rename_and_sort_fasta "./assemblies/nextDenovo.result.fasta" "$tmp_renamed" "nextDenovo" && mv -f "$tmp_renamed" "./assemblies/nextDenovo.result.fasta"

      cp ./peregrine-2021/asm_ctgs_m_p.fa ./assemblies/peregrine.result.fasta
      tmp_renamed="assemblies/.peregrine.renamed.tmp.fasta"
      rename_and_sort_fasta "./assemblies/peregrine.result.fasta" "$tmp_renamed" "peregrine" && mv -f "$tmp_renamed" "./assemblies/peregrine.result.fasta"

      cp ./ipa/assembly-results/final.p_ctg.fasta ./assemblies/ipa.result.fasta
      tmp_renamed="assemblies/.ipa.renamed.tmp.fasta"
      rename_and_sort_fasta "./assemblies/ipa.result.fasta" "$tmp_renamed" "ipa" && mv -f "$tmp_renamed" "./assemblies/ipa.result.fasta"

      cp ./flye/assembly.fasta ./assemblies/flye.result.fasta
      tmp_renamed="assemblies/.flye.renamed.tmp.fasta"
      rename_and_sort_fasta "./assemblies/flye.result.fasta" "$tmp_renamed" "flye" && mv -f "$tmp_renamed" "./assemblies/flye.result.fasta"

      cp ./hifiasm/hifiasm.fasta ./assemblies/hifiasm.result.fasta
      tmp_renamed="assemblies/.hifiasm.renamed.tmp.fasta"
      rename_and_sort_fasta "./assemblies/hifiasm.result.fasta" "$tmp_renamed" "hifiasm" && mv -f "$tmp_renamed" "./assemblies/hifiasm.result.fasta"

      if [[ -n "$external_fasta" ]]; then
        cp "$external_fasta" "./assemblies/external.result.fasta"
      fi
      ;;
    8)
      echo "Step 8 - Run BUSCO on all assembled genomes (including external)"

      shopt -s nullglob
      a=(assemblies/*.result.fasta)
      if [[ ${#a[@]} -eq 0 ]]; then
        echo "[error] No assemblies in ./assemblies for BUSCO."
        exit 1
      fi

      cols=("canu" "external" "flye" "ipa" "nextDenovo" "peregrine" "hifiasm")
      lineage=${busco_lineage:-"ascomycota_odb10"}
      threads=${threads:-8}

      [[ -d busco ]] || mkdir -p busco

      has_busco_metrics() {
        local base="$1"
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

      # Build assemblies/assembly.busco.csv
      pyfile="$(mktemp -t busco_matrix.XXXXXX.py)"
      cat > "$pyfile" <<'PY'
import os, re, glob, csv

OUT_CSV = os.path.join("assemblies", "assembly.busco.csv")
DESIRED = ["canu","external","flye","ipa","nextDenovo","peregrine","hifiasm"]

def newest(paths):
    return max(paths, key=os.path.getmtime) if paths else None

def find_run_root(base):
    candidates = []
    p1 = os.path.join("busco", base)
    p2 = os.path.join("busco", f"run_{base}")
    if os.path.isdir(p1): candidates.append(p1)
    if os.path.isdir(p2): candidates.append(p2)
    if not candidates:
        candidates = sorted(glob.glob(os.path.join("busco", f"{base}*"))) + \
                     sorted(glob.glob(os.path.join("busco", f"run_{base}*")))
    return newest(candidates)

def read_counts_from_anywhere(run_root):
    if not run_root: return None
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
    sums = glob.glob(os.path.join(run_root, "**", "short_summary*.txt"), recursive=True)
    if sums:
        txt = open(newest(sums), "r", errors="ignore").read()
        m = re.search(r"C:(\d+(?:\.\d+)?)%.*?S:(\d+(?:\.\d+)?)%.*?D:(\d+(?:\.\d+)?)%.*?F:(\d+(?:\.\d+)?)%.*?M:(\d+(?:\.\d+)?)%.*?n:(\d+)", txt, re.S)
        if not m: return None
        Cpct, Spct, Dpct, Fpct, Mpct = map(float, m.groups()[:5])
        n = int(m.group(6))
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
      mkdir -p assemblies
      shopt -s nullglob
      existing_assemblies=(assemblies/*.result.fasta)
      if [[ ${#existing_assemblies[@]} -eq 0 ]]; then
        [[ -s ./hicanu/canu.contigs.fasta ]] && cp ./hicanu/canu.contigs.fasta ./assemblies/canu.result.fasta
        [[ -s ./NextDenovo/03.ctg_graph/nd.asm.fasta ]] && cp ./NextDenovo/03.ctg_graph/nd.asm.fasta ./assemblies/nextDenovo.result.fasta
        [[ -s ./peregrine-2021/asm_ctgs_m_p.fa ]] && cp ./peregrine-2021/asm_ctgs_m_p.fa ./assemblies/peregrine.result.fasta
        [[ -s ./ipa/assembly-results/final.p_ctg.fasta ]] && cp ./ipa/assembly-results/final.p_ctg.fasta ./assemblies/ipa.result.fasta
        [[ -s ./flye/assembly.fasta ]] && cp ./flye/assembly.fasta ./assemblies/flye.result.fasta
        [[ -s ./hifiasm/hifiasm.fasta ]] && cp ./hifiasm/hifiasm.fasta ./assemblies/hifiasm.result.fasta
        if [[ -n "$external_fasta" && -s "$external_fasta" ]]; then
          cp "$external_fasta" ./assemblies/external.result.fasta
        fi
      fi

      existing_assemblies=(assemblies/*.result.fasta)
      if [[ ${#existing_assemblies[@]} -eq 0 ]]; then
        echo "[error] No assemblies found in ./assemblies. Run step 7 first or supply --fasta." >&2
        exit 1
      fi

      cols=("canu" "external" "flye" "ipa" "nextDenovo" "peregrine" "hifiasm")
      declare -A tdouble tsingle

      for fasta in assemblies/*.result.fasta; do
        asm="${fasta##*/}"; asm="${asm%.result.fasta}"; asm="${asm%.fasta}"
        list="${fasta%.result.fasta}.telo.list"
        out="${fasta%.result.fasta}.telo.fasta"

        log_version "seqtk" "seqtk"
        seqtk telo -s 1 -m "$motif" "$fasta" > "$list"
        check_command

        awk '{print $1}' "$list" | sed 's/[[:space:]].*$//' | tr -d $'\r' | sort -u > "${list}.ids"
        if [[ -s "${list}.ids" ]]; then
          seqtk subseq "$fasta" "${list}.ids" > "$out"
          [[ -s "$out" ]] && echo "[ok] Wrote $(basename "$out")" || echo "[warn] ${out} is empty"
        else
          echo "[warn] No telomere-containing contigs for $(basename "$fasta"); writing empty ${out}"
          : > "$out"
        fi

        double=$(awk 'NF>=4 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $4~/^[0-9]+$/ {c=$1; sub(/[ \t].*$/,"",c); if($2==0 && $3==$4) print c}' "$list" | sort -u | wc -l)
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

      matrix_csv="assemblies/assembly.telo.csv"
      {
        printf "Metric"
        for c in "${cols[@]}"; do printf ",%s" "$c"; done
        printf "\n"
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
      echo "Step 10 - Build optimized telomere contig pool"

      shopt -s nullglob

      fasta_files=(assemblies/*.telo.fasta)
      if [[ ${#fasta_files[@]} -eq 0 ]]; then
        echo "[info] No per-assembly *.telo.fasta files found; generating per-assembly telomere-contig FASTAs from ./assemblies/*.result.fasta using motif '${motif:-TTAGGG}'."

        for f in assemblies/*.result.fasta; do
          [[ -s "$f" ]] || continue
          base=$(basename "$f")
          base="${base%.result.fasta}"
          base="${base%.fasta}"

          M="${motif:-TTAGGG}"
          RC="$(echo "$M" | tr 'ACGTacgtnN' 'TGCAtgcanN' | rev)"

          awk -v W="${TELO_WINDOW:-200}" -v R="${TELO_MIN_REPEATS:-3}" -v M="$M" -v RC="$RC" '
            BEGIN { pat="(" M "|" RC "){" R ",}" }
            /^>/ {
              if (seq != "") {
                left=substr(seq,1,W)
                right=substr(seq,length(seq)-W+1,W)
                if ((left ~ pat) || (right ~ pat)) {
                  print header
                  print seq
                }
              }
              header=$0
              seq=""
              next
            }
            {
              gsub(/[ 	
]/,"")
              seq=seq $0
            }
            END {
              if (seq != "") {
                left=substr(seq,1,W)
                right=substr(seq,length(seq)-W+1,W)
                if ((left ~ pat) || (right ~ pat)) {
                  print header
                  print seq
                }
              }
            }
          ' "$f" > "assemblies/${base}.telo.fasta"

          [[ -s "assemblies/${base}.telo.fasta" ]] || rm -f "assemblies/${base}.telo.fasta"

          if [[ -s "assemblies/${base}.telo.fasta" ]]; then
            require_cmd seqtk
            seqtk telo -s 1 -m "$motif" "assemblies/${base}.telo.fasta" > "assemblies/${base}.telo.list" || true
          fi
        done

        fasta_files=(assemblies/*.telo.fasta)
        if [[ ${#fasta_files[@]} -eq 0 ]]; then
          echo "[error] Still no *.telo.fasta files after generation. Check Step 7 and motif '$motif'."
          exit 1
        fi
      fi

      fasta_files=(assemblies/*.telo.fasta)

      if [[ ${#fasta_files[@]} -eq 1 ]]; then
        echo "[info] Only one telo FASTA found: ${fasta_files[0]}"
        cp -f "${fasta_files[0]}" allmerged_telo.fasta
      else
        rm -f merged_*.fasta

        for ((i=0; i<${#fasta_files[@]}; i++)); do
          for ((j=i+1; j<${#fasta_files[@]}; j++)); do
            file1="${fasta_files[i]}"
            file2="${fasta_files[j]}"
            base1=$(basename "$file1" .telo.fasta)
            base2=$(basename "$file2" .telo.fasta)

            require_cmd merge_wrapper.py
            log_version "merge_wrapper.py" "merge_wrapper.py" 2>/dev/null || true
            merge_wrapper.py -l 1000000 "$file1" "$file2" --prefix merged_"$base1"_"$base2"
            check_command
          done
        done

        merged_list=(merged_*.fasta)
        if [[ ${#merged_list[@]} -eq 0 ]]; then
          echo "[warn] No merged_* files produced; falling back to concatenating input telo FASTAs."
          cat "${fasta_files[@]}" > allmerged_telo.fasta
        else
          cat "${merged_list[@]}" > allmerged_telo.fasta
        fi
      fi

      require_cmd funannotate
      funannotate sort -i allmerged_telo.fasta -b contig -o allmerged_telo_sort.fasta --minlen 500
      check_command

      require_cmd seqtk
      seqtk telo -s 1 -m "$motif" allmerged_telo_sort.fasta > allmerged.telo.list
      check_command

      t2t_list_inline -i allmerged.telo.list -o t2t.list -s single_tel.list -u telomere_supported.list

      extract_by_list t2t.list allmerged_telo_sort.fasta t2t.fasta
      extract_by_list single_tel.list allmerged_telo_sort.fasta single_tel.fasta
      extract_by_list telomere_supported.list allmerged_telo_sort.fasta telomere_supported.fasta

      if [[ -s single_tel.fasta ]] && command -v minimap2 >/dev/null 2>&1; then
        echo "[INFO] Optimizing single-end telomeric pool by redundancy reduction"
        log_version "minimap2" "minimap2" 2>/dev/null || true

        minimap2 -x asm20 -D -P -t "${threads:-8}" single_tel.fasta single_tel.fasta > single_tel.self.paf
        check_command

        python3 - <<'PY'
import sys
from pathlib import Path

fa = Path("single_tel.fasta")
paf = Path("single_tel.self.paf")
out_ids = Path("single_tel_best.ids")
out_tsv = Path("telomere_cluster_summary.tsv")

def read_fa(path):
    seqs = {}
    name = None
    buf = []
    with open(path) as f:
        for ln in f:
            if ln.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = ln[1:].strip().split()[0]
                buf = []
            else:
                buf.append(ln.strip())
        if name is not None:
            seqs[name] = "".join(buf)
    return seqs

seqs = read_fa(fa)
names = sorted(seqs)
parent = {n: n for n in names}

def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

def union(a, b):
    ra = find(a)
    rb = find(b)
    if ra != rb:
        parent[rb] = ra

for ln in open(paf):
    if not ln.strip() or ln.startswith("#"):
        continue
    p = ln.rstrip().split("	")
    if len(p) < 12:
        continue

    qname = p[0]
    qlen = int(p[1])
    qs = int(p[2])
    qe = int(p[3])
    tname = p[5]
    matches = int(p[9])
    alnlen = int(p[10])

    if qname == tname:
        continue
    if qname not in seqs or tname not in seqs or alnlen <= 0:
        continue

    ident = matches / max(1, alnlen)
    qcov = (qe - qs) / max(1, qlen)

    if ident >= 0.95 and qcov >= 0.85:
        union(qname, tname)

clusters = {}
for n in names:
    clusters.setdefault(find(n), []).append(n)

best = []
with open(out_tsv, "w") as o:
    o.write("cluster_id	representative	members
")
    for cid, members in sorted(clusters.items()):
        members = sorted(members, key=lambda x: (len(seqs[x]), x), reverse=True)
        rep = members[0]
        best.append(rep)
        o.write(f"{cid}	{rep}	{','.join(members)}
")

with open(out_ids, "w") as o:
    for n in sorted(best):
        o.write(n + "
")
PY
        check_command

        extract_by_list single_tel_best.ids single_tel.fasta single_tel_best.fasta
      else
        cp -f single_tel.fasta single_tel_best.fasta 2>/dev/null || : > single_tel_best.fasta
      fi

      : > telomere_supported_best.fasta
      if [[ -s t2t.fasta ]]; then
        cat t2t.fasta >> telomere_supported_best.fasta
      fi
      if [[ -s single_tel_best.fasta ]]; then
        cat single_tel_best.fasta >> telomere_supported_best.fasta
      fi
      if [[ ! -s telomere_supported_best.fasta && -s telomere_supported.fasta ]]; then
        cp -f telomere_supported.fasta telomere_supported_best.fasta
      fi

      require_cmd funannotate

      if [[ -s t2t.fasta ]]; then
        funannotate clean -i t2t.fasta -p 30 -o t2t_clean.fasta --exhaustive
        check_command
      else
        : > t2t_clean.fasta
      fi

      if [[ -s single_tel_best.fasta ]]; then
        funannotate clean -i single_tel_best.fasta -p 30 -o single_tel_best_clean.fasta --exhaustive
        check_command
      else
        : > single_tel_best_clean.fasta
      fi

      if [[ -s telomere_supported_best.fasta ]]; then
        funannotate clean -i telomere_supported_best.fasta -p 30 -o telomere_supported_best_clean.fasta --exhaustive
        check_command
      else
        : > telomere_supported_best_clean.fasta
      fi

      cp -f single_tel_best.fasta single_tel.fasta 2>/dev/null || true
      cp -f single_tel_best_clean.fasta single_tel_clean.fasta 2>/dev/null || true
      cp -f telomere_supported_best.fasta telomere_supported.fasta 2>/dev/null || true
      cp -f telomere_supported_best_clean.fasta telomere_supported_clean.fasta 2>/dev/null || true

      strict_t2t_n=$(grep -c '^>' t2t.fasta 2>/dev/null || true)
      single_tel_n=$(grep -c '^>' single_tel_best.fasta 2>/dev/null || true)
      tel_supported_n=$(grep -c '^>' telomere_supported_best.fasta 2>/dev/null || true)

      {
        echo "strict_t2t_contigs,${strict_t2t_n:-0}"
        echo "single_telomere_best_contigs,${single_tel_n:-0}"
        echo "telomere_supported_best_contigs,${tel_supported_n:-0}"
      } > telomere_support_summary.csv

      echo "[INFO] Telomere support summary:"
      cat telomere_support_summary.csv

      if [[ -s t2t_clean.fasta ]]; then
        cp -f t2t_clean.fasta protected_telomere_contigs.fasta
        echo "strict_t2t" > protected_telomere_mode.txt
        echo "[INFO] Protected contigs mode: strict_t2t"
      elif [[ -s single_tel_best_clean.fasta ]]; then
        cp -f single_tel_best_clean.fasta protected_telomere_contigs.fasta
        echo "single_tel_best" > protected_telomere_mode.txt
        echo "[INFO] Protected contigs mode: single_tel_best"
      elif [[ -s telomere_supported_best_clean.fasta ]]; then
        cp -f telomere_supported_best_clean.fasta protected_telomere_contigs.fasta
        echo "telomere_supported_best" > protected_telomere_mode.txt
        echo "[INFO] Protected contigs mode: telomere_supported_best"
      else
        : > protected_telomere_contigs.fasta
        echo "none" > protected_telomere_mode.txt
        echo "[WARN] No telomere-supported contigs found"
      fi
      ;;

    11)
      echo "Step 11 - QUAST metrics for all assemblies"

      shopt -s nullglob
      a=(assemblies/*.result.fasta)
      if [[ ${#a[@]} -eq 0 ]]; then
        echo "[error] No assemblies to run QUAST." >&2
        exit 1
      fi

      mkdir -p quast_out assemblies

      QUAST_BIN="$(command -v quast.py || command -v quast || true)"
      if [[ -n "$QUAST_BIN" ]]; then
        "$QUAST_BIN" "${a[@]}" --threads "${threads:-8}" -o quast_out 2>&1 | tee quast.log || true
      else
        echo "[warn] quast.py/quast not found; will use existing quast_out/report files if present." >&2
      fi

      python3 - <<'PY'
import os, csv, sys

OUT = os.path.join("assemblies", "assembly.quast.csv")
TREPORT = os.path.join("quast_out", "transposed_report.tsv")
REPORT  = os.path.join("quast_out", "report.tsv")

DESIRED = ["canu","external","flye","ipa","nextDenovo","peregrine","hifiasm"]

def norm(s):
    return (s or "").lower().replace("-", "").replace("_", "").replace(".result", "").replace(".fasta","").strip()

def write_csv(rows, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        csv.writer(f).writerows(rows)

def build_from_report(pth):
    with open(pth, newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    if not rows: return None
    header = rows[0]
    asm_headers = header[1:]
    n_asm = len(asm_headers)

    mapcol = {}
    used = set()
    for d in DESIRED:
        nd = norm(d)
        idx = None
        for j, h in enumerate(asm_headers):
            if j in used: continue
            if nd in (norm(h),) or nd == norm(h) or norm(h).startswith(nd) or nd in norm(h):
                idx = j; break
        if idx is None and nd == "hifiasm":
            for j, h in enumerate(asm_headers):
                nh = norm(h)
                if "hifiasm" in nh:
                    idx = j; break
        mapcol[d] = idx
        if idx is not None: used.add(idx)

    out = [["Metric"] + DESIRED]
    for r in rows[1:]:
        if not r: continue
        metric = r[0]
        r = r + [""] * (1 + n_asm - len(r))
        vals = []
        for d in DESIRED:
            j = mapcol[d]
            vals.append(r[1 + j] if j is not None else "")
        out.append([metric] + vals)
    return out

def build_from_transposed(pth):
    with open(pth, newline="") as f:
        rows = list(csv.reader(f, delimiter="\t"))
    if not rows: return None
    header = rows[0]
    metrics = header[1:]
    maprow = {d: None for d in DESIRED}
    for i in range(1, len(rows)):
        asm_name = rows[i][0] if rows[i] else ""
        na = norm(asm_name)
        for d in DESIRED:
            nd = norm(d)
            if maprow[d] is not None: continue
            if nd == na or nd in na or na in nd:
                maprow[d] = i
            elif nd == "hifiasm" and ("hifiasm" in na):
                maprow[d] = i

    out = [["Metric"] + DESIRED]
    for j, m in enumerate(metrics, start=1):
        line = [m]
        for d in DESIRED:
            i = maprow[d]
            if i is None:
                line.append("")
            else:
                row = rows[i]
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
      echo "[ok] Wrote assemblies/assembly.quast.csv"
      ;;
    12)
      echo "Step 12 - Final assembly refinement with optimized telomere-end replacement"

      assembler="${assembler:-}"
      CHOOSE_FLAG="${CHOOSE_FLAG:-0}"
      AUTO_MODE="${auto_mode:-smart}"
      MERQURY_ENABLE=0
      [[ "${merqury_enable:-false}" == true ]] && MERQURY_ENABLE=1
      MERQURY_DB="${merqury_db:-}"

      echo "[info] Auto-selection mode: $AUTO_MODE"
      echo "[info] Merqury enabled: $MERQURY_ENABLE"

      mkdir -p merqury assemblies

      # ------------------------------------------------------------
      # 12A. Optional Merqury pre-selection
      # ------------------------------------------------------------
      if [[ "$MERQURY_ENABLE" -eq 1 ]]; then
        if [[ -z "$MERQURY_DB" ]]; then
          for cand in reads.meryl meryl/reads.meryl merqury/reads.meryl *.meryl; do
            [[ -e "$cand" ]] || continue
            if [[ -d "$cand" ]]; then
              MERQURY_DB="$cand"
              break
            fi
          done
        fi

        if command -v merqury.sh >/dev/null 2>&1 && [[ -n "$MERQURY_DB" && -d "$MERQURY_DB" ]]; then
          echo "[INFO] Running Merqury pre-selection using database: $MERQURY_DB"
          log_version "merqury.sh" "merqury.sh" 2>/dev/null || true

          declare -A asm_paths=(
            [canu]="assemblies/canu.result.fasta"
            [external]="assemblies/external.result.fasta"
            [flye]="assemblies/flye.result.fasta"
            [ipa]="assemblies/ipa.result.fasta"
            [nextDenovo]="assemblies/nextDenovo.result.fasta"
            [peregrine]="assemblies/peregrine.result.fasta"
            [hifiasm]="assemblies/hifiasm.result.fasta"
          )

          for asm_name in canu external flye ipa nextDenovo peregrine hifiasm; do
            asm_fa_i="${asm_paths[$asm_name]}"
            [[ -s "$asm_fa_i" ]] || continue
            if [[ ! -f "merqury/${asm_name}.qv" || ! -f "merqury/${asm_name}.completeness.stats" ]]; then
              merqury.sh "$MERQURY_DB" "$asm_fa_i" "merqury/${asm_name}" || \
                echo "[warn] Merqury failed for ${asm_name}" >&2
            fi
          done
        else
          echo "[warn] Merqury requested but merqury.sh or a valid .meryl database was not found; skipping Merqury." >&2
        fi
      else
        echo "[INFO] Merqury disabled for this run"
      fi

      # Always write Merqury summary CSV
      python3 - <<'PY'
import os, csv, re

assemblers = ["canu","external","flye","ipa","nextDenovo","peregrine","hifiasm"]
rows = [
    ["Metric"] + assemblers,
    ["Merqury QV"],
    ["Merqury completeness (%)"],
]

def parse_first_float(path):
    if not os.path.exists(path):
        return ""
    txt = open(path, "r", errors="ignore").read()
    patterns = [
        r'(?i)\bqv\b[^0-9]*([0-9]+(?:\.[0-9]+)?)',
        r'(?i)\bcompleteness\b[^0-9]*([0-9]+(?:\.[0-9]+)?)',
        r'([0-9]+(?:\.[0-9]+)?)'
    ]
    for pat in patterns:
        m = re.search(pat, txt)
        if m:
            return m.group(1)
    return ""

for asm in assemblers:
    qv_file = os.path.join("merqury", f"{asm}.qv")
    comp_file = os.path.join("merqury", f"{asm}.completeness.stats")
    rows[1].append(parse_first_float(qv_file))
    rows[2].append(parse_first_float(comp_file))

with open("assemblies/assembly.merqury.csv", "w", newline="") as f:
    csv.writer(f).writerows(rows)

print("Wrote assemblies/assembly.merqury.csv")
PY
      check_command

      build_assembly_info_v2

      # ------------------------------------------------------------
      # 12B. Auto-select backbone assembler
      # ------------------------------------------------------------
      if [[ $CHOOSE_FLAG -eq 0 && -z "$assembler" ]]; then
        if [[ ! -s assemblies/assembly_info.csv ]]; then
          echo "[warn] assemblies/assembly_info.csv missing or empty; cannot auto-select assembler." >&2
        else
          echo "[info] Selection criteria: BUSCO + telomere + Merqury + contiguity + N50"

          assembler="$(
AUTO_MODE_ENV="$AUTO_MODE" python3 - <<'PY'
import csv, math, os, sys
from pathlib import Path

mode = os.environ.get("AUTO_MODE_ENV", "smart").strip().lower()
info = Path("assemblies") / "assembly_info.csv"
debug_tsv = Path("assemblies") / "selection_debug.tsv"
decision_txt = Path("assemblies") / "selection_decision.txt"

try:
    with info.open(newline="") as f:
        rows = list(csv.reader(f))
except Exception:
    sys.exit(0)

if not rows:
    sys.exit(0)

header = [h.strip() for h in rows[0]]
rows = rows[1:]

busco_row = None
contig_row = None
n50_row = None
t2t_row = None
single_row = None
merqury_qv_row = None
merqury_comp_row = None

for r in rows:
    if not r:
        continue
    name = r[0].strip().lower()
    if "busco c (%)" in name:
        busco_row = r
    elif name == "# contigs":
        contig_row = r
    elif name == "n50":
        n50_row = r
    elif "telomere double-end contigs" in name:
        t2t_row = r
    elif "telomere single-end contigs" in name:
        single_row = r
    elif "merqury qv" in name:
        merqury_qv_row = r
    elif "merqury completeness (%)" in name:
        merqury_comp_row = r

best_name = None
best_score = None
records = []

for idx, asm in enumerate(header[1:], start=1):
    asm = asm.strip()
    if not asm:
        continue

    fa = Path("assemblies") / f"{asm}.result.fasta"
    if not fa.is_file() or fa.stat().st_size == 0:
        continue

    def get(row, default=0.0):
        try:
            return float(row[idx])
        except Exception:
            return default

    busco = get(busco_row, 0.0)
    contigs = get(contig_row, 0.0)
    n50 = get(n50_row, 0.0)
    t2t = get(t2t_row, 0.0)
    single = get(single_row, 0.0)
    merqury_qv = get(merqury_qv_row, 0.0)
    merqury_comp = get(merqury_comp_row, 0.0)

    if mode == "n50":
        if n50 <= 0:
            continue
        score = n50
    else:
        if contigs <= 0 or n50 <= 0:
            continue
        score = (
            busco * 1000
            + t2t * 600
            + single * 250
            + merqury_comp * 200
            + merqury_qv * 20
            - contigs * 10
            + math.log10(n50) * 100
        )

    records.append({
        "assembler": asm,
        "busco_c": busco,
        "t2t": t2t,
        "single_tel": single,
        "merqury_qv": merqury_qv,
        "merqury_comp": merqury_comp,
        "contigs": contigs,
        "n50": n50,
        "score": score,
    })

    print(
        f"[DEBUG] {asm}: BUSCO={busco} T2T={t2t} single={single} "
        f"MerquryQV={merqury_qv} MerquryComp={merqury_comp} "
        f"contigs={contigs} N50={n50} score={score}",
        file=sys.stderr
    )

    if best_score is None or score > best_score:
        best_name = asm
        best_score = score

with debug_tsv.open("w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow([
        "assembler","busco_c","t2t","single_tel",
        "merqury_qv","merqury_comp","contigs","n50","score"
    ])
    for r in records:
        w.writerow([
            r["assembler"], r["busco_c"], r["t2t"], r["single_tel"],
            r["merqury_qv"], r["merqury_comp"], r["contigs"], r["n50"], r["score"]
        ])

with decision_txt.open("w") as f:
    f.write(f"auto_mode\t{mode}\n")
    f.write(f"selected_assembler\t{best_name or ''}\n")
    f.write(f"selected_score\t{best_score if best_score is not None else ''}\n")
    f.write("score_formula\tBUSCO*1000 + T2T*600 + single*250 + MerquryComp*200 + MerquryQV*20 - contigs*10 + log10(N50)*100\n")

if best_name:
    print(best_name)
PY
          )"

          if [[ -n "$assembler" ]]; then
            echo "[info] Auto-selected assembler: $assembler"
          else
            echo "[warn] Auto-selection failed; fallback required." >&2
          fi
        fi
      fi

      valid_assemblers="canu external flye ipa nextDenovo peregrine hifiasm"
      if [[ -z "${assembler:-}" ]]; then
        if ! prompt_from_tty assembler "Enter the assembler to use for the final refinement (e.g., ${valid_assemblers// /, }): "; then
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

      asm_fa="assemblies/${assembler}.result.fasta"
      if [[ ! -s "$asm_fa" ]]; then
        echo "[error] Selected assembler '$assembler' does not have '$asm_fa'." >&2
        exit 1
      fi

      protected_fasta=""
      protected_mode="none"

      if [[ -s protected_telomere_contigs.fasta ]]; then
        protected_fasta="protected_telomere_contigs.fasta"
        if [[ -s protected_telomere_mode.txt ]]; then
          protected_mode="$(cat protected_telomere_mode.txt)"
        else
          protected_mode="telomere_supported"
        fi
      fi

      echo "[INFO] Protected mode for final refinement: $protected_mode"
      echo "[INFO] Selected backbone assembly: $asm_fa"

      # ------------------------------------------------------------
      # 12C. Prepare cleaned backbone
      # ------------------------------------------------------------
      python3 - "$asm_fa" <<'PY'
import sys
from pathlib import Path

inp = Path(sys.argv[1])
out = Path("assemblies/backbone.clean.fa")

seen = {}
with inp.open() as f, out.open("w") as o:
    name = None
    seq = []

    def flush():
        if name is None:
            return
        base = name.split()[0]
        base = "".join(c if c.isalnum() or c in "._-" else "_" for c in base)
        if base in seen:
            seen[base] += 1
            new = f"{base}_{seen[base]}"
        else:
            seen[base] = 1
            new = base
        s = "".join(seq)
        if not s:
            return
        o.write(f">{new}\n")
        for i in range(0, len(s), 60):
            o.write(s[i:i+60] + "\n")

    for ln in f:
        if ln.startswith(">"):
            flush()
            name = ln[1:].strip()
            seq = []
        else:
            seq.append(ln.strip())
    flush()
PY
      check_command

      backbone_fa="assemblies/backbone.clean.fa"

      # ------------------------------------------------------------
      # 12D. Remove backbone contigs redundant to protected telomere pool
      # ------------------------------------------------------------
      if [[ -n "$protected_fasta" && -s "$protected_fasta" ]]; then
        echo "[INFO] Using protected telomere contigs from $protected_fasta"
        cp -f "$protected_fasta" assemblies/protected.telomere.fa

        if command -v minimap2 >/dev/null 2>&1; then
          log_version "minimap2" "minimap2" 2>/dev/null || true
          paf="assemblies/backbone_vs_protected.paf"
          minimap2 -x asm20 -t "${threads:-8}" assemblies/protected.telomere.fa "$backbone_fa" > "$paf"
          check_command

          python3 - "$paf" "$backbone_fa" > assemblies/backbone.keep.ids <<'PY'
import os, sys

cov_thr = float(os.getenv("PROTECT_COV", "0.95"))
id_thr = float(os.getenv("PROTECT_ID", "0.95"))

best = {}
with open(sys.argv[1]) as f:
    for ln in f:
        if not ln.strip() or ln.startswith("#"):
            continue
        p = ln.rstrip().split("\t")
        if len(p) < 12:
            continue

        qname = p[0]
        qlen = int(p[1])
        qs = int(p[2])
        qe = int(p[3])
        matches = int(p[9])
        alnlen = int(p[10])

        if qlen <= 0 or alnlen <= 0:
            continue

        ident = matches / max(1, alnlen)
        cov = (qe - qs) / max(1, qlen)

        cur = best.get(qname, (0.0, 0.0))
        if cov > cur[0] or (abs(cov - cur[0]) < 1e-12 and ident > cur[1]):
            best[qname] = (cov, ident)

drop = {q for q, (cov, ident) in best.items() if cov >= cov_thr and ident >= id_thr}

keep = []
with open(sys.argv[2]) as f:
    name = None
    for ln in f:
        if ln.startswith(">"):
            if name is not None and name not in drop:
                keep.append(name)
            name = ln[1:].strip().split()[0]
    if name is not None and name not in drop:
        keep.append(name)

print("\n".join(keep))
PY
          check_command

          python3 - assemblies/backbone.keep.ids "$backbone_fa" assemblies/backbone.filtered.fa <<'PY'
import sys

ids_file, inp, outp = sys.argv[1:4]

with open(ids_file) as f:
    ids = {ln.strip() for ln in f if ln.strip()}

with open(inp) as f, open(outp, "w") as o:
    name = None
    seq = []

    def flush():
        if name is not None and name in ids:
            o.write(f">{name}\n")
            s = "".join(seq)
            for i in range(0, len(s), 60):
                o.write(s[i:i+60] + "\n")

    for ln in f:
        if ln.startswith(">"):
            flush()
            name = ln[1:].strip().split()[0]
            seq = []
        else:
            seq.append(ln.strip())
    flush()
PY
          check_command
        else
          echo "[warn] minimap2 not found; cannot filter redundant backbone contigs. Keeping all selected-assembly contigs." >&2
          cp -f "$backbone_fa" assemblies/backbone.filtered.fa
        fi

        awk '
        BEGIN {
          while ((getline < PFA) > 0) {
            if ($0 ~ /^>/) {
              h=substr($0,2)
              split(h,a,/[\t ]/)
              protected[a[1]]=1
            }
          }
        }
        /^>/ {
          h=substr($0,2)
          split(h,a,/[\t ]/)
          keep=!(a[1] in protected)
        }
        keep
        ' PFA="assemblies/protected.telomere.fa" assemblies/backbone.filtered.fa > assemblies/backbone.filtered.nodup.fa
      else
        echo "[WARN] No protected telomere contigs found; keeping backbone for telomere rescue only"
        : > assemblies/protected.telomere.fa
        cp -f "$backbone_fa" assemblies/backbone.filtered.nodup.fa
      fi

      # ------------------------------------------------------------
      # 12E. Telomere rescue using optimized best single-end pool
      # ------------------------------------------------------------
      single_tel_src=""
      if [[ -s single_tel_best_clean.fasta ]]; then
        single_tel_src="single_tel_best_clean.fasta"
      elif [[ -s single_tel_clean.fasta ]]; then
        single_tel_src="single_tel_clean.fasta"
      elif [[ -s single_tel.fasta ]]; then
        single_tel_src="single_tel.fasta"
      fi

      if [[ -n "$single_tel_src" ]]; then
        echo "[INFO] Running telomere rescue using $single_tel_src"
        cp -f assemblies/backbone.filtered.nodup.fa assemblies/backbone.telomere_rescued.fa

        if command -v minimap2 >/dev/null 2>&1; then
          rescue_paf="assemblies/single_tel_vs_backbone.paf"
          minimap2 -x asm20 -t "${threads:-8}" assemblies/backbone.telomere_rescued.fa "$single_tel_src" > "$rescue_paf"
          check_command

          python3 - "$single_tel_src" assemblies/backbone.telomere_rescued.fa "$rescue_paf" assemblies/single_tel.replaced.ids assemblies/backbone.telomere_rescued.next.fa <<'PY'
import sys
from pathlib import Path

single_fa = Path(sys.argv[1])
backbone_fa = Path(sys.argv[2])
paf = Path(sys.argv[3])
out_ids = Path(sys.argv[4])
out_fa = Path(sys.argv[5])

min_ident = 0.90
min_cov = 0.80
min_ext = 500

def read_fa(path):
    seqs = {}
    name = None
    buf = []
    with open(path) as f:
        for ln in f:
            if ln.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = ln[1:].strip().split()[0]
                buf = []
            else:
                buf.append(ln.strip())
        if name is not None:
            seqs[name] = "".join(buf)
    return seqs

single = read_fa(single_fa)
backbone = read_fa(backbone_fa)

best = {}
for ln in open(paf):
    if not ln.strip() or ln.startswith("#"):
        continue
    p = ln.rstrip().split("\t")
    if len(p) < 12:
        continue

    qname = p[0]
    qlen = int(p[1])
    qs = int(p[2])
    qe = int(p[3])
    tname = p[5]
    matches = int(p[9])
    alnlen = int(p[10])

    if qname not in single or tname not in backbone:
        continue
    if alnlen <= 0 or qlen <= 0:
        continue

    ident = matches / max(1, alnlen)
    cov_q = (qe - qs) / max(1, qlen)
    left_overhang = qs
    right_overhang = qlen - qe
    ext = max(left_overhang, right_overhang)

    if ident < min_ident or cov_q < min_cov or ext < min_ext:
        continue

    cur = best.get(tname)
    score = (ext, ident, cov_q, len(single[qname]), qname)
    if cur is None or score > cur[0]:
        best[tname] = (score, qname)

replaced = set()
with out_ids.open("w") as ids:
    for tname, (_, qname) in sorted(best.items()):
        replaced.add(tname)
        ids.write(f"{tname}\t{qname}\n")

used_single = set(qname for _, qname in best.values())

with out_fa.open("w") as out:
    for qname in sorted(used_single):
        seq = single[qname]
        out.write(f">{qname}\n")
        for i in range(0, len(seq), 60):
            out.write(seq[i:i+60] + "\n")

    for tname, seq in backbone.items():
        if tname in replaced:
            continue
        out.write(f">{tname}\n")
        for i in range(0, len(seq), 60):
            out.write(seq[i:i+60] + "\n")
PY
          check_command

          if [[ -s assemblies/backbone.telomere_rescued.next.fa ]]; then
            mv -f assemblies/backbone.telomere_rescued.next.fa assemblies/backbone.telomere_rescued.fa
          fi
        else
          echo "[warn] minimap2 not found; skipping telomere rescue." >&2
        fi
      else
        echo "[INFO] No single-end telomeric contigs available for telomere rescue"
        cp -f assemblies/backbone.filtered.nodup.fa assemblies/backbone.telomere_rescued.fa
      fi

      # ------------------------------------------------------------
      # 12F. Remove rescued contigs redundant to protected set
      # ------------------------------------------------------------
      if [[ -s assemblies/protected.telomere.fa && -s assemblies/backbone.telomere_rescued.fa ]] && command -v minimap2 >/dev/null 2>&1; then
        dedup_paf="assemblies/rescued_vs_protected.paf"
        minimap2 -x asm20 -t "${threads:-8}" assemblies/protected.telomere.fa assemblies/backbone.telomere_rescued.fa > "$dedup_paf"
        check_command

        python3 - "$dedup_paf" assemblies/backbone.telomere_rescued.fa > assemblies/backbone.telomere_rescued.keep.ids <<'PY'
import sys

cov_thr = 0.95
id_thr = 0.95
best = {}

with open(sys.argv[1]) as f:
    for ln in f:
        if not ln.strip() or ln.startswith("#"):
            continue
        p = ln.rstrip().split("\t")
        if len(p) < 12:
            continue

        qname = p[0]
        qlen = int(p[1])
        qs = int(p[2])
        qe = int(p[3])
        matches = int(p[9])
        alnlen = int(p[10])

        if qlen <= 0 or alnlen <= 0:
            continue

        ident = matches / max(1, alnlen)
        cov = (qe - qs) / max(1, qlen)

        cur = best.get(qname, (0.0, 0.0))
        if cov > cur[0] or (abs(cov-cur[0]) < 1e-12 and ident > cur[1]):
            best[qname] = (cov, ident)

drop = {q for q, (cov, ident) in best.items() if cov >= cov_thr and ident >= id_thr}

keep = []
with open(sys.argv[2]) as f:
    name = None
    for ln in f:
        if ln.startswith(">"):
            if name is not None and name not in drop:
                keep.append(name)
            name = ln[1:].strip().split()[0]
    if name is not None and name not in drop:
        keep.append(name)

print("\n".join(keep))
PY
        check_command

        python3 - assemblies/backbone.telomere_rescued.keep.ids assemblies/backbone.telomere_rescued.fa assemblies/backbone.telomere_rescued.dedup.fa <<'PY'
import sys

ids_file, inp, outp = sys.argv[1:4]

with open(ids_file) as f:
    ids = {ln.strip() for ln in f if ln.strip()}

with open(inp) as f, open(outp, "w") as o:
    name = None
    seq = []

    def flush():
        if name is not None and name in ids:
            o.write(f">{name}\n")
            s = "".join(seq)
            for i in range(0, len(s), 60):
                o.write(s[i:i+60] + "\n")

    for ln in f:
        if ln.startswith(">"):
            flush()
            name = ln[1:].strip().split()[0]
            seq = []
        else:
            seq.append(ln.strip())
    flush()
PY
        check_command
      else
        cp -f assemblies/backbone.telomere_rescued.fa assemblies/backbone.telomere_rescued.dedup.fa
      fi

      # ------------------------------------------------------------
      # 12G. Final combine
      # ------------------------------------------------------------
      if [[ -s assemblies/protected.telomere.fa && -s assemblies/backbone.telomere_rescued.dedup.fa ]]; then
        cat assemblies/protected.telomere.fa assemblies/backbone.telomere_rescued.dedup.fa > assemblies/final_merge.raw.fasta
      elif [[ -s assemblies/protected.telomere.fa ]]; then
        cp -f assemblies/protected.telomere.fa assemblies/final_merge.raw.fasta
      elif [[ -s assemblies/backbone.telomere_rescued.dedup.fa ]]; then
        cp -f assemblies/backbone.telomere_rescued.dedup.fa assemblies/final_merge.raw.fasta
      else
        cp -f "$asm_fa" assemblies/final_merge.raw.fasta
      fi

      echo "[ok] Built assemblies/final_merge.raw.fasta (mode: $protected_mode)"

      require_cmd funannotate
      log_version "funannotate" "funannotate" 2>/dev/null || true
      funannotate sort -i assemblies/final_merge.raw.fasta -b contig -o "merged_${assembler}_sort.fa" --minlen 500
      check_command

      [[ -d assemblies ]] || { echo "[error] 'assemblies' folder not found (expected to exist)." >&2; exit 1; }
      cp -f "merged_${assembler}_sort.fa" "assemblies/final.merged.fasta"
      echo "[ok] Wrote assemblies/final.merged.fasta"
      ;;
    13)
      echo "Step 13 - BUSCO analysis"

      threads=${threads:-8}
      lineage=${busco_lineage:-"ascomycota_odb10"}

      final_fa="assemblies/final.merged.fasta"
      if [[ ! -s "$final_fa" ]]; then
        echo "[error] Final merged FASTA '$final_fa' not found." >&2
        exit 1
      fi
      [[ -d assemblies ]] || { echo "[error] 'assemblies' folder not found."; exit 1; }
      [[ -d busco ]] || mkdir busco

      has_final_metrics() {
        for d in "busco/final" "busco/run_final"; do
          [[ -d "$d" ]] || continue
          if find "$d" -type f \( -name 'short_summary*.json' -o -name 'short_summary*.txt' -o -name 'full_table*.tsv' \) -print -quit | grep -q .; then
            return 0
          fi
        done
        return 1
      }

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

      pyfile="$(mktemp -t final_busco.XXXXXX.py)"
      cat > "$pyfile" <<'PY'
import os, re, glob, json, csv

OUT_CSV = os.path.join("assemblies", "merged.busco.csv")

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
    ["Metric", "merged"],
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
      echo "[ok] Wrote assemblies/merged.busco.csv (lineage: ${lineage})"
      ;;
    14)
      echo "Step 14 - Telomere analysis (final assembly)"

      [[ -d assemblies ]] || { echo "[error] 'assemblies' folder not found."; exit 1; }
      final_fa="assemblies/final.merged.fasta"
      if [[ ! -s "$final_fa" ]]; then
        echo "[error] Final merged FASTA '$final_fa' not found." >&2
        exit 1
      fi
      if ! command -v seqtk >/dev/null 2>&1; then
        echo "[error] seqtk not found in the active conda environment." >&2
        exit 127
      fi

      list="assemblies/final.telo.list"
      ids="assemblies/final.telo.list.ids"
      out="assemblies/final.telo.fasta"
      csv="assemblies/merged.telo.csv"
      tel_window="${TEL_WINDOW:-100}"

      seqtk telo -s 1 -m "${motif:-TTAGGG}" "$final_fa" > "$list"
      check_command

      awk -v w="$tel_window" '
        NF>=4 { gsub(/
/, "") }
        NF>=4 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $4~/^[0-9]+$/ {
          c=$1; sub(/[ 	].*$/, "", c)
          left=$2; right=$3; len=$4
          if (left <= w || (len - right) <= w) print c
        }
      ' "$list" | sort -u > "$ids"

      if [[ -s "$ids" ]]; then
        seqtk subseq "$final_fa" "$ids" > "$out"
        [[ -s "$out" ]] && echo "[ok] Wrote $(basename "$out")" || echo "[warn] ${out} is empty"
      else
        echo "[warn] No telomere-containing contigs for final assembly; writing empty ${out}"
        : > "$out"
      fi

      double=$(awk -v w="$tel_window" '
        NF>=4 { gsub(/
/, "") }
        NF>=4 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $4~/^[0-9]+$/ {
          c=$1; sub(/[ 	].*$/, "", c)
          left=$2; right=$3; len=$4
          has_left=(left <= w); has_right=((len - right) <= w)
          if (has_left && has_right) print c
        }
      ' "$list" | sort -u | wc -l)

      single=$(awk -v w="$tel_window" '
        NF>=4 { gsub(/
/, "") }
        NF>=4 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $4~/^[0-9]+$/ {
          c=$1; sub(/[ 	].*$/, "", c)
          left=$2; right=$3; len=$4
          has_left=(left <= w); has_right=((len - right) <= w)
          if ((has_left + has_right) == 1) print c
        }
      ' "$list" | sort -u | wc -l)

      total_supported=$(( ${double:-0} + ${single:-0} ))

      protected_mode="none"
      [[ -s protected_telomere_mode.txt ]] && protected_mode="$(cat protected_telomere_mode.txt)"

      strict_pool_n=$(grep -c '^>' t2t_clean.fasta 2>/dev/null || echo 0)
      single_best_pool_n=$(grep -c '^>' single_tel_best_clean.fasta 2>/dev/null || echo 0)
      tel_pool_n=$(grep -c '^>' telomere_supported_best_clean.fasta 2>/dev/null || echo 0)

      rescue_n=0
      if [[ -s assemblies/single_tel.replaced.ids ]]; then
        rescue_n=$(wc -l < assemblies/single_tel.replaced.ids 2>/dev/null || echo 0)
      fi

      {
        echo "Metric,merged"
        echo "Telomere double-end contigs,${double:-0}"
        echo "Telomere single-end contigs,${single:-0}"
        echo "Telomere-supported contigs,${total_supported:-0}"
        echo "Protected telomere mode,${protected_mode}"
        echo "Step10 strict T2T pool contigs,${strict_pool_n:-0}"
        echo "Step10 best single-end telomere pool contigs,${single_best_pool_n:-0}"
        echo "Step10 optimized telomere-supported pool contigs,${tel_pool_n:-0}"
        echo "Step12 rescued telomere replacements,${rescue_n:-0}"
      } > "$csv"

      echo "[ok] Wrote $(basename "$csv")"
      ;;

    15)
      echo "Step 15 - QUAST metrics for final assembly"

      final_fa="assemblies/final.merged.fasta"
      if [[ ! -s "$final_fa" ]]; then
        echo "[error] Final merged FASTA '$final_fa' not found." >&2
        exit 1
      fi
      [[ -d assemblies ]] || { echo "[error] 'assemblies' folder not found."; exit 1; }

      mkdir -p quast_final
      QUAST_BIN="$(command -v quast.py || command -v quast || true)"
      if [[ -n "$QUAST_BIN" ]]; then
        "$QUAST_BIN" "$final_fa" --labels final --threads "${threads:-8}" -o quast_final 2>&1 | tee quast_final.log || true
      else
        echo "[warn] quast.py/quast not found; will try to use any existing quast_final/*report*.tsv" >&2
      fi

      python3 - <<'PY'
import os, csv, sys

out_dir_csv = os.path.join("assemblies")
if not os.path.isdir(out_dir_csv):
    sys.stderr.write("[error] 'assemblies' folder not found.\n")
    sys.exit(1)

report  = os.path.join("quast_final", "report.tsv")
treport = os.path.join("quast_final", "transposed_report.tsv")

final_report_tsv = os.path.join("assemblies", "merged-quast.tsv")
final_csv        = os.path.join("assemblies", "merged.quast.csv")

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
    out = [["Metric","merged"]]
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
    csv_rows = [["Metric","merged"]]
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
      echo "[ok] Wrote assemblies/merged-quast.tsv and assemblies/merged.quast.csv"
      ;;
    16)
      echo "Step 16 - Final assembly comparison"

      MERQURY_ENABLE=0
      [[ "${merqury_enable:-false}" == true ]] && MERQURY_ENABLE=1
      MERQURY_DB="${merqury_db:-}"

      if [[ "$MERQURY_ENABLE" -eq 1 ]]; then
        if [[ -z "$MERQURY_DB" ]]; then
          for cand in reads.meryl meryl/reads.meryl merqury/reads.meryl *.meryl; do
            [[ -e "$cand" ]] || continue
            if [[ -d "$cand" ]]; then
              MERQURY_DB="$cand"
              break
            fi
          done
        fi

        if command -v merqury.sh >/dev/null 2>&1 && [[ -n "$MERQURY_DB" && -d "$MERQURY_DB" ]]; then
          echo "[INFO] Running Merqury on final assembly"
          mkdir -p merqury
          log_version "merqury.sh" "merqury.sh" 2>/dev/null || true
          merqury.sh "$MERQURY_DB" assemblies/final.merged.fasta "merqury/final" || \
            echo "[warn] Merqury failed on final assembly" >&2
        else
          echo "[warn] Merqury requested for final assembly but merqury.sh or .meryl database was not found; skipping." >&2
        fi
      else
        echo "[INFO] Merqury disabled for final assembly comparison"
      fi

      mkdir -p final_results

      python3 - <<'PY'
import os, csv, re
qv_file = os.path.join("merqury", "final.qv")
comp_file = os.path.join("merqury", "final.completeness.stats")
def parse_first_float(path):
    if not os.path.exists(path):
        return ""
    txt = open(path, "r", errors="ignore").read()
    patterns = [
        r'(?i)qv[^0-9]*([0-9]+(?:\.[0-9]+)?)',
        r'(?i)completeness[^0-9]*([0-9]+(?:\.[0-9]+)?)',
        r'([0-9]+(?:\.[0-9]+)?)'
    ]
    for pat in patterns:
        m = re.search(pat, txt)
        if m:
            return m.group(1)
    return ""
qv = parse_first_float(qv_file); comp = parse_first_float(comp_file)
with open("assemblies/merged.merqury.csv", "w", newline="") as f:
    csv.writer(f).writerows([["Metric","merged"],["Merqury QV",qv],["Merqury completeness (%)",comp]])
print("Wrote assemblies/merged.merqury.csv")
PY

      python3 - <<'PY'
import os, csv
from collections import OrderedDict
TEL_WINDOW = int(os.getenv("TEL_WINDOW", "100"))
FINAL_FILES=[("telo","assemblies/merged.telo.csv"),("busco","assemblies/merged.busco.csv"),("quast","assemblies/merged.quast.csv"),("merqury","assemblies/merged.merqury.csv")]
ASM_INFO="assemblies/assembly_info.csv"
DECISION_FILE="assemblies/selection_decision.txt"
OUT="final_results/final_result.csv"
def read_metric_merged(path):
    rows=[]
    if not os.path.exists(path): return rows
    with open(path,newline="") as f:
        r=csv.reader(f); hdr=next(r,None)
        if not hdr or len(hdr)<2: return rows
        idx_metric=0; idx_merged=1
        for i,h in enumerate(hdr):
            name=(h or "").strip().lower()
            if name=="metric": idx_metric=i
            if name=="merged": idx_merged=i
        for row in r:
            if not row: continue
            m=(row[idx_metric] if idx_metric < len(row) else "").strip()
            v=(row[idx_merged] if idx_merged < len(row) else "").strip()
            if m: rows.append((m,v))
    return rows
def read_assembly_info(path):
    if not os.path.exists(path): return None,None
    with open(path,newline="") as f:
        r=csv.reader(f); hdr=next(r,None)
        if not hdr: return None,None
        hdr=list(hdr); hdr[0]="Metric"; body=[row for row in r if row]
    return hdr,body
def compute_telo_counts(list_path, tel_window=100):
    try:
        s={}; e={}
        with open(list_path,'r',newline='') as f:
            for line in f:
                line=line.strip().replace('
','')
                if not line: continue
                parts=line.split()
                if len(parts)<4: continue
                c=parts[0]
                try:
                    a=int(parts[1]); b=int(parts[2]); L=int(parts[3])
                except ValueError:
                    continue
                has_left=(a <= tel_window); has_right=((L-b) <= tel_window)
                s[c]=1 if has_left else 0; e[c]=1 if has_right else 0
        keys=set(s)|set(e); double=sum(1 for k in keys if s.get(k,0)>0 and e.get(k,0)>0); single=sum(1 for k in keys if (s.get(k,0)+e.get(k,0))==1); total=double+single
        return str(double),str(single),str(total)
    except Exception:
        return None
def count_fasta(path):
    try:
        with open(path) as f: return str(sum(1 for ln in f if ln.startswith(">")))
    except Exception: return ""
final_map=OrderedDict(); order_seen=[]
for _,p in FINAL_FILES:
    for m,v in read_metric_merged(p):
        if m not in final_map: order_seen.append(m)
        final_map[m]=v
_tcounts=compute_telo_counts("assemblies/final.telo.list", TEL_WINDOW)
if _tcounts is not None:
    _d,_s,_t=_tcounts
    final_map["Telomere double-end contigs"]=_d
    final_map["Telomere single-end contigs"]=_s
    final_map["Telomere-supported contigs"]=_t
if os.path.exists(DECISION_FILE):
    with open(DECISION_FILE) as f:
        for line in f:
            line=line.rstrip("
")
            if not line: continue
            k,*rest=line.split("	",1); v=rest[0] if rest else ""
            if k=="selected_score": final_map["Selection score"]=v
            elif k=="selected_assembler": final_map["Selected assembler"]=v
            elif k=="auto_mode": final_map["Auto-selection mode"]=v
            elif k=="score_formula": final_map["Score formula"]=v
if os.path.exists("protected_telomere_mode.txt"):
    final_map["Protected telomere mode"]=open("protected_telomere_mode.txt").read().strip()
final_map["Step10 strict T2T pool contigs"]=count_fasta("t2t_clean.fasta")
final_map["Step10 best single-end telomere pool contigs"]=count_fasta("single_tel_best_clean.fasta")
final_map["Step10 optimized telomere-supported pool contigs"]=count_fasta("telomere_supported_best_clean.fasta")
try:
    with open("assemblies/single_tel.replaced.ids") as f:
        final_map["Step12 rescued telomere replacements"]=str(sum(1 for _ in f))
except Exception:
    final_map["Step12 rescued telomere replacements"]=""
hdr,body=read_assembly_info(ASM_INFO)
if hdr is not None:
    hdr_out=list(hdr)+(["merged"] if "merged" not in [h.lower() for h in hdr] else [])
    out_rows=[]
    for row in body:
        if len(row) < len(hdr_out)-1: row = row + [""]*(len(hdr_out)-1-len(row))
        metric=(row[0] if row else "").strip()
        out_rows.append(row + [final_map.get(metric,"")])
    present_metrics=set(r[0].strip() for r in body if r and r[0].strip())
    for m in final_map.keys():
        if m not in present_metrics:
            r=[""]*len(hdr_out); r[0]=m; r[-1]=final_map.get(m,""); out_rows.append(r)
    with open(OUT,"w",newline="") as f:
        w=csv.writer(f); w.writerow(hdr_out); w.writerows(out_rows)
else:
    with open(OUT,"w",newline="") as f:
        w=csv.writer(f); w.writerow(["Metric","merged"])
        for m,v in final_map.items(): w.writerow([m,v])
print(f"Wrote {OUT}")
PY
      echo "[ok] Wrote final_results/final_result.csv"
      ;;

    17)
      echo "Step 17 - Cleanup temporary files"

      cleanup_temp() {
        echo "Cleanup: organizing temporary outputs..."
        shopt -s nullglob dotglob extglob nocaseglob

        mkdir -p temp/merge temp/merge/fasta temp/merge/param temp/busco temp/log final_results

        find . -maxdepth 1 -type f -name 'aln_summary_merged*.tsv' -exec mv -f -- {} temp/merge/ \; 2>/dev/null || true
        find . -maxdepth 1 -type f -name 'anchor_summary_merged_*.txt' -exec mv -f -- {} temp/merge/ \; 2>/dev/null || true

        find . -maxdepth 1 -type f \(           -name '.merged_*' -o -name 'merged_*.delta' -o -name 'merged_*.coords' -o -name 'merged_*.snps' -o           -name 'merged_*.delta.*' -o -name 'merged_*.crunch' -o -name 'merged_*.filter' -o           -name 'merged_*.qdiff' -o -name 'merged_*.rdiff' -o -name 'merged_*.mcoords'         \) -exec mv -f -- {} temp/merge/ \; 2>/dev/null || true

        if [[ -d busco ]]; then
          find busco -maxdepth 2 -type f -name '*.log' -exec mv -f -- {} temp/busco/ \; 2>/dev/null || true
        fi

        find . -maxdepth 1 -type f \( -name 'busco_*.log' -o -name '*busco*.log' \) -exec mv -f -- {} temp/busco/ \; 2>/dev/null || true
        find . -maxdepth 1 -type f \( -name 'merged_*.fasta' -o -name 'merged_*.fa' \) -exec mv -f -- {} temp/merge/fasta/ \; 2>/dev/null || true
        find . -maxdepth 1 -type f -name 'param_summary_merged_*.txt' -exec mv -f -- {} temp/merge/param/ \; 2>/dev/null || true

        for f in *.log; do
          [[ -f "$f" ]] || continue
          [[ "$f" == *busco* ]] && continue
          mv -f -- "$f" temp/log/
        done

        mv -f -- assemblies/final.merged.fasta final_results/ 2>/dev/null || true
        mv -f -- assemblies/selection_debug.tsv final_results/ 2>/dev/null || true
        mv -f -- assemblies/selection_decision.txt final_results/ 2>/dev/null || true
        mv -f -- assemblies/merged.merqury.csv final_results/ 2>/dev/null || true
        mv -f -- assemblies/merged.busco.csv final_results/ 2>/dev/null || true
        mv -f -- assemblies/merged.quast.csv final_results/ 2>/dev/null || true
        mv -f -- assemblies/merged.telo.csv final_results/ 2>/dev/null || true
        mv -f -- telomere_cluster_summary.tsv final_results/ 2>/dev/null || true
        mv -f -- telomere_support_summary.csv final_results/ 2>/dev/null || true
        mv -f -- protected_telomere_mode.txt final_results/ 2>/dev/null || true

        echo "[ok] Cleanup complete."
      }

      cleanup_temp
      ;;
    18)
      echo "Step 18 - Assembly-only comparison summary"

      mkdir -p assemblies merqury final_results

      MERQURY_ENABLE=0
      [[ "${merqury_enable:-false}" == true ]] && MERQURY_ENABLE=1
      MERQURY_DB="${merqury_db:-}"

      if [[ "$MERQURY_ENABLE" -eq 1 ]]; then
        if [[ -z "$MERQURY_DB" ]]; then
          for cand in reads.meryl meryl/reads.meryl merqury/reads.meryl *.meryl; do
            [[ -e "$cand" ]] || continue
            if [[ -d "$cand" ]]; then
              MERQURY_DB="$cand"
              break
            fi
          done
        fi

        if command -v merqury.sh >/dev/null 2>&1 && [[ -n "$MERQURY_DB" && -d "$MERQURY_DB" ]]; then
          echo "[INFO] Running Merqury for assembly-only comparison using database: $MERQURY_DB"
          log_version "merqury.sh" "merqury.sh" 2>/dev/null || true

          declare -A asm_paths=(
            [canu]="assemblies/canu.result.fasta"
            [external]="assemblies/external.result.fasta"
            [flye]="assemblies/flye.result.fasta"
            [ipa]="assemblies/ipa.result.fasta"
            [nextDenovo]="assemblies/nextDenovo.result.fasta"
            [peregrine]="assemblies/peregrine.result.fasta"
            [hifiasm]="assemblies/hifiasm.result.fasta"
          )

          for asm_name in canu external flye ipa nextDenovo peregrine hifiasm; do
            asm_fa_i="${asm_paths[$asm_name]}"
            [[ -s "$asm_fa_i" ]] || continue
            if [[ ! -f "merqury/${asm_name}.qv" || ! -f "merqury/${asm_name}.completeness.stats" ]]; then
              merqury.sh "$MERQURY_DB" "$asm_fa_i" "merqury/${asm_name}" || \
                echo "[warn] Merqury failed for ${asm_name}" >&2
            fi
          done
        else
          echo "[warn] Merqury requested but merqury.sh or a valid .meryl database was not found; skipping Merqury." >&2
        fi
      else
        echo "[INFO] Merqury disabled for assembly-only comparison"
      fi

      python3 - <<'PY'
import os, csv, re
assemblers = ["canu","external","flye","ipa","nextDenovo","peregrine","hifiasm"]
rows = [["Metric"] + assemblers, ["Merqury QV"], ["Merqury completeness (%)"]]
def parse_first_float(path):
    if not os.path.exists(path): return ""
    txt = open(path, "r", errors="ignore").read()
    for pat in [r'(?i)qv[^0-9]*([0-9]+(?:\.[0-9]+)?)', r'(?i)completeness[^0-9]*([0-9]+(?:\.[0-9]+)?)', r'([0-9]+(?:\.[0-9]+)?)']:
        m = re.search(pat, txt)
        if m: return m.group(1)
    return ""
for asm in assemblers:
    rows[1].append(parse_first_float(os.path.join("merqury", f"{asm}.qv")))
    rows[2].append(parse_first_float(os.path.join("merqury", f"{asm}.completeness.stats")))
with open("assemblies/assembly.merqury.csv", "w", newline="") as f:
    csv.writer(f).writerows(rows)
print("Wrote assemblies/assembly.merqury.csv")
PY
      build_assembly_info_v2
      check_command
      cp -f assemblies/assembly_info.csv final_results/assembly_only_result.csv 2>/dev/null || true
      echo "[ok] Wrote assemblies/assembly_info.csv"
      echo "[ok] Wrote final_results/assembly_only_result.csv"
      ;;
    esac
    echo "===== [$(timestamp)] STEP ${step} END ====="
  } 2>&1 | awk "$LOG_AWK" | tee -a "$STEP_LOG"
  step_status=${PIPESTATUS[0]}
  step_finished_epoch=$(date +%s)
  append_step_benchmark "$step" "$step_started_epoch" "$step_finished_epoch" "$step_status" "$STEP_LOG"
  if [[ "$step_status" -ne 0 ]]; then
    echo "[error] Step ${step} failed. See ${STEP_LOG}" >&2
    write_benchmark_summary
    exit "$step_status"
  fi
done
write_benchmark_summary
