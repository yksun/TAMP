#!/bin/bash

# Default values
genomesize=""
threads=20
fastq=""
steps=()
motif="AACCCT"
assembler="peregrine" # Default assembler

# Function to display usage
usage() {
  echo "Usage: $0 -g <genomesize> -t <threads> (--fastq <HiFi.fastq> | --fasta <assembly.fasta>) -m <telomere_motif> [-s <steps>] [--busco <lineage>] [--choose]"
  echo ""
  echo "Options:"
  echo "  -g    Genome size (required), e.g. 2g, 500m, 100k"
  echo "  -t    Threads (required), e.g. 32"
  echo "  --fastq    Path to PacBio HiFi FASTQ reads (required unless you pass --fasta)"
  echo "  --fasta    Path to a pre-assembled genome FASTA (optional; can be used alone or with --fastq)"
  echo "  -m    Telomere repeat motif (required), e.g. AACCCT"
  echo "  -s    Steps to run (optional). Comma/range allowed, e.g. 1,3-5 (default: all applicable steps)"
  echo "        You can specify individual steps or ranges (e.g., 1,2,3 or 1-6 or 2-4)."
  echo "  --choose  Prompt to choose an assembler for the final merge (optional, default: peregrine).
  --busco   Lineage for BUSCO to evaluate each assembly after step 7 (default: ascomycota_odb10)"
  echo ""
  echo "Steps:"
  echo "  1. Assembly of the genome using HiCanu"
  echo "  2. Assembly of the genome using NextDenovo"
  echo "  3. Assembly of the genome using Peregrine"
  echo "  4. Assembly of the genome using IPA"
  echo "  5. Assembly of the genome using Flye"
  echo "  6. Assembly of the genome using Hifiasm"
  echo "  7. Copy all assemblies (adds rows to assembly_info.csv; runs BUSCO too if --busco)"
  echo "  8. Extract telomere-containing contigs"
  echo "  9. Merge all assemblies"
  echo "  10. Run quast.py for all assembler results"
  echo "  11. Final merge using selected assembler"
  echo "  12. BUSCO analysis"
  echo "  13. Telomere analysis"
  echo ""
  echo "Example:"
  echo "  $0 -g 2g -t 16 --fastq /path/to/reads.fastq -m AACCCT -s 1,3-5 --choose"
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

busco_lineage="ascomycota_odb10"
busco_run=false

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
         fasta) fasta="${!OPTIND}"; OPTIND=$((OPTIND + 1))
         ;;
         choose) 
           echo "Please enter the assembler you want to use for the final merge (canu, nextDenovo, peregrine, ipa, flye, RAFT-hifiasm, external):"
           read assembler
           ;;
         *) echo "Unknown option --${OPTARG}" >&2; usage
         ;;
         busco)
           busco_run=true
           busco_lineage="${!OPTIND}"; OPTIND=$((OPTIND + 1))
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
if [[ -z "$genomesize" || -z "$threads" || ( -z "$fastq" && -z "$fasta" ) || -z "$motif" ]]; then
  usage
fi

if [[ -n "$fastq" ]]; then
  project=$(basename "$fastq" .fastq)
else
  # derive from fasta name, stripping common extensions
  base=$(basename "$fasta")
  project=${base%%.*}
fi

# Modify run.cfg with genomesize and threads
if [[ -n "$fastq" ]]; then
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
fi

# Function to check the success of the previous command
check_command() {
  if [ $? -ne 0 ]; then
    echo "Error: Command failed. Exiting."
    exit 1
  fi
}

echo "Activating assembly environment..."
eval "$(conda shell.bash hook)"
conda activate pacbiohifi

# Collect assembler versions (best-effort)
declare -A ASM_VER
get_ver() { cmd="$1"; eval "$cmd" 2>/dev/null | head -n1 | tr -d '\r' || true; }
ASM_VER["canu"]="$(get_ver 'canu --version')"
ASM_VER["nextDenovo"]="$(get_ver 'nextDenovo --version || nextDenovo -V')"
ASM_VER["peregrine"]="$(get_ver 'pg_asm --version')"
ASM_VER["ipa"]="$(get_ver 'ipa --version')"
ASM_VER["flye"]="$(get_ver 'flye --version')"
ASM_VER["RAFT-hifiasm"]="$(get_ver 'hifiasm --version'); RAFT_VER=$(get_ver '"'"'~/opt/RAFT/raft --version'"'"'); [[ -n "$RAFT_VER" ]] && ASM_VER["RAFT-hifiasm"]="${ASM_VER[RAFT-hifiasm]} / RAFT: $RAFT_VER"'
ASM_VER["external"]="(user-provided FASTA)"
# Ensure CSV has a Version column
csv_init
if ! head -n1 "$assembly_csv" | grep -q ",Version"; then
  tmp="${assembly_csv}.tmp"
  head -n1 "$assembly_csv" | awk '{print $0",Version"}' > "$tmp"
  tail -n +2 "$assembly_csv" | awk '{print $0","}' >> "$tmp"
  mv "$tmp" "$assembly_csv"
fi

check_command

# CSV helpers (simple append-or-replace)

# ---------- Pretty print & Markdown export helpers ----------
# csv_pretty_and_markdown preferred_cols...
# Dynamically appends all Quast_* metrics not explicitly listed at the end.
csv_pretty_and_markdown() {
  csv_init
  local md_out="assembly_info.md"
  local pref_list="$*"

  awk -v PREF="$pref_list" -v MDOUT="$md_out" -F, '
    BEGIN{
      OFS=",";
      # Build set of preferred columns
      nPref=split(PREF, P, /[ \t]+/);
      for(i=1;i<=nPref;i++){ if(P[i]!="") pref[P[i]]=1; }
    }
    NR==1{
      # Map header names to indices
      N=split($0, H);
      for(i=1;i<=N;i++){ idx[H[i]]=i; }
      # Build final order: preferred (if present), then all other Quast_* not already included, then remaining not included
      for(i=1;i<=nPref;i++){ if(P[i] != "" && idx[P[i]]>0){ order[++M]=P[i]; seen[P[i]]=1; } }
      for(i=1;i<=N;i++){
        if (H[i] ~ /^Quast_/ && !seen[H[i]]) { order[++M]=H[i]; seen[H[i]]=1; }
      }
      for(i=1;i<=N;i++){ if(!seen[H[i]]){ order[++M]=H[i]; seen[H[i]]=1; } }
      # Print header in selected order to tmp CSV (stdout) and prepare Markdown
      for(i=1;i<=M;i++){
        printf("%s%s", order[i], (i<M?OFS:RS));
      }
      # Prepare markdown header row
      mdHeader="|"; mdSep="|";
      for(i=1;i<=M;i++){ mdHeader=mdHeader""order[i]"|"; mdSep=mdSep"---|"; }
      mdRows=mdHeader"\n"mdSep"\n";
      next;
    }
    NR>1{
      # Emit re-ordered CSV row
      for(i=1;i<=M;i++){
        j=idx[order[i]]; val=(j>0? $(j) : "");
        printf("%s%s", val, (i<M?OFS:RS));
      }
      # Build Markdown row (escape pipes minimally)
      mdRow="|";
      for(i=1;i<=M;i++){
        j=idx[order[i]]; val=(j>0? $(j) : "");
        gsub(/\|/,"\\|",val);
        mdRow=mdRow""val"|";
      }
      mdRows=mdRows""mdRow"\n";
    }
    END{
      # Write markdown file
      if (MDOUT!="") {
        f=MDOUT; print mdRows > f;
      }
    }
  ' "$assembly_csv" | column -s, -t || cat

  echo "Markdown summary written to ${md_out}"
}

assembly_csv="assembly_info.csv"
csv_init() {
  if [[ ! -f "$assembly_csv" ]]; then
    echo "Name,Lineage,Timestamp,C_count,C_pct,S_count,S_pct,D_count,D_pct,F_count,F_pct,M_count,M_pct,Total_BUSCOs,DoubleEndTelo,SingleEndTelo,Quast_TotalLen,Quast_Contigs,Quast_N50,Quast_ReportAdded" > "$assembly_csv"
  fi
}
# replace_row <Name> <csv_line_without_newline>
replace_row() {
  local name="$1"; shift
  local newline="$*"
  csv_init
  if grep -q "^${name}," "$assembly_csv"; then
    tmp="${assembly_csv}.tmp"
    awk -F, -v OFS=, -v target="$name" -v new="$newline" 'NR==1{print; next} $1==target{print new; next} $1!=target{print}' "$assembly_csv" > "$tmp" && mv "$tmp" "$assembly_csv"
  else
    echo "$newline" >> "$assembly_csv"
  fi
}
# merge_row <Name> key=value ... (only sets keys provided, keeps other existing values)
merge_row() {
  local name="$1"; shift
  csv_init
  # Read header into array
  IFS=, read -r -a H < <(head -n1 "$assembly_csv")
  # Build associative array of current values
  local current=""
  if grep -q "^${name}," "$assembly_csv"; then
    line=$(grep -m1 "^${name}," "$assembly_csv")
    IFS=, read -r -a V <<< "$line"
  else
    V=()
    for ((i=0;i<${#H[@]};i++)); do V[i]=""; done
    V[0]="$name"
  fi
  # Apply provided key=val pairs
  for kv in "$@"; do
    key="${kv%%=*}"; val="${kv#*=}"
    # find column
    for ((i=0;i<${#H[@]};i++)); do
      if [[ "${H[i]}" == "$key" ]]; then
        V[i]="$val"
        break
      fi
    done
  done
  # Reconstruct row
  out="${V[0]}"; for ((i=1;i<${#H[@]};i++)); do out+=",""${V[i]}"; done
  replace_row "$name" "$out"
}


# If no specific steps provided, default to running all steps (1-13)
if [ ${#steps[@]} -eq 0 ]; then
  if [[ -n "$fastq" ]]; then
    steps=(1 2 3 4 5 6 7 8 9 10 11 12 13)
  else
    steps=(7 8 9 10 11 12 13)
  fi
fi

# Execute the specified steps
for step in "${steps[@]}"; do
  case $step in
    1)
      echo "Step 1 - Assembly of the genome using HiCanu"
      canu -p canu -d hicanu genomeSize=$genomesize maxThreads=$threads -pacbio-hifi $fastq
      check_command
      ;;
    2)
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
      echo "Step 4 - Assembly of the genome using IPA"
      ipa local --nthreads $threads --njobs 1 --run-dir ipa -i $fastq
      check_command
      ;;
    5)
      echo "Step 5 - Assembly of the genome using Flye"
      flye --pacbio-hifi $fastq --out-dir flye --threads $threads
      check_command
      ;;
    6)
      echo "Step 6 - Assembly of the genome using Hifiasm"
      mkdir hifiasm
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
      ;;
    7)
      if [[ "$busco_run" == true ]]; then echo "Step 7 - Copy all assemblies + BUSCO"; else echo "Step 7 - Copy all assemblies"; fi
      if [[ -n "$fasta" ]]; then
        cp "$fasta" ./external.result.fasta
      fi
      [[ -d ./hicanu ]] && cp ./hicanu/canu.contigs.fasta ./canu.result.fasta
      [[ -d ./NextDenovo/03.ctg_graph ]] && cp ./NextDenovo/03.ctg_graph/nd.asm.fasta ./nextDenovo.result.fasta
      [[ -d ./peregrine-2021 ]] && cp ./peregrine-2021/asm_ctgs_m_p.fa ./peregrine.result.fasta
      [[ -d ./ipa/assembly-results ]] && cp ./ipa/assembly-results/final.p_ctg.fasta ./ipa.result.fasta
      [[ -d ./flye ]] && cp ./flye/assembly.fasta ./flye.result.fasta
      [[ -d ./hifiasm ]] && cp ./hifiasm/RAFT-hifiasm.fasta ./RAFT-hifiasm.result.fasta

      # Initialize/mark assemblies in CSV with timestamp
      csv_init
      run_ts="$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
      for f in *.result.fasta; do
        base="${f%.result.fasta}"
        ver="${ASM_VER[$base]}"; merge_row "$base" Timestamp="$run_ts" Version="${ver}"
      done
      ;;

      # Optional BUSCO run after copying assemblies
      if [[ "$busco_run" == true ]]; then
        echo "Running BUSCO on assemblies with lineage: ${busco_lineage}"
        eval "$(conda shell.bash hook)"
        conda deactivate
        conda activate busco
        # Initialize CSV
        busco_csv="busco_summary.csv"
        echo "Name,Lineage,Timestamp,C_count,C_pct,S_count,S_pct,D_count,D_pct,F_count,F_pct,M_count,M_pct,Total_BUSCOs,DoubleEndTelo,SingleEndTelo" > "$busco_csv"
        ts="$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
        shopt -s nullglob
        for f in *.result.fasta; do
          base="${f%.result.fasta}"
          outdir="busco_${base}"
          # Run BUSCO
          busco -i "$f" -l "$busco_lineage" -m genome -o "$outdir" --cpu "$threads"
          # Find short summary file
          summary_path="$(find \"$outdir\" -maxdepth 3 -type f -name 'short_summary*.txt' | head -n1)"
          if [[ -z "$summary_path" ]]; then summary_path="$(find \"$outdir\" -maxdepth 5 -type f -name 'short_summary*.txt' | head -n1)"; fi
          Cc=""; Sc=""; Dc=""; Fc=""; Mc=""; Tot=""; Cpct=""; Spct=""; Dpct=""; Fpct=""; Mpct=""
          if [[ -f "$summary_path" ]]; then
            Cc=$(grep -F "Complete BUSCOs (C)" "$summary_path" | awk '{print $1}')
            Sc=$(grep -F "Complete and single-copy BUSCOs (S)" "$summary_path" | awk '{print $1}')
            Dc=$(grep -F "Complete and duplicated BUSCOs (D)" "$summary_path" | awk '{print $1}')
            Fc=$(grep -F "Fragmented BUSCOs (F)" "$summary_path" | awk '{print $1}')
            Mc=$(grep -F "Missing BUSCOs (M)" "$summary_path" | awk '{print $1}')
            Tot=$(grep -F "Total BUSCO groups searched" "$summary_path" | awk '{print $1}')
            line=$(grep -m1 "C:" "$summary_path")
            Cpct=$(echo "$line" | sed -n 's/.*C:\([0-9\.]+%\).*/\1/p')
            Spct=$(echo "$line" | sed -n 's/.*S:\([0-9\.]+%\).*/\1/p' | cut -d',' -f1)
            Dpct=$(echo "$line" | sed -n 's/.*D:\([0-9\.]+%\).*/\1/p')
            Fpct=$(echo "$line" | sed -n 's/.*F:\([0-9\.]+%\).*/\1/p')
            Mpct=$(echo "$line" | sed -n 's/.*M:\([0-9\.]+%\).*/\1/p')
          fi
          # Telomere stats
          telo_list="${base}.telo.list"
          if [[ ! -f "$telo_list" && -n "$motif" ]]; then
            seqtk telo -s 1 -m "$motif" "$f" > "$telo_list" || true
          fi
          double=0; single=0
          if [[ -f "$telo_list" ]]; then
            awk '{count[$1]++} END{for (k in count){if(count[k]==2) d++; else if(count[k]==1) s++}; print d","s}' "$telo_list" |             awk -F',' '{print $1,$2}' | { read d s; double=${d:-0}; single=${s:-0}; echo >/dev/null; }
          fi
          merge_row "$base" Lineage="${busco_lineage}" Timestamp="${ts}" C_count="${Cc}" C_pct="${Cpct}" S_count="${Sc}" S_pct="${Spct}" D_count="${Dc}" D_pct="${Dpct}" F_count="${Fc}" F_pct="${Fpct}" M_count="${Mc}" M_pct="${Mpct}" Total_BUSCOs="${Tot}"
        done
        conda deactivate
      fi

    8)
      echo "Step 8 - Extract telomere-containing contigs"
      for fasta in *.result.fasta; do
        seqtk telo -s 1 -m "$motif" $fasta > "${fasta%.result.fasta}.telo.list"
        bash ./extract_contig_T_V3.sh -i "${fasta%.result.fasta}.telo.list" -f $fasta -o "${fasta%.result.fasta}.telo.fasta"
        check_command
      done
      ;;
    9)
      echo "Step 9 - Merge all telo assemblies"
      fasta_files=(*.telo.fasta)
    for file1 in "${fasta_files[@]}"; do
        for file2 in "${fasta_files[@]}"; do
          base1=$(basename "$file1" .telo.fasta)
          base2=$(basename "$file2" .telo.fasta)
          merge_wrapper.py -l 1000000 "$file1" "$file2" --prefix merged_"$base1"_"$base2"
          check_command
        done
      done
      cat merged_*.fasta > allmerged_telo.fasta
      eval "$(conda shell.bash hook)"
      conda deactivate
      conda activate funannotate
      funannotate sort -i allmerged_telo.fasta -b contig -o allmerged_telo_sort.fasta --minlen 500
      seqtk telo -s 1 -m "$motif" allmerged_telo_sort.fasta > allmerged.telo.list
      bash ./t2t_list.sh -i allmerged.telo.list -o t2t.list
      ~/opt/scripts/faSomeRecords allmerged_telo_sort.fasta t2t.list t2t.fasta
      funannotate clean -i  t2t.fasta -p 30 -o  t2t_clean.fasta --exhaustive
      check_command
      ;;
    10)
      echo "Step 10 - Run quast.py for all assembler results"
      echo "Activating assembly environment..."
      eval "$(conda shell.bash hook)"
      conda activate pacbiohifi
      quast.py *.result.fasta --threads $threads
      check_command
      
      # Update CSV with QUAST results (all metrics)
      report_tsv="quast_results/latest/report.tsv"
      report_txt="quast_results/latest/report.txt"
      if [[ -f "$report_tsv" ]]; then
        header=$(head -n1 "$report_tsv")
        IFS=$'	' read -r -a H <<< "$header"
        # Build list of Quast_ columns (skip "Assembly")
        declare -a NEWCOLS=()
        for i in "${!H[@]}"; do
          col="${H[$i]}"
          [[ "$col" == "Assembly" ]] && continue
          sc=$(echo "Quast_${col}" | sed -e 's/[ %()\/-]/_/g' -e 's/__/_/g')
          NEWCOLS+=("$sc")
        done
        # Ensure all columns exist
        if ((${#NEWCOLS[@]}>0)); then
          header0=$(head -n1 "$assembly_csv")
          need=0; addcols=""
          for col in "${NEWCOLS[@]}"; do
            if ! echo "$header0" | grep -q ",${col}\(,\|$\)"; then
              need=1; addcols+=",${col}"
            fi
          done
          if [[ $need -eq 1 ]]; then
            tmp="${assembly_csv}.tmp"
            echo "${header0}${addcols}" > "$tmp"
            tail -n +2 "$assembly_csv" | while IFS= read -r line; do
              pads=$(echo "$addcols" | sed 's/[^,]//g' | wc -c); pads=$((pads-1))
              # pads equals number of commas in addcols => number of new fields
              empties=""; for ((k=0;k<pads;k++)); do empties+=",\"\""; done
              echo "${line}${empties}"
            done >> "$tmp"
            mv "$tmp" "$assembly_csv"
          fi
        fi
        # Map column name to index
        declare -A IDX
        for i in "${!H[@]}"; do IDX["${H[$i]}"]=$i; done
        tail -n +2 "$report_tsv" | while IFS=$'	' read -r -a R; do
          asm="${R[${IDX[Assembly]}]}"
          base="${asm%.result.fasta}"
          # Build merge_row args dynamically
          args=()
          for key in "${!IDX[@]}"; do
            [[ "$key" == "Assembly" ]] && continue
            sc=$(echo "Quast_${key}" | sed -e 's/[ %()\/-]/_/g' -e 's/__/_/g')
            val="${R[${IDX[$key]}]}"
            args+=("${sc}=${val}")
          done
          merge_row "$base" "${args[@]}"
        done
      elif [[ -f "$report_txt" ]]; then
        # Fallback parser: assemble key-value pairs from report.txt lines "key	value	..."
        # First, ensure columns exist by scanning keys
        keys=$(awk -F'	' 'NR>1 && $1!=""{print $1}' "$report_txt" | sed 's/[ %()\/-]/_/g' | sort -u)
        if [[ -n "$keys" ]]; then
          header0=$(head -n1 "$assembly_csv")
          addcols=""; need=0
          while IFS= read -r k; do
            col="Quast_${k}"
            if ! echo "$header0" | grep -q ",${col}\(,\|$\)"; then
              need=1; addcols+=",${col}"
            fi
          done <<< "$keys"
          if [[ $need -eq 1 ]]; then
            tmp="${assembly_csv}.tmp"
            echo "${header0}${addcols}" > "$tmp"
            tail -n +2 "$assembly_csv" | while IFS= read -r line; do
              pads=$(echo "$addcols" | sed 's/[^,]//g' | wc -c); pads=$((pads-1))
              empties=""; for ((k=0;k<pads;k++)); do empties+=",\"\""; done
              echo "${line}${empties}"
            done >> "$tmp"
            mv "$tmp" "$assembly_csv"
          fi
        fi
        # Assign values per assembly: report.txt groups per assembly; detect lines that begin with assembly name or include it
        # We’ll map block-wise by "Assembly" column in the header line above the table
        # Parse using quast's simple format: "Assembly	<asm>
Metric	Value"
        current=""
        while IFS=$'	' read -r k v rest; do
          if [[ "$k" == "Assembly" ]]; then
            current="$v"
            base="${current%.result.fasta}"
            continue
          fi
          [[ -z "$current" || -z "$k" ]] && continue
          scol=$(echo "Quast_${k}" | sed -e 's/[ %()\/-]/_/g' -e 's/__/_/g')
          merge_row "$base" "${scol}=${v}"
        done < "$report_txt"
      else
        # No QUAST report found at all
        :
      fi

      ;;
    11)
      # Step 11 - Final merge using selected assembler and redundancy reduction using redundans
      echo "Step 11 - Final merge"

      # Preview QUAST and BUSCO results before choosing assembler

      if [[ -f "assembly_info.csv" ]]; then
        echo "================ Assembly summary (assembly_info.csv) =========================="
        csv_pretty_and_markdown Name Version Lineage Timestamp C_count C_pct S_count S_pct D_count D_pct F_count F_pct M_count M_pct Total_BUSCOs DoubleEndTelo SingleEndTelo Quast_TotalLen Quast_Contigs Quast_N50
        echo "==============================================================================="
      fi

      if [[ -f "quast_results/latest/report.txt" ]]; then
        echo "================ QUAST report (quast_results/latest/report.txt) ================"
        sed -n '1,200p' quast_results/latest/report.txt
        echo "==============================================================================="
      else
        echo "QUAST report not found at quast_results/latest/report.txt (run step 10 to generate)."
      fi
      if [[ -f "busco_summary.csv" ]]; then
        echo "================ BUSCO summary (busco_summary.csv) ============================="
        column -s, -t busco_summary.csv | sed -n '1,200p'
        echo "==============================================================================="
      else
        echo "BUSCO summary CSV not found. (Use --busco to generate after step 7.)"
      fi

      # Map external choice to external.result.fasta if provided
      if [[ "$assembler" == "external" && ! -f external.result.fasta ]]; then
        if [[ -n "$fasta" ]]; then
          cp "$fasta" external.result.fasta
        else
          echo "Error: external assembly requested but no --fasta provided."
          exit 1
        fi
      fi
      # Add assembler choice logic
      if [[ "$1" == "--choose" ]]; then
          echo "Please enter the assembler you want to use for the final merge (canu, nextDenovo, peregrine, ipa, flye, RAFT-hifiasm, external):"
          read assembler
      if [[ ! "$assembler" =~ ^(canu|nextDenovo|peregrine|ipa|flye|RAFT-hifiasm|external)$ ]]; then
        echo "Invalid assembler selected. Exiting."
        exit 1
      fi
      merge_wrapper.py -l 1000000 ${assembler}.result.fasta t2t_clean.fasta --prefix "$assembler"
      check_command
      eval "$(conda shell.bash hook)"
      conda deactivate
      conda activate redundans
      redundans.py --noscaffolding --nogapclosing -t $threads -f merged_${assembler}.fasta --identity 0.50 --overlap 0.80 --log redundans.log
      eval "$(conda shell.bash hook)"
      conda deactivate
      conda activate redundans
      check_command
      funannotate sort -i redundans/scaffolds.reduced.fa -b contig -o merged_${assembler}_sort.fa --minlen 500
      check_command
      ;;
    12)
      # BUSCO analysis
      echo "Step 12 - BUSCO analysis"
      eval "$(conda shell.bash hook)"
      conda deactivate
      conda activate busco
      busco -o merged -i merged_${assembler}_sort.fa -l ascomycota_odb10 --cpu $threads -m genome
      check_command
      ;;
    13)
      # Telomere analysis
       echo "Step 13 - Telomere analysis"
      eval "$(conda shell.bash hook)"
      conda deactivate
      conda activate quartet
      ~/opt/telomere_analysis.sh merged_${assembler}_sort.fa
      ;;
    *)
      echo "Invalid step: $step"
      exit 1
      ;;
  esac
done