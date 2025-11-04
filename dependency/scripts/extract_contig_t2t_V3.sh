#!/bin/bash
# extract_contig_t2t_V3
# Extract contigs that have telomere at BOTH ends from a FASTA,
# using a seqtk 'telo' result file (columns: contig start end length ...)

set -euo pipefail

# Defaults
output="output.fasta"
fasta_file=""
seqtk_file=""

usage() {
  echo "Usage: $0 -f <fasta_file> -i <seqtk_telo_result> [-o <output_fasta>]"
  echo "  -f  Input FASTA file"
  echo "  -i  seqtk telo result file (expects columns: contig start end length ...)"
  echo "  -o  Output FASTA file (default: output.fasta)"
  exit 1
}

# Parse args
while getopts "f:i:o:" opt; do
  case "$opt" in
    f) fasta_file="$OPTARG" ;;
    i) seqtk_file="$OPTARG" ;;
    o) output="$OPTARG" ;;
    *) usage ;;
  endsw
done

# Checks
[[ -z "${fasta_file}" || -z "${seqtk_file}" ]] && usage
[[ -s "${fasta_file}" ]] || { echo "[error] FASTA not found or empty: ${fasta_file}" >&2; exit 1; }
[[ -s "${seqtk_file}" ]] || { echo "[error] seqtk telo result not found or empty: ${seqtk_file}" >&2; exit 1; }
command -v seqtk >/dev/null 2>&1 || { echo "[error] seqtk not found in PATH" >&2; exit 1; }

# Prepare
: > "${output}"
ids="$(mktemp "${TMPDIR:-/tmp}/t2t_ids.XXXXXX")"
trap 'rm -f "${ids}"' EXIT

# Build ID list of contigs that have BOTH: (at least one start==0) AND (at least one end==length)
# We normalize to the FIRST TOKEN of the contig header (before any whitespace),
# so it matches seqtk's header-ID handling even if FASTA headers have descriptions.
awk '
  {
    gsub(/\r$/, "", $0);               # strip CR if present
    if (NF < 4) next;                  # need contig, start, end, length
    if ($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $4 !~ /^[0-9]+$/) next;

    contig=$1; sub(/[ \t].*$/, "", contig);
    start=$2+0; end=$3+0; len=$4+0;

    if (start==0) has_start[contig]=1;
    if (end==len) has_end[contig]=1;
  }
  END {
    for (c in has_start) if (c in has_end) print c;
  }
' "${seqtk_file}" | sort -u > "${ids}"

if [[ ! -s "${ids}" ]]; then
  echo "[warn] No contigs have telomeres at BOTH ends; ${output} will be empty."
  exit 0
fi

echo "Extracting sequences for $(wc -l < "${ids}") contig(s) with telomeres at BOTH ends..."
seqtk subseq "${fasta_file}" "${ids}" > "${output}" || {
  echo "[error] seqtk subseq failed" >&2
  exit 1
}

if [[ ! -s "${output}" ]]; then
  echo "[warn] Output FASTA is empty. Possible header/ID mismatch."
  echo "  Checks:"
  echo "    head -n1 ${fasta_file} | sed 's/^>//; s/ .*//'"
  echo "    head -n1 ${ids}"
else
  echo "Extraction complete. Output written to ${output}"
fi