#!/bin/bash
set -euo pipefail

usage() {
  echo "Usage: $0 -f <input_fasta> -i <input_txt> -o <output_fasta>"
  echo "  -f  Input FASTA file"
  echo "  -i  Input text file with contig details (e.g., from 'seqtk telo')"
  echo "  -o  Output FASTA file with extracted contigs"
  exit 1
}

fasta=""
input_txt=""
output_fasta=""

while getopts ":f:i:o:" opt; do
  case ${opt} in
    f ) fasta=$OPTARG ;;
    i ) input_txt=$OPTARG ;;
    o ) output_fasta=$OPTARG ;;
    * ) usage ;;
  esac
done

[[ -z "${fasta}" || -z "${input_txt}" || -z "${output_fasta}" ]] && usage

if [[ ! -s "${fasta}" ]]; then
  echo "[error] FASTA not found or empty: ${fasta}" >&2
  exit 1
fi
if [[ ! -s "${input_txt}" ]]; then
  echo "[error] Input list not found or empty: ${input_txt}" >&2
  exit 1
fi
if ! command -v seqtk >/dev/null 2>&1; then
  echo "[error] seqtk not found in PATH" >&2
  exit 1
fi

# Create a temp file for IDs
ids_file="$(mktemp "${TMPDIR:-/tmp}/telo_ids.XXXXXX")"
trap 'rm -f "${ids_file}"' EXIT

# The input is typically: contig start end motif strand ...
# Some pipelines want only contigs where telomere is at the very start or very end.
# We'll detect if columns 2,3,4 are numeric; if yes, apply (start==0 || end==length) filter.
# Otherwise, just use column 1.
awk '
  BEGIN { OFS="\n" }
  {
    # Normalize DOS newlines and collapse multiple spaces/tabs
    gsub(/\r$/, "", $0)
    # Require at least one field
    if (NF < 1) next

    # Check if we can interpret start/end/length as integers (columns 2,3,4)
    has_numeric = (NF >= 4 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $4 ~ /^[0-9]+$/)

    contig = $1
    # Trim contig to the first token (up to whitespace)
    sub(/[ \t].*$/, "", contig)

    if (has_numeric) {
      start=$2; end=$3; len=$4
      if (start == 0 || end == len) {
        print contig
      }
    } else {
      print contig
    }
  }
' "${input_txt}" \
| awk 'length>0' \
| sort -u > "${ids_file}"

if [[ ! -s "${ids_file}" ]]; then
  echo "[warn] No contig IDs found to extract (check ${input_txt})" >&2
  # Create an empty output to be explicit
  : > "${output_fasta}"
  exit 0
fi

echo "Extracting sequences..."
# seqtk matches the first token of the header (before any whitespace), which is what we prepared.
seqtk subseq "${fasta}" "${ids_file}" > "${output_fasta}" || {
  echo "[error] seqtk subseq failed" >&2
  exit 1
}

if [[ ! -s "${output_fasta}" ]]; then
  echo "[warn] Output FASTA is empty. Possible header/ID mismatch between ${fasta} and ${input_txt}." >&2
  echo "       Inspect a header:   head -n1 ${fasta} | sed 's/^>//; s/ .*//'"
  echo "       Inspect an ID line: head -n1 ${ids_file}"
else
  echo "Sequences written to ${output_fasta}"
fi