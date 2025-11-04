#!/bin/bash

# Default output file
output="t2t.list"

# Parse command-line options
while getopts "i:o:" opt; do
    case "$opt" in
        i) seqtk_file="$OPTARG" ;;  # Seqtk telo result file
        o) output="$OPTARG" ;;      # Output file for the list
        *) echo "Usage: $0 -i <seqtk_telo_result> [-o <output_list>]"; exit 1 ;;
    esac
done

# Check for the required input argument
if [ -z "$seqtk_file" ]; then
    echo "Seqtk telo result file (-i) is required."
    exit 1
fi

# Initialize the output file
> "$output"  # Clear any existing content in the output file

# Variables to track the previous line
prev_contig=""
prev_start=-1
prev_end=-1
prev_length=-1

# Read the seqtk telo file line by line
while read -r line; do
    # Parse fields from the current line
    contig=$(echo "$line" | awk '{print $1}')
    start=$(echo "$line" | awk '{print $2}')
    end=$(echo "$line" | awk '{print $3}')
    length=$(echo "$line" | awk '{print $4}')

    # Skip invalid or incomplete lines
    if [ -z "$start" ] || [ -z "$end" ] || [ -z "$length" ] || ! [[ "$start" =~ ^[0-9]+$ ]] || ! [[ "$end" =~ ^[0-9]+$ ]] || ! [[ "$length" =~ ^[0-9]+$ ]]; then
        continue
    fi

    # Check conditions
    if [ "$start" -eq 0 ]; then
        # Start a new potential match
        prev_contig="$contig"
        prev_start="$start"
        prev_end="$end"
        prev_length="$length"
    elif [ "$contig" == "$prev_contig" ] && [ "$end" -eq "$length" ] && [ "$prev_start" -eq 0 ]; then
        # Match found: write the contig to the output list
        echo "$contig" >> "$output"
        # Reset the tracking variables for a fresh match
        prev_contig=""
        prev_start=-1
        prev_end=-1
        prev_length=-1
    else
        # Reset the variables when conditions are not met
        prev_contig=""
        prev_start=-1
        prev_end=-1
        prev_length=-1
    fi
done < "$seqtk_file"

echo "List of matching contigs written to $output."
