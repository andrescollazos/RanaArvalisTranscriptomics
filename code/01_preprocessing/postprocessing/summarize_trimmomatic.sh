#!/bin/bash

source ../../../.env

# Output file
output="$RESULTS/01_preprocessing/02_trimming/00_trimmomatic_summary.tsv"
echo -e "Sample_id\tInput_read_pairs\tBoth_surviving\tForward_only_surviving\tReverse_only_surviving\tDropped" > "$output"

# Iterate over all trim_P32262_*.er files
for file in $RESULTS/trim_P32262_*.er; do
    # Extract sample id (remove 'trim_' and everything after first dot)
    sample_id=$(basename "$file" | sed -E 's/^trim_([^.]*)\..*/\1/')

    line=$(sed -n '6p' "$file")

    input=$(echo "$line" | grep -oP 'Input Read Pairs:\s*\K[0-9]+')
    both=$(echo "$line" | grep -oP 'Both Surviving:\s*\K[0-9]+\s*\([^)]+\)')
    forward=$(echo "$line" | grep -oP 'Forward Only Surviving:\s*\K[0-9]+\s*\([^)]+\)')
    reverse=$(echo "$line" | grep -oP 'Reverse Only Surviving:\s*\K[0-9]+\s*\([^)]+\)')
    dropped=$(echo "$line" | grep -oP 'Dropped:\s*\K[0-9]+\s*\([^)]+\)')

    echo -e "${sample_id}\t${input}\t${both}\t${forward}\t${reverse}\t${dropped}" >> "$output"
done
