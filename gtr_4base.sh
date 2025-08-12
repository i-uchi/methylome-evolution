#!/bin/bash

# Set paths
current_path=`pwd`
alignment_path="$current_path/orig_align"
output_dir="gtr8_out"
files=($(ls "$alignment_path"))

run_raxml () {
    echo "=========================================================================="
    echo "Analysis of $1"
    echo "=========================================================================="

    name=$(echo "$1" | sed -e 's/[^0-9_]//g')
    mkdir -p "$output_dir/$name"
    cd "$output_dir/$name" || exit 1

    raxml-ng --msa "$alignment_path/$1" \
             --model GTR+G+I \
             --prefix 4b \
             --threads 1 \
             --msa-format FASTA
}

# Run in parallel (batch of 100)
for ((i = 0; i < ${#files[@]}; i++)); do
    run_raxml "${files[i]}" &
    if (( i % 100 == 99 )); then
        wait
    fi
done

# Final wait to catch any remaining background jobs
wait
