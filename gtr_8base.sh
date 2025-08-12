#!/bin/bash

# Define input/output paths (genericized)
current_path=`pwd`
alignment_path="$current_path/ali10b"
output_dir="gtr8_out"
files=($(ls "$alignment_path"))

run_raxml_evaluate () {
    echo "=========================================================================="
    echo "Analysis of $1"
    echo "=========================================================================="

    name=$(echo "$1" | sed -e 's/[^0-9_]//g')
    mkdir -p "$output_dir/$name"
    cd "$output_dir/$name" || exit 1

    raxml-ng --evaluate \
             --msa "$alignment_path/$1" \
             --model MULTI8_GTR+G+I+M{ACEFGIJT}{N-} \
             --prefix 8b \
             --tree 4b.raxml.bestTree \
             --threads 1 \
             --msa-format FASTA
}

# Run in parallel (batch of 100)
for ((i = 0; i < ${#files[@]}; i++)); do
    run_raxml_evaluate "${files[i]}" &
    if (( i % 100 == 99 )); then
        wait
    fi
done

# Final wait to ensure all jobs complete
wait
