#!/bin/bash

set -euo pipefail

# === Usage ===
usage() {
  cat <<EOF
Usage: $(basename "$0") -i INPUT_DIR -c COLUMN_DIR -o OUTPUT_DIR [-s SCRIPT_PATH]

Applies column selections (produced by trimAl with -colnumbering) to a batch of
aligned FASTA files, writing trimmed alignments to the output directory.

Required:
  -i INPUT_DIR    Directory of input aligned FASTA files
  -c COLUMN_DIR   Directory of column index files (one per input, stem-matched)
  -o OUTPUT_DIR   Directory to write trimmed FASTA files (created if missing)

Optional:
  -s SCRIPT_PATH  Path to selectColsFromFasta.py (default: ./src/selectColsFromFasta.py)
  -h              Show this help and exit

Notes:
  - For an input file like INPUT_DIR/1234a.fa, this script looks for COLUMN_DIR/1234a
    and writes OUTPUT_DIR/1234a.fa
EOF
}

# === Parse arguments ===
CURRENT_PATH=$(pwd)
INPUT_DIR=""
COLUMN_DIR=""
OUTPUT_DIR=""
SCRIPT_PATH="${CURRENT_PATH}/src/selectColsFromFasta.py"

while getopts ":i:c:o:s:h" opt; do
  case $opt in
    i) INPUT_DIR="$OPTARG" ;;
    c) COLUMN_DIR="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    s) SCRIPT_PATH="$OPTARG" ;;
    h) usage; exit 0 ;;
    :) echo "ERROR: Option -$OPTARG requires an argument" >&2; usage; exit 2 ;;
    \?) echo "ERROR: Invalid option -$OPTARG" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "$INPUT_DIR" || -z "$COLUMN_DIR" || -z "$OUTPUT_DIR" ]]; then
  echo "ERROR: -i, -c, and -o are required" >&2
  usage
  exit 2
fi

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "ERROR: INPUT_DIR not found: $INPUT_DIR" >&2
  exit 2
fi
if [[ ! -d "$COLUMN_DIR" ]]; then
  echo "ERROR: COLUMN_DIR not found: $COLUMN_DIR" >&2
  exit 2
fi
mkdir -p "$OUTPUT_DIR"

if [[ ! -f "$SCRIPT_PATH" ]]; then
  echo "ERROR: SCRIPT_PATH not found: $SCRIPT_PATH" >&2
  exit 2
fi

# Get list of input files (names only)
mapfile -t FILES < <(cd "$INPUT_DIR" && ls -1)

# === Trimming function ===
run_trimal () {
  local fasta_file=$1
  echo "================================================"
  echo " Trimming file: $fasta_file"
  echo "================================================"

  # Extract stem (e.g., 1234a.fa -> 1234a)
  local stem
  stem=$(echo "$fasta_file" | sed 's/\.[a-zA-Z0-9]*$//')
  local col_file="${COLUMN_DIR}/${stem}"

  local input_path="${INPUT_DIR}/${fasta_file}"
  local output_path="${OUTPUT_DIR}/${fasta_file}"

  # Skip non-regular files
  if [[ ! -f "$input_path" ]]; then
    echo "Skipping non-file: $input_path" >&2
    return 0
  fi

  # Run trimming script
  if [[ -f "$col_file" ]]; then
    python "$SCRIPT_PATH" "$input_path" "$col_file" "$output_path"
  else
    echo "WARNING: Column file not found: $col_file" >&2
  fi
}

# === Main loop ===
for file in "${FILES[@]}"; do
  run_trimal "$file"
done
