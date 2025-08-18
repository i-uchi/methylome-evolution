**Purpose**
Trim poorly aligned columns on the 4-base alignment, then apply the same column selection to the extended-base (8 base) alignment to keep them in sync.

**Requirements**
- **`trimal`**: https://github.com/scapella/trimal
- **Python 3**: used by `src/selectColsFromFasta.py`

**Batch Script Usage**
```
./trimal_batch.sh -i INPUT_DIR -c COLUMN_DIR -o OUTPUT_DIR [-s SCRIPT_PATH]
```
- `-i`: directory containing input aligned FASTA files (e.g., `ali10b`)
- `-c`: directory containing kept-column lists from trimAl (one file per input; stem must match)
- `-o`: output directory to write trimmed FASTA files
- `-s`: optional path to `src/selectColsFromFasta.py` (defaults to `./src/selectColsFromFasta.py`)


**Step 1 — Produce column lists with trimAl (on 4-base alignments)**
For each 4-base alignment run trimAl with your preferred criteria.


Notes:
- `-gt`, `-st`, `-cons`, and other trimAl options can be used depending on stringency.
- `-colnumbering` outputs 1-based indices, typically whitespace-separated.


**Step 2 — Apply the same trimming to extended-base alignments**
- The batch script reads `COLUMN_DIR/<gene_id>` files and applies those columns to the corresponding aligned FASTA in `INPUT_DIR`.
  - The filename stem must match between `INPUT_DIR/<gene_id>.*` and `COLUMN_DIR/<gene_id>`.

**Tips**
- Ensure your 4-base and extended-base alignments represent the same ortholog groups and have matching sequence ordering/length after conversion.
- Consider testing with a single gene to validate your trimAl parameters before batch runs.
