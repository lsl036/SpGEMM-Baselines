#!/bin/bash
# Run HASpGEMM on all datasets in dataset.txt, record each "HASpGEMM,..." result line.

DATA_PATH="${DATA_PATH:-/data/suitesparse_collection}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATASET_FILE="${SCRIPT_DIR}/dataset.txt"
# DATASET_FILE="${SCRIPT_DIR}/20data.txt"
RESULT_FILE="${SCRIPT_DIR}/haspgemm_results0402.txt"
BIN="${SCRIPT_DIR}/haspgemm12"

if [ ! -f "$DATASET_FILE" ]; then
  echo "Error: dataset file not found: $DATASET_FILE"
  exit 1
fi

if [ ! -x "$BIN" ]; then
  echo "Error: haspgemm12 not found or not executable: $BIN"
  exit 1
fi

# Start result file with CSV header (match the printf order in main.c)
echo "HASpGEMM,mat_name,nnzA,mA,avg_nnzCub,avg_gflops,max_gflops,avg_time,analyze_time" > "$RESULT_FILE"

total=$(wc -l < "$DATASET_FILE")
count=0

while IFS= read -r name || [ -n "$name" ]; do
  name=$(echo "$name" | tr -d '\r')
  [ -z "$name" ] && continue

  count=$((count + 1))
  mtx_path="${DATA_PATH}/${name}/${name}.mtx"

  if [ ! -f "$mtx_path" ]; then
    echo "[$count/$total] Skip (file not found): $name"
    echo "HASpGEMM,${name},SKIP,,,,,," >> "$RESULT_FILE"
    continue
  fi

  echo "[$count/$total] Running: $name"
  output=$("$BIN" "$mtx_path" "$mtx_path" 2>&1)
  exitcode=$?

  if [ $exitcode -ne 0 ]; then
    echo "  -> run failed (exit $exitcode)"
    echo "HASpGEMM,${name},ERROR,,,,,," >> "$RESULT_FILE"
    continue
  fi

  line=$(echo "$output" | grep "^HASpGEMM," | tail -n 1)
  if [ -n "$line" ]; then
    echo "$line" >> "$RESULT_FILE"
    echo "  -> recorded"
  else
    echo "  -> no HASpGEMM line in output"
    echo "HASpGEMM,${name},NO_OUTPUT,,,,,," >> "$RESULT_FILE"
  fi
done < "$DATASET_FILE"

echo "Done. Results written to: $RESULT_FILE"
