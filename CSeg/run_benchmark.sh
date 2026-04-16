#!/bin/bash
# Run SpGEMM benchmark on all datasets in dataset.txt, average 5 runs, save [mtx_name runtime] to result file.

DATA_PATH="${DATA_PATH:-/data/suitesparse_collection}"
# DATASET_FILE="$(dirname "$0")/casedataset.txt"
DATASET_FILE="$(dirname "$0")/600dataset.txt"
RESULT_FILE="$(dirname "$0")/spgemm_runtime_results0401.txt"
MAIN_BIN="$(dirname "$0")/main"

if [ ! -f "$DATASET_FILE" ]; then
  echo "Error: dataset file not found: $DATASET_FILE"
  exit 1
fi

if [ ! -x "$MAIN_BIN" ]; then
  echo "Error: main executable not found or not executable: $MAIN_BIN"
  exit 1
fi

# Overwrite result file with header
echo "mtx_name	runtime_sec" > "$RESULT_FILE"

total=$(wc -l < "$DATASET_FILE")
count=0

while IFS= read -r name || [ -n "$name" ]; do
  name=$(echo "$name" | tr -d '\r')
  [ -z "$name" ] && continue

  count=$((count + 1))
  mtx_path="${DATA_PATH}/${name}/${name}.mtx"

  if [ ! -f "$mtx_path" ]; then
    echo "[$count/$total] Skip (file not found): $name"
    echo "${name}	SKIP" >> "$RESULT_FILE"
    continue
  fi

  echo "[$count/$total] Running: $name"
  output=$("$MAIN_BIN" -A "$mtx_path" -B "$mtx_path" -w 1 -e 5 2>&1)

  # Extract "Total runtime (secs) : X.XXXXXX" lines (perf runs only; exclude warmup line)
  times=($(echo "$output" | grep "Total runtime (secs) :" | sed 's/.*Total runtime (secs) : *//'))

  if [ ${#times[@]} -lt 5 ]; then
    echo "  Warning: got ${#times[@]} run(s), expected 5. Check output for errors."
    avg="ERROR"
  else
    # Average the 5 runs (use awk for float)
    avg=$(printf '%s\n' "${times[@]}" | awk '{ sum += $1; n++ } END { if (n>0) printf "%.6f", sum/n; else print "ERROR" }')
  fi

  echo "${name}	${avg}" >> "$RESULT_FILE"
  echo "  -> ${avg} sec"
done < "$DATASET_FILE"

echo "Done. Results written to: $RESULT_FILE"
