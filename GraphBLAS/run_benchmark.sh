#!/bin/bash
# Run GraphBLAS spgemm_demo on all datasets in dataset.txt.
# Records results as: "mtx_name avg_time"

set -euo pipefail

DATA_PATH="${DATA_PATH:-/data/suitesparse_collection}"
NTRIALS="${NTRIALS:-5}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATASET_FILE="${SCRIPT_DIR}/600dataset.txt"
# DATASET_FILE="${SCRIPT_DIR}/casedataset.txt"
RESULT_FILE="${RESULT_FILE:-${SCRIPT_DIR}/graphblas_spgemm_results0402.txt}"
BIN="${BIN:-${SCRIPT_DIR}/build/spgemm_demo}"

if [ ! -f "$DATASET_FILE" ]; then
  echo "Error: dataset file not found: $DATASET_FILE" >&2
  exit 1
fi

if [ ! -x "$BIN" ]; then
  echo "Error: spgemm_demo not found or not executable: $BIN" >&2
  echo "Hint: build it first with: (cd \"$SCRIPT_DIR\" && make demos)" >&2
  exit 1
fi

echo "# mtx_name avg_time_sec" > "$RESULT_FILE"

total=$(wc -l < "$DATASET_FILE" | tr -d ' ')
count=0

run_cmd() {
  # Use `timeout` if TIMEOUT is set and the command exists.
  if [ -n "${TIMEOUT:-}" ] && command -v timeout >/dev/null 2>&1; then
    timeout "${TIMEOUT}" "$@"
  else
    "$@"
  fi
}

while IFS= read -r name || [ -n "$name" ]; do
  name=$(echo "$name" | tr -d '\r')
  [ -z "$name" ] && continue

  count=$((count + 1))
  mtx_path="${DATA_PATH}/${name}/${name}.mtx"

  if [ ! -f "$mtx_path" ]; then
    echo "[$count/$total] Skip (file not found): $name" >&2
    echo "${name} SKIP" >> "$RESULT_FILE"
    continue
  fi

  echo "[$count/$total] Running: $name" >&2

  set +e
  output=$(run_cmd "$BIN" "$mtx_path" "$mtx_path" "$NTRIALS" 2>&1)
  exitcode=$?
  set -e

  if [ $exitcode -ne 0 ]; then
    echo "  -> run failed (exit $exitcode)" >&2
    echo "${name} ERROR" >> "$RESULT_FILE"
    continue
  fi

  # Parse line like: "avg  time: 0.123 sec"
  avg_time=$(echo "$output" | awk '
    /^avg[[:space:]]+time:/ {v=$3}
    END { if (v != "") print v; }
  ')

  if [ -n "$avg_time" ]; then
    echo "${name} ${avg_time}" >> "$RESULT_FILE"
    echo "  -> avg_time_sec = $avg_time" >&2
  else
    echo "  -> no avg time in output" >&2
    echo "${name} NO_OUTPUT" >> "$RESULT_FILE"
  fi
done < "$DATASET_FILE"

echo "Done. Results written to: $RESULT_FILE" >&2

