#!/bin/bash
# Run HASpGEMM on all datasets and record peak memory via /usr/bin/time -v.

set -u -o pipefail

DATA_PATH="${DATA_PATH:-/data/suitesparse_collection}"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATASET_FILE="${DATASET_FILE:-${SCRIPT_DIR}/20data.txt}"
# DATASET_FILE="${SCRIPT_DIR}/20data.txt"
RESULT_FILE="${RESULT_FILE:-${SCRIPT_DIR}/haspgemm_results_with_mem.txt}"
BIN="${BIN:-${SCRIPT_DIR}/haspgemm12}"

if [ ! -f "$DATASET_FILE" ]; then
  echo "Error: dataset file not found: $DATASET_FILE"
  exit 1
fi

if [ ! -x "$BIN" ]; then
  echo "Error: haspgemm12 not found or not executable: $BIN"
  exit 1
fi

# Keep original columns and append MaxMem_GB.
echo "HASpGEMM,mat_name,nnzA,mA,avg_nnzCub,avg_gflops,max_gflops,avg_time,analyze_time,MaxMem_GB" > "$RESULT_FILE"

total=$(wc -l < "$DATASET_FILE")
count=0

while IFS= read -r name || [ -n "$name" ]; do
  name="${name%$'\r'}"
  name="${name#$'\xef\xbb\xbf'}"
  [ -z "$name" ] && continue
  case "$name" in
    \#*) continue ;;
  esac

  count=$((count + 1))
  mtx_path="${DATA_PATH}/${name}/${name}.mtx"

  if [ ! -f "$mtx_path" ]; then
    echo "[$count/$total] Skip (file not found): $name"
    echo "HASpGEMM,${name},SKIP,,,,,,,NaN" >> "$RESULT_FILE"
    continue
  fi

  echo "[$count/$total] Running: $name"
  timed_out="$(
    /usr/bin/time -v "$BIN" "$mtx_path" "$mtx_path" 2>&1
  )"
  exitcode=$?

  max_rss_kb="$(printf "%s\n" "$timed_out" | awk -F': *' '/Maximum resident set size \(kbytes\)/{print $2; exit}')"
  if [ -n "$max_rss_kb" ]; then
    max_mem_gb="$(awk -v kb="$max_rss_kb" 'BEGIN { printf "%.6f", kb / 1024 / 1024 }')"
  else
    max_mem_gb="NaN"
  fi

  if [ $exitcode -ne 0 ]; then
    echo "  -> run failed (exit $exitcode)"
    echo "HASpGEMM,${name},ERROR,,,,,,,${max_mem_gb}" >> "$RESULT_FILE"
    continue
  fi

  line="$(printf "%s\n" "$timed_out" | awk '/^HASpGEMM,/{last=$0} END{print last}')"
  if [ -n "$line" ]; then
    echo "${line},${max_mem_gb}" >> "$RESULT_FILE"
    echo "  -> recorded (MaxMem=${max_mem_gb} GB)"
  else
    echo "  -> no HASpGEMM line in output"
    echo "HASpGEMM,${name},NO_OUTPUT,,,,,,,${max_mem_gb}" >> "$RESULT_FILE"
  fi
done < "$DATASET_FILE"

echo "Done. Results written to: $RESULT_FILE"

