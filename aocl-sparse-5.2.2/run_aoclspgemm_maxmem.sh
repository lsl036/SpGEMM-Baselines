#!/usr/bin/env bash
# Batch AOCL SpGEMM test with peak memory collection via /usr/bin/time -v
# Output columns:
# mtx_name, nnz_count_ms, finalize_ms, total_ms, MaxMem_GB

set -u -o pipefail
shopt -s nullglob

DATASET_LIST="${DATASET_LIST:-casedataset.txt}"
MTX_ROOT="${MTX_ROOT:-/data/suitesparse_collection}"
BENCH="${BENCH:-./build-amd-gcc/tests/examples/sample_spgemm_bench}"
OUT_FILE="${OUT_FILE:-AOCL_SpGEMM_with_mem.txt}"

if [ ! -x "$BENCH" ]; then
  echo "Error: benchmark executable not found or not executable: $BENCH" >&2
  exit 1
fi

if [ ! -f "$DATASET_LIST" ]; then
  echo "Error: dataset list not found: $DATASET_LIST" >&2
  exit 1
fi

# CSV header
echo "mtx_name, nnz_count_ms, finalize_ms, total_ms, MaxMem_GB" > "$OUT_FILE"

total=0
failed=0

while IFS= read -r name; do
  # Normalize line endings / BOM to avoid hidden-char filename issues.
  name="${name%$'\r'}"
  name="${name#$'\xef\xbb\xbf'}"
  [ -z "$name" ] && continue
  case "$name" in
    \#*) continue ;;
  esac

  total=$((total + 1))
  resolved_name="$name"
  MTX_FILE="${MTX_ROOT}/${resolved_name}/${resolved_name}.mtx"

  # Fallback: if exact path not found, try unique fuzzy directory match.
  if [ ! -f "$MTX_FILE" ]; then
    candidates=( "${MTX_ROOT}"/*"${name}"* )
    if [ "${#candidates[@]}" -eq 1 ] && [ -d "${candidates[0]}" ]; then
      cand_base="$(basename "${candidates[0]}")"
      cand_mtx="${candidates[0]}/${cand_base}.mtx"
      if [ -f "$cand_mtx" ]; then
        resolved_name="$cand_base"
        MTX_FILE="$cand_mtx"
        echo "Info: auto-resolved dataset '${name}' -> '${resolved_name}'" >&2
      fi
    fi
  fi

  if [ ! -f "$MTX_FILE" ]; then
    echo "Warning: matrix file not found: ${MTX_FILE}" >&2
    echo "${name}, NaN, NaN, NaN, NaN" >> "$OUT_FILE"
    failed=$((failed + 1))
    continue
  fi

  # Capture stdout (bench result line) and stderr (/usr/bin/time -v report)
  timed_out="$(
    /usr/bin/time -v "$BENCH" "$MTX_FILE" 2>&1
  )"
  rc=$?

  if [ $rc -ne 0 ]; then
    echo "ERROR: benchmark failed on ${name} (exit=${rc})" >&2
    echo "${name}, NaN, NaN, NaN, NaN" >> "$OUT_FILE"
    failed=$((failed + 1))
    continue
  fi

  # Benchmark line is expected like:
  # mtx_name, nnz_count_ms, finalize_ms, total_ms
  bench_line="$(printf "%s\n" "$timed_out" | awk -F',' '/,/{print $0; exit}')"
  max_rss_kb="$(printf "%s\n" "$timed_out" | awk -F': *' '/Maximum resident set size \(kbytes\)/{print $2; exit}')"

  if [ -z "$bench_line" ] || [ -z "$max_rss_kb" ]; then
    echo "ERROR: parse failed on ${name}" >&2
    echo "${name}, NaN, NaN, NaN, NaN" >> "$OUT_FILE"
    failed=$((failed + 1))
    continue
  fi

  max_mem_gb="$(awk -v kb="$max_rss_kb" 'BEGIN { printf "%.6f", kb / 1024 / 1024 }')"
  echo "${bench_line}, ${max_mem_gb}" >> "$OUT_FILE"
  echo "OK: ${name} MaxMem=${max_mem_gb} GB"
done < "$DATASET_LIST"

echo "Done: $((total - failed))/${total} passed, ${failed} failed. Output: ${OUT_FILE}" >&2

