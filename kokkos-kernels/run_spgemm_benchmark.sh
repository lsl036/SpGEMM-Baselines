#!/usr/bin/env bash
# Kokkos-Kernels SpGEMM benchmark: run on datasets from a list file,
# parse mm_time (seconds), average and convert to ms, save as [mtx_name mm_time(ms)].
#
# === Performance: if results are slow, check ===
# 1) Build type: cmake with -DCMAKE_BUILD_TYPE=Release (and build Kokkos with Release).
# 2) OpenMP: script uses --openmp N; set OMP_PROC_BIND=spread and OMP_PLACES=threads (done below).
# 3) Algorithm: default is KKSPGEMM; optionally pass --algorithm KKSPGEMM or KKDENSE for tuning.
# 4) Do NOT redirect stderr so you see "Running on OpenMP backend." and thread topology.

set -e

# Defaults (override with env or args)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="${SCRIPT_DIR}/perf_test/sparse/sparse_spgemm"
DATASET_LIST="${SCRIPT_DIR}/600dataset.txt"
MATRIX_PREFIX="${SPGEMM_DATASET_PREFIX:-/data/suitesparse_collection}"
REPEAT="${REPEAT:-5}"
OPENMP_THREADS="${OPENMP_THREADS:-64}"
# Optional: KKSPGEMM (default), KKDENSE, KKMEM, KKLP. Empty = use program default.
SPGEMM_ALGORITHM="${SPGEMM_ALGORITHM:-}"
OUTPUT_FILE=""

usage() {
  echo "Usage: $0 [OPTIONS]"
  echo "  -d FILE    Dataset list file (one matrix name per line). Default: casedataset.txt"
  echo "  -p DIR     Matrix path prefix. Path = DIR/name/name.mtx. Default: \$SPGEMM_DATASET_PREFIX or /data/suitesparse_collection"
  echo "  -r N       Repeat count per matrix. Default: 5"
  echo "  -t N       OpenMP threads. Default: 28"
  echo "  -a ALG     SpGEMM algorithm: KKSPGEMM (default), KKDENSE, KKMEM, KKLP. Default: program default"
  echo "  -o FILE    Output file for results [mtx_name mm_time(ms)]. Default: spgemm_results_<date>.txt"
  echo "  -h         Print this help"
  exit 0
}

while getopts "d:p:r:t:a:o:h" opt; do
  case "$opt" in
    d) DATASET_LIST="$OPTARG" ;;
    p) MATRIX_PREFIX="$OPTARG" ;;
    r) REPEAT="$OPTARG" ;;
    t) OPENMP_THREADS="$OPTARG" ;;
    a) SPGEMM_ALGORITHM="$OPTARG" ;;
    o) OUTPUT_FILE="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

if [[ ! -x "$EXE" ]]; then
  echo "Error: executable not found or not runnable: $EXE" >&2
  exit 1
fi

if [[ ! -f "$DATASET_LIST" ]]; then
  echo "Error: dataset list not found: $DATASET_LIST" >&2
  exit 1
fi

if [[ -z "$OUTPUT_FILE" ]]; then
  OUTPUT_FILE="spgemm_results_$(date +%Y%m%d_%H%M%S).txt"
fi

# Recommended for OpenMP performance (Kokkos prints a warning if unset)
export OMP_PROC_BIND="${OMP_PROC_BIND:-spread}"
export OMP_PLACES="${OMP_PLACES:-threads}"

echo "SpGEMM benchmark"
echo "  executable:    $EXE"
echo "  dataset list:  $DATASET_LIST"
echo "  matrix prefix: $MATRIX_PREFIX"
echo "  repeat:        $REPEAT"
echo "  OpenMP threads: $OPENMP_THREADS"
[[ -n "$SPGEMM_ALGORITHM" ]] && echo "  algorithm:     $SPGEMM_ALGORITHM"
echo "  OMP_PROC_BIND: $OMP_PROC_BIND  OMP_PLACES: $OMP_PLACES"
echo "  output file:   $OUTPUT_FILE"
echo ""

# Clear or create output file with header
printf "%-40s %s\n" "mtx_name" "mm_time(ms)" > "$OUTPUT_FILE"

run_one() {
  local name="$1"
  local mtx_path="${MATRIX_PREFIX}/${name}/${name}.mtx"
  if [[ ! -f "$mtx_path" ]]; then
    echo "  SKIP (file not found): $name"
    printf "%-40s %s\n" "$name" "SKIP" >> "$OUTPUT_FILE"
    return 0
  fi

  local raw cmd
  cmd=("$EXE" --amtx "$mtx_path" --repeat "$REPEAT" --openmp "$OPENMP_THREADS")
  [[ -n "$SPGEMM_ALGORITHM" ]] && cmd+=(--algorithm "$SPGEMM_ALGORITHM")
  raw=$("${cmd[@]}" 2>/dev/null || true)
  if [[ -z "$raw" ]]; then
    echo "  FAIL (no output): $name"
    printf "%-40s %s\n" "$name" "FAIL" >> "$OUTPUT_FILE"
    return 0
  fi

  # Extract mm_time values (seconds); lines like "mm_time:0.0750248 ..."
  local avg_ms count
  avg_ms=$(echo "$raw" | grep 'mm_time:' | sed 's/.*mm_time:\([0-9.]*\).*/\1/' | awk '{sum+=$1; n++} END {if(n>0) printf "%.4f", (sum/n)*1000; else print "FAIL"}')
  count=$(echo "$raw" | grep -c 'mm_time:' || echo 0)

  if [[ "$avg_ms" == "FAIL" ]] || [[ "$count" -eq 0 ]]; then
    echo "  FAIL (no mm_time): $name"
    printf "%-40s %s\n" "$name" "FAIL" >> "$OUTPUT_FILE"
    return 0
  fi
  echo "  $name  avg mm_time = ${avg_ms} ms (from $count samples)"
  printf "%-40s %s\n" "$name" "$avg_ms" >> "$OUTPUT_FILE"
}

total=0
while IFS= read -r name || [[ -n "$name" ]]; do
  name=$(echo "$name" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
  [[ -z "$name" ]] && continue
  [[ "$name" =~ ^# ]] && continue
  ((total++)) || true
  echo "[$total] $name"
  run_one "$name"
done < "$DATASET_LIST"

echo ""
echo "Results written to: $OUTPUT_FILE"
