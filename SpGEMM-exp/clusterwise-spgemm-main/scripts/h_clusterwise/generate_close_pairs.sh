#!/bin/bash
# Generate close pairs for mypaper datasets with per-matrix timing log.

export OMP_PROC_BIND=spread
export OMP_PLACES=cores

mypaper_dataset=(
  "cant" "pwtk" "offshore" "cage12" "scircuit" "wiki-Vote"
  "poisson3Da" "pdb1HYS" "rma10" "shipsec1" "consph" "filter3D"
  "mac_econ_fwd500" "af_shell10" "hood" "case39" "gupta3"
  "TSOPF_FS_b300_c2" "com-LiveJournal" "wikipedia-20070206"
)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CLUSTERWISE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
BIN_PATH="${CLUSTERWISE_ROOT}/bin/GenerateCandidatePairs_hw"
RESULTS_FILE="${RESULTS_FILE:-${SCRIPT_DIR}/close_pairs_timing.log}"

if [ -z "${DATA_PATH}" ]; then
  echo "Error: DATA_PATH is not set."
  exit 1
fi

if [ -z "${CLOSE_PAIR_DATA_PATH}" ]; then
  echo "Error: CLOSE_PAIR_DATA_PATH is not set."
  exit 1
fi

extract_generate_time() {
  local output="$1"
  printf '%s\n' "${output}" | awk '
    /Generate A and B in/ { print $(NF-1); exit }
  '
}

extract_preprocess_time() {
  local output="$1"
  printf '%s\n' "${output}" | awk '
    /Preprocess A and B in/ { print $(NF-1); exit }
  '
}

extract_aat_time() {
  local output="$1"
  printf '%s\n' "${output}" | awk '
    /done C = A \* A_T in/ { print $(NF-1); exit }
  '
}

extract_max_rss_kb() {
  local output="$1"
  printf '%s\n' "${output}" | awk -F': *' '
    /Maximum resident set size \(kbytes\)/ { print $2; exit }
  '
}

record_timing() {
  local dataset="$1"
  local output="$2"
  local generate_time preprocess_time aat_time max_rss_kb line

  generate_time="$(extract_generate_time "${output}")"
  preprocess_time="$(extract_preprocess_time "${output}")"
  aat_time="$(extract_aat_time "${output}")"
  max_rss_kb="$(extract_max_rss_kb "${output}")"

  generate_time="${generate_time:-NA}"
  preprocess_time="${preprocess_time:-NA}"
  aat_time="${aat_time:-NA}"
  max_rss_kb="${max_rss_kb:-NA}"

  line=$(printf '%s %s %s %s %s' \
    "${dataset}" "${generate_time}" "${preprocess_time}" "${aat_time}" "${max_rss_kb}")
  echo "${line}" | tee -a "${RESULTS_FILE}"
}

run_on_datasets() {
  local datasets=("$@")
  local output

  mkdir -p "${CLOSE_PAIR_DATA_PATH}"
  if [[ ! -f "${RESULTS_FILE}" ]]; then
    printf '%s %s %s %s %s\n' \
      "dataset" "GenerateAB(s)" "PreprocessAB(s)" "AAT(s)" "MaxRSS(kB)" > "${RESULTS_FILE}"
  fi

  for dataset in "${datasets[@]}"; do
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~<DATASET: ${dataset}>~~~~~~~~~~~~~~~~~~~~~~~~~"
    output=$(
      /usr/bin/time -v "${BIN_PATH}" text \
        "${DATA_PATH}/${dataset}/${dataset}.mtx" \
        "${DATA_PATH}/${dataset}/${dataset}.mtx" \
        7 2>&1
    )
    printf '%s\n' "${output}"
    record_timing "${dataset}" "${output}"
    echo ""
  done
}

echo "Generate close pairs ..."
echo "Timing results will be appended to: ${RESULTS_FILE}"
run_on_datasets "${mypaper_dataset[@]}"
