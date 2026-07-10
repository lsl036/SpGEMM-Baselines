#!/bin/bash
# PaToH hypergraph partition for mypaper datasets only.
# Run from this directory (Matrix-Partitioning-Utility/PaToH). Requires DATA_PATH.

export OMP_PROC_BIND=spread
export OMP_PLACES=cores

# Paper dataset: same names as mypaper_dataset in 1_generate_close_pairs.sh
mypaper_dataset=(
  "cant" "pwtk" "offshore" "cage12" "scircuit" "wiki-Vote"
  "poisson3Da" "pdb1HYS" "rma10" "shipsec1" "consph" "filter3D"
  "mac_econ_fwd500" "af_shell10" "hood" "case39" "gupta3"
  "TSOPF_FS_b300_c2" "com-LiveJournal" "wikipedia-20070206"
)

if [ -z "$DATA_PATH" ]; then
  echo "Error: DATA_PATH is not set."
  exit 1
fi

TIME_LOG_FILE="${TIME_LOG_FILE:-hp_partition_mypaper_timev.log}"
RESULTS_FILE="${RESULTS_FILE:-hp_partition_mypaper_timing.log}"

extract_partition_time() {
  local output="$1"
  printf "%s\n" "${output}" | awk '
    /PartitionTime:/ {
      for (i = 1; i <= NF; i++) {
        if ($i ~ /^[0-9]+(\.[0-9]+)?$/) {
          print $i
          exit
        }
      }
    }
  '
}

extract_max_rss_kb() {
  local output="$1"
  printf "%s\n" "${output}" | awk -F': *' '
    /Maximum resident set size \(kbytes\)/ { print $2; exit }
  '
}

record_timing() {
  local dataset="$1"
  local output="$2"
  local rc="$3"
  local partition_time max_rss_kb line

  partition_time="$(extract_partition_time "${output}")"
  max_rss_kb="$(extract_max_rss_kb "${output}")"

  partition_time="${partition_time:-NA}"
  max_rss_kb="${max_rss_kb:-NA}"

  line=$(printf '%s %s %s' "${dataset}" "${partition_time}" "${max_rss_kb}")
  echo "${line}" | tee -a "${RESULTS_FILE}"
}

if [[ ! -f "${RESULTS_FILE}" ]]; then
  printf '%s %s %s\n' "dataset" "PartitionTime(s)" "MaxRSS(kB)" > "${RESULTS_FILE}"
fi

echo "PaToH partition ..."
echo "Timing results will be appended to: ${RESULTS_FILE}"

for dataset in "${mypaper_dataset[@]}"; do
  echo "~~~~~~~~~~~~~~~~~~~~~~~~~<DATASET: ${dataset}>~~~~~~~~~~~~~~~~~~~~~~~~~"
  timed_out="$(
    /usr/bin/time -v ./a.out "${DATA_PATH}/${dataset}/${dataset}.mtx" cutpart 64 quality 1 1 2>&1
  )"
  rc=$?

  printf "%s\n" "${timed_out}"

  {
    echo "===== ${dataset} ====="
    printf "%s\n" "${timed_out}"
    echo "Exit code: ${rc}"
    echo
  } >> "${TIME_LOG_FILE}"

  record_timing "${dataset}" "${timed_out}" "${rc}"
  echo ""
done
