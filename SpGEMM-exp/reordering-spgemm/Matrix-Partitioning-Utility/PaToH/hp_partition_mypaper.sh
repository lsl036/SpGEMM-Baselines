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
SUMMARY_FILE="${SUMMARY_FILE:-hp_partition_mypaper_mem_summary.csv}"
echo "dataset,exit_code,MaxMem_GB" > "${SUMMARY_FILE}"

for dataset in "${mypaper_dataset[@]}"; do
  echo "~~~~~~~~~~~~~~~~~~~~~~~~~<DATASET: ${dataset}>~~~~~~~~~~~~~~~~~~~~~~~~~"
  timed_out="$(
    /usr/bin/time -v ./a.out "${DATA_PATH}/${dataset}/${dataset}.mtx" cutpart 56 quality 1 1 2>&1
  )"
  rc=$?

  # Keep command output visible on terminal.
  printf "%s\n" "${timed_out}"

  # Persist full /usr/bin/time -v output.
  {
    echo "===== ${dataset} ====="
    printf "%s\n" "${timed_out}"
    echo "Exit code: ${rc}"
    echo
  } >> "${TIME_LOG_FILE}"

  max_rss_kb="$(printf "%s\n" "${timed_out}" | awk -F': *' '/Maximum resident set size \(kbytes\)/{print $2; exit}')"
  if [ -n "${max_rss_kb}" ]; then
    max_mem_gb="$(awk -v kb="${max_rss_kb}" 'BEGIN { printf "%.6f", kb / 1024 / 1024 }')"
  else
    max_mem_gb="NaN"
  fi
  echo "${dataset},${rc},${max_mem_gb}" >> "${SUMMARY_FILE}"

  echo ""
done
