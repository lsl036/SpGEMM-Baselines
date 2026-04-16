#!/bin/bash

# OpenMP settings
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

# Common path prefixes
CLUSTERWISE_ROOT="/data/lsl/SpGEMM/SpGEMM-exp/clusterwise-spgemm-main"
BIN_PATH="${CLUSTERWISE_ROOT}/bin/GenerateCandidatePairs_hw"

# DATA_PATH 可以在外部环境中提前 export，
# 这里提供一个默认值，按需修改
: "${DATA_PATH:=/data/suitesparse_collection}"

# 第一个参数：数据集列表文件（默认使用当前目录下的 dataset.txt）
DATASET_FILE="${1:-${CLUSTERWISE_ROOT}/scripts/h_clusterwise/casedata.txt}"

# 第二个参数：模式，"run" 表示只运行，"save" 表示带 -s 选项保存（默认 run）
MODE="${2:-run}"

if [[ ! -f "${DATASET_FILE}" ]]; then
  echo "Dataset file not found: ${DATASET_FILE}"
  exit 1
fi

echo "Using DATA_PATH = ${DATA_PATH}"
echo "Reading dataset list from ${DATASET_FILE}"
echo "Mode = ${MODE} (run=normal, save=with -s)"
echo ""

if [[ "${MODE}" == "run" ]]; then
  TIME_LOG_FILE="${TIME_LOG_FILE:-hclusterwise_run_timev.log}"
  SUMMARY_FILE="${SUMMARY_FILE:-hclusterwise_run_mem_summary.csv}"
  echo "dataset,exit_code,MaxMem_GB" > "${SUMMARY_FILE}"
fi

while IFS= read -r dataset || [[ -n "$dataset" ]]; do
  # 跳过空行和以 # 开头的注释
  if [[ -z "${dataset}" || "${dataset}" =~ ^# ]]; then
    continue
  fi

  echo "~~~~~~~~~~~~~~~~~~~~~~~~~<DATASET: ${dataset}>~~~~~~~~~~~~~~~~~~~~~~~~~"

  if [[ "${MODE}" == "save" ]]; then
    "${BIN_PATH}" text \
      "${DATA_PATH}/${dataset}/${dataset}.mtx" \
      "${DATA_PATH}/${dataset}/${dataset}.mtx" \
      7 -s
  else
    timed_out="$(
      /usr/bin/time -v "${BIN_PATH}" text \
        "${DATA_PATH}/${dataset}/${dataset}.mtx" \
        "${DATA_PATH}/${dataset}/${dataset}.mtx" \
        7 2>&1
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
    if [[ -n "${max_rss_kb}" ]]; then
      max_mem_gb="$(awk -v kb="${max_rss_kb}" 'BEGIN { printf "%.6f", kb / 1024 / 1024 }')"
    else
      max_mem_gb="NaN"
    fi
    echo "${dataset},${rc},${max_mem_gb}" >> "${SUMMARY_FILE}"
  fi

  echo ""
done < "${DATASET_FILE}"

