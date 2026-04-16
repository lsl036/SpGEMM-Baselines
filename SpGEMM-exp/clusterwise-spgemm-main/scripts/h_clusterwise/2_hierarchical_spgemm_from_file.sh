#!/bin/bash

# OpenMP settings
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

# Common path prefixes (与 2_hierarchical_spgemm.sh 一致)
CLUSTERWISE_ROOT="/data/lsl/SpGEMM/SpGEMM-exp/clusterwise-spgemm-main"
BIN_PATH="${CLUSTERWISE_ROOT}/bin/HierarchicalClusterSpGEMM_hw"

# 矩阵数据根目录、close_pairs 输出目录，可在外部 export 覆盖
: "${DATA_PATH:=/data/suitesparse_collection}"
: "${CLOSE_PAIR_DATA_PATH:=/data/linshengle_data/SpGEMM-Reordering/close_pairs}"
# GP candidate
# : "${CLOSE_PAIR_DATA_PATH:=/data2/linshengle_data/SpGEMM-Reordering/gp_order}"
# hp candidate
# CLOSE_PAIR_DATA_PATH="/data2/linshengle_data/SpGEMM-Reordering/hp_order"

# 第一个参数：数据集列表文件（默认使用当前目录下的 dataset.txt）
DATASET_FILE="${1:-${CLUSTERWISE_ROOT}/scripts/h_clusterwise/casedata.txt}"

if [[ ! -f "${DATASET_FILE}" ]]; then
  echo "Dataset file not found: ${DATASET_FILE}"
  exit 1
fi

echo "Using DATA_PATH = ${DATA_PATH}"
echo "Using CLOSE_PAIR_DATA_PATH = ${CLOSE_PAIR_DATA_PATH}"
echo "Reading dataset list from ${DATASET_FILE}"
echo ""

while IFS= read -r dataset || [[ -n "$dataset" ]]; do
  # 跳过空行和以 # 开头的注释
  if [[ -z "${dataset}" || "${dataset}" =~ ^# ]]; then
    continue
  fi

  echo "~~~~~~~~~~~~~~~~~~~~~~~~~<DATASET: ${dataset}>~~~~~~~~~~~~~~~~~~~~~~~~~"
  "${BIN_PATH}" text \
    "${DATA_PATH}/${dataset}/${dataset}.mtx" \
    "${DATA_PATH}/${dataset}/${dataset}.mtx" \
    "${CLOSE_PAIR_DATA_PATH}/${dataset}.mtx"
  echo ""
done < "${DATASET_FILE}"
