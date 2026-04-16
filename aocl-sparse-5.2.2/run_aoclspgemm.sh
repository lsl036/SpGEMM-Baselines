#!/usr/bin/env bash
# 不使用 -e，避免某个矩阵出错时整个脚本提前退出
set -u pipefail

# 根路径按你的说明
DATASET_LIST="600dataset.txt"
MTX_ROOT="/data/suitesparse_collection"
BENCH="./build-amd-gcc/tests/examples/sample_spgemm_bench"

if [ ! -x "$BENCH" ]; then
  echo "Error: benchmark executable not found or not executable: $BENCH" >&2
  exit 1
fi

if [ ! -f "$DATASET_LIST" ]; then
  echo "Error: dataset list not found: $DATASET_LIST" >&2
  exit 1
fi

# 可选：先输出表头
echo "mtx_name, nnz_count_ms, finalize_ms, total_ms"

while IFS= read -r name; do
  # 跳过空行或注释
  [ -z "$name" ] && continue
  case "$name" in
    \#*) continue ;;
  esac

  MTX_DIR="${MTX_ROOT}/${name}"
  MTX_FILE="${MTX_DIR}/${name}.mtx"

  if [ ! -f "$MTX_FILE" ]; then
    echo "Warning: matrix file not found: ${MTX_FILE}" >&2
    continue
  fi

  # 调用 benchmark（它本身输出一行：mtx_name, nnz_count, finalize, total）
  # 若单个矩阵出错（例如段错误），记录错误并继续后续数据集
  if ! "$BENCH" "$MTX_FILE"; then
    echo "ERROR: benchmark failed on ${name}" >&2
    continue
  fi
done < "$DATASET_LIST"