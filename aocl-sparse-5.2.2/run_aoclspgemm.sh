#!/usr/bin/env bash
# 不使用 -e，避免某个矩阵出错时整个脚本提前退出
set -u pipefail

# 根路径按你的说明
# DATASET_LIST="600dataset.txt"
DATASET_LIST="casedataset.txt"
MTX_ROOT="${MTX_ROOT:-/data/suitesparse_collection}"
BENCH="./build-amd-gcc/tests/examples/sample_spgemm_bench"

# 脚本所在目录（用于定位编译产物与共享库）
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# 库的实际安装位置：build/library/libaoclsparse.so
LIB_DIR="${SCRIPT_DIR}/build/library"
export LD_LIBRARY_PATH="${LIB_DIR}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

# 结果文件（CSV），放在脚本同目录下
RESULT_FILE="${RESULT_FILE:-${SCRIPT_DIR}/aoclspgemm_results.txt}"

if [ ! -x "$BENCH" ]; then
  echo "Error: benchmark executable not found or not executable: $BENCH" >&2
  exit 1
fi

if [ ! -f "$DATASET_LIST" ]; then
  echo "Error: dataset list not found: $DATASET_LIST" >&2
  exit 1
fi

# 写入 CSV 表头（与 benchmark 单行输出顺序一致）
echo "mtx_name,nnz_count_ms,finalize_ms,total_ms" > "$RESULT_FILE"

total=$(wc -l < "$DATASET_LIST")
count=0

while IFS= read -r name || [ -n "$name" ]; do
  # 跳过空行或注释
  [ -z "$name" ] && continue
  case "$name" in
    \#*) continue ;;
  esac

  count=$((count + 1))

  MTX_DIR="${MTX_ROOT}/${name}"
  MTX_FILE="${MTX_DIR}/${name}.mtx"

  if [ ! -f "$MTX_FILE" ]; then
    echo "[$count/$total] Skip (file not found): ${name}"
    echo "${name},SKIP,,," >> "$RESULT_FILE"
    continue
  fi

  echo "[$count/$total] Running: ${name}"
  # 捕获 benchmark 输出（它本身输出一行：mtx_name, nnz_count, finalize, total）
  output=$("$BENCH" "$MTX_FILE" 2>&1)
  exitcode=$?

  if [ $exitcode -ne 0 ]; then
    echo "  -> run failed (exit $exitcode)"
    echo "${name},ERROR,,," >> "$RESULT_FILE"
    continue
  fi

  # 取 benchmark 输出的最后一行数据（跳过可能的 banner/提示）
  line=$(echo "$output" | grep -E '^[^[:space:]]+,.*[0-9]' | tail -n 1)
  if [ -n "$line" ]; then
    echo "$line" >> "$RESULT_FILE"
    echo "  -> recorded"
  else
    echo "  -> no data line in output"
    echo "${name},NO_OUTPUT,,," >> "$RESULT_FILE"
  fi
done < "$DATASET_LIST"

echo "Done. Results written to: $RESULT_FILE"