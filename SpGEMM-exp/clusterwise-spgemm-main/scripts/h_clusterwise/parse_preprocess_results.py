#!/usr/bin/env python3
"""解析 1_generate_cp_MemTime.txt 预处理结果，提取关键数据并输出到 CSV。"""

import re
import csv

# INPUT_FILE = "1_generate_cp_MemTime.txt"
INPUT_FILE = "generate_candidate_pairs_all.txt"
OUTPUT_CSV = "preprocess_600dataset_results.csv"

# 表头：与用户要求的顺序一致
HEADERS = [
    "DATA_name",
    "mtx Rows",
    "mtx nnz",
    "sizeofCSC",
    "GenAandBtime",
    "PreprocessAandBtime",
    "Floating-pointOperations",
    "C.nnz",
    "CTotalSize",
    "BinSize",
    "AA_Ttime",
]

#  整理cluster-spgemm 的预处理开销
def parse_file(path: str) -> list[dict]:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        content = f.read()

    # 按数据集块分割：每个块以 ~~~~<DATASET: xxx>~~~~ 开始
    blocks = re.split(r"~+<DATASET:\s*([^>]+)>~+", content)
    # blocks[0] 可能是文件开头的空或无关内容，从 blocks[1] 开始是 (name, block_content) 交替
    records = []
    i = 1
    while i + 1 < len(blocks):
        data_name = blocks[i].strip()
        block = blocks[i + 1]
        i += 2

        row = {"DATA_name": data_name}

        # mtx Rows
        m = re.search(r"Input Matrix:\s*Rows\s*=\s*(\d+)", block)
        row["mtx Rows"] = m.group(1) if m else ""

        # mtx nnz
        m = re.search(r"nnz\s*=\s*(\d+)", block)
        row["mtx nnz"] = m.group(1) if m else ""

        # sizeofCSC：取第一个 "Size of CSC matrix:"（在 Converting to csc 之后那个，即主矩阵）
        m = re.search(r"Size of CSC matrix:\s*([\d.e+-]+)\s*GB", block)
        row["sizeofCSC"] = m.group(1) if m else ""

        # GenAandBtime
        m = re.search(r"Generate A and B in ([\d.]+) \[seconds\]", block)
        row["GenAandBtime"] = m.group(1) if m else ""

        # PreprocessAandBtime
        m = re.search(r"Preprocess A and B in ([\d.]+) \[seconds\]", block)
        row["PreprocessAandBtime"] = m.group(1) if m else ""

        # Floating-pointOperations
        m = re.search(
            r"Total number of floating-point operations including addition and multiplication in SpGEMM \(A \* B\):\s*(\d+)",
            block,
        )
        row["Floating-pointOperations"] = m.group(1) if m else ""

        # C.nnz（hash_symbolic_topK 后面的数字）
        m = re.search(r"\[DONE\] hash_symbolic_topK\s+(\d+)", block)
        row["C.nnz"] = m.group(1) if m else ""

        # CTotalSize
        m = re.search(r"C total size:\s*([\d.]+)\s*GB", block)
        row["CTotalSize"] = m.group(1) if m else ""

        # BinSize
        m = re.search(r"Total size of BIN class instance:\s*([\d.]+)\s*GB", block)
        row["BinSize"] = m.group(1) if m else ""

        # AA_Ttime
        m = re.search(
            rf"Dataset:\s*{re.escape(data_name)}\s+done C = A \* A_T in ([\d.]+) \[seconds\]",
            block,
        )
        if not m:
            m = re.search(r"done C = A \* A_T in ([\d.]+) \[seconds\]", block)
        row["AA_Ttime"] = m.group(1) if m else ""

        records.append(row)

    return records


def main():
    records = parse_file(INPUT_FILE)
    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=HEADERS, extrasaction="ignore")
        w.writeheader()
        w.writerows(records)
    print(f"已解析 {len(records)} 条记录，已写入 {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
