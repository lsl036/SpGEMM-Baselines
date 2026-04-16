#!/usr/bin/env python3
"""解析 2_hierarchical_spgemm_results.txt 数值计算结果，提取 DATA_name, threads, Time, MFLOPS 并输出到 CSV。"""

import re
import csv

# INPUT_FILE = "2_hierarchical_spgemm_results.txt"
# INPUT_FILE = "2_56cores_Performance.txt"
INPUT_FILE = "runspgemm_candidate_pairs_all.txt"
# OUTPUT_CSV = "spgemm_results.csv"
OUTPUT_CSV = "spgemm_results_600datasets.csv"

HEADERS = ["DATA_name", "threads", "Time", "MFLOPS"]

# HashSpGEMMVLCluster with  64 threads computes C = A * B in 70.661354 [milli seconds] (7627.549077 [MFLOPS])
PATTERN = re.compile(
    r"HashSpGEMMVLCluster with\s+(\d+)\s+threads computes C = A \* B in ([\d.]+) \[milli seconds\] \(([\d.]+) \[MFLOPS\]\)"
)

#  处理 SC'25 clutser-spgemm 的实验结果
def parse_file(path: str) -> list[dict]:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        content = f.read()

    blocks = re.split(r"~+<DATASET:\s*([^>]+)>~+", content)
    records = []
    i = 1
    while i + 1 < len(blocks):
        data_name = blocks[i].strip()
        block = blocks[i + 1]
        i += 2

        m = PATTERN.search(block)
        if m:
            records.append({
                "DATA_name": data_name,
                "threads": m.group(1),
                "Time": m.group(2),
                "MFLOPS": m.group(3),
            })
        else:
            records.append({
                "DATA_name": data_name,
                "threads": "",
                "Time": "",
                "MFLOPS": "",
            })

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
