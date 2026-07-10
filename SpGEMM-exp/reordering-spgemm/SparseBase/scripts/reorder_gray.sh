#!/bin/bash

fig10_datasesets=("cant" "pwtk" "offshore" "cage12" "scircuit" "wiki-Vote" "poisson3Da"
"pdb1HYS" "rma10" "shipsec1" "consph" "filter3D" "mac_econ_fwd500" "af_shell10"
"hood" "case39" "gupta3" "TSOPF_FS_b300_c2" "com-LiveJournal" "wikipedia-20070206")

GRAY_ORDER_BIN="${SPARSEBASE_HOME}/build/examples/gray_order/gray_order"
RESULTS_FILE="${GRAY_DATA_PATH}/gray_timing.log"

extract_value() {
  local pattern="$1"
  local field="$2"
  local output="$3"
  echo "$output" | awk -v pattern="$pattern" -v field="$field" '
    index($0, pattern) == 1 { print $field; exit }
  '
}

extract_max_rss_kb() {
  local output="$1"
  printf '%s\n' "$output" | awk -F': *' '
    /Maximum resident set size \(kbytes\)/ { print $2; exit }
  '
}

record_timing() {
  local dataset="$1"
  local output="$2"
  local vertices edges reorder_time max_rss_kb line

  vertices=$(extract_value "Number of vertices:" 4 "$output")
  edges=$(extract_value "Number of edges:" 4 "$output")
  reorder_time=$(extract_value "Reordering takes:" 3 "$output")
  max_rss_kb=$(extract_max_rss_kb "$output")

  vertices=${vertices:-NA}
  edges=${edges:-NA}
  reorder_time=${reorder_time:-NA}
  max_rss_kb=${max_rss_kb:-NA}

  line=$(printf '%s %s %s %s %s' \
    "$dataset" "$vertices" "$edges" "$reorder_time" "$max_rss_kb")
  echo "$line" | tee -a "$RESULTS_FILE"
}

run_gray_on_datasets() {
  local data_path="$1"; shift
  local datasets=("$@")
  local output

  mkdir -p "${GRAY_DATA_PATH}"
  if [[ ! -f "$RESULTS_FILE" ]]; then
    printf '%s %s %s %s %s\n' \
      "dataset" "vertices" "edges" "Time(s)" "MaxRSS(kB)" > "$RESULTS_FILE"
  fi

  for dataset in "${datasets[@]}"; do
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~<DATASET: ${dataset}>~~~~~~~~~~~~~~~~~~~~~~~~~"
    output=$(
      /usr/bin/time -v "${GRAY_ORDER_BIN}" \
        "${data_path}/${dataset}/${dataset}.mtx" \
        "${GRAY_DATA_PATH}/${dataset}.grayorder" 1 2>&1
    )
    printf '%s\n' "$output"
    record_timing "$dataset" "$output"
    echo ""
  done
}

echo "Gray reordering ..."
echo "Timing results will be appended to: ${RESULTS_FILE}"
run_gray_on_datasets "${DATA_PATH}" "${fig10_datasesets[@]}"
