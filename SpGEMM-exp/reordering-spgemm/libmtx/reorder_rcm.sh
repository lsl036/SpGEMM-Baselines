#!/bin/bash

# Dataset lists as strings
initial_datasets=(
  "2cubes_sphere" "cage12" "cop20k_A" "m133-b3" "mac_econ_fwd500" "mario002"
  "mc2depi" "poisson3Da" "pwtk" "webbase-1M" "cage15" "cant" "conf5_4-8x8-05"
  "consph" "delaunay_n24" "filter3D" "hood" "majorbasis" "mono_500Hz" "offshore"
  "patents_main" "pdb1HYS" "rma10" "scircuit" "shipsec1" "wb-edu"
)

ispass_datasets=(
  "333SP" "adaptive" "af_shell10" "AS365" "as-Skitter" "channel-500x100x100-b050"
  "cit-Patents" "com-LiveJournal" "CurlCurl_4" "europe_osm" "G3_circuit" "GAP-road"
  "hugebubbles-00020" "hugetrace-00020" "hugetric-00010" "HV15R" "kkt_power" "M6"
  "nlpkkt240" "NLR" "packing-500x100x100-b050" "patents" "Queen_4147" "rajat31"
  "rgg_n_2_24_s0" "soc-LiveJournal1" "stokes" "sx-stackoverflow" "venturiLevel3"
  "wikipedia-20070206" "wiki-Talk" "wiki-topcats"
)

new_datasets=(
  "mip1" "mouse_gene" "BenElechi1" "atmosmodl" "pkustk14" "great-britain_osm"
  "hugebubbles-00000" "hugebubbles-00010" "hugetrace-00000" "hugetrace-00010"
  "hugetrace-00010" "hugetric-00000" "hugetric-00020" "italy_osm"
  "kron_g500-logn17" "kron_g500-logn18" "kron_g500-logn19" "kron_g500-logn20"
  "rgg_n_2_20_s0" "rgg_n_2_21_s0" "rgg_n_2_22_s0" "rgg_n_2_23_s0" "road_central"
  "road_usa" "asia_osm" "coPapersCiteseer" "coPapersDBLP" "delaunay_n21"
  "delaunay_n22" "delaunay_n23" "germany_osm" "halfb" "dielFilterV3real"
  "kmer_V1r" "bmw3_2" "audikw_1" "Hardesty3" "JP" "msdoor" "relat9" "F1"
  "uk-2002" "fem_hifreq_circuit" "bundle_adj" "mycielskian18" "nd24k" "torso1"
  "bone010" "Ga41As41H72" "gearbox" "thermal2" "higgs-twitter"
  "TSOPF_RS_b2383_c1"
)

# res_datasets=(
#   "HV15R" "kmer_V1r" "bmw3_2" "audikw_1" "Hardesty3" "JP" "msdoor" "relat9" "F1"
#   "uk-2002" "fem_hifreq_circuit" "bundle_adj" "mycielskian18" "nd24k" "torso1"
#   "bone010" "Ga41As41H72" "gearbox" "thermal2" "higgs-twitter"
#   "TSOPF_RS_b2383_c1"
# )

fig10_datasesets=("cant" "pwtk" "offshore" "cage12" "scircuit" "wiki-Vote" "poisson3Da" 
"pdb1HYS" "rma10" "shipsec1" "consph" "filter3D" "mac_econ_fwd500" "af_shell10" 
"hood" "case39" "gupta3" "TSOPF_FS_b300_c2" "com-LiveJournal" "wikipedia-20070206")

BIN="mtxreorder --verbose"
RESULTS_FILE="${RCM_DATA_PATH}/rcm_timing.log"

extract_time() {
  local label="$1"
  local output="$2"
  echo "$output" | awk -v label="$label" '
    $1 == label ":" { print $2; exit }
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
  local read_time reorder_time fwrite_time write_time max_rss_kb
  local line

  read_time=$(extract_time "mtxfile_read" "$output")
  reorder_time=$(extract_time "mtxfile_reorder" "$output")
  fwrite_time=$(extract_time "mtxfile_fwrite" "$output")
  write_time=$(extract_time "mtxfile_write" "$output")
  max_rss_kb=$(extract_max_rss_kb "$output")

  read_time=${read_time:-NA}
  reorder_time=${reorder_time:-NA}
  fwrite_time=${fwrite_time:-NA}
  write_time=${write_time:-NA}
  max_rss_kb=${max_rss_kb:-NA}

  line=$(printf '%s %s %s %s %s %s' \
    "$dataset" "$read_time" "$reorder_time" "$fwrite_time" "$write_time" "$max_rss_kb")
  echo "$line" | tee -a "$RESULTS_FILE"
}

# Run RCM reordering following the requested pattern:
# echo "~~~~~~~~~~~~~~~~~~~~~~~~~<DATASET: <name>-RCM>~~~~~~~~~~~~~~~~~~~~~~~~~"
# mtxreorder --verbose <DATA_PATH>/<name>/<name>.mtx --ordering=rcm \
#   --rowperm-path=<DATA_PATH>/reordering/rcm_order/<name>.rcmorder \
#   > <DATA_PATH>/reordering/rcm_order/<name>_rcm.mtx
run_rcm_on_datasets() {
  local data_path="$1"; shift
  local datasets=("$@")
  local output

  mkdir -p "${RCM_DATA_PATH}"
  if [[ ! -f "$RESULTS_FILE" ]]; then
    printf '%s %s %s %s %s %s\n' \
      "dataset" "read" "reorder" "fwrite" "write" "MaxRSS(kB)" > "$RESULTS_FILE"
  fi

  for dataset in "${datasets[@]}"; do
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~<DATASET: ${dataset}-RCM>~~~~~~~~~~~~~~~~~~~~~~~~~"
    output=$(
      /usr/bin/time -v ${BIN} "${data_path}/${dataset}/${dataset}.mtx" \
        --ordering=rcm \
        --rowperm-path="${RCM_DATA_PATH}/${dataset}.rcmorder" \
        2>&1 1> "${RCM_DATA_PATH}/${dataset}_rcm.mtx"
    )
    printf '%s\n' "$output"
    record_timing "$dataset" "$output"
    echo ""
  done
}

echo "RCM reordering ..."
echo "Timing results will be appended to: ${RESULTS_FILE}"
# run_rcm_on_datasets "${DATA_PATH}" "${initial_datasets[@]}"
# run_rcm_on_datasets "${DATA_PATH}" "${ispass_datasets[@]}"
# run_rcm_on_datasets "${DATA_PATH}" "${new_datasets[@]}"
run_rcm_on_datasets "${DATA_PATH}" "${fig10_datasesets[@]}"