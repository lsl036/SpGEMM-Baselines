#!/bin/bash

#export LD_LIBRARY_PATH="${RABBIT_HOME}/gperftools-2.9.1/lib":$LD_LIBRARY_PATH

# 转换工具路径（如果不在当前目录，请修改此路径）
MTX_TO_EL="${RABBIT_HOME:-.}/converter/mtx_to_el"
if [ ! -f "$MTX_TO_EL" ] && [ -f "./converter/mtx_to_el" ]; then
    MTX_TO_EL="./converter/mtx_to_el"
fi

# 函数：检查并转换 .mtx 到 .el（如果需要）
# 参数：数据集名称（例如 "2cubes_sphere"）
convert_mtx_to_el_if_needed() {
    local dataset="$1"
    local el_file="${DATA_PATH}/${dataset}/${dataset}.el"
    local mtx_file="${DATA_PATH}/${dataset}/${dataset}.mtx"
    
    # 如果 .el 文件已存在，直接返回
    if [ -f "$el_file" ]; then
        return 0
    fi
    
    # 如果 .el 不存在但 .mtx 存在，进行转换
    if [ -f "$mtx_file" ]; then
        echo "  [转换] ${dataset}.mtx -> ${dataset}.el"
        if [ -f "$MTX_TO_EL" ]; then
            "$MTX_TO_EL" "$mtx_file" "$el_file"
            if [ $? -eq 0 ]; then
                echo "  [完成] 转换成功"
                return 0
            else
                echo "  [错误] 转换失败"
                return 1
            fi
        else
            echo "  [错误] 找不到转换工具 mtx_to_el，请先编译 converter/mtx_to_el.cc"
            echo "  编译命令: g++ converter/mtx_to_el.cc -o converter/mtx_to_el -std=c++11"
            return 1
        fi
    else
        echo "  [错误] 找不到 ${dataset}.mtx 或 ${dataset}.el 文件"
        return 1
    fi
}

# 函数：处理单个数据集的重排序
# 参数：数据集名称
run_reorder() {
    local dataset="$1"
    convert_mtx_to_el_if_needed "$dataset" && \
    ./reorder "${DATA_PATH}/${dataset}/${dataset}.el" \
              "${RABBIT_DATA_PATH}/${dataset}.rabbitorder" \
              "${RABBIT_DATA_PATH}/${dataset}.off" -s
}

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<2cubes_sphere>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "2cubes_sphere"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<cage12>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "cage12"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<cage15>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "cage15"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<cant>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "cant"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<com-Amazon>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "com-Amazon"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<conf5_4-8x8-05>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "conf5_4-8x8-05"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<consph>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "consph"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<cop20k_A>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "cop20k_A"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<delaunay_n24>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "delaunay_n24"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<filter3D>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "filter3D"

#echo "~~~~~~~~~~~~~~~~~~~~~~~~~<Freescale2>~~~~~~~~~~~~~~~~~~~~~~~~~"
#./reorder "${DATA_PATH}/Freescale2/Freescale2.el" "${RABBIT_DATA_PATH}/Freescale2.rabbitorder" "${RABBIT_DATA_PATH}/Freescale2.off"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<hood>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "hood"

#echo "~~~~~~~~~~~~~~~~~~~~~~~~~<kmer_V1r>~~~~~~~~~~~~~~~~~~~~~~~~~"
#./reorder "${DATA_PATH}/kmer_V1r/kmer_V1r.el" "${RABBIT_DATA_PATH}/kmer_V1r.rabbitorder" "${RABBIT_DATA_PATH}/kmer_V1r.off"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<m133-b3>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "m133-b3"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<mac_econ_fwd500>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "mac_econ_fwd500"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<majorbasis>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "majorbasis"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<mario002>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "mario002"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<mc2depi>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "mc2depi"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<mono_500Hz>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "mono_500Hz"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<offshore>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "offshore"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<patents_main>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "patents_main"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<pdb1HYS>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "pdb1HYS"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<poisson3Da>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "poisson3Da"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<pwtk>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "pwtk"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<rma10>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "rma10"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<scircuit>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "scircuit"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<shipsec1>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "shipsec1"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<wb-edu>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "wb-edu"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<webbase-1M>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "webbase-1M"


###ISPASS-2023 datasets
echo "~~~~~~~~~~~~~~~~~~~~~~~~~<333SP>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "333SP"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<adaptive>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "adaptive"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<af_shell10>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "af_shell10"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<AS365>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "AS365"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<as-Skitter>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "as-Skitter"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<channel-500x100x100-b050>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "channel-500x100x100-b050"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<circuit5M>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "circuit5M"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<cit-Patents>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "cit-Patents"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<com-LiveJournal>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "com-LiveJournal"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<com-Orkut>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "com-Orkut"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<CurlCurl_4>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "CurlCurl_4"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<europe_osm>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "europe_osm"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<G3_circuit>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "G3_circuit"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<GAP-road>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "GAP-road"

#echo "~~~~~~~~~~~~~~~~~~~~~~~~~<GAP-twitter>~~~~~~~~~~~~~~~~~~~~~~~~~"
#./reorder "${DATA_PATH}/GAP-twitter/GAP-twitter.el" "${RABBIT_DATA_PATH}/GAP-twitter.rabbitorder" "${RABBIT_DATA_PATH}/GAP-twitter.off"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<hugebubbles-00020>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "hugebubbles-00020"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<hugetrace-00020>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "hugetrace-00020"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<hugetric-00010>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "hugetric-00010"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<HV15R>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "HV15R"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<kkt_power>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "kkt_power"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<kron_g500-logn21>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "kron_g500-logn21"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<M6>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "M6"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<mawi_201512020330>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "mawi_201512020330"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<nlpkkt240>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "nlpkkt240"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<NLR>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "NLR"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<packing-500x100x100-b050>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "packing-500x100x100-b050"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<patents>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "patents"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<Queen_4147>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "Queen_4147"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<rajat31>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "rajat31"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<rgg_n_2_24_s0>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "rgg_n_2_24_s0"

#echo "~~~~~~~~~~~~~~~~~~~~~~~~~<sk-2005>~~~~~~~~~~~~~~~~~~~~~~~~~"
#./reorder "${DATA_PATH}/sk-2005/sk-2005.el" "${RABBIT_DATA_PATH}/sk-2005.rabbitorder" "${RABBIT_DATA_PATH}/sk-2005.off"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<soc-LiveJournal1>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "soc-LiveJournal1"

#echo "~~~~~~~~~~~~~~~~~~~~~~~~~<soc-Pokec>~~~~~~~~~~~~~~~~~~~~~~~~~"
#./reorder "${DATA_PATH}/soc-Pokec/soc-Pokec.el" "${RABBIT_DATA_PATH}/soc-Pokec.rabbitorder" "${RABBIT_DATA_PATH}/soc-Pokec.off"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<stokes>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "stokes"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<sx-stackoverflow>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "sx-stackoverflow"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<venturiLevel3>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "venturiLevel3"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<wikipedia-20070206>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "wikipedia-20070206"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<wiki-Talk>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "wiki-Talk"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<wiki-topcats>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "wiki-topcats"


### new datasets
# dataset from the following criteria:
#   - Square matrices with more than 8 Million nnz and less than 10 Billion nnz
#   - This suites with our L2 cache size of 64MB and memory size of 512GB
echo "~~~~~~~~~~~~~~~~~~~~~~~~~<mip1>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "mip1"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<mouse_gene>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "mouse_gene"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<BenElechi1>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "BenElechi1"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<atmosmodl>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "atmosmodl"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<12month1>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "12month1"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<pkustk14>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "pkustk14"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<great-britain_osm>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "great-britain_osm"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<hugebubbles-00000>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "hugebubbles-00000"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<hugebubbles-00010>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "hugebubbles-00010"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<hugetrace-00000>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "hugetrace-00000"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<hugetrace-00010>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "hugetrace-00010"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<hugetrace-00010>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "hugetrace-00010"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<hugetric-00000>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "hugetric-00000"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<hugetric-00020>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "hugetric-00020"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<italy_osm>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "italy_osm"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<kron_g500-logn17>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "kron_g500-logn17"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<kron_g500-logn18>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "kron_g500-logn18"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<kron_g500-logn19>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "kron_g500-logn19"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<kron_g500-logn20>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "kron_g500-logn20"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<rgg_n_2_20_s0>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "rgg_n_2_20_s0"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<rgg_n_2_21_s0>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "rgg_n_2_21_s0"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<rgg_n_2_22_s0>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "rgg_n_2_22_s0"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<rgg_n_2_23_s0>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "rgg_n_2_23_s0"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<road_central>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "road_central"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<road_usa>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "road_usa"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<asia_osm>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "asia_osm"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<coPapersCiteseer>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "coPapersCiteseer"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<coPapersDBLP>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "coPapersDBLP"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<delaunay_n21>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "delaunay_n21"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<delaunay_n22>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "delaunay_n22"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<delaunay_n23>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "delaunay_n23"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<germany_osm>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "germany_osm"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<halfb>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "halfb"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<dielFilterV3real>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "dielFilterV3real"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<kmer_V1r>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "kmer_V1r"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<bmw3_2>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "bmw3_2"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<audikw_1>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "audikw_1"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<Hardesty3>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "Hardesty3"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<JP>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "JP"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<msdoor>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "msdoor"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<bibd_22_8>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "bibd_22_8"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<GL7d19>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "GL7d19"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<relat9>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "relat9"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<kim2>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "kim2"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<F1>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "F1"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<uk-2002>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "uk-2002"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<fem_hifreq_circuit>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "fem_hifreq_circuit"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<bundle_adj>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "bundle_adj"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<tp-6>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "tp-6"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<spal_004>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "spal_004"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<mycielskian18>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "mycielskian18"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<mycielskian20>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "mycielskian20"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<nd24k>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "nd24k"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<torso1>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "torso1"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<bone010>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "bone010"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<Ga41As41H72>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "Ga41As41H72"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<gearbox>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "gearbox"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<thermal2>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "thermal2"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<3Dspectralwave>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "3Dspectralwave"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<higgs-twitter>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "higgs-twitter"

#echo "~~~~~~~~~~~~~~~~~~~~~~~~~<twitter7>~~~~~~~~~~~~~~~~~~~~~~~~~"
#./reorder "${DATA_PATH}/twitter7/twitter7.el" "${RABBIT_DATA_PATH}/twitter7.rabbitorder" "${RABBIT_DATA_PATH}/twitter7.off"
#
#echo "~~~~~~~~~~~~~~~~~~~~~~~~~<com-Friendster>~~~~~~~~~~~~~~~~~~~~~~~~~"
#./reorder "${DATA_PATH}/com-Friendster/com-Friendster.el" "${RABBIT_DATA_PATH}/com-Friendster.rabbitorder" "${RABBIT_DATA_PATH}/com-Friendster.off"
#
#echo "~~~~~~~~~~~~~~~~~~~~~~~~~<MOLIERE_2016>~~~~~~~~~~~~~~~~~~~~~~~~~"
#./reorder "${DATA_PATH}/MOLIERE_2016/MOLIERE_2016.el" "${RABBIT_DATA_PATH}/MOLIERE_2016.rabbitorder" "${RABBIT_DATA_PATH}/MOLIERE_2016.off"

echo "~~~~~~~~~~~~~~~~~~~~~~~~~<TSOPF_RS_b2383_c1>~~~~~~~~~~~~~~~~~~~~~~~~~"
run_reorder "TSOPF_RS_b2383_c1"
