#!/bin/bash

export OMP_PROC_BIND=spread
export OMP_PLACES=cores

#res 1
./a.out "${DATA_PATH}/webbase-1M/webbase-1M.mtx" cutpart 56 quality 1 1

#3
./a.out "${DATA_PATH}/mip1/mip1.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/mouse_gene/mouse_gene.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/BenElechi1/BenElechi1.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/atmosmodl/atmosmodl.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/12month1/12month1.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/pkustk14/pkustk14.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/great-britain_osm/great-britain_osm.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/hugebubbles-00000/hugebubbles-00000.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/hugebubbles-00010/hugebubbles-00010.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/hugetrace-00000/hugetrace-00000.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/hugetrace-00010/hugetrace-00010.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/hugetrace-00010/hugetrace-00010.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/hugetric-00000/hugetric-00000.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/hugetric-00020/hugetric-00020.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/italy_osm/italy_osm.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/kron_g500-logn17/kron_g500-logn17.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/kron_g500-logn18/kron_g500-logn18.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/kron_g500-logn19/kron_g500-logn19.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/kron_g500-logn20/kron_g500-logn20.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/rgg_n_2_20_s0/rgg_n_2_20_s0.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/rgg_n_2_21_s0/rgg_n_2_21_s0.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/rgg_n_2_22_s0/rgg_n_2_22_s0.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/rgg_n_2_23_s0/rgg_n_2_23_s0.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/road_central/road_central.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/road_usa/road_usa.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/asia_osm/asia_osm.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/coPapersCiteseer/coPapersCiteseer.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/coPapersDBLP/coPapersDBLP.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/delaunay_n21/delaunay_n21.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/delaunay_n22/delaunay_n22.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/delaunay_n23/delaunay_n23.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/germany_osm/germany_osm.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/halfb/halfb.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/dielFilterV3real/dielFilterV3real.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/kmer_V1r/kmer_V1r.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/bmw3_2/bmw3_2.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/audikw_1/audikw_1.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/Hardesty3/Hardesty3.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/JP/JP.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/msdoor/msdoor.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/bibd_22_8/bibd_22_8.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/GL7d19/GL7d19.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/relat9/relat9.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/kim2/kim2.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/F1/F1.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/uk-2002/uk-2002.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/fem_hifreq_circuit/fem_hifreq_circuit.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/bundle_adj/bundle_adj.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/tp-6/tp-6.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/spal_004/spal_004.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/mycielskian18/mycielskian18.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/mycielskian20/mycielskian20.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/nd24k/nd24k.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/torso1/torso1.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/bone010/bone010.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/Ga41As41H72/Ga41As41H72.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/gearbox/gearbox.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/thermal2/thermal2.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/3Dspectralwave/3Dspectralwave.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/higgs-twitter/higgs-twitter.mtx" cutpart 56 quality 1 1
# ./a.out "${DATA_PATH}/twitter7/twitter7.mtx" cutpart 56 quality 1 1
# ./a.out "${DATA_PATH}/com-Friendster/com-Friendster.mtx" cutpart 56 quality 1 1
# ./a.out "${DATA_PATH}/MOLIERE_2016/MOLIERE_2016.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/TSOPF_RS_b2383_c1/TSOPF_RS_b2383_c1.mtx" cutpart 56 quality 1 1