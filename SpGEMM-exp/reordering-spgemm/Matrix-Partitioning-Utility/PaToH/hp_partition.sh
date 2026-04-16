#!/bin/bash

export OMP_PROC_BIND=spread
export OMP_PLACES=cores

#1
./a.out "${DATA_PATH}/2cubes_sphere/2cubes_sphere.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/cage12/cage12.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/cage15/cage15.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/cant/cant.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/com-Amazon/com-Amazon.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/conf5_4-8x8-05/conf5_4-8x8-05.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/consph/consph.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/cop20k_A/cop20k_A.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/delaunay_n24/delaunay_n24.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/filter3D/filter3D.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/Freescale2/Freescale2.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/hood/hood.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/kmer_V1r/kmer_V1r.mtx" cutpart 56 quality 1 1 #卡很久
./a.out "${DATA_PATH}/m133-b3/m133-b3.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/mac_econ_fwd500/mac_econ_fwd500.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/majorbasis/majorbasis.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/mario002/mario002.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/mc2depi/mc2depi.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/mono_500Hz/mono_500Hz.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/offshore/offshore.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/patents_main/patents_main.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/pdb1HYS/pdb1HYS.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/poisson3Da/poisson3Da.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/pwtk/pwtk.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/rma10/rma10.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/scircuit/scircuit.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/shipsec1/shipsec1.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/wb-edu/wb-edu.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/webbase-1M/webbase-1M.mtx" cutpart 56 quality 1 1

#2
./a.out "${DATA_PATH}/333SP/333SP.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/adaptive/adaptive.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/af_shell10/af_shell10.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/AS365/AS365.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/as-Skitter/as-Skitter.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/channel-500x100x100-b050/channel-500x100x100-b050.mtx" cutpart 56 quality 1 1
#./a.out "${DATA_PATH}/circuit5M/circuit5M.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/cit-Patents/cit-Patents.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/com-LiveJournal/com-LiveJournal.mtx" cutpart 56 quality 1 1
#./a.out "${DATA_PATH}/com-Orkut/com-Orkut.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/CurlCurl_4/CurlCurl_4.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/europe_osm/europe_osm.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/G3_circuit/G3_circuit.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/GAP-road/GAP-road.mtx" cutpart 56 quality 1 1
#./a.out "${DATA_PATH}/GAP-twitter/GAP-twitter.mtx" cutpart 56 quality 1 1 # killed
./a.out "${DATA_PATH}/hugebubbles-00020/hugebubbles-00020.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/hugetrace-00020/hugetrace-00020.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/hugetric-00010/hugetric-00010.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/HV15R/HV15R.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/kkt_power/kkt_power.mtx" cutpart 56 quality 1 1
#./a.out "${DATA_PATH}/kron_g500-logn21/kron_g500-logn21.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/M6/M6.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/mawi_201512020330/mawi_201512020330.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/nlpkkt240/nlpkkt240.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/NLR/NLR.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/packing-500x100x100-b050/packing-500x100x100-b050.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/patents/patents.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/Queen_4147/Queen_4147.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/rajat31/rajat31.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/rgg_n_2_24_s0/rgg_n_2_24_s0.mtx" cutpart 56 quality 1 1
#./a.out "${DATA_PATH}/sk-2005/sk-2005.mtx" cutpart 56 quality 1 1 # killed
./a.out "${DATA_PATH}/soc-LiveJournal1/soc-LiveJournal1.mtx" cutpart 56 quality 1 1
#./a.out "${DATA_PATH}/soc-Pokec/soc-Pokec.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/stokes/stokes.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/sx-stackoverflow/sx-stackoverflow.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/venturiLevel3/venturiLevel3.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/wikipedia-20070206/wikipedia-20070206.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/wiki-Talk/wiki-Talk.mtx" cutpart 56 quality 1 1
./a.out "${DATA_PATH}/wiki-topcats/wiki-topcats.mtx" cutpart 56 quality 1 1

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
# ./a.out "${DATA_PATH}/kmer_V1r/kmer_V1r.mtx" cutpart 56 quality 1 1
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
