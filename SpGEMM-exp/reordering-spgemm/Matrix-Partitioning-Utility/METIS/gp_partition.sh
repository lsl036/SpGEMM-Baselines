#!/bin/bash

export OMP_PROC_BIND=spread
export OMP_PLACES=cores

#1
./run -i "${DATA_PATH}/2cubes_sphere/2cubes_sphere.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/cage12/cage12.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/cage15/cage15.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/cant/cant.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/com-Amazon/com-Amazon.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/conf5_4-8x8-05/conf5_4-8x8-05.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/consph/consph.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/cop20k_A/cop20k_A.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/delaunay_n24/delaunay_n24.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/filter3D/filter3D.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/Freescale2/Freescale2.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/hood/hood.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/kmer_V1r/kmer_V1r.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/m133-b3/m133-b3.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/mac_econ_fwd500/mac_econ_fwd500.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/majorbasis/majorbasis.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/mario002/mario002.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/mc2depi/mc2depi.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/mono_500Hz/mono_500Hz.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/offshore/offshore.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/patents_main/patents_main.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/pdb1HYS/pdb1HYS.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/poisson3Da/poisson3Da.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/pwtk/pwtk.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/rma10/rma10.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/scircuit/scircuit.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/shipsec1/shipsec1.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/wb-edu/wb-edu.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/webbase-1M/webbase-1M.mtx" -k 56 -o edge-cut -t -s

#2
./run -i "${DATA_PATH}/333SP/333SP.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/adaptive/adaptive.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/af_shell10/af_shell10.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/AS365/AS365.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/as-Skitter/as-Skitter.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/channel-500x100x100-b050/channel-500x100x100-b050.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/circuit5M/circuit5M.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/cit-Patents/cit-Patents.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/com-LiveJournal/com-LiveJournal.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/com-Orkut/com-Orkut.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/CurlCurl_4/CurlCurl_4.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/europe_osm/europe_osm.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/G3_circuit/G3_circuit.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/GAP-road/GAP-road.mtx" -k 56 -o edge-cut -t -s
# ./run -i "${DATA_PATH}/GAP-twitter/GAP-twitter.mtx" -k 56 -o edge-cut -t -s # killed
./run -i "${DATA_PATH}/hugebubbles-00020/hugebubbles-00020.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/hugetrace-00020/hugetrace-00020.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/hugetric-00010/hugetric-00010.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/HV15R/HV15R.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/kkt_power/kkt_power.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/kron_g500-logn21/kron_g500-logn21.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/M6/M6.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/mawi_201512020330/mawi_201512020330.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/nlpkkt240/nlpkkt240.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/NLR/NLR.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/packing-500x100x100-b050/packing-500x100x100-b050.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/patents/patents.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/Queen_4147/Queen_4147.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/rajat31/rajat31.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/rgg_n_2_24_s0/rgg_n_2_24_s0.mtx" -k 56 -o edge-cut -t -s
#./run -i "${DATA_PATH}/sk-2005/sk-2005.mtx" -k 56 -o edge-cut -t -s # killed
./run -i "${DATA_PATH}/soc-LiveJournal1/soc-LiveJournal1.mtx" -k 56 -o edge-cut -t -s
#./run -i "${DATA_PATH}/soc-Pokec/soc-Pokec.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/stokes/stokes.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/sx-stackoverflow/sx-stackoverflow.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/venturiLevel3/venturiLevel3.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/wikipedia-20070206/wikipedia-20070206.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/wiki-Talk/wiki-Talk.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/wiki-topcats/wiki-topcats.mtx" -k 56 -o edge-cut -t -s

#3
./run -i "${DATA_PATH}/mip1/mip1.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/mouse_gene/mouse_gene.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/BenElechi1/BenElechi1.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/atmosmodl/atmosmodl.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/12month1/12month1.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/pkustk14/pkustk14.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/great-britain_osm/great-britain_osm.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/hugebubbles-00000/hugebubbles-00000.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/hugebubbles-00010/hugebubbles-00010.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/hugetrace-00000/hugetrace-00000.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/hugetrace-00010/hugetrace-00010.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/hugetrace-00010/hugetrace-00010.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/hugetric-00000/hugetric-00000.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/hugetric-00020/hugetric-00020.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/italy_osm/italy_osm.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/kron_g500-logn17/kron_g500-logn17.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/kron_g500-logn18/kron_g500-logn18.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/kron_g500-logn19/kron_g500-logn19.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/kron_g500-logn20/kron_g500-logn20.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/rgg_n_2_20_s0/rgg_n_2_20_s0.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/rgg_n_2_21_s0/rgg_n_2_21_s0.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/rgg_n_2_22_s0/rgg_n_2_22_s0.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/rgg_n_2_23_s0/rgg_n_2_23_s0.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/road_central/road_central.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/road_usa/road_usa.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/asia_osm/asia_osm.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/coPapersCiteseer/coPapersCiteseer.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/coPapersDBLP/coPapersDBLP.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/delaunay_n21/delaunay_n21.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/delaunay_n22/delaunay_n22.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/delaunay_n23/delaunay_n23.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/germany_osm/germany_osm.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/halfb/halfb.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/dielFilterV3real/dielFilterV3real.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/kmer_V1r/kmer_V1r.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/bmw3_2/bmw3_2.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/audikw_1/audikw_1.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/Hardesty3/Hardesty3.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/JP/JP.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/msdoor/msdoor.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/bibd_22_8/bibd_22_8.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/GL7d19/GL7d19.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/relat9/relat9.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/kim2/kim2.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/F1/F1.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/uk-2002/uk-2002.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/fem_hifreq_circuit/fem_hifreq_circuit.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/bundle_adj/bundle_adj.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/tp-6/tp-6.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/spal_004/spal_004.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/mycielskian18/mycielskian18.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/mycielskian20/mycielskian20.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/nd24k/nd24k.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/torso1/torso1.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/bone010/bone010.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/Ga41As41H72/Ga41As41H72.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/gearbox/gearbox.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/thermal2/thermal2.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/3Dspectralwave/3Dspectralwave.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/higgs-twitter/higgs-twitter.mtx" -k 56 -o edge-cut -t -s
# ./run -i "${DATA_PATH}/twitter7/twitter7.mtx" -k 56 -o edge-cut -t -s
# ./run -i "${DATA_PATH}/com-Friendster/com-Friendster.mtx" -k 56 -o edge-cut -t -s
# ./run -i "${DATA_PATH}/MOLIERE_2016/MOLIERE_2016.mtx" -k 56 -o edge-cut -t -s
./run -i "${DATA_PATH}/TSOPF_RS_b2383_c1/TSOPF_RS_b2383_c1.mtx" -k 56 -o edge-cut -t -s
