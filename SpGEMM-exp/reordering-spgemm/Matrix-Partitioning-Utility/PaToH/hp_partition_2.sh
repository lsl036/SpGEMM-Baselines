#!/bin/bash

export OMP_PROC_BIND=spread
export OMP_PLACES=cores

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