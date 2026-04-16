#!/bin/bash

export CLUSTERWISE_ROOT=/data/lsl/SpGEMM/SpGEMM-exp/clusterwise-spgemm-main

export DATA_PATH=/data/suitesparse_collection

export SHUFFLED_DATA_PATH=/data2/linshengle_data/SpGEMM-Reordering/shuffle_order

export ND_DATA_PATH=/data2/linshengle_data/SpGEMM-Reordering/nd_order

export RCM_DATA_PATH=/data2/linshengle_data/SpGEMM-Reordering/rcm_order

export AMD_DATA_PATH=/data2/linshengle_data/SpGEMM-Reordering/amd_order

export GRAY_DATA_PATH=/data2/linshengle_data/SpGEMM-Reordering/gray_order

export DEGREE_DATA_PATH=/data2/linshengle_data/SpGEMM-Reordering/degree_order

export SLASHBURN_DATA_PATH=/data2/linshengle_data/SpGEMM-Reordering/slashburn_order

export GP_DATA_PATH=/data2/linshengle_data/SpGEMM-Reordering/gp_order

export HP_DATA_PATH=/data2/linshengle_data/SpGEMM-Reordering/hp_order

export RABBIT_DATA_PATH=/data2/linshengle_data/SpGEMM-Reordering/rabbit_order


# bash rowwise/1_rowwise_original.sh > rowwise/1_rowwise_original.out
bash rowwise/2_rowwise_random.sh > rowwise/2_rowwise_random.out
bash rowwise/3_rowwise_rabbit.sh > rowwise/3_rowwise_rabbit.out
bash rowwise/4_rowwise_amd.sh > rowwise/4_rowwise_amd.out
bash rowwise/5_rowwise_rcm.sh > rowwise/5_rowwise_rcm.out
bash rowwise/6_rowwise_nd.sh > rowwise/6_rowwise_nd.out
bash rowwise/7_rowwise_gp.sh > rowwise/7_rowwise_gp.out
bash rowwise/8_rowwise_hp.sh > rowwise/8_rowwise_hp.out
bash rowwise/9_rowwise_gray.sh > rowwise/9_rowwise_gray.out
bash rowwise/10_rowwise_degree.sh > rowwise/10_rowwise_degree.out
bash rowwise/11_rowwise_slashburn.sh > rowwise/11_rowwise_slashburn.out