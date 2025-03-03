#!/bin/bash

exp_name="GWAS_NMR_MDD_Subtype_Clustering"
f_names=("1_VS_3" "1_VS_2" "3_VS_2" "0_VS_123" "0_VS_1" "0_VS_2" "0_VS_3")

for file_name in "${f_names[@]}"; do
    echo "${file_name}"
    ./filter_individual.sh  "${exp_name}" "${file_name}"
    ./run_GWAS_QC.sh "${exp_name}" "${file_name}"
done

