#!/bin/bash

# exp_name="GWAS_NMR_MDD_Subtype_Clustering_bgen"
#f_names=("1_VS_3" "1_VS_2" "3_VS_2" "0_VS_123" "0_VS_1" "0_VS_2" "0_VS_3")
exp_name=$1
file_name=$2

# for file_name in "${f_names[@]}"; do
echo "${exp_name}/${file_name}"
cp  -f "./manhattan_QQ.py"   "./${exp_name}/${file_name}/"
cd "./${exp_name}/${file_name}/"
python "manhattan_QQ.py"
cd "../../"
# done