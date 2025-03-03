#!/bin/bash
numbers=(2 8 11 18)
for number in "${numbers[@]}"; do
#for((i=7;i<=21;i++)); do
  plink2 --bgen "raw_data/ukb22828_c${number}_b0_v3.bgen" ref-first \
       --sample raw_data/ukb22828_c1_b0_v3_s487150.sample \
       --extract "final_result/filter_rs_txt/CHR_${number}.txt" \
       --make-bed --out "final_result/filter_rs_bed/chr_${number}"
done