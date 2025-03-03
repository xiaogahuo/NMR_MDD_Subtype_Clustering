#!/bin/bash

#i=$1

 for((i=1;i<=22;i++)); do
  plink2 --bgen "raw_data/ukb22828_c${i}_b0_v3.bgen" ref-first \
       --sample raw_data/ukb22828_c1_b0_v3_s487150.sample \
       --extract "SNP_QC_raw_data/keep_snp_chr${i}.txt" \
       --make-bed --out "SNP_QC_raw_data/c${i}_b0_v3"
done