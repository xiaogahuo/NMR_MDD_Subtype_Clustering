#!/bin/bash

exp_name="GWAS_NMR_MDD_Subtype_Clustering_bgen"
f_names=("1_VS_3" "1_VS_2" "3_VS_2" "0_VS_123" "0_VS_1" "0_VS_2" "0_VS_3")
#f_names=("0_VS_123")

for file_name in "${f_names[@]}"; do
    echo "${exp_name}/${file_name}"
    echo "Start filter Individuals on each Chrome"
    cd "/data/UKBioBank/bgen_Data/"
    for ((chr=1; chr<=22; chr++))
    do
      awk '{print $2}' "./SNP_QC_raw_data/c${chr}_b0_v3.bim" | sort | uniq -d > "./${exp_name}/${file_name}/c${chr}_dup_snp.txt"
      plink2 --bfile "./SNP_QC_raw_data/c${chr}_b0_v3" --exclude "./${exp_name}/${file_name}/c${chr}_dup_snp.txt" --keep "./${exp_name}/${file_name}/${file_name}.csv" --make-bed --out "./${exp_name}/${file_name}/selIndivi_c${chr}"
    done
    cd "/data/UKBioBank/bgen_Data/${exp_name}/${file_name}/"
    plink2 --pmerge-list /data/UKBioBank/bgen_Data/allchromes.txt  --make-bed  --out "${file_name}"
    echo "Done: filter Individuals on each Chrome and Merge 22 Chromes"
done

