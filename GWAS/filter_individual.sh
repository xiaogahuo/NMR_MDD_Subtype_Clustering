#!/bin/bash
# conda activate GWAS
# ./filter_individual.sh {exp_name} {file_name}
    exp_name=$1
    file_name=$2

    cd ./SNP_QC_raw_data
    echo "Start filter Individuals on each Chrome"
    for ((chr=1; chr<=22; chr++))
    do
      plink --bfile "c${chr}_b0_v3"  --keep "../${exp_name}/${file_name}/${file_name}.csv" --make-bed --out "../${exp_name}/${file_name}/selIndivi_c${chr}"
      plink --bfile "../${exp_name}/${file_name}/selIndivi_c${chr}" --exclude /data/UKBioBank/bed_Data/exclude_snps.txt --make-bed --out "../${exp_name}/${file_name}/chr${chr}"
    done
    echo "Done: filter Individuals on each Chrome"
    echo "start merge all chromes on ${file_name}"
    cd "../${exp_name}/${file_name}/"
    plink --merge-list /data/UKBioBank/bed_Data/allchromes.txt --make-bed --out "${file_name}_1"
    rm -rf *.log *.nosex selIndivi* chr*.bed chr*.bim chr*.fam
    echo "Done: merge chromes all chromes on ${exp_name}/${file_name}/${file_name}_1"
    exit 0
