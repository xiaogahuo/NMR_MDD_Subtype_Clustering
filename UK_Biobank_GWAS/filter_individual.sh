#!/bin/bash
# conda activate GWAS
# nohup ./filter_individual.sh  {exp_name} > {exp_name}_filter_individual.log &
    exp_name=$1
    if [ -d "./${exp_name}" ]; then
       rm -rf "./${exp_name}"
       echo "删除文件夹 ${exp_name} ！"
    else
       echo "文件夹 ${exp_name} 不存在！"
    fi
    mkdir -p "./${exp_name}"
    echo "文件夹 ${exp_name} 创建成功！"

    cd ./raw_data
    echo "Start filter Individuals on each Chrome"
    for ((chr=1; chr<=22; chr++))
    do
      plink --bfile "c${chr}"  --keep "../${exp_name}.csv" --make-bed --out "selIndivi_c${chr}"
    done
    plink --bfile "cX"  --keep "../${exp_name}.csv" --make-bed --out "selIndivi_cX"
    echo "Done: filter Individuals on each Chrome"

    mv -f selIndivi* "../${exp_name}"
    echo "Done: move all selIndivi results to ${exp_name}"
    
    echo "copy allchromes.txt to ${exp_name}"
    cp -f "../allchromes.txt" "../${exp_name}"

    echo "start merge all chromes on ${exp_name}"
    cd "../${exp_name}"
    plink --merge-list allchromes.txt --make-bed --out "${exp_name}"
    rm -rf *.log
    rm -rf *.nosex
    echo "Done: merge chromes all chromes on ${exp_name}"

    if [ -d "./GWAS_QC" ]; then
       rm -rf "./GWAS_QC"
       echo "删除文件夹 GWAS_QC ！"
    else
       mkdir -p "./GWAS_QC"
       echo "文件夹 GWAS_QC 创建成功！"
    fi
    mv -f "${exp_name}".*  "./GWAS_QC"
    # mv -f "./GWAS_QC/${exp_name}.fam" "./GWAS_QC/${exp_name}.fam.bak"
    # cp -f "../${exp_name}.fam" "./GWAS_QC/${exp_name}.fam"
    exit 0
