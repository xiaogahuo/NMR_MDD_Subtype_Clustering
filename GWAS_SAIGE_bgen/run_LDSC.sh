#!/bin/bash
# https://www.cnblogs.com/chenwenyan/p/11321272.html
# https://cloud.tencent.com/developer/article/1670338
# https://zenodo.org/records/10515792

#exp_name="GWAS_NMR_MDD_Subtype_Clustering_bgen"
#f_names=("1_VS_3" "1_VS_2" "3_VS_2" "0_VS_123" "0_VS_1" "0_VS_2" "0_VS_3")
exp_name=$1
file_name=$2
Num=$3

echo "${exp_name}/${file_name}"

/home/user/ldsc/munge_sumstats.py \
    --sumstats "/data/UKBioBank/bgen_Data/${exp_name}/${file_name}/SAIGE_results.ldsc.txt" \
    --N ${Num} \
    --merge-alleles "/data/UKBioBank/LDScore_EUR/w_hm3.snplist" \
    --out "/data/UKBioBank/bgen_Data/${exp_name}/${file_name}/scz"

/home/user/ldsc/ldsc.py \
    --h2 "/data/UKBioBank/bgen_Data/${exp_name}/${file_name}/scz.sumstats.gz" \
    --ref-ld-chr "/data/UKBioBank/LDScore_EUR/LDscore/" \
    --w-ld-chr "/data/UKBioBank/LDScore_EUR/LDscore/" \
    --out "/data/UKBioBank/bgen_Data/${exp_name}/${file_name}/scz_h2"


sumstats_107_file_name="/data/UKBioBank/LDScore_EUR/sumstats_107/traits_indep107.txt"
# 构建路径
out_path="/data/UKBioBank/bgen_Data/${exp_name}/${file_name}/LDGC/"
# 判断路径是否存在，如果不存在则创建
if [ ! -d "${out_path}" ]; then
    echo "路径 ${out_path} 不存在，正在创建..."
    mkdir -p "${out_path}"
else
    echo "路径 ${out_path} 已存在。"
fi

# 按行读取文件
while IFS= read -r line; do
    echo "${line}"
    /home/user/ldsc/ldsc.py \
        --rg "/data/UKBioBank/LDScore_EUR/sumstats_107/${line}.sumstats.gz","/data/UKBioBank/bgen_Data/${exp_name}/${file_name}/scz.sumstats.gz" \
        --ref-ld-chr "/data/UKBioBank/LDScore_EUR/LDscore/" \
        --w-ld-chr "/data/UKBioBank/LDScore_EUR/LDscore/" \
        --out "${out_path}${line}"
done < "${sumstats_107_file_name}"