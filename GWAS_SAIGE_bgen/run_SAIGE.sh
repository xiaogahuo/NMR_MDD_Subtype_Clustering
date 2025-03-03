#!/bin/bash

#exp_name="GWAS_NMR_MDD_Subtype_Clustering_bgen"
#f_names=("1_VS_3" "1_VS_2" "3_VS_2" "0_VS_123" "0_VS_1" "0_VS_2" "0_VS_3")
exp_name=$1
file_name=$2

echo "${exp_name}/${file_name}"

Rscript /pro2/SAIGE-main/extdata/step1_fitNULLGLMM.R     \
      --plinkFile="/pro2/${exp_name}/${file_name}/${file_name}"  \
      --phenoFile="/pro2/${exp_name}/${file_name}/${file_name}_pheno.txt" \
      --phenoCol=final_cluster \
      --covarColList=Sex_cate,Age,genotyping_array,PCA1,PCA2,PCA3,PCA4,PCA5 \
      --qCovarColList=Sex_cate,genotyping_array \
      --sampleIDColinphenoFile=eid \
      --sexCol=Sex_cate \
      --MaleCode=1 \
      --FemaleCode=0 \
      --traitType=binary        \
      --outputPrefix="/pro2/${exp_name}/${file_name}/Step1_fitNULLGLMM" \
      --nThreads=80	\
      --IsOverwriteVarianceRatioFile=TRUE

Rscript /pro2/SAIGE-main/extdata/step2_SPAtests.R        \
      --bedFile="/pro2/${exp_name}/${file_name}/${file_name}.bed"       \
      --bimFile="/pro2/${exp_name}/${file_name}/${file_name}.bim"       \
      --famFile="/pro2/${exp_name}/${file_name}/${file_name}.fam"       \
      --AlleleOrder=alt-first \
      --SAIGEOutputFile="/pro2/${exp_name}/${file_name}/Step2_SPAtests.txt" \
      --minMAF=0 \
      --minMAC=20 \
      --GMMATmodelFile="/pro2/${exp_name}/${file_name}/Step1_fitNULLGLMM.rda" \
      --varianceRatioFile="/pro2/${exp_name}/${file_name}/Step1_fitNULLGLMM.varianceRatio.txt"   \
      --LOCO=FALSE \
      --is_output_moreDetails=TRUE

# for((i=1;i<=22;i++)); do
#  Rscript /pro2/SAIGE-main/extdata/step1_fitNULLGLMM.R     \
#      --plinkFile="/pro2/${exp_name}/${file_name}/selIndivi_c${i}"  \
#      --phenoFile="/pro2/${exp_name}/${file_name}/${file_name}_pheno.txt" \
#      --phenoCol=final_cluster \
#      --covarColList=Sex_cate,Age,genotyping_array,genotyping_batch,PCA1,PCA2,PCA3,PCA4,PCA5 \
#      --qCovarColList=Sex_cate,genotyping_array,genotyping_batch \
#      --sampleIDColinphenoFile=eid \
#      --sexCol=Sex_cate \
#      --MaleCode=1 \
#      --FemaleCode=0 \
#      --traitType=binary        \
#      --outputPrefix="/pro2/${exp_name}/${file_name}/Step1_fitNULLGLMM_c${i}" \
#      --nThreads=80	\
#      --IsOverwriteVarianceRatioFile=TRUE

#  Rscript /pro2/SAIGE-main/extdata/step2_SPAtests.R        \
#      --bedFile="/pro2/${exp_name}/${file_name}/selIndivi_c${i}.bed"       \
#      --bimFile="/pro2/${exp_name}/${file_name}/selIndivi_c${i}.bim"       \
#      --famFile="/pro2/${exp_name}/${file_name}/selIndivi_c${i}.fam"       \
#      --AlleleOrder=alt-first \
#      --SAIGEOutputFile="/pro2/${exp_name}/${file_name}/Step2_SPAtests_c${i}.txt" \
#      --chrom=${i}       \
#      --minMAF=0 \
#      --minMAC=20 \
#      --GMMATmodelFile="/pro2/${exp_name}/${file_name}/Step1_fitNULLGLMM_c${i}.rda" \
#      --varianceRatioFile="/pro2/${exp_name}/${file_name}/Step1_fitNULLGLMM_c${i}.varianceRatio.txt"   \
#      --LOCO=FALSE \
#      --is_output_moreDetails=TRUE
# done