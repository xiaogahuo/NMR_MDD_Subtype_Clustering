#!/bin/bash
# conda activate GWAS
# ./run_GWAS_QC.sh {exp_name} {file_name}
   exp_name=$1
   file_name=$2

#cp  -f "./GWAS_Script/hist_miss.R" "./${exp_name}/${file_name}/"
#cp  -f "./GWAS_Script/gender_check.R" "./${exp_name}/${file_name}/"
#cp  -f "./GWAS_Script/MAF_check.R" "./${exp_name}/${file_name}/"
#cp  -f "./GWAS_Script/hwe.R" "./${exp_name}/${file_name}/"
#cp  -f "./GWAS_Script/check_heterozygosity_rate.R" "./${exp_name}/${file_name}/"
#cp  -f "./GWAS_Script/heterozygosity_outliers_list.R" "./${exp_name}/${file_name}/"
#cp  -f "./GWAS_Script/Relatedness.R" "./${exp_name}/${file_name}/"
#cp  -f "./GWAS_Script/inversion.txt" "./${exp_name}/${file_name}/"


cd "./${exp_name}/${file_name}/"

## Investigate missingness per individual and per SNP and make histograms.
#plink --bfile "${file_name}_1" --missing
#Rscript --no-save hist_miss.R

## Delete SNPs with missingness >0.2.
#plink --bfile "${file_name}_1" --geno 0.2 --make-bed --out "${file_name}_2"
## Delete individuals with missingness >0.2.
#plink --bfile "${file_name}_2" --mind 0.2 --make-bed --out "${file_name}_3"

# Delete SNPs with missingness >0.02.
plink --bfile "${file_name}_1" --geno 0.02 --make-bed --out "${file_name}_4"
## Delete individuals with missingness >0.02.
#plink --bfile "${file_name}_4" --mind 0.02 --make-bed --out "${file_name}_5"

## Generate a plot of the MAF distribution.
#plink --bfile "${file_name}_5" --freq --out MAF_check
#Rscript --no-save MAF_check.R

# Remove SNPs with a low MAF frequency.
plink --bfile "${file_name}_4" --maf 0.05 --make-bed --out "${file_name}_8"
mv -f "${file_name}_8.fam" "${file_name}_8.fam.backup"

