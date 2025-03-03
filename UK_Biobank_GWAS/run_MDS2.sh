    file_name=$1
################ Explanation of the main script ##########################

cd "./${file_name}/MDS/"

# Extract these individuals in HapMap data.
plink --bfile "${file_name}_11" --keep EUR_MDS_merge2 --make-bed --out "${file_name}_12"
# Note, since our HapMap data did include any ethnic outliers, no individuls were removed at this step. However, if our data would have included individuals outside of the thresholds we set, then these individuals would have been removed.

## Create covariates based on MDS.
# Perform an MDS ONLY on HapMap data without ethnic outliers. The values of the 10 MDS dimensions are subsequently used as covariates in the association analysis in the third tutorial.
#plink --bfile MDS_merge2 --extract indepSNP.prune.in --genome --out MDS_merge2 --memory 480000
#plink --bfile MDS_merge2 --read-genome MDS_merge2.genome --cluster --mds-plot 10 --out MDS_merge2
#plink2 --bfile MDS_merge2 --pca 10 --out MDS_merge2
#plink --bfile "${file_name}_12" --extract indepSNP.prune.in --genome --out "${file_name}_12"
#plink --bfile "${file_name}_12" --read-genome "${file_name}_12.genome" --cluster --mds-plot 10 --out "${file_name}_12_mds"
plink2 --bfile "${file_name}_12" --pca 10 --out "${file_name}_12_mds"  --memory 480000

# Change the format of the .mds file into a plink covariate file.
awk '{print$1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' "${file_name}_12_mds.eigenvec" > covar_mds.txt

# The values in covar_mds.txt will be used as covariates, to adjust for remaining population stratification, in the third tutorial where we will perform a genome-wide association analysis.

##########################################################################################################################################################################

## CONGRATULATIONS you have succesfully controlled your data for population stratification!

# For the next tutorial you need the following files:
# - "${file_name}_12" (the bfile, i.e., "${file_name}_12".bed,"${file_name}_12".bim,and "${file_name}_12".fam
# - covar_mds.txt





