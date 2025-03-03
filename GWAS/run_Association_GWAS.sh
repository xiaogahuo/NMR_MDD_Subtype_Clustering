#!/bin/bash
# conda activate GWAS
# ./run_script_association_GWAS.sh  GWAS_sample_cluster1
  exp_name=$1
  file_name=$2

################ Explanation of the main script ##########################


###########################################################
### Association analyses ###
cp -f "./GWAS_Script/Manhattan_plot.R" "./${exp_name}/${file_name}/"
cp -f "./GWAS_Script/QQ_plot.R" "./${exp_name}/${file_name}/"

# For binary traits.
cd "./${exp_name}/${file_name}/"

# assoc
plink --bfile "${file_name}_8" --assoc --out assoc_results
# Note, the --assoc option does not allow to correct covariates such as principal components (PC's)/ MDS components, which makes it less suited for association analyses.

# logistic 
# We will be using 10 principal components as covariates in this logistic analysis. We use the MDS components calculated from the previous tutorial: covar_mds.txt.
 plink --bfile "${file_name}_8" --covar "${file_name}_covr.txt" --logistic --hide-covar --out logistic_results
# plink --bfile "${file_name}_8" --covar "${file_name}_covr.txt" --linear --hide-covar --out assoc_covar_results
# Note, we use the option -�hide-covar to only show the additive results of the SNPs in the output file.

# Remove NA values, those might give problems generating plots in later steps.
awk '!/'NA'/' logistic_results.assoc.logistic > logistic_results.assoc.logistic2

# The results obtained from these GWAS analyses will be visualized in the last step. This will also show if the data set contains any genome-wide significant SNPs.

# Note, in case of a quantitative outcome measure the option --logistic should be replaced by --linear. The use of the --assoc option is also possible for quantitative outcome measures (as metioned previously, this option does not allow the use of covariates).

#################################################################

# Multiple testing
# There are various way to deal with multiple testing outside of the conventional genome-wide significance threshold of 5.0E-8, below we present a couple. 

#adjust
plink --bfile "${file_name}_8" -assoc --adjust --out adjusted_assoc_results
# This file gives a Bonferroni corrected p-value, along with FDR and others.

## Permutation
# This is a computational intensive step. Further pros and cons of this method, which can be used for association and dealing with multiple testing, are described in our article corresponding to this tutorial (https://www.ncbi.nlm.nih.gov/pubmed/29484742).
# The reduce computational time we only perform this test on a subset of the SNPs from chromosome 22.
# The EMP2 collumn provides the for multiple testing corrected p-value.

## Generate subset of SNPs
#awk '{ if ($4 >= 21595000 && $4 <= 21605000) print $2 }' "${file_name}_12.bim" > subset_snp_chr_22.txt
## Filter your bfile based on the subset of SNPs generated in the step above.
#plink --bfile "${file_name}_12" --extract subset_snp_chr_22.txt --make-bed --out HapMap_subset_for_perm
## Perform 1000000 perrmutations.
#plink --bfile HapMap_subset_for_perm --assoc --mperm 1000000 --out subset_1M_perm_result
#
## Order your data, from lowest to highest p-value.
#sort -gk 4 subset_1M_perm_result.assoc.mperm > sorted_subset.txt
## Check ordered permutation results
#head sorted_subset.txt

#####################################################################

# Generate Manhattan and QQ plots.

# These scripts assume R >= 3.0.0.
# If you changed the name of the .assoc file or to the assoc.logistic file, please assign those names also to the Rscripts for the Manhattan and QQ plot, otherwise the scripts will not run.

# The following Rscripts require the R package qqman, the scripts provided will automatically download this R package and install it in /home/{user}/ . Additionally, the scripts load the qqman library and can therefore, similar to all other Rscript on this GitHub page, be executed from the command line.
# This location can be changed to your desired directory

Rscript --no-save Manhattan_plot.R
Rscript --no-save QQ_plot.R

# Please read below when you encountered an error:
# Note, the mirror used to download the package qqman can no longer by active, which will result in an error (so below for specific error).
# If you encounter this error, please contact me at a.t.marees@amc.uva.nl.
# This error can simply be resolved by changing the addresses in the scripts: Manhattan_plot.R and QQ_plot.R.
# Simply change the address (http://cran...) in line below, in both Rscripts (Manhattan_plot.R and QQ_plot.R) with (for example) https://cran.univ-paris1.fr/ .
# install.packages("qqman",repos="http://cran.cnr.berkeley.edu/",lib="~" ) # location of installation can be changed but has to correspond with the library location

# Example error:
# "Error in library("qqman", lib.loc = "-") :
#	there is no package called 'qqman'
#Execution halted"

#######################################################

## CONGRATULATIONS you have succesfully conducted a GWAS analyses!!

# If you are also interested in learning how to conduct a polygenic risk score (PRS) analysis please see our fourth tutorial.
# The tutorial explaining PRS is independent from the previous tutorials.
