if (!requireNamespace("cutoff", quietly = TRUE)) {
  install.packages("cutoff")
}
if (!requireNamespace("stats", quietly = TRUE)) {
  install.packages("stats")
}
# install.packages("readxl")
library(readr)
library(readxl)
library(rms)
library(rmda)
library(dplyr)
library(plyr)
library(Cairo)#保存图片的包
library(ggplot2)

setwd("/data/UKBioBank/FileHandler_And_Dataset/NMR_MDD_Subtype_Clustering/dataset/")

trans_cate <- function(df){
  df$Sex_cate <- as.factor(df$Sex_cate)
  df$final_cluster <- as.factor(df$final_cluster)
  return(df)
}

df <- read_csv("PRS_SCORE_Dataset.csv")
df <- trans_cate(df)

# colnames(df)
X_cols <- c('final_cluster','Sex_cate','Age')
Y_cols <-c('5e08', '1e07', '5e07', '1e06', '5e06', '1e05', '5e05', '1e04', '5e04', '1e03', '5e03', '1e02', '5e02')
Y_cols <- lapply(Y_cols, function(x) paste0('SCORE_', x))

# 打印每个字段值范围，包含离散变量
for (col_name in X_cols) {
  if (is.factor(df[[col_name]])) {
    cat("分类变量：", col_name, "类别：", levels(df[[col_name]]), "\n")
  } else {
    cat("连续变量：", col_name, "范围：", min(df[[col_name]]), "到", max(df[[col_name]]), "\n")
  }
}


############ 线性回归   ########################
var_str <- paste(X_cols, collapse=' + ')

lm_four_cates <- function(header_str, file_name){
  result_lst <-list()
  result_lst <- append(result_lst, header_str)
  for (i in 1:length(Y_cols)) {
    formula_tmp <- paste(Y_cols[i], var_str, sep=" ~ ") # quadratic model
    print(formula_tmp)
    model <- lm(formula_tmp, data = df)
    summary_result <- summary(model)
    print(summary_result)

    esti_2 <- summary_result$coefficients[2,1]
    std_2  <- summary_result$coefficients[2,2]
    tval_2 <- summary_result$coefficients[2,3]
    pval_2 <- summary_result$coefficients[2,4]

    esti_3 <- summary_result$coefficients[3,1]
    std_3  <- summary_result$coefficients[3,2]
    tval_3 <- summary_result$coefficients[3,3]
    pval_3 <- summary_result$coefficients[3,4]

    esti_4 <- summary_result$coefficients[4,1]
    std_4  <- summary_result$coefficients[4,2]
    tval_4 <- summary_result$coefficients[4,3]
    pval_4 <- summary_result$coefficients[4,4]

    Fstatistic <- summary_result$fstatistic[[1]]
    result_lst <- append(result_lst, paste(Y_cols[i],esti_2,std_2,tval_2,pval_2,esti_3,std_3,tval_3,pval_3,esti_4,std_4,tval_4,pval_4,Fstatistic, sep=","))
  }
  writeLines(unlist(result_lst), paste(file_name, "PRS_Linear_Regression_Result.csv", sep="_"))
#   df_result <- read.csv(file=paste(file_name, "all.csv", sep="_"), sep=",")
#   df_result$P2_fdr <- p.adjust(df_result$pval_2, method = "fdr")
#   df_result$P3_fdr <- p.adjust(df_result$pval_3, method = "fdr")
#   df_result$P4_fdr <- p.adjust(df_result$pval_4, method = "fdr")
#   write.csv(df_result, file = paste(file_name, "all_fdr.csv", sep="_"), row.names = FALSE)
#   df_P2_fdr <- subset(df_result, P2_fdr<=0.05 )
#   write.csv(df_P2_fdr, file = paste(file_name, "only_fdr_P2.csv", sep="_"), row.names = FALSE)
#   df_P3_fdr <- subset(df_result, P3_fdr<=0.05 )
#   write.csv(df_P3_fdr, file = paste(file_name, "only_fdr_P3.csv", sep="_"), row.names = FALSE)
#   df_P4_fdr <- subset(df_result, P4_fdr<=0.05 )
#   write.csv(df_P4_fdr, file = paste(file_name, "only_fdr_P4.csv", sep="_"), row.names = FALSE)
}


header_str <- "PRS_Col,esti_1,std_1,tval_1,pval_1,esti_2,std_2,tval_2,pval_2,esti_3,std_3,tval_3,pval_3,Fstatistic"
lm_four_cates(header_str, "MDD")


