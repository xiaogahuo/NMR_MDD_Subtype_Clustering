library("IsingSampler")
library("IsingFit")
library("bootnet")
library("NetworkComparisonTest")
library("readr")
library("dplyr")

trans_cate <- function(data){
  data$Sex_cate <- factor(data$Sex_cate, levels = c("L1_Female", "L2_Male"), labels = c(1, 2))
  data$BMI_cate <- factor(data$BMI_cate, levels = c('BMI1', 'BMI2', 'BMI3', 'BMI4', 'BMI5_Unknown'), labels = c(1, 2, 3, 4, 5))
  # 将因子变量转换为数值
  data$Sex_cate <- as.numeric(data$Sex_cate)
  data$BMI_cate <- as.numeric(data$BMI_cate)
  return(data)
}

compare_Subtypes <- function (est1, est2, exp_name){
  ### Compare networks of data sets using NCT ###
  # with gamma = 0.
  # Iterations (it) set to 10 to save time.
  # Low number of iterations can give unreliable results. Should be 1000 at least.
  # Testing the three aspects that are validated (network invariance, global strength, edge weight)
  # 2 edges are tested here: between variable 1 and 2,
  # and between 3 and 6 (can be list(c(2,1),c(6,3)) as well)
  Res <- NCT(est1, est2, gamma=0, it=1000, test.edges=TRUE, edges="all")
  ## Plotting of NCT results
  ## See the help file of plot.NCT for more information about the plotting function and its arguments
  # Plot results of the network structure invariance test (not reliable with only 10 permutations!):
  head <- sprintf("**************************** %s **************************************", exp_name)
  print(head)
  # print(Res)
  summary_res <- capture.output(summary(Res))
  writeLines(summary_res, con = paste(exp_name, "NCT.txt", sep = "_"))
}

setwd("/data/UKBioBank/FileHandler_And_Dataset/NMR_MDD_Subtype_Clustering/")
Sys.setenv(LANGUAGE = "en")

df_MDD <- read_csv("dataset/Diag_MDD.csv")
df_NC <- read_csv("dataset/Diag_NC.csv")
df_NC$final_cluster <- 0

df_Y <- read_csv("result/NMF/2_components_cluster_rst_3.csv", col_select = c("eid", "final_cluster"))
# 替换 final_cluster 列的值
df_Y <- df_Y %>%
  mutate(final_cluster = case_when(
    final_cluster == 2 ~ 3,
    final_cluster == 3 ~ 2,
    TRUE ~ final_cluster  # 对于其他值保持不变
  ))

# merge df_MDD 和 df_Y
df_MDD <- merge(df_MDD, df_Y, by = "eid", all.x = TRUE)
# 合并 df_MDD 和 df_NC
ALL_Data <- bind_rows(df_MDD, df_NC)
ALL_Data <- trans_cate(ALL_Data)

TG_cols <- c('Total_TG')
U_cols <- c('Unsaturation')
V_cols <- c('VLDL_FC','VLDL_L','VLDL_P','VLDL_PL','XXL_VLDL_CE','XL_VLDL_FC','XL_VLDL_P','L_VLDL_C','L_VLDL_FC','L_VLDL_P','L_VLDL_PL','M_VLDL_L', 'M_VLDL_P', 'S_VLDL_L', 'S_VLDL_P', 'S_VLDL_TG','XS_VLDL_TG')
I_cols <- c('IDL_TG')
L_cols <- c('LDL_TG','L_LDL_TG','M_LDL_TG','M_LDL_FC(%)', 'S_LDL_TG')
H_cols <- c('L_HDL_TG(%)','S_HDL_TG')
FA_cols <- c('DHA(%)', 'MUFA', 'MUFA(%)', 'Omega_6_by_Omega_3','PUFA(%)','PUFA_by_MUFA')
G_cols <- c('GlycA')
cols <- c(TG_cols, U_cols, V_cols, I_cols, L_cols, H_cols, FA_cols, G_cols)
covar_cols <- c('Age', 'chronic_num', 'Sex_cate', 'BMI_cate')   # sex: Fmale Male; BMI: BMI1,BMI2,BMI3,BMI4,BMI5
ALL_cols <- c(cols, covar_cols)

# 设置全局选项以屏蔽Note信息
options(show.error.messages =  TRUE) # TRUE FALSE
options(warn = 0) # 0 -1
# c(0), c(1), c(2), c(3), c(1, 2, 3)
df1 <- ALL_Data %>% filter(final_cluster %in% c(1)) %>% select(all_of(cols))
df2 <- ALL_Data %>% filter(final_cluster %in% c(2)) %>% select(all_of(cols))
df3 <- ALL_Data %>% filter(final_cluster %in% c(3)) %>% select(all_of(cols))
df123 <- ALL_Data %>% filter(final_cluster %in% c(1, 2, 3)) %>% select(all_of(cols))
df0 <- ALL_Data %>% filter(final_cluster %in% c(0)) %>% select(all_of(cols))
est1 <- estimateNetwork(df1, default = "EBICglasso")
est2 <- estimateNetwork(df2, default = "EBICglasso")
est3 <- estimateNetwork(df3, default = "EBICglasso")
est123 <- estimateNetwork(df123, default = "EBICglasso")
est0 <- estimateNetwork(df0, default = "EBICglasso")
compare_Subtypes(est1, est2, "1_VS_2")
compare_Subtypes(est1, est3, "1_VS_3")
compare_Subtypes(est2, est3, "2_VS_3")
compare_Subtypes(est1, est123, "1_VS_123")
compare_Subtypes(est2, est123, "2_VS_123")
compare_Subtypes(est3, est123, "3_VS_123")
compare_Subtypes(est1, est0, "1_VS_0")
compare_Subtypes(est2, est0, "2_VS_0")
compare_Subtypes(est3, est0, "3_VS_0")
compare_Subtypes(est123, est0, "123_VS_0")

df1_covar <- ALL_Data %>% filter(final_cluster %in% c(1)) %>% select(all_of(ALL_cols))
df2_covar <- ALL_Data %>% filter(final_cluster %in% c(2)) %>% select(all_of(ALL_cols))
df3_covar <- ALL_Data %>% filter(final_cluster %in% c(3)) %>% select(all_of(ALL_cols))
df123_covar <- ALL_Data %>% filter(final_cluster %in% c(1, 2, 3)) %>% select(all_of(ALL_cols))
df0_covar <- ALL_Data %>% filter(final_cluster %in% c(0)) %>% select(all_of(ALL_cols))
est1_covar <- estimateNetwork(df1_covar, default = "EBICglasso")
est2_covar <- estimateNetwork(df2_covar, default = "EBICglasso")
est3_covar <- estimateNetwork(df3_covar, default = "EBICglasso")
est123_covar <- estimateNetwork(df123_covar, default = "EBICglasso")
est0_covar <- estimateNetwork(df0_covar, default = "EBICglasso")
compare_Subtypes(est1_covar, est2_covar, "1_VS_2_covar")
compare_Subtypes(est1_covar, est3_covar, "1_VS_3_covar")
compare_Subtypes(est2_covar, est3_covar, "2_VS_3_covar")
compare_Subtypes(est1_covar, est123_covar, "1_VS_123_covar")
compare_Subtypes(est2_covar, est123_covar, "2_VS_123_covar")
compare_Subtypes(est3_covar, est123_covar, "3_VS_123_covar")
compare_Subtypes(est1_covar, est0_covar, "1_VS_0_covar")
compare_Subtypes(est2_covar, est0_covar, "2_VS_0_covar")
compare_Subtypes(est3_covar, est0_covar, "3_VS_0_covar")
compare_Subtypes(est123_covar, est0_covar, "123_VS_0_covar")





