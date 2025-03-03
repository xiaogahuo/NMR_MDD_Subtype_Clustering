# Note: In RStudio change all dev.new() calls to x11(). Also, some
# figures might look different depending on the computer used to
# compile these. This is because the Fruchterman-Reingold algorithm
# behaves chaotically and floating point errors might cause the layout
# to look very different. The general structure should be the same.

###################################################
### Load library and data
###################################################
# library("psych")
library("readr")
library("qgraph")
library("openxlsx")
library("bootnet")
library("mgm")
library("dplyr")


setwd("/data/UKBioBank/FileHandler_And_Dataset/NMR_MDD_Subtype_Clustering/")
Sys.setenv(LANGUAGE = "en")
# dataTotal <- openxlsx::read.xlsx("D:/PythonWorkspace/5.university/xiu武汉理工二师数据合并.xlsx", sheet = 1) # 1-all 2-Male 3-Female  4-ALL-Suicide_C 5-Male-Suicide 6-Female-Suicide
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

TG_cols <- c('Total_TG')
U_cols <- c('Unsaturation')
V_cols <- c('VLDL_FC','VLDL_L','VLDL_P','VLDL_PL','XXL_VLDL_CE','XL_VLDL_FC','XL_VLDL_P','L_VLDL_C','L_VLDL_FC','L_VLDL_P','L_VLDL_PL','M_VLDL_L', 'M_VLDL_P', 'S_VLDL_L', 'S_VLDL_P', 'S_VLDL_TG','XS_VLDL_TG')
I_cols <- c('IDL_TG')
L_cols <- c('LDL_TG','L_LDL_TG','M_LDL_TG','M_LDL_FC(%)', 'S_LDL_TG')
H_cols <- c('L_HDL_TG(%)','S_HDL_TG')
FA_cols <- c('DHA(%)', 'MUFA', 'MUFA(%)', 'Omega_6_by_Omega_3','PUFA(%)','PUFA_by_MUFA')
G_cols <- c('GlycA')
cols <- c(TG_cols, U_cols, V_cols, I_cols, L_cols, H_cols, FA_cols, G_cols)
TG_Names <- c('TG')
U_Names <- c('U')
V_Names <- c('V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9', 'V10', 'V11', 'V12', 'V13', 'V14', 'V15', 'V16', 'V17')
I_Names <- c('I')
L_Names <- c('L1', 'L2', 'L3', 'L4', 'L5')
H_Names <- c('H1', 'H2')
FA_Names <- c('FA1', 'FA2', 'FA3', 'FA4', 'FA5', 'FA6')
G_Names <- c('G')
cols_Names <- c(TG_Names, U_Names, V_Names, I_Names, L_Names, H_Names, FA_Names, G_Names)
type <- rep('g', 34)  # "g" for Gaussian, "p" for Poisson, "c" for categorical.
lev <- rep(1, 34)  # vector indicating the number of categories of each variable. For continuous variables set to 1.
GroupLabels  <- list("Total_TG" = c(1), "Unsaturation" = c(2), "VLDL" = c(3:19), "IDL" = c(20),
                "LDL" = c(21:25), "HDL" = c(26:27), "FA" = c(28:33), "GlycA" = c(34))

covar_cols <- c('Age', 'chronic_num', 'Sex_cate', 'BMI_cate')   # sex: Fmale Male; BMI: BMI1,BMI2,BMI3,BMI4,BMI5
covar_Names <- c('C1', 'C2', 'C3', 'C4')
type_covar <- c(rep('g', 36), rep('c', 2))   # "g" for Gaussian, "p" for Poisson, "c" for categorical.
lev_covar <- c(rep(1, 36), c(2, 5))  # vector indicating the number of categories of each variable. For continuous variables set to 1.
GroupLabelsWithCovar  <- list("Total_TG" = c(1), "Unsaturation" = c(2), "VLDL" = c(3:19), "IDL" = c(20),
                "LDL" = c(21:25), "HDL" = c(26:27), "FA" = c(28:33), "GlycA" = c(34), "Covar"=c(35:38))
ALL_cols <- c(cols, covar_cols)
ALL_Names <- c(cols_Names, covar_Names)

trans_cate <- function(data){
  data$Sex_cate <- factor(data$Sex_cate, levels = c("L1_Female", "L2_Male"), labels = c(1, 2))
  data$BMI_cate <- factor(data$BMI_cate, levels = c('BMI1', 'BMI2', 'BMI3', 'BMI4', 'BMI5_Unknown'), labels = c(1, 2, 3, 4, 5))
  # 将因子变量转换为数值
  data$Sex_cate <- as.numeric(data$Sex_cate)
  data$BMI_cate <- as.numeric(data$BMI_cate)
  return(data)
}
ALL_Data <- trans_cate(ALL_Data)

# 打印每个字段值范围，包含离散变量
for (col_name in ALL_cols) {
  if (is.factor(ALL_Data[[col_name]])) {
    cat("分类变量：", col_name, "类别：", levels(ALL_Data[[col_name]]), "\n")
  } else {
    cat("连续变量：", col_name, "范围：", min(ALL_Data[[col_name]]), "到", max(ALL_Data[[col_name]]), "\n")
  }
}

df_Total <- ALL_Data %>% filter(final_cluster %in% c(0)) # c(0), c(1), c(2), c(3), c(1, 2, 3)
# 自定义颜色调色板
my_colors <- c("red", "darkblue", "darkgreen", "purple", "orange", "pink", "brown", "darkcyan", "black", "yellow")

# -------------------------------------------------------------------------
# ------------------------- (34 variables) --------------------------------
# -------------------------------------------------------------------------
df <- df_Total[, cols]
Network <- estimateNetwork(df, default = "EBICglasso")
# plot(Network, layout="spring", labels=TRUE)

centralityPlot(Network, labels= cols, scale = c("z-scores"), include=c("Strength", "Closeness", "Betweenness"), orderBy = "Strength")

Results <- bootnet(Network, nBoots = 1000, nCores = 8)

# Plot bootstrapped edge CIs:
CIresult <- plot(Results, labels = FALSE, order = "sample")
print(CIresult)
print(CIresult$data$id)
# print(CIresult$data$CIlower)
# print(CIresult$data$CIupper)

# differences
# plot(Results, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")
plot(Results, "strength")

# Compute CS-coefficients:
CSresults <- bootnet(Network, nBoots = 1000, nCores = 8, type="case")
plot(CSresults)
corStability(CSresults)

# ---------- Gamma = 0 model ----------
# Fit model
gam0 <- mgm(data = as.matrix(df),
            type = type,
            level = lev,
            labels = cols,
            lambdaSel = "EBIC",
            lambdaGam = 0)

# Compute Predictability
Pred <- predict(gam0, df)
Pred$errors
# pie11 <- as.numeric(as.character(Pred_phq9$errors[1:2, 3]))
# pie12 <- as.numeric(as.character(Pred_phq9$errors[3:4, 5]))
# pie13 <- as.numeric(as.character(Pred_phq9$errors[5:13, 3]))#predictability estimates
# pie<- c(pie11, pie12, pie13 ) #predictability estimates as one piece , pie13, pie14
# View(pie1a)
# mean(pie)
pie <- as.numeric(as.character(Pred$errors[1:34, 2]))

# Plot Network
N1amgm0 <- qgraph(gam0$pairwise$wadj, layout = "spring", cut=0,
                  maximum = 0.4,
                  groups=GroupLabels , palette = "colorblind",
                  nodeNames=cols_Names , labels = cols_Names , vsize=4.5,
                  label.cex=1.3, legend.cex=.38, GLratio = 1.8,
                  label.color = "white", label.prop = 0.65, color = my_colors[1:8],
                  pie = pie, pieBorder = 0.25, threshold= 0.1)
N1amgm0$Edgelist
gam0$pairwise$wadj
cent<- centrality_auto(N1amgm0, weighted = TRUE, signed = TRUE)
cent$node.centrality
# cent$edge.betweenness.centrality


# -------------------------------------------------------------------------
# -------------------- (34 variables + 4 covariates) ----------------------
# -------------------------------------------------------------------------
# Get data
df <- df_Total[, ALL_cols]
# View(data)

Network <- estimateNetwork(df, default = "EBICglasso")

# plot(Network, layout="spring", labels=TRUE)

centralityPlot(Network, labels = ALL_cols, scale = c("z-scores"), include=c("Strength", "Closeness", "Betweenness"), orderBy = "Strength")

Results <- bootnet(Network, nBoots = 1000, nCores = 8)
# Plot bootstrapped edge CIs:
CIresult <- plot(Results, labels = FALSE, order = "sample")
print(CIresult)
print(CIresult$data$id)

# differences
# plot(Results, "edge", plot = "difference", onlyNonZero = TRUE, order = "sample")
plot(Results, "strength")

# Compute CS-coefficients:
CSresults <- bootnet(Network, nBoots = 1000, nCores = 8, type="case")
plot(CSresults)
corStability(CSresults)

for (col_name in ALL_cols) {
  if (is.factor(df[[col_name]])) {
    cat("分类变量：", col_name, "类别：", levels(df[[col_name]]), "\n")
  } else {
    cat("连续变量：", col_name, "范围：", min(df[[col_name]]), "到", max(df[[col_name]]), "\n")
  }
}

# cluster 为3 时 BMI_cate == 2 只有一个样本，无法执行mgm分析，剔除该样本并将BIM cate Number 由 5 变为 4
# df <- df %>%  filter(BMI_cate != 2)
# lev_covar <- c(rep(1, 36), c(2, 4)) 
# ---------- Gamma = 0 model ----------
# Fit model
gam0 <- mgm(data = as.matrix(df),
            type = type_covar,
            level = lev_covar,
            lambdaSel = "EBIC",
            lambdaGam = 0)


# Compute Predictability
Pred <- predict(gam0, df)
Pred$errors
# pie11 <- as.numeric(as.character(Pred_phq9$errors[1:2, 3]))
# pie12 <- as.numeric(as.character(Pred_phq9$errors[3:4, 5]))
# pie13 <- as.numeric(as.character(Pred_phq9$errors[5:13, 3]))#predictability estimates
# pie<- c(pie11, pie12, pie13 )#predictability estimates as one piece , pie13, pie14
# View(pie1a)
# mean(pie)
pie <- as.numeric(as.character(Pred$errors[1:36, 2])) # 两个分类变量值为NA， 设置为0.01
pie <- c(pie, rep(0.01, 2))
print(pie)

# Plot Network
N1amgm0 <- qgraph(gam0$pairwise$wadj, layout = "spring", cut=0,
                  maximum = 0.4,
                  groups=GroupLabelsWithCovar , palette = "colorblind",
                  nodeNames=ALL_Names , labels = ALL_Names , vsize=4.5,
                  label.cex=1.3, legend.cex=.38, GLratio = 1.8,
                  label.color = "white", label.prop = 0.65,color = my_colors[1:9],
                  pie = pie, pieBorder = 0.25, threshold= 0.1)
N1amgm0$Edgelist
gam0$pairwise$wadj
cent<- centrality_auto(N1amgm0, weighted = TRUE, signed = TRUE)
cent$node.centrality
# cent$edge.betweenness.centrality

