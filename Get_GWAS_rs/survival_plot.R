library(survival)  # Survival Analysis
library(survminer) # Drawing Survival Curves using 'ggplot2'
library(pcoxtime) # Penalized Cox Proportional Hazard Model for Time-Dependent Covariates
library(ggplot2)
library(plyr)
library(aqp) # 输出结果到文件

setwd("/data/UKBioBank/bgen_Data/final_result/img/")
Sys.setenv(LANGUAGE = "en")
options(max.print=1000000) 

data <- read.csv(file='/data/UKBioBank/bgen_Data/final_result/filter_rs_csv/sum_rs964184.csv', sep=",")
sex <- data$Sex_cate
data$Sex_cate<- as.numeric(mapvalues(sex, c("L1_Female", "L2_Male"),c(0, 1)))
head(data)
summary(data)

##########################################################################################################################
# https://wenku.baidu.com/view/e5a5c3044b2fb4daa58da0116c175f0e7dd11977.html
# Reference: https://www.bbsmax.com/A/1O5E7QRaz7/

fit <- survfit(Surv(survival_time, MDD) ~ factor(rs964184), data=data)
# 绘制发病曲线
ggsurvplot(fit, data=data, conf.int=TRUE,risk.table=TRUE, palette="lancet",
          pval=TRUE, pval.coord = c(0.5,0.03), pval.size =4,
          title="Genotype-Depression Curve",
          risk.table.height =0.25,
          xlab="Time(Day)", ylab="Incidence Rate of Depression",
          break.x.by=180, break.y.by=0.01, ylim=c(0, 0.18),
     #     legend.labs=c("0","1","2"), 
          legend.title="A2 Allele Num",
          axes.offset=FALSE,
          fun = "cumhaz")
sink("/data/UKBioBank/bgen_Data/final_result/img/Incidence_rate_depression.txt")#保存为txt文档（定向到txt文件）
summary(fit)
sink()#结束重定向，不能少，相当于close()

# conf.int=TRUE 显示生存率的95% CI
# risk.table=TRUE 显示风险表
# palette="lancet" 柳叶刀配色
# pval = TRUE 显示p值    pval.coord = c(0,0.2) p值坐标位置; pval.size =4 p值字体大小
# title 大标题
# risk.table.height =0.25 风险表的高度比例
# break.x.by=12 X轴刻度间距    break.y.by=0.01 Y轴刻度间距   ylim=c(0, 0.1) y轴范围
# axes.offset=FALSE, 曲线坐标轴从原点开始
# fun ="cumhaz" # 累计发病曲线

# 绘制生存曲线
ggsurvplot(fit, data=data, conf.int=TRUE,risk.table=TRUE, palette="lancet",
          pval=TRUE, pval.coord = c(0.5,0.91), pval.size =4,
          title="Genotype-Depression Curve",
          risk.table.height =0.25,
          xlab="Time(Day)",
          break.x.by=180, break.y.by=0.01, ylim=c(0.82, 1.0),
          legend.labs=c("0","1","2"), legend.title="A2 Allele Num",
          axes.offset=FALSE)




