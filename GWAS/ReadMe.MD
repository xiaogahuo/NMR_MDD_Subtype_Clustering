##### 1. 过滤 Sample
    1.1 执行GWAS_data_preproc.py 的 Sample_QC()  生成用于过滤 sample 的csv文件
    conda activate GWAS
    cd /data/UKBioBank/bed_Data
    1.2 执行 nohup ./filter_individual.sh {exp_name} {file_name} > filter_individual.log & 生成用于执行GWAS QC 的数据集到 /data/UKBioBank/bed_Data/{exp_name}/{file_name} 子目录
        tail -f filter_individual.log
    1.3 执行 nohup ./run_GWAS_QC.sh  {exp_name} {file_name} > GWAS_QC.log & 进行质控
        tail -f GWAS_QC.log        

##### 2. 执行 GWAS 关联分析
    2.1 执行GWAS_data_preproc.py 的 split_cluster() 生成用于执行 GWAS 关联分析的 fam文件及covr协变量文件 
    conda activate GWAS
    cd /data/UKBioBank/bed_Data
    2.2 执行 nohup ./run_Association_GWAS.sh {exp_name} {file_name} > GWAS_Association.log & 执行 GWAS 关联分析
        tail -f GWAS_Association.log
    2.3 执行GWAS_data_preproc.py 的 cnt_sig_pvalue() 将 A2 等位基因信息合并到 GWAS summary 的 logistic_results.assoc.logistic3 文件中, 用于后续的LD分析

##### 3. 基于第2步的GWAS summary结果计算LD SCore、基因多效性、遗传度、遗传相关性
[//]: # (https://www.cnblogs.com/chenwenyan/p/11321272.html)
    虚拟环境名称为ldsc (python 2.7 环境): 
    conda activate ldsc 
    cd /data/UKBioBank/bed_Data
    3.1 执行 nohup ./run_LDSC.sh ${exp_name} ${file_name} > /data/UKBioBank/bed_Data/${exp_name}/${file_name}/ldsc.log & 用于生成LD 分析结果

    nohup ./run_LDSC.sh GWAS_NMR_MDD_Subtype_Clustering 0_VS_1 > /data/UKBioBank/bed_Data/GWAS_NMR_MDD_Subtype_Clustering/0_VS_1/ldsc.log &
    nohup ./run_LDSC.sh GWAS_NMR_MDD_Subtype_Clustering 0_VS_2 > /data/UKBioBank/bed_Data/GWAS_NMR_MDD_Subtype_Clustering/0_VS_2/ldsc.log &
    nohup ./run_LDSC.sh GWAS_NMR_MDD_Subtype_Clustering 0_VS_3 > /data/UKBioBank/bed_Data/GWAS_NMR_MDD_Subtype_Clustering/0_VS_3/ldsc.log &
    nohup ./run_LDSC.sh GWAS_NMR_MDD_Subtype_Clustering 0_VS_123 > /data/UKBioBank/bed_Data/GWAS_NMR_MDD_Subtype_Clustering/0_VS_123/ldsc.log &
    nohup ./run_LDSC.sh GWAS_NMR_MDD_Subtype_Clustering 1_VS_2 > /data/UKBioBank/bed_Data/GWAS_NMR_MDD_Subtype_Clustering/1_VS_2/ldsc.log &
    nohup ./run_LDSC.sh GWAS_NMR_MDD_Subtype_Clustering 1_VS_3 > /data/UKBioBank/bed_Data/GWAS_NMR_MDD_Subtype_Clustering/1_VS_3/ldsc.log &
    nohup ./run_LDSC.sh GWAS_NMR_MDD_Subtype_Clustering 2_VS_3 > /data/UKBioBank/bed_Data/GWAS_NMR_MDD_Subtype_Clustering/2_VS_3/ldsc.log &
    
    3.2 执行GWAS_data_preproc.py 的 format_LDGC_rst() 将LDGC 结果写入LDGC.csv 文件中
    
#### 注意： 注意执行 df -h 后查看 /data/ 目录是否已满
执行完实验将 QQ-Plot_assoc.jpeg、 assoc_manhattan.jpeg、logistic_results.1e-5.txt、 {file_name}_gwas_summary.txt.gz、LDGC文件夹、LDGC.csv、ldsc.log结果分类拷贝到本地

