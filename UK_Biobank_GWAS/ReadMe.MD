1. 执行GWAS_data_preproc.py 的 Sample_QC() 进行初步的 Sample 过滤
2. cd /data/UKBioBank/bed_Data
   conda activate GWAS
3. 执行 ./filter_individual.sh {exp_name} 生成用于执行GWAS QC 的数据集到 {exp_name}/GWAS_QC 子目录
4. 执行 ./run_GWAS_QC.sh {exp_name} 质控

5. GWAS_data_preproc.py 的 split_cluster() 生成用于执行GWAS 的 Sample
6. 执行./run_script_association_GWAS.sh  {exp_name} {file_name} 生成 GWAS result
7. 执行GWAS_data_preproc.py 的 cnt_sig_pvalue() 获取 significant snp
保存 GWAS_QC  GWAS  两个子目录下的GWAS分析结果