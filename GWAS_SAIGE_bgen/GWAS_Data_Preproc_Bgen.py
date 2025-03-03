import numpy as np
import pandas as pd
import gwaslab as gl
from SUtils import *

exp_name = "GWAS_NMR_MDD_Subtype_Clustering_bgen"

# minor allele frequency >= 0.01 and imputation INFO score >= 0.9
def SNP_QC():
    # df_w_hm3 = pd.read_csv("/data/UKBioBank/LDScore_EUR/w_hm3.snplist", sep=r'[ \t]', engine='python', usecols=['SNP'])
    # df_w_hm3.rename(columns={'SNP': 'RS_id'}, inplace=True)
    total_snp = 0
    for chr in range(1, 23):
        df_snp = pd.read_csv(bgen_path + "raw_data/ukb_mfi_chr{}_v3.txt".format(chr), sep=r'[ \t]', engine='python',
            names=['Alternate_id', 'RS_id', 'Position', 'Allele1', 'Allele2', 'MAF', 'Minor_Allele', 'Info_score'],
            dtype={'MAF': float, 'Info_score': float})
        print("before chr{}: {}".format(chr, df_snp.shape[0]))
        df_snp = df_snp[(df_snp['MAF'] >= 0.05) & (df_snp['Info_score'] >= 0.9)]
        # df_snp = df_snp.merge(df_w_hm3, on='RS_id', how='inner')
        # print(df_snp.head())
        # # 找到重复的 RS_ID 并打印出来
        # duplicate_rs_ids = df_snp[df_snp.duplicated('RS_id', keep=False)]
        # print(duplicate_rs_ids.shape[0])
        # print(duplicate_rs_ids)
        df_snp = df_snp.drop_duplicates(subset='RS_id')
        print("after chr{}: {}".format(chr, df_snp.shape[0]))
        total_snp += df_snp.shape[0]
        df_snp.to_csv(bgen_path + "SNP_QC_raw_data/ukb_mfi_chr{}_v3.txt".format(chr), sep='\t', header=None, index=False)
        df_snp.to_csv(bgen_path + "SNP_QC_raw_data/keep_snp_chr{}.txt".format(chr), columns=['RS_id'], sep='\t', header=None, index=False)
    print("Total_snp: {}".format(total_snp))

def trans_genotyping_array(x):
    if -11 <= x <= -1:
        return 'BiLEVE'
    elif 1 <= x <= 95:
        return 'Axiom'
    else:
        return None  # 或者根据需要返回其他值


# ukb22828_c1_b0_v3_s487150.sample是储存样本信息的文件，为有标题行的文本文件，每个样本占一行，共有以下4列：
# 1. ID_1：样本的第一列ID，通常是家系 ID（Family ID, FID），在UK Biobank中，这列通常与ID_2相同，UK Biobank样本大多是单独个体，而不是家系数据。
# 2. ID_2：通常是个体ID（Individual ID, IID），UK Biobank使用的ID_2是样本的唯一标识符，通常与参与者的eid（UK Biobank 分配给参与者的唯一编码）相对应
# 3. missing：指示样本是否有缺失数据。这个字段的值通常为 0 或 1
#    - 0 表示样本数据可用，没有缺失。
#    - 1 表示样本有缺失数据，或该样本被标记为不包含在分析中
# 4. sex：样本的性别，通常用以下编码表示
#    - 1：男性（Male）
#    - 2：女性（Female）
#    - 0 或 NA：性别未知或未指定。
# 第二行 (0 0 0 D) 作为占位符，常用于确保文件格式的一致性。它不包含实际的样本数据，也不用于后续的分析。
# 这种结构是为了保持 .sample 文件的标准化格式，而实际的样本信息从第三行开始（即真正的样本数据行）。

# 22006 白种人
# 22005 高缺失值 98 %
# 22021 亲缘关系 -1 0 1 10 只要0
# 22001 0女性1男性  与自我报告的性别进行匹配
# 22000 基因分型阵列/批次
#       1000 BiLEVE   -11 到 -1
#       2000 Axiom     1  到 95
# 22009 取前10个主成分
def Sample_QC():
    df_bgen_sample = pd.read_csv(bgen_path + "raw_data/ukb22828_c1_b0_v3_s487150.sample", sep=r'[ \t]', header=0, engine='python')
    df_bgen_sample.rename(columns={'ID_2': 'eid'}, inplace=True)
    print(df_bgen_sample.columns.tolist())
    print(df_bgen_sample.shape[0])
    df_MDS = pd.read_csv(root_path + "Sample_QC.csv")
    print(df_MDS.columns.tolist())
    df = pd.read_csv(data_path + "ALL.csv", usecols=['eid', 'Age', 'Sex_cate'], low_memory=False)  # 'L1_Female', 'L2_Male'
    print(df.columns.tolist())
    df['Sex_cate'] = df['Sex_cate'].replace({'L1_Female': 0, 'L2_Male': 1})
    df_MDS = reduce(lambda left, right: pd.merge(left, right, on=['eid'], how='inner'), [df_MDS, df, df_bgen_sample])
    df_MDS = df_MDS[df_MDS["22006-0.0"] == 1]  # 22006-0.0 代表高加索人群   22020-0.0 Used in genetic principal components
    print("filter after 22006: {}".format(df_MDS.shape[0]))
    df_MDS = df_MDS[df_MDS["22005-0.0"] <= 0.02]
    print("filter after 22005: {}".format(df_MDS.shape[0]))
    df_MDS = df_MDS[df_MDS["22021-0.0"] == 0]
    print("filter after 22021: {}".format(df_MDS.shape[0]))
    # 比较两列22001和Sex_cate并赋值给新列 'Sex_match', 相等的值被赋值为1，不相等的值被赋值为0
    df_MDS['Sex_match'] = (df_MDS['22001-0.0'] == df_MDS['Sex_cate']).astype(int)
    df_MDS = df_MDS[df_MDS["Sex_match"] == 1]
    print("filter after Sex Match Compare: {}".format(df_MDS.shape[0]))
    PCA_original_cols, PCA_new_cols = ["22009-0."+ str(i) for i in range(1, 41)], ["PCA"+str(i) for i in range(1, 41)]
    df_MDS.rename(columns=dict(zip(PCA_original_cols, PCA_new_cols)), inplace=True)
    df_MDS['genotyping_array'] = df_MDS['22000-0.0'].apply(trans_genotyping_array)
    df_MDS.rename(columns={'22000-0.0': 'genotyping_batch'}, inplace=True)
    print(df_MDS.head())
    df_MDS.to_csv(path + "NMR_Sample_QC_filter.csv", sep=',', index=False)


# fam文件包括以下信息：
# （1）第一列：家系编号(“FID”)；
# （2）第二列：个体编号(“IID”; 不能是”0”)；
# （3）第三列：父系编号 (“0”表示父系信息缺失)；
# （4）第四列：母系编号(“0”表示母系信息缺失)；
# （5）第五列：性别编号(1= 男, 2=女, 0=性别未知)；
# （6）第六列：表型值 (1=对照,2=病例, -9/0/表示表型缺失)
def split_cluster(l1, l2, f_name):
    df_MDS = pd.read_csv(path + "NMR_Sample_QC_filter.csv", low_memory=False)
    # print("Sample_MDS_QC: {}".format(df_MDS.shape[0]))

    df_NC = pd.read_csv(data_path + "Diag_NC.csv",  usecols=['eid'], low_memory=False)
    df_NC['final_cluster'] = 0
    rst_file = result_path + "{}/{}_components_cluster_rst_{}.csv".format("NMF", 2, 3)
    df_Y = pd.read_csv(rst_file, low_memory=False, usecols=['eid', 'final_cluster'])
    print(df_Y.groupby('final_cluster').size())
    df = pd.concat([df_NC, df_Y])
    df = reduce(lambda left, right: pd.merge(left, right, on=['eid'], how='inner'), [df, df_MDS])
    print("df after merge: {}".format(df.shape[0]))

    # 0 vs 1  ----- 0为对照组，1为病例组
    df_0, df_1 = df[df['final_cluster'].isin(l1)], df[df['final_cluster'].isin(l2)]
    df_0['final_cluster'], df_1['final_cluster'] = 0, 1
    df = pd.concat([df_0, df_1])
    print("f_name: {}".format(f_name))
    print("df final: {}".format(df.shape[0]))
    print(df.groupby('final_cluster').size())

    if not os.path.exists(bgen_path + "{}/{}/".format(exp_name, f_name)):
        os.makedirs(bgen_path + "{}/{}/".format(exp_name, f_name))
    # 根据年龄、性别、基因分型批次、基因分型阵列和10个主成分进行调整，以纠正人口分层
    df.to_csv(bgen_path + "{}/{}/{}_pheno.txt".format(exp_name, f_name, f_name),
              columns=['eid', 'Sex_cate', 'Age', 'genotyping_array', 'genotyping_batch', 'final_cluster']+["PCA"+str(i) for i in range(1, 41)], sep='\t', header=True, index=False)
    df.to_csv(bgen_path + "{}/{}/{}.csv".format(exp_name, f_name, f_name), columns=['eid', 'eid'], sep='\t',
              header=None, index=False)

    # # fam 文件生成, 用于后续LDSC计算
    # df['father_ID'] = 0
    # df['mother_ID'] = 0
    # df['Sex_cate'] = df['Sex_cate'].replace({0: 2})   # 1 男 0 女  转为-->  1 男 2 女
    # df['final_cluster'] = df['final_cluster'].replace({0: 1, 1: 2})    # 0为对照组，1为病例组  -->  1为对照组，2为病例组
    # df[['eid', 'eid', 'father_ID', 'mother_ID', 'Sex_cate', 'final_cluster']].to_csv(bgen_path + "{}/{}/{}.fam".format(exp_name, f_name, f_name), sep='\t', header=False, index=False)


def plot_manhattan_QQ(df, file_name, cut=20, sig_level=5e-8):
    df_m = df[["CHR", "POS", "MarkerID", "Allele1", "Allele2", "AF_Allele2", "BETA", "SE", "p.value"]]
    print("file_name: {}, snp number: {}".format(file_name, df_m.shape[0]))
    df_m = gl.Sumstats(df_m, snpid="MarkerID", chrom="CHR", pos="POS",
             ea="Allele2", nea="Allele1", eaf="AF_Allele2",
             beta="BETA", se="SE", p="p.value")
    df_m.plot_mqq(skip=3, cut=cut, anno=True, sig_level=sig_level, sig_level_lead=sig_level,
        colors=["#ff0000", "#fc4444", "#fc6404", "#fcd444", "#8cc43c", "#029658", "#1abc9c", "#5bc0de", "#6454ac", "#fc8c84"],
        fontsize=8, anno_fontsize=15, title_fontsize=13, marker_size=(5, 25), figargs={"figsize": (15, 5), "dpi":300},
        save=bgen_path + "{}/{}/{}manhattan_QQ_plots.png".format(exp_name, file_name, file_name))


def cnt_sig_pvalue(f_s, cuts, sig_levels):
    for file_name, cut, sig_level in zip(f_s, cuts, sig_levels):
        # df = []
        # for chr in range(1, 23):
        #     print("Reading chr{}: ".format(chr))
        #     df_chr = pd.read_csv(bgen_path + "{}/{}/Step2_SPAtests_c{}.txt".format(exp_name, file_name, chr),
        #         sep=r'\s+', engine='python')
        #     df.append(df_chr)
        # df = pd.concat(df)
        df= pd.read_csv(bgen_path + "{}/{}/Step2_SPAtests.txt".format(exp_name, file_name), sep=r'\s+', engine='python')
        df['OR'] = df['BETA'].apply(lambda x: np.exp(x))
        print(df.head())
        # df.to_csv(bgen_path + "{}/{}/SAIGE_results.txt".format(exp_name, file_name), sep='\t', index=False)
        plot_manhattan_QQ(df, file_name, cut, sig_level)
        df.rename(columns={'MarkerID': 'snpid', 'CHR': 'hg18chr', 'POS': 'bp', 'Allele1': 'a2', 'Allele2': 'a1',
                           'OR': 'or', 'p.value': 'pval'}, inplace=True)
        df_ldsc = df[['hg18chr', 'bp', 'snpid', 'a1', 'a2', 'pval', 'or']]
        df_ldsc.to_csv(bgen_path + "{}/{}/SAIGE_results.ldsc.txt".format(exp_name, file_name), sep='\t', index=False)
        df = df[df['pval'] <= 1e-5]
        df = df.sort_values(by="pval")
        df.to_csv(bgen_path + "{}/{}/SAIGE_results.1e-5.csv".format(exp_name, file_name), sep=',', index=False)


if __name__ == '__main__':
    SNP_QC()
    # Sample_QC()

    # # before GWAS
    # split_cluster([0], [1, 2, 3], '0_VS_123')
    # split_cluster([0], [1], '0_VS_1')
    # split_cluster([0], [2], '0_VS_2')
    # split_cluster([0], [3], '0_VS_3')
    # split_cluster([1], [3], '1_VS_3')
    # split_cluster([1], [2], '1_VS_2')
    # split_cluster([3], [2], '3_VS_2')

    #after GWAS
    # f_s, manhattan_cuts, sig_levels = ["1_VS_3", "1_VS_2", "3_VS_2", "0_VS_2", "0_VS_1", "0_VS_3", "0_VS_123"], [8, 12, 8, 8, 8, 8, 10], [5e-6, 5e-8, 5e-6, 5e-6, 5e-6, 5e-6, 5e-8]
    # f_s, manhattan_cuts, sig_levels = ["1_VS_3", "1_VS_2", "3_VS_2"], [8, 8, 8], [5e-6, 5e-6, 5e-6]
    # f_s, manhattan_cuts, sig_levels = ["0_VS_2"], [8], [5e-6]
    # cnt_sig_pvalue(f_s, manhattan_cuts, sig_levels)

    # for file_name, cut, sig_level in zip(f_s, manhattan_cuts, sig_levels):
    #     df = pd.read_csv(bgen_path + "{}/{}/Step2_SPAtests.txt".format(exp_name, file_name), sep=r'[ \t]',
    #                      engine="python", usecols=["CHR", "POS", "MarkerID", "Allele1", "Allele2", "AF_Allele2",
    #                                                "BETA", "SE", "p.value"])
    #     plot_manhattan_QQ(df, file_name, cut, sig_level)

    # for file_name in f_s:
    #     df = pd.read_csv(bgen_path + "{}/{}/SAIGE_results.txt".format(exp_name, file_name), sep=r'[ \t]',
    #                      engine="python", usecols=["CHR", "POS", "MarkerID", "Allele1", "Allele2", "AF_Allele2",
    #                                                "BETA", "SE", "p.value", "N_case", "N_ctrl"])
    #     df.rename(columns={'POS': 'BP', 'MarkerID': 'SNP', 'Allele1': 'A2', 'Allele2': 'A1', "AF_Allele2": "FREQ",
    #         'p.value': 'P'}, inplace=True)
    #     df["N"] = df["N_case"] + df["N_ctrl"]
    #     df[["CHR", "BP", "SNP", "A1", "A2", "FREQ", "BETA", "SE", "P", "N"]].to_csv(bgen_path + "{}/{}/{}_SAIGE_LDGC.txt".format(exp_name, file_name, file_name), sep='\t', index=False)

    # for file_name in f_s:
    #     df = pd.read_csv(bgen_path + "{}/{}/SAIGE_results.txt".format(exp_name, file_name), sep=r'[ \t]',
    #                      engine="python")
    #     df.rename(columns={'MarkerID': 'snpid', 'CHR': 'hg18chr', 'POS': 'bp', 'Allele1': 'a1', 'Allele2': 'a2',
    #                        'OR': 'or', 'p.value': 'pval'}, inplace=True)
    #     df_ldsc = df[['hg18chr', 'bp', 'snpid', 'a1', 'a2', 'pval', 'or']]
    #     df_ldsc.to_csv(bgen_path + "{}/{}/SAIGE_results.ldsc.txt".format(exp_name, file_name), sep='\t', index=False)

    # df = pd.read_csv(bed_path + "assoc_results.sig.csv", low_memory=False)
    # print(df.columns.tolist())
    # df_ldsc = df[['CHR', 'BP', 'SNP', 'A1', 'A2', 'P', 'OR']]
    # df_ldsc.to_csv(bed_path + "bed_gwas_summary.txt", sep='\t', index=False)

