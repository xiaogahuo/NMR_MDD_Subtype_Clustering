import pandas as pd
from SUtils import *


exp_name = "GWAS_NMR_MDD_Subtype_Clustering"

# ukbconv ukb669272.enc_ukb csv -ifields_id_Sample_QC.txt -oSample_QC    获取 Sample_QC.csv
def Sample_QC():
    df = pd.read_csv(root_path + "Sample_QC.csv")
    print("Sample_QC: {}".format(df.shape[0]))
    df = df[(df["22006-0.0"] == 1) & (df["22020-0.0"] == 1)]   # 22006-0.0 代表高加索人群   22020-0.0 Used in genetic principal components
    print("Sample_QC after filter: {}".format(df.shape[0]))

    df_NC = pd.read_csv(data_path + "Diag_NC.csv", low_memory=False, usecols=['eid'])
    df_MDD = pd.read_csv(data_path + "Diag_MDD.csv", low_memory=False, usecols=['eid'])
    df_Prote = pd.concat([df_NC, df_MDD])
    print("NMR NUM: {}".format(df_Prote.shape[0]))

    df_bed_sample = pd.read_csv(bed_path + "raw_data/c1.fam", sep=r'[ \t]', header=None, engine='python')
    df_bed_sample.rename(columns={0: 'fid', 1: 'eid', 2: 'pid', 3: 'mid', 4: 'sid'}, inplace=True)
    print("bed NUM: {}".format(df_bed_sample.shape[0]))

    df = pd.merge(df, df_Prote, on='eid', how='inner')
    df = reduce(lambda left, right: pd.merge(left, right, on=['eid'], how='inner'), [df, df_Prote, df_bed_sample])
    print("df: {}".format(df.shape[0]))

    df = df.sort_values(by='eid')  # 匹配filter_individual.sh 脚本执行后sample样本的顺序
    df.to_csv(bed_path + "{}.csv".format(exp_name), columns=['fid', 'eid'], sep='\t', header=None, index=False)


# .fam文件是储存样本信息的文件，为没有标题行的文本文件，每个样本占一行，共有以下6列：
# 1. 家系ID, Family ID (FID)
# 2. 个体ID, Within-family ID (IID;不能为0)
# 3. 父本个体ID(如果父本不在样本中则为0)
# 4. 母本个体ID(如果母本不在样本中则为0)
# 5. 性别(1为雄性,2为雌性, 0表示未知)
# 6. 表型值 (1为对照组,2为病例组,-9、0和非数值在二分类性状时表示缺失值)。如果存在除-9，0，1，2以外的任何数字表型值，则表型为数量性状，而不是病例/对照的二分类性状。在这种情况下，-9仍表示表型缺失值。
covar_num = 11
def split_cluster(l1, l2, f_name):
    rst_file = result_path + "{}/{}_components_cluster_rst_{}.csv".format("NMF", 2, 3)
    df_Y = pd.read_csv(rst_file, low_memory=False, usecols=['eid', 'final_cluster'])
    print(df_Y.groupby('final_cluster').size())
    df_NC = pd.read_csv(data_path + "Diag_NC.csv", low_memory=False, usecols=['eid'])
    df_NC['final_cluster'] = 0
    df = pd.concat([df_NC, df_Y])

    #  after QC bed file
    df_bed_sample = pd.read_csv(bed_path + "{}/GWAS_QC/{}_8.fam".format(exp_name, exp_name), sep=r'[ \t]', header=None, engine='python')
    df_bed_sample.rename(columns={0: 'fid', 1: 'eid', 2: 'pid', 3: 'mid', 4: 'sid'}, inplace=True)
    print("df_bed_sample: {}".format(df_bed_sample.shape[0]))

    df_covar = pd.read_csv(root_path + "Sample_QC.csv")
    key, value = ['22009-0.'+str(c) for c in range(1, covar_num)], ['PC'+str(c) for c in range(1, covar_num)]
    df_covar.rename(columns=dict(zip(key, value)), inplace=True)
    df_covar = df_covar[['eid']+value]
    print("df_covar: {}".format(df_covar.shape[0]))

    df = reduce(lambda left, right: pd.merge(left, right, on=['eid'], how='inner'), [df, df_bed_sample, df_covar])
    print("df after merge: {}".format(df.shape[0]))

    # 1 vs 2  ----- 1为对照组，2为病例组
    df_1, df_2 = df[df['final_cluster'].isin(l1)], df[df['final_cluster'].isin(l2)]
    df_1['final_cluster'], df_2['final_cluster'] = 1, 2
    weight = max(1, int(df_1.shape[0]/df_2.shape[0]))
    print("weight: {}".format(weight))
    df_1['weight'], df_2['weight'] = 1, weight
    df = pd.concat([df_1, df_2])
    print(df.groupby('final_cluster').size())
    print(df.head())

    df.rename(columns={'fid': '#FID', 'eid': 'IID'}, inplace=True)
    df = df.sort_values(by='IID')  # 匹配sample样本的顺序,方便后续进行GWAS分析
    df.to_csv(bed_path+"{}/{}_{}.csv".format(exp_name, exp_name, f_name), columns=['#FID', 'IID'], sep='\t', header=None, index=False)
    df.to_csv(bed_path+"{}/{}_{}_weight.txt".format(exp_name, exp_name, f_name), columns=['#FID', 'IID', 'weight'], sep='\t', header=None, index=False)
    df.to_csv(bed_path+"{}/{}_{}.fam".format(exp_name, exp_name, f_name), columns=['#FID', 'IID', 'pid', 'mid', 'sid', 'final_cluster'], sep='\t', header=None, index=False)
    df.to_csv(bed_path+"{}/{}_{}_covar.txt".format(exp_name, exp_name, f_name), columns=['#FID', 'IID']+value, sep='\t', index=False)


def cnt_sig_pvalue():
    # rst_l = ["GWAS_sample_cluster_1", "GWAS_sample_cluster_2", "GWAS_sample_cluster_3", "GWAS_sample_cluster_ALL", "GWAS_sample_cluster_1vs2", "GWAS_sample_cluster_1vs3", "GWAS_sample_cluster_3vs2"]
    f_s = ['0_VS_123', '0_VS_1', '0_VS_2', "0_VS_3", "1_VS_3", "1_VS_2", "3_VS_2"]
    for file_name in f_s:
        # df = pd.read_csv(bed_path + "{}/GWAS_{}/logistic_results.assoc_2.logistic".format(exp_name, file_name), sep=r'\s+', engine='python')
        df = pd.read_csv(bed_path + "{}/GWAS_{}/assoc_results.assoc".format(exp_name, file_name), sep=r'\s+', engine='python')
        df = df.sort_values(by="P")
        # df.to_csv(bed_path + "{}/GWAS_{}/logistic_results.sig.csv".format(exp_name, file_name), sep=',', index=False) # sep='\t'
        df.to_csv(bed_path + "{}/GWAS_{}/assoc_results.sig.csv".format(exp_name, file_name), sep=',', index=False)


if __name__ == '__main__':
    # 1. before QC
    # Sample_QC()

    # 2. before GWAS
    split_cluster([0], [1, 2, 3], '0_VS_123')
    split_cluster([0], [1], '0_VS_1')
    split_cluster([0], [2], '0_VS_2')
    split_cluster([0], [3], '0_VS_3')
    split_cluster([1], [3], '1_VS_3')
    split_cluster([1], [2], '1_VS_2')
    split_cluster([3], [2], '3_VS_2')

    # # 3. after GWAS
    # cnt_sig_pvalue()

