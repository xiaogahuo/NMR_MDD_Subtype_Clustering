from functools import reduce
import pandas as pd
from sklearn.model_selection import train_test_split, StratifiedKFold
from SUtils import *
import os


def cross_valid_exp(df, exp_Name, clu_list, hc_list):
    df_new_MDD, df_new_NC = df[df['final_cluster'].isin(clu_list)], df[df['final_cluster'].isin(hc_list)]
    df_new_MDD['final_cluster'], df_new_NC['final_cluster'] = 1, 0

    df_ALL = pd.concat([df_new_MDD, df_new_NC], ignore_index=True)
    print(df_ALL.groupby('final_cluster').size())
    # 重置索引
    df_ALL.reset_index(drop=True, inplace=True)
    # # 使用train_test_split分割DataFrame，确保'target'字段的比例在训练集和验证集中相似
    # train_df, test_df = train_test_split(df, test_size=0.3, stratify=df['final_cluster'], random_state=0)

    # 创建 StratifiedKFold 对象，设置折数为 5
    skf = StratifiedKFold(n_splits=5, random_state=42, shuffle=True)
    X, y = df_ALL[NMR_cols], df_ALL['final_cluster']
    # 使用 StratifiedKFold 对象的 split 方法进行数据集划分
    for fold, (train_index, test_index) in enumerate(skf.split(X, y)):
        if not os.path.exists(data_path + "cross_valid/{}/{}/".format(fold, exp_Name)):
            os.makedirs(data_path + "cross_valid/{}/{}/".format(fold, exp_Name))
        # print("TRAIN:", train_index, "TEST:", test_index)
        X_train, X_test = X.loc[train_index], X.loc[test_index]
        y_train, y_test = y.loc[train_index], y.loc[test_index]
        X_train.to_csv(data_path + "cross_valid/{}/{}/X_train.csv".format(fold, exp_Name), index=None)
        y_train.to_csv(data_path + "cross_valid/{}/{}/y_train.csv".format(fold, exp_Name), index=None)
        X_test.to_csv(data_path + "cross_valid/{}/{}/X_test.csv".format(fold, exp_Name), index=None)
        y_test.to_csv(data_path + "cross_valid/{}/{}/y_test.csv".format(fold, exp_Name), index=None)


if __name__ == '__main__':
    rst_file = result_path + "{}/{}_components_cluster_rst_{}.csv".format("NMF", 2, 3)
    df_Y = pd.read_csv(rst_file, low_memory=False, usecols=['eid', 'final_cluster'])
    df_MDD = pd.read_csv(data_path + "Diag_MDD.csv", low_memory=False)
    df_MDD = pd.merge(df_MDD, df_Y, on='eid', how='inner')

    df_NC = pd.read_csv(data_path + "Diag_NC.csv", low_memory=False)
    df_NC['final_cluster'] = 0
    df = pd.concat([df_MDD, df_NC], ignore_index=True)

    # # cross_valid_exp(df.copy(), 'cluster1_vs_0', [1], [0])
    # # cross_valid_exp(df.copy(), 'cluster2_vs_0',[2], [0])
    # # cross_valid_exp(df.copy(), 'cluster3_vs_0', [3], [0])
    # cross_valid_exp(df.copy(), 'cluster123_vs_0', [1,2,3], [0])
    # cross_valid_exp(df.copy(), 'cluster1_vs_other', [1], [0,2,3])
    # cross_valid_exp(df.copy(), 'cluster2_vs_other', [2], [0,1,3])
    # cross_valid_exp(df.copy(), 'cluster3_vs_other', [3], [0,1,2])
    #
    # df_NO_Chronic = df[df['chronic_num'] == 0]
    # cross_valid_exp(df_NO_Chronic.copy(), 'NO_Chronic_cluster123_vs_0', [1,2,3], [0])
    # cross_valid_exp(df_NO_Chronic.copy(), 'NO_Chronic_cluster1_vs_other', [1], [0, 2, 3])
    # cross_valid_exp(df_NO_Chronic.copy(), 'NO_Chronic_cluster2_vs_other', [2], [0, 1, 3])
    # cross_valid_exp(df_NO_Chronic.copy(), 'NO_Chronic_cluster3_vs_other', [3], [0, 1, 2])

    # df_Male = df[df['Sex_cate'] == 'L2_Male']
    # cross_valid_exp(df_Male.copy(), 'Male_cluster123_vs_0', [1, 2, 3], [0])
    # cross_valid_exp(df_Male.copy(), 'Male_cluster1_vs_other', [1], [0, 2, 3])
    # cross_valid_exp(df_Male.copy(), 'Male_cluster2_vs_other', [2], [0, 1, 3])
    # cross_valid_exp(df_Male.copy(), 'Male_cluster3_vs_other', [3], [0, 1, 2])
    #
    # df_Female = df[df['Sex_cate'] == 'L1_Female']
    # cross_valid_exp(df_Female.copy(), 'Female_cluster123_vs_0', [1, 2, 3], [0])
    # cross_valid_exp(df_Female.copy(), 'Female_cluster1_vs_other', [1], [0, 2, 3])
    # cross_valid_exp(df_Female.copy(), 'Female_cluster2_vs_other', [2], [0, 1, 3])
    # cross_valid_exp(df_Female.copy(), 'Female_cluster3_vs_other', [3], [0, 1, 2])

    df_BMI = df[df['BMI_cate'].isin(['BMI1', 'BMI3'])]
    cross_valid_exp(df_BMI.copy(), 'BMI_cluster123_vs_0', [1, 2, 3], [0])
    cross_valid_exp(df_BMI.copy(), 'BMI_cluster1_vs_other', [1], [0, 2, 3])
    cross_valid_exp(df_BMI.copy(), 'BMI_cluster2_vs_other', [2], [0, 1, 3])
    cross_valid_exp(df_BMI.copy(), 'BMI_cluster3_vs_other', [3], [0, 1, 2])
