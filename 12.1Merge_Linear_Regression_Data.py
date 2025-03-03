import pandas as pd

from SUtils import *
from functools import reduce
from scipy.stats.mstats import kruskal

def merge_data():
    df_MDD, df_NC = pd.read_csv(data_path + "Diag_MDD.csv", low_memory=False), pd.read_csv(data_path + "Diag_NC.csv", low_memory=False)
    df_NC['final_cluster'] = 0
    df_Y = pd.read_csv(result_path + "{}/{}_components_cluster_rst_{}.csv".format("NMF", 2, 3), low_memory=False,
                       usecols=['eid', 'final_cluster'])
    df_MDD = pd.merge(df_MDD, df_Y, how='inner', on='eid', sort=True, suffixes=('_x', '_y'), copy=True, indicator=False)
    df = pd.concat([df_MDD, df_NC], ignore_index=True)

    # cols_MRI = pd.read_excel(path + "Cate192_cols.xlsx")['FieldID'].tolist()
    # print(cols_MRI)
    cols_MRI_192 = [str(c) + '-2.0' for c in cols_MRI_c192_Area+cols_MRI_c192_Thickness+cols_MRI_c192_Volume]
    df_MRI_192 = pd.read_csv(root_path + "Cate_192.csv", low_memory=False, usecols=['eid']+cols_MRI_192)
    df_MRI_192.dropna(thresh=2, inplace=True)
    print("df_MRI_192: {}".format(df_MRI_192.shape[0]))

    cols_MRI_135 = [str(c) + '-2.0' for c in cols_MRI_c135_MD+cols_MRI_c135_FA]
    df_MRI_135 = pd.read_csv(root_path + "Cate_135.csv", low_memory=False, usecols=['eid']+cols_MRI_135)
    df_MRI_135.dropna(thresh=2, inplace=True)
    print("df_MRI_135: {}".format(df_MRI_135.shape[0]))

    cols_MRI_110 = [str(c) + '-2.0' for c in cols_MRI_c1102_subcortical_volume]
    df_MRI_110 = pd.read_csv(root_path + "Cate_110.csv", low_memory=False, usecols=['eid'] + cols_MRI_110)
    df_MRI_110.dropna(thresh=2, inplace=True)
    print("df_MRI_110: {}".format(df_MRI_110.shape[0]))

    # df_rsfMRI = pd.read_csv(path + "25753_2_0.csv", low_memory=False, usecols=['eid']+cols_MRI_field_25753)
    # print("df_rsfMRI: {}".format(df_rsfMRI.shape[0]))
    # grouped = df_rsfMRI.groupby('final_cluster').size()
    # print(grouped)

    df_MRI = reduce(lambda left, right: pd.merge(left, right, on=['eid'], how='inner'), [df[['eid', 'final_cluster']+covr_cols], df_MRI_192, df_MRI_135, df_MRI_110])
    print("df_MRI: {}".format(df_MRI.shape[0]))
    for col in cols_MRI_192+cols_MRI_135+cols_MRI_110:
        df_MRI[col] = (df_MRI[col] - df_MRI[col].mean())/df_MRI[col].std()
    df_MRI.fillna(0, inplace=True)  # 空值补零- 即均值
    cols_MRI = cols_MRI_192+cols_MRI_135+cols_MRI_110
    new_cols = ['f' + c[:-4] for c in cols_MRI]
    print(new_cols)
    df_MRI.rename(columns=dict(zip(cols_MRI, new_cols)), inplace=True)
    print(df_MRI.groupby('final_cluster').size())
    df_MRI.to_csv(data_path + "Linear_Regression_MRI_Dataset.csv", index=False)


def social_stats(exps):
    for exp in exps:
        df = pd.read_csv(data_path + exp + '.csv')
        print('---------------exp: {}, number: {}------------------------'.format(exp, df.shape[0]))
        for col in cat_cols+['final_cluster']:
            count_1feature(col, df)
        for col in conti_cols:
            count_mean_std(col, df)

        for col in cat_cols:
            count_2features(col, 'final_cluster', df)
            count_Chi2(col, 'final_cluster', df)
        for col in conti_cols:
            count_mean_std_group_by_div_feature(col, 'final_cluster', df)
            kwtest(col, 'final_cluster', df)


def kw_compare_MRI():
    df = pd.read_csv(data_path + "Linear_Regression_MRI_Dataset.csv", low_memory=False)
    f = open(path + "result/compare_kw_result_MRI.csv", "w")
    header_str = "col,statistic,pvalue,M_all,STD_all"
    div_feature_value_list = sorted(
        df.drop_duplicates(subset=['final_cluster'], keep='first')['final_cluster'].tolist())
    for c in div_feature_value_list:
        header_str = header_str + ',M{},STD{}'.format(c, c)
    print(header_str, file=f)
    print(header_str)
    cols_MRI_192 = [str(c) + '-2.0' for c in cols_MRI_c192_Area + cols_MRI_c192_Thickness + cols_MRI_c192_Volume]
    cols_MRI_135 = [str(c) + '-2.0' for c in cols_MRI_c135_MD + cols_MRI_c135_FA]
    cols_MRI_110 = [str(c) + '-2.0' for c in cols_MRI_c1102_subcortical_volume]
    cols_MRI = cols_MRI_192 + cols_MRI_135 + cols_MRI_110
    cols_MRI = ['f' + c[:-4] for c in cols_MRI]
    for i, col in zip(range(len(cols_MRI)), cols_MRI):
        print("MRI, index:{}, col:{}".format(i, col))
        # 按final_cluster列的值对col列进行分组
        grp_list = []
        for c in div_feature_value_list:
            grp_list.append(df[df['final_cluster'] == c][col].to_list())
        statistic, pvalue = kruskal(*grp_list)
        M_all, STD_all = df[col].mean(), df[col].std(ddof=1)
        val_str = col + ',' + str(statistic) + ',' + str(pvalue) + ',' + str(M_all) + ',' + str(STD_all)
        for c in div_feature_value_list:
            M, STD = df[df['final_cluster'] == c][col].mean(), df[df['final_cluster'] == c][col].std(ddof=1)
            val_str = val_str + ',' + str(M) + ',' + str(STD)
        print(val_str, file=f)


if __name__ == '__main__':
    # merge_data()

    # 4. 社会人口学数据统计
    exps = ['Linear_Regression_MRI_Dataset']
    social_stats(exps)

    kw_compare_MRI()



