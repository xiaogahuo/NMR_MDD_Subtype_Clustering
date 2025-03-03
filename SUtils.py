import pandas as pd
import os
import numpy as np
from scipy.stats import chi2_contingency, zscore
from scipy.stats.mstats import kruskalwallis
from sklearn.feature_selection import SelectKBest, f_regression
from scipy.stats import ttest_ind
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import norm
from functools import reduce


root_path = "/data/UKBioBank/FileHandler_And_Dataset/"
path = root_path + "NMR_MDD_Subtype_Clustering/"
data_path = path + "dataset/"
PRS_path = root_path + "PRS_MDD_Scores/"
img_path = path +"image/"
result_path = path +"result/"
bgen_path = "/data/UKBioBank/bgen_Data/"
bed_path = "/data/UKBioBank/bed_Data/"

NC, MDD, UNNORMAL = 0, 1, -1
cate100042_cols = ['2188-0.0','2178-0.0','2296-0.0','2306-0.0']
cate100057_cols = ['1160-0.0','1170-0.0','1180-0.0','1190-0.0','1200-0.0','1210-0.0','1220-0.0']
cate100060_cols = [c+'-0.0' for c in ['20123','20124','20125','4609','4620','5375','5386','1920','1930','1940','1950','1960','1970','1980','1990','2000','2010','2020','2030','2040','2090','2100','4598','4631','4642','4653','4526','4548','4559','4570','4581','4537','2050','2060','2070','2080','6156','5663','5674','6145','10721','20122','20126','20127']]
cate100061_cols = [c+'-0.0' for c in ['1031','6160','2110','10740']]
cate100065_cols = [c+'-0.0' for c in ['21000']]
cate17518_cols = [c+'-0.0' for c in ['30600','30610','30620','30630','30640','30650','30660','30670','30680','30690','30700','30710','30720','30730','30740','30750','30760','30770','30780','30790','30800','30810','30820','30830','30840','30850','30860','30870','30880','30890']]
cate17518_new_cols = ['ALB','ALP','ALT','APOA1','APOB','AST','DBIL','UREA','CALC','CHOL','CREA','CRP','CYSC','GGT','GLU','HbA1c','HDL','IGF-1','LDL','LP-a','OEST','PHOS','RF','SHBG','TBIL','TEST','TP','TRIG','UA','VITD']

cate145_cols = [c+'-0.0' for c in ['20487','20488','20489','20490','20491','20521','20522','20523','20524','20525','20526','20527','20528','20529','20530','20531','20494','20495','20496','20497','20498']]

cat_cols = ['Sex_cate', 'Townsend_index_cate', 'Qualifications_cate', 'Smoking_status_cate', 'Alcohol_status_cate', 'BMI_cate', 'IPAQ_cate', 'family_history', 'Drug']
conti_cols = ['Age', 'chronic_num']
covr_cols = cat_cols + conti_cols

cols_MRI_c192_Area = [26721, 26822, 26722, 26823, 26723, 26824, 26724, 26825, 26725, 26826, 26726, 26827, 26752, 26853, 26727, 26828, 26728, 26829, 26729, 26830, 26754, 26855, 26730, 26831, 26731, 26832, 26732, 26833, 26733, 26834, 26734, 26835, 26735, 26836, 26737, 26838, 26736, 26837, 26738, 26839, 26739, 26840, 26740, 26841, 26741, 26842, 26742, 26843, 26743, 26844, 26744, 26845, 26745, 26846, 26746, 26847, 26747, 26848, 26748, 26849, 26749, 26850, 26750, 26851, 26751, 26852, 26753, 26854]
cols_MRI_c192_Thickness = [26755, 26856, 26756, 26857, 26757, 26858, 26758, 26859, 26759, 26860, 26760, 26861, 26786, 26887, 26761, 26862, 26762, 26863, 26763, 26864, 26788, 26889, 26764, 26865, 26765, 26866, 26766, 26867, 26767, 26868, 26768, 26869, 26769, 26870, 26771, 26872, 26770, 26871, 26772, 26873, 26773, 26874, 26774, 26875, 26775, 26876, 26776, 26877, 26777, 26878, 26778, 26879, 26779, 26880, 26780, 26881, 26781, 26882, 26782, 26883, 26783, 26884, 26784, 26885, 26785, 26886, 26787, 26888]
cols_MRI_c192_Volume = [26789, 26890, 26790, 26891, 26791, 26892, 26792, 26893, 26793, 26894, 26819, 26920, 26794, 26895, 26795, 26896, 26796, 26897, 26821, 26922, 26797, 26898, 26798, 26899, 26799, 26900, 26800, 26901, 26801, 26902, 26802, 26903, 26804, 26905, 26803, 26904, 26805, 26906, 26806, 26907, 26807, 26908, 26808, 26909, 26809, 26910, 26810, 26911, 26811, 26912, 26812, 26913, 26813, 26914, 26814, 26915, 26815, 26916, 26816, 26917, 26817, 26918, 26818, 26919, 26820, 26921]

cols_MRI_c135_MD = [25515,25516,25517,25518,25519,25520,25523,25524,25525,25526,25527,25528,25529,25530,25532,25533,25531,25521,25522,25534,25535,25536,25537,25538,25539,25540,25541]
cols_MRI_c135_FA = [25488,25489,25490,25491,25492,25493,25496,25497,25498,25499,25500,25501,25502,25503,25505,25506,25504,25494,25495,25507,25508,25509,25510,25511,25512,25513,25514]
cols_MRI_c1102_subcortical_volume = [25023, 25024, 25021, 25022, 25013, 25014, 25019, 25020, 25017, 25018, 25015, 25016, 25011, 25012]
cols_MRI_field_25753 = ["col"+str(i) for i in range(0, 1485)]
# print(cols_MRI_field_25753)

def get_NMR_names():
    df = pd.read_excel(path+'Feature_Name.xlsx', sheet_name=0, header=0)
    keys, values = [str(val)+'-0.0' for val in df.iloc[:, 0].tolist()], df.iloc[:, 2].tolist()
    return keys, values
raw_NMR_cols, NMR_cols = get_NMR_names()

TRL = ['VLDL_C', 'VLDL_TG', 'VLDL_PL', 'VLDL_CE', 'VLDL_FC', 'VLDL_L', 'VLDL_P', 'VLDL_size', 'IDL_P', 'IDL_L', 'IDL_PL', 'IDL_C', 'IDL_CE', 'IDL_FC', 'IDL_TG']
VLDL = ['VLDL_C', 'VLDL_TG', 'VLDL_PL', 'VLDL_CE', 'VLDL_FC', 'VLDL_L', 'VLDL_P', 'VLDL_size']
TRL_TG = ['VLDL_TG','XXL_VLDL_TG','XL_VLDL_TG','L_VLDL_TG','M_VLDL_TG','S_VLDL_TG','XS_VLDL_TG','IDL_TG']
TRL_C = ['VLDL_C','XXL_VLDL_C','XL_VLDL_C','L_VLDL_C','M_VLDL_C','S_VLDL_C','XS_VLDL_C','IDL_C']
RC = ['Total_C','LDL_C','HDL_C', 'Remnant_C', 'non_HDL_C']
TRL_P =['VLDL_P','XXL_VLDL_P','XL_VLDL_P','L_VLDL_P','M_VLDL_P','S_VLDL_P','XS_VLDL_P','IDL_P']
TRL_PL = ['VLDL_PL','XXL_VLDL_PL','XL_VLDL_PL','L_VLDL_PL','M_VLDL_PL','S_VLDL_PL','XS_VLDL_PL','IDL_PL']
TRL_FC = ['VLDL_FC','XXL_VLDL_FC','XL_VLDL_FC','L_VLDL_FC','M_VLDL_FC','S_VLDL_FC','XS_VLDL_FC','IDL_FC']
TRL_CE = ['VLDL_CE','XXL_VLDL_CE','XL_VLDL_CE','L_VLDL_CE','M_VLDL_CE','S_VLDL_CE','XS_VLDL_CE','IDL_CE']
TRL_L = ['VLDL_L','XXL_VLDL_L','XL_VLDL_L','L_VLDL_L','M_VLDL_L','S_VLDL_L','XS_VLDL_L','IDL_L']
HDL = ['HDL_C', 'HDL_TG', 'HDL_PL', 'HDL_CE', 'HDL_FC', 'HDL_L', 'HDL_P', 'HDL_size']
LDL = ['LDL_C', 'LDL_TG', 'LDL_PL', 'LDL_CE', 'LDL_FC', 'LDL_L', 'LDL_P', 'LDL_size']
IDL = ['IDL_C', 'IDL_TG', 'IDL_PL', 'IDL_CE', 'IDL_FC', 'IDL_L', 'IDL_P']

C = ['HDL_C', 'LDL_C', 'IDL_C', 'VLDL_C', 'Total_C']
TG = ['HDL_TG', 'LDL_TG', 'IDL_TG', 'VLDL_TG', 'Total_TG']
PL = ['HDL_PL', 'LDL_PL', 'IDL_PL', 'VLDL_PL', 'Total_PL']
CE = ['HDL_CE', 'LDL_CE', 'IDL_CE', 'VLDL_CE', 'Total_CE']
FC = ['HDL_FC', 'LDL_FC', 'IDL_FC', 'VLDL_FC', 'Total_FC']
L = ['HDL_L', 'LDL_L', 'IDL_L', 'VLDL_L', 'Total_L']
P = ['HDL_P', 'LDL_P', 'IDL_P', 'VLDL_P', 'Total_P']
Size = ['HDL_size', 'LDL_size', 'VLDL_size']
Total = ['Total_C', 'Total_TG', 'Total_PL', 'Total_CE', 'Total_FC', 'Total_L', 'Total_P']
FA = ['Total_FA','Unsaturation','Omega_3','Omega_6','PUFA','MUFA','SFA','LA','DHA','Omega_3(%)','Omega_6(%)','PUFA(%)','MUFA(%)','SFA(%)','LA(%)','DHA(%)','PUFA_by_MUFA','Omega_6_by_Omega_3']

# 将labels转换为numpy数组，并使用广播和向量化操作来计算相似度矩阵。这样可以避免使用双重循环，提高代码性能。
def compute_similarity_matrix(labels):
    labels = np.array(labels)
    similarity_matrix = (labels[:, None] == labels).astype(int)
    return similarity_matrix


def annot_trans(pval):
    if pval < 0.001:
        return '***'
    elif pval < 0.01:
        return '**'
    elif pval < 0.05:
        return '*'
    else:
        return ' '


def social_stats(exps, cmp_col):  # 'MDD'
    for exp in exps:
        df = pd.read_csv(data_path + exp + '.csv')
        print('---------------exp: {}, number: {}------------------------'.format(exp, df.shape[0]))
        for col in cat_cols+[cmp_col]:
            count_1feature(col, df)
        for col in conti_cols:
            count_mean_std(col, df)

        for col in cat_cols:
            count_2features(col, cmp_col, df)
            count_Chi2(col, cmp_col, df)
        for col in conti_cols:
            count_mean_std_group_by_div_feature(col, cmp_col, df)
            kwtest(col, cmp_col, df)


def count_1feature(feature_name, df):
    print('----------打印字段: <' + feature_name + "> 的人数分组统计情况----------")
    target_feature_value_list = df.drop_duplicates(subset=[feature_name], keep='first')[feature_name].tolist()
    target_feature_value_list.sort()
    for value in target_feature_value_list:
        print("字段< %s > 类别 < %s > 人数统计: %d" % (feature_name, str(value), df[df[feature_name] == value].shape[0]))


def count_2features(feature_name1, feature_name2, df):
    print('----------打印字段: <' + feature_name1 + "> 和 字段< " + feature_name2 + " > 的分组统计情况----------")
    target_feature_value_list1 = df.drop_duplicates(subset=[feature_name1], keep='first')[feature_name1].tolist()
    target_feature_value_list1.sort()
    target_feature_value_list2 = df.drop_duplicates(subset=[feature_name2], keep='first')[feature_name2].tolist()
    target_feature_value_list2.sort()

    for value1 in target_feature_value_list1:
        for value2 in target_feature_value_list2:
            print("字段< %s > 类别 < %s > 字段< %s > 类别 < %s > 人数统计: %d" % (
            feature_name1, str(value1), feature_name2, str(value2),
            df[(df[feature_name1] == value1) & (df[feature_name2] == value2)].shape[0]))


def count_Chi2(cnt_feature, div_feature, df):
    print('卡方检验结果: <' + cnt_feature + "> on <" + div_feature + ">")
    # print(df.drop_duplicates(subset=[target_feature], keep='first')[target_feature])
    div_feature_value_list = sorted(df.drop_duplicates(subset=[div_feature], keep='first')[div_feature].tolist())
    cnt_feature_value_list = sorted(df.drop_duplicates(subset=[cnt_feature], keep='first')[cnt_feature].tolist())
    print("{}: {}".format(div_feature, div_feature_value_list))
    print("{}: {}".format(cnt_feature, cnt_feature_value_list))

    kf_data = []
    for div in div_feature_value_list:
        c_g =[]
        for value in cnt_feature_value_list:
            count = df[(df[div_feature] == div) & (df[cnt_feature] == value)].shape[0]
            c_g.append(count)
        kf_data.append(c_g)
    print(kf_data)
    kf = chi2_contingency(kf_data)
    print('chisq-statistic=%.4f, p-value=%.4f, df=%i expected_frep=%s' % kf)
    print('-'*25)
    return kf, kf_data


def kwtest(feature_name, target_feature, df):
    print('kwtest结果: <' + feature_name + "> on <" + target_feature + ">")
    # print(df.drop_duplicates(subset=[target_feature], keep='first')[target_feature])
    target_feature_value_list = df.drop_duplicates(subset=[target_feature], keep='first')[target_feature].tolist()

    kw_data = []
    for value in target_feature_value_list:
        feature_name_value_list = df[df[target_feature] == value][feature_name].tolist()
        kw_data.append(feature_name_value_list)
    statistic, pvalue = kruskalwallis(*kw_data) #  将kw_data中的值作为参数传入kruskalwallis中, 第1个值是第1个参数， 第2个值是第2个参数
    print("column: {}, statistic: {}, pvalue: {}".format(feature_name, statistic, pvalue))
    print('-'*25)


def count_mean_std_group_by_div_feature(feature_name, div_feature_name, df):
    print(
        '----------打印字段: <' + feature_name + "> 根据字段< " + div_feature_name + " > 的分组 均值(Mean) 标准差(Std) 统计情况----------")
    print('字段< %s > 总体Mean: %f' % (feature_name, df[feature_name].mean()))
    print('字段< %s > 总体Std: %f' % (feature_name, df[feature_name].std(ddof=1)))
    print('-' * 25)
    div_feature_value_list = df.drop_duplicates(subset=[div_feature_name], keep='first')[div_feature_name].tolist()
    div_feature_value_list.sort()

    for div_value in div_feature_value_list:
        print("字段< %s > 按字段< %s > 类别 < %s > 统计的均值 (Mean): %f" % (
            feature_name, div_feature_name, str(div_value), df[df[div_feature_name] == div_value][feature_name].mean()))
        print("字段< %s > 按字段< %s > 类别 < %s > 统计的标准差 (Std): %f" % (
            feature_name, div_feature_name, str(div_value),
            df[df[div_feature_name] == div_value][feature_name].std(ddof=1)))


def count_mean_std(feature_name, df):
    print('----------打印字段: <' + feature_name + "> 的均值(Mean) 标准差(Std) 统计情况----------")
    print('字段< %s > 总体Mean: %f' % (feature_name, df[feature_name].mean()))
    print('字段< %s > 总体Std: %f' % (feature_name, df[feature_name].std(ddof=1)))
    print('-' * 25)


def count_T(feature_name, target_feature, df):
    count_mean_std_group_by_div_feature(feature_name, target_feature, df)
    print('T检验结果: <' + feature_name + "> on <" + target_feature + ">")
    # print(df.drop_duplicates(subset=[target_feature], keep='first')[target_feature])
    target_feature_value_list = df.drop_duplicates(subset=[target_feature], keep='first')[target_feature].tolist()
    c1 = target_feature_value_list[0]
    c2 = target_feature_value_list[1]
    v1_list = df[df[target_feature] == c1][feature_name].tolist()
    v2_list = df[df[target_feature] == c2][feature_name].tolist()
    v1_list = [v for v in v1_list if not pd.isnull(v)]
    v2_list = [v for v in v2_list if not pd.isnull(v)]
    print("size:{}".format(len(v1_list)))
    print("size:{}".format(len(v2_list)))

    # Run a two sample t-test to compare the two samples
    tstat, pval = ttest_ind(a=v1_list, b=v2_list, alternative="two-sided")

    # Display results
    print("t-stat: {:.2f}   pval: {:.4f}".format(tstat, pval))
    return tstat, pval

# "#D7191C" "#FDAE61" "#FFFFBF" "#ABDDA4" "#2B83BA"
palette = {"HC": "#2B83BA", "Subtpye1": "#ABDDA4", "Subtpye2": "#FDAE61", "Subtpye3": "#D7191C", "MDD": "#FFFFBF"}
def plot_bar(df, file_name, feature_lst, width, high, hue_order, y_label):
    df_bar = []
    for feature in feature_lst:
        df_sub = df[[feature, 'final_cluster']].copy()
        df_sub['Features'] = feature
        df_sub.rename(columns={feature: y_label}, inplace=True)
        # print(df_sub.head())
        df_bar.append(df_sub)
    df_bar = pd.concat(df_bar, ignore_index=True)
    # print("before: {}".format(df_bar.shape))
    print(df_bar.head())
    # df_bar = df_bar[(df_bar['z-score']<=5) & (df_bar['z-score']>=-5)]
    print(" after: {}".format(df_bar.shape))

    plt.figure(figsize=(width, high))
    ax = sns.barplot(hue="final_cluster", hue_order=hue_order, y=y_label, x="Features", data=df_bar,
                     capsize=.2, errorbar=("ci", 95), err_kws={"color": ".9", "linewidth": 0.6},
                palette=palette, width=.5, gap=0.1) # errorbar=None
    ax.tick_params(axis='x', labelrotation=0, labelsize=6)
    ax.tick_params(axis='y', labelrotation=0, labelsize=6)
    plt.legend(loc='upper center', bbox_to_anchor=(1.07, 1.02))  # 使用这行代码把图例放到图外。
    plt.savefig(img_path + "Bar_{}.png".format(file_name), dpi=720)
    plt.clf()


def plot_violin(df, file_name, feature_lst, width, high, hue_order, y_label):
    df_violin = []
    for feature in feature_lst:
        df_sub = df[[feature, 'final_cluster']].copy()
        df_sub['Features'] = feature
        df_sub.rename(columns={feature: y_label}, inplace=True)
        # print(df_sub.head())
        df_violin.append(df_sub)
    df_violin = pd.concat(df_violin, ignore_index=True)
    # print("before: {}".format(df_violin.shape))
    print(df_violin.head())

    plt.figure(figsize=(width, high))
    ax = sns.violinplot(hue="final_cluster", hue_order=hue_order, y=y_label, x="Features", data=df_violin,
                palette=palette, width=.9, gap=0.1, legend=False) # errorbar=None
    ax.tick_params(axis='x', labelrotation=0, labelsize=8)
    ax.tick_params(axis='y', labelrotation=0, labelsize=8)
    plt.legend(loc='upper center', bbox_to_anchor=(1.07, 1.02))  # 使用这行代码把图例放到图外。
    plt.savefig(img_path + "Violin_{}.png".format(file_name), dpi=720)
    plt.clf()


def plot_box(df, file_name, feature_lst, width, high, hue_order, y_label):
    df_box = []
    for feature in feature_lst:
        df_sub = df[[feature, 'final_cluster']].copy()
        df_sub['Features'] = feature
        df_sub.rename(columns={feature: y_label}, inplace=True)
        # print(df_sub.head())
        df_box.append(df_sub)
    df_box = pd.concat(df_box, ignore_index=True)
    # print("before: {}".format(df_box.shape))
    # print(df_box.head())
    # df_box = df_box[(df_box['z-score']<=3) & (df_box['z-score']>=-3)]
    # df_box.to_csv("{}{}_{}.csv".format(data_path, file_name, sub), index=None)  # Dataset 2
    # print(" after: {}".format(df_box.shape))

    # 绘制箱型图
    plt.figure(figsize=(width, high))
    ax = sns.boxplot(hue="final_cluster", hue_order=hue_order, y=y_label, x="Features", data=df_box,
                     width=0.6, linewidth=0.5, gap=0.1, legend=False, fliersize=0.5)
    # ax = sns.barplot(data=df_box, y="z-score", x="Features", hue="MDD", hue_order=['No', 'Yes'] )
    ax.tick_params(axis='x', labelrotation=15, labelsize=6)
    ax.tick_params(axis='y', labelrotation=0, labelsize=6)
    plt.legend(loc='upper center', bbox_to_anchor=(1.47, 1.02))  # 使用这行代码把图例放到图外。
    plt.savefig(img_path + "Box_{}.png".format(file_name), dpi=720)
    plt.clf()  # CLEAR FIGURE


def cnt_OR_CI95(f_list):
    for f_name in f_list:
        z = norm.ppf(0.975)  # 对应于95%置信水平的z值
        df = pd.read_csv(result_path+'{}.csv'.format(f_name))
        # df['clust1_OR'] = df.apply(lambda row: np.exp(row['clust1_Estimate']), axis=1)
        df['clust1_lower95'] = df.apply(lambda row: row['clust1_Estimate']-z*row['clust1_StdError'], axis=1)
        df['clust1_upper95'] = df.apply(lambda row: row['clust1_Estimate']+z*row['clust1_StdError'], axis=1)
        # df['clust1'] = 'clust1'

        # df['clust2_OR'] = df.apply(lambda row: np.exp(row['clust2_Estimate']), axis=1)
        df['clust2_lower95'] = df.apply(lambda row: row['clust2_Estimate']-z*row['clust2_StdError'], axis=1)
        df['clust2_upper95'] = df.apply(lambda row: row['clust2_Estimate']+z*row['clust2_StdError'], axis=1)
        # df['clust2'] = 'clust2'

        # df_rst = pd.DataFrame({'Feature': df['Col'].tolist()+df['Col'].tolist(),
        #                        'OR': df['clust1_OR'].tolist()+df['clust2_OR'].tolist(),
        #                        'CI_lower': df['clust1_lower95'].tolist()+df['clust2_lower95'].tolist(),
        #                        'CI_upper': df['clust1_upper95'].tolist()+df['clust2_upper95'].tolist(),
        #                        'Clust': df['clust1'].tolist()+df['clust2'].tolist()})
        df.to_csv(result_path+'{}.csv'.format(f_name), index=False)


if __name__ == '__main__':
    # 创建一个示例dataframe
    # rst = set(pd.read_csv(path + "proteomics.csv").columns.tolist())-set(Prote_cols_all)
    # print(rst)

    # df1 = pd.read_excel(path+'protein2923.xlsx', sheet_name='Sheet1', usecols=['coding', 'meaning'])
    # df2 = pd.read_excel(path+'protein2923.xlsx', sheet_name='Sheet2', usecols=['Assay Target', 'Protein panel'])
    # df = pd.merge(df1, df2, left_on='coding', right_on='Assay Target', how='inner')
    # df.to_csv(path + "protein2923.csv", index=None)
    # print(df.shape[0])

    pass
