from scipy.stats.mstats import kruskal
import pandas as pd
from SUtils import *
import math


def social_clust_stats(df, cmp_col):  # 'final_cluster'
    print('---------------exp: cluster, number: {}------------------------'.format(df.shape[0]))
    # for col in cat_cols+[cmp_col]:
    #     count_1feature(col, df)
    # # for col in conti_cols+NMR_cols:
    for col in conti_cols:
        count_mean_std(col, df)

    # for col in cat_cols:
    #     count_2features(col, cmp_col, df)
    #     count_Chi2(col, cmp_col, df)
    # # for col in conti_cols+NMR_cols:
    for col in conti_cols:
        count_mean_std_group_by_div_feature(col, cmp_col, df)
        kwtest(col, cmp_col, df)


def compare_KW(df, cols, group_by_col, fname):
    f = open(result_path + "compare_kw_result_{}.csv".format(fname), "w")
    header_str = "col,statistic,pvalue,M_all,STD_all"
    div_feature_value_list = sorted(df.drop_duplicates(subset=[group_by_col], keep='first')[group_by_col].tolist())
    for c in div_feature_value_list:
        header_str = header_str + ',M{},STD{}'.format(c,c)
    print(header_str, file=f)
    print(header_str)
    for i, col in zip(range(len(cols)), cols):
        print("clust_N:{}, index:{}, col:{}".format(fname, i, col))
        # 按final_cluster列的值对col列进行分组
        grp_list =[]
        for c in div_feature_value_list:
            grp_list.append(df[df[group_by_col]==c][col].to_list())
        statistic, pvalue = kruskal(*grp_list)
        M_all, STD_all = df[col].mean(), df[col].std(ddof=1)
        val_str = col + ',' + str(statistic) + ',' + str(pvalue) + ',' + str(M_all) + ',' + str(STD_all)
        for c in div_feature_value_list:
            M, STD = df[df[group_by_col]==c][col].mean(), df[df[group_by_col]==c][col].std(ddof=1)
            val_str = val_str + ',' + str(M) + ',' + str(STD)
        print(val_str, file=f)

        
def compare_Chi2(df, fname):
    f = open(result_path + "compare_kafang_result_{}.csv".format(fname), "w")
    header_str = "col,statistic,pvalue,degree_freedom,cnt_feature_values,chi2_data"
    print(header_str, file=f)
    cate_100060=[c+'-0.0' for c in ['20123','20124','20125','1920','1930','1940','1950','1960','1970','1980','1990','2000','2010','2020','2030','2040','2090','2100','4598','4631','4642','4653','4526','4548','4559','4570','4581','4537','6156','5663','5674','6145','10721','20122','20126']]
    cols =cat_cols+ cate100042_cols + cate100057_cols[1:] + cate_100060 + cate100061_cols + cate100065_cols + cate145_cols
    for i, col in zip(range(len(cols)), cols):
        print("--------clust_N:{}, index:{}, col:{}--------".format(fname, i, col))
        null_ratio = df[col].isnull().mean()
        print("{}: null_ratio: {}".format(col, null_ratio))
        if null_ratio >= 0.5:
            continue
        print('卡方检验结果: <' + col + "> on < final_cluster >")
        div_feature_value_list = sorted(df.drop_duplicates(subset=['final_cluster'], keep='first')['final_cluster'].tolist())
        cnt_feature_value_list = df.drop_duplicates(subset=[col], keep='first')[col].tolist()
        cnt_feature_value_series = pd.Series(cnt_feature_value_list)
        cnt_feature_value_series = cnt_feature_value_series.dropna()
        cnt_feature_value_list = sorted(cnt_feature_value_series.to_list())
        print("{}: {}".format('final_cluster', div_feature_value_list))
        print("{}: {}".format(col, cnt_feature_value_list))

        kf_data = []
        for div in div_feature_value_list:
            c_g =[]
            for value in cnt_feature_value_list:
                count = df[(df['final_cluster'] == div) & (df[col] == value)].shape[0]
                c_g.append(count)
            kf_data.append(c_g)
        print(kf_data)
        statistic, pvalue, dof, expected_freq = chi2_contingency(kf_data)
        print("{},{},{},{},{},{}".format(col, statistic, pvalue, dof, cnt_feature_value_list, kf_data), file=f)


def bar_biochem_plot(df):
    # 使用 replace 方法实现 cluster2 和 cluster3 的互换
    df['final_cluster'] = df['final_cluster'].replace({2: 3, 3: 2})
    df['final_cluster'] = df['final_cluster'].replace({0: 'HC', 1: 'Subtpye1', 2: 'Subtpye2', 3: 'Subtpye3'})
    df_biochem = df[cate17518_cols+['final_cluster']]
    df_biochem.rename(columns=dict(zip(cate17518_cols, cate17518_new_cols)), inplace=True)

    for col in cate17518_new_cols:
        df_biochem[col] = (df_biochem[col] - df_biochem[col].mean())/df_biochem[col].std()
    plot_bar(df_biochem, 'biochem_1', ['IGF-1','ALP','ALT','APOA1','AST','CRP','CYSC','GGT','HbA1c','HDL','SHBG','TRIG','UA','VITD'], 12, 4, ['HC', 'Subtpye3', 'Subtpye1', 'Subtpye2'], 'Z-score')
    plot_bar(df_biochem, 'biochem_2', ['ALB','APOB','DBIL','UREA','CALC','CHOL','CREA','GLU','LDL','LP-a','PHOS','TBIL','TEST','TP'], 12, 4, ['HC', 'Subtpye3', 'Subtpye1', 'Subtpye2'], 'Z-score')


def PHQ_plot(df):
    df['final_cluster'] = df['final_cluster'].replace({0: 'HC', 1: 'Subtpye1', 2: 'Subtpye2', 3: 'Subtpye3'})
    df_PHQ = df[['20127-0.0', '2050-0.0', '2060-0.0', '2070-0.0', '2080-0.0', 'PHQ', 'PHQ_Total', 'final_cluster']]
    for col in ['20127-0.0', '2050-0.0', '2060-0.0', '2070-0.0', '2080-0.0', 'PHQ', 'PHQ_Total']:
        plot_violin(df_PHQ, col, [col], 2.8, 4, ['HC', 'Subtpye1', 'Subtpye2', 'Subtpye3'], 'Score')

    df_GlycA = df[['GlycA', 'final_cluster']]
    df_GlycA = df_GlycA[(df_GlycA['GlycA'] > -4) & (df_GlycA['GlycA'] < 5)]
    plot_violin(df_GlycA, 'GlycA', ['GlycA'], 2.8, 4, ['HC', 'Subtpye1', 'Subtpye2', 'Subtpye3'], 'Score')

    df_CRP = df[['30710-0.0', 'final_cluster']]
    df_CRP = df_CRP[df_CRP['30710-0.0'] < 10]
    plot_violin(df_CRP, '30710-0.0', ['30710-0.0'], 2.8, 4, ['HC', 'Subtpye1', 'Subtpye2', 'Subtpye3'], 'Score')

    df_1160 = df[['1160-0.0', 'final_cluster']]
    df_1160 = df_1160[(df_1160['1160-0.0'] > -4) & (df_1160['1160-0.0'] < 14)]
    plot_violin(df_1160, '1160-0.0', ['1160-0.0'], 2.8, 4, ['HC', 'Subtpye1', 'Subtpye2', 'Subtpye3'], 'Score')

if __name__ == '__main__':
    df_MDD, df_NC = pd.read_csv(data_path + "Diag_MDD.csv", low_memory=False), pd.read_csv(data_path + "Diag_NC.csv", low_memory=False)
    df_NC['final_cluster'] = 0
    df_Y = pd.read_csv(result_path + "{}/{}_components_cluster_rst_{}.csv".format("NMF", 2, 3), low_memory=False, usecols=['eid', 'final_cluster'])
    df_Y['final_cluster'] = df_Y['final_cluster'].replace({2: 3, 3: 2})
    df_MDD = pd.merge(df_MDD, df_Y, how='inner', on='eid', sort=True, suffixes=('_x', '_y'), copy=True, indicator=False)
    df = pd.concat([df_MDD, df_NC])
    print(df.columns.tolist())
    # compare_KW(df, NMR_cols, 'MDD', 'MDD_NMR')
    social_clust_stats(df_MDD, 'final_cluster')
    # bar_biochem_plot(df)
    # PHQ_plot(df)

    # cols = NMR_cols + [c + '-0.0' for c in ['1160', '4609', '4620', '5375', '5386', '20127', '2050', '2060', '2070', '2080']] + cate17518_cols + conti_cols + ['PHQ_Total', 'PHQ']
    # df01 = df[(df['final_cluster'] == 0) | (df['final_cluster'] == 1)]
    # compare_KW(df01, cols, 'final_cluster', '0_VS_1')
    # compare_Chi2(df01, '0_VS_1')
    #
    # df02 = df[(df['final_cluster'] == 0) | (df['final_cluster'] == 2)]
    # compare_KW(df02, cols, 'final_cluster', '0_VS_2')
    # compare_Chi2(df02, '0_VS_2')
    #
    # df03 = df[(df['final_cluster'] == 0) | (df['final_cluster'] == 3)]
    # compare_KW(df03, cols, 'final_cluster', '0_VS_3')
    # compare_Chi2(df03, '0_VS_3')
    #
    # df12 = df[(df['final_cluster'] == 1) | (df['final_cluster'] == 2)]
    # compare_KW(df12, cols, 'final_cluster', '1_VS_2')
    # compare_Chi2(df12, '1_VS_2')
    #
    # df13 = df[(df['final_cluster'] == 1) | (df['final_cluster'] == 3)]
    # compare_KW(df13, cols, 'final_cluster',  '1_VS_3')
    # compare_Chi2(df13, '1_VS_3')
    #
    # df23 = df[(df['final_cluster'] == 2) | (df['final_cluster'] == 3)]
    # compare_KW(df23, cols, 'final_cluster', , '2_VS_3')
    # compare_Chi2(df23, '2_VS_3')



