from SUtils import *


def heatmap_plot_bio_chem():
    df_MDD = pd.read_csv(data_path + "Diag_MDD.csv", low_memory=False)
    df_Y = pd.read_csv(result_path + "{}/{}_components_cluster_rst_{}.csv".format("NMF", 2, 3), low_memory=False,
                       usecols=['eid', 'final_cluster'])
    # 使用 replace 方法实现 cluster2 和 cluster3 的互换
    df_Y['final_cluster'] = df_Y['final_cluster'].replace({2: 3, 3: 2})
    df_MDD = pd.merge(df_MDD, df_Y, how='inner', on='eid', sort=True, suffixes=('_x', '_y'), copy=True, indicator=False)
    df_MDD.rename(columns=dict(zip(cate17518_cols, cate17518_new_cols)), inplace=True)
    # bio_chem_lst = [item for item in cate17518_new_cols if item not in ['RF', 'OEST']]
    df = df_MDD[['final_cluster']+cate17518_new_cols]
    print(df.shape[0])
    df.dropna(thresh=20, inplace=True)
    print(df.shape[0])
    null_ratios = df.isnull().sum() / df.shape[0]
    # 筛选出空值比例小于等于20%的列
    columns_to_keep = null_ratios[null_ratios <= 0.1].index
    keep_cols = columns_to_keep.tolist()
    print("{}, {}".format(len(keep_cols), keep_cols))
    del_col = [item for item in cate17518_new_cols if item not in keep_cols]
    print(del_col)

    df = df[keep_cols].sort_values(by='final_cluster')
    for col in keep_cols:
        df[col] = (df[col] - df[col].mean()) / df[col].std()

    # 设置索引并选择值
    heatmap_data = df.T
    plt.figure(figsize=(12, 6))  # 设置宽度为30，高度为6
    sns.heatmap(heatmap_data, annot=False, cmap='RdBu_r', center=0, vmin=-3, vmax=3)
    plt.xticks([], [])
    plt.title('Heatmap by Feature Order')
    plt.ylabel('Features')
    plt.savefig(img_path + "Heatmap_bio_chem.png", dpi=720)
    plt.show()
    plt.clf()


def heatmap_plot_NMR():
    df_MDD = pd.read_csv(data_path + "Diag_MDD.csv", low_memory=False)
    df_Y = pd.read_csv(result_path + "{}/{}_components_cluster_rst_{}.csv".format("NMF", 2, 3), low_memory=False,
                       usecols=['eid', 'final_cluster'])
    # 使用 replace 方法实现 cluster2 和 cluster3 的互换
    df_Y['final_cluster'] = df_Y['final_cluster'].replace({2: 3, 3: 2})
    df_MDD = pd.merge(df_MDD, df_Y, how='inner', on='eid', sort=True, suffixes=('_x', '_y'), copy=True, indicator=False)

    NMR_all_cols = []
    # for sheet_name, high in zip(['Total', 'HDL', 'LDL', 'VLDL', '%', 'OTHER'], [11, 6, 5, 10, 16, 10]):
    for sheet_name in ['Total', 'HDL', 'LDL', 'VLDL', '%', 'OTHER']:
        NMR_cols = pd.read_excel(path+'NMR_Heatmap_cols.xlsx', sheet_name=sheet_name, header=None).iloc[:, 0].tolist()
        NMR_all_cols = NMR_all_cols + NMR_cols

    print("{}, {}".format(len(NMR_all_cols), NMR_all_cols))
    df = df_MDD[['final_cluster']+NMR_all_cols]
    for col in df.columns.tolist():
        df[col] = (df[col] - df[col].mean()) / df[col].std()
    df = df.sort_values(by='final_cluster')
    # 设置索引并选择值
    heatmap_data = df.T

    plt.figure(figsize=(12, 58))  # 设置宽度、高度为
    sns.heatmap(heatmap_data, annot=False, cmap='RdBu_r', center=0, vmin=-3, vmax=3, cbar=False)
    plt.xticks([], [])
    plt.title('Heatmap by Feature Order')
    plt.savefig(img_path + "Heatmap_NMR.png", dpi=720)
    plt.show()
    plt.clf()

if __name__ == '__main__':
    # heatmap_plot_bio_chem()
    heatmap_plot_NMR()

