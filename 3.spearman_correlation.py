from SUtils import *
from scipy.stats import spearmanr


def transform_value(x):
    return 0 if -np.log10(x) < -np.log10(0.05) else -np.log10(x)


def corr(df, f_name):
    # 计算Spearman相关性
    spearman_corr, p_values = spearmanr(df, axis=0)
    print(df.columns.tolist())
    # 将相关性和p值转换为DataFrame
    corr_df = pd.DataFrame(spearman_corr, columns=df.columns, index=df.columns)
    p_values_df = pd.DataFrame(p_values, columns=df.columns, index=df.columns)
    corr_df.to_csv(img_path + f_name + "_NMR_spearman_corr.csv")
    p_values_df.to_csv(img_path + f_name + "_NMR_spearman_p_values.csv")

    # 创建一个掩码，保留左下半角
    mask = np.triu(np.ones_like(corr_df, dtype=bool))
    plt.figure(figsize=(60, 60))
    sns.heatmap(corr_df, mask=mask, cmap='RdBu_r', annot=False, square=True)
    plt.title('Spearman Correlation Heatmap (Lower Triangle)')
    plt.savefig(img_path + f_name + "_NMR_spearman_corr.png")
    plt.show()
    plt.clf()

    # 创建一个掩码，保留右上半角
    p_values_df = p_values_df.map(transform_value)
    mask = np.tril(np.ones_like(p_values_df, dtype=bool))
    plt.figure(figsize=(60, 60))
    sns.heatmap(p_values_df, mask=mask, cmap='Reds', vmax=20, annot=False, square=True, xticklabels=False, yticklabels=False)
    plt.title('Spearman P_value Heatmap (Upper Triangle)')
    plt.savefig(img_path + f_name + "_NMR_spearman_p_values.png")
    plt.show()
    plt.clf()


if __name__ == '__main__':
    NMR_order_cols = pd.read_excel(path+'热图NMR顺序.xlsx', sheet_name=0, header=None).iloc[:, 0]
    df_MDD = pd.read_csv(data_path + "Diag_MDD.csv", usecols=NMR_order_cols, low_memory=False)
    df_MDD = df_MDD.reindex(columns=NMR_order_cols)
    corr(df_MDD, "MDD")
    df_NC = pd.read_csv(data_path + "Diag_NC.csv", usecols=NMR_order_cols, low_memory=False)
    df_NC = df_NC.reindex(columns=NMR_order_cols)
    corr(df_NC, "NC")
    df = pd.concat([df_MDD, df_NC])
    corr(df, "MDD&NC")

