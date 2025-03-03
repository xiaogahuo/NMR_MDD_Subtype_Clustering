from SUtils import *
from scipy.stats import spearmanr

def corr(df, clu_list):
    df_corr = df[df['final_cluster'].isin(clu_list)]

    rst_df_corr = pd.DataFrame(columns=NMR_shap_Cols)
    rst_df_pval = pd.DataFrame(columns=NMR_shap_Cols)

    for prs_c in PRS_Cols:
        lst_corr, lst_pval = [], []
        for nmr_c in NMR_shap_Cols:
            corr, pval = spearmanr(df_corr[prs_c], df_corr[nmr_c])
            lst_corr.append(corr)
            lst_pval.append(pval)
        rst_df_corr = pd.concat([rst_df_corr, pd.DataFrame([lst_corr], columns=rst_df_corr.columns)], ignore_index=True)
        rst_df_pval = pd.concat([rst_df_pval, pd.DataFrame([lst_pval], columns=rst_df_pval.columns)], ignore_index=True)

    rst_df_corr.index, rst_df_pval.index = PRS_Cols, PRS_Cols
    rst_df_corr.to_csv(result_path + "PRS_NMR_Spearmanr_Corr_clust_{}.csv".format(','.join(map(str, clu_list))), ",")
    rst_df_pval.to_csv(result_path + "PRS_NMR_Spearmanr_Pval_clust_{}.csv".format(','.join(map(str, clu_list))), ",")


if __name__ == '__main__':
    NMR_shap_Cols = pd.read_excel(path+'TOP20.xlsx', sheet_name=0, header=None).iloc[:, 0].tolist()
    print(NMR_shap_Cols)
    PRS_Cols = ['SCORE_5e08', 'SCORE_1e07', 'SCORE_5e07', 'SCORE_1e06', 'SCORE_5e06', 'SCORE_1e05', 'SCORE_5e05', 'SCORE_1e04', 'SCORE_5e04', 'SCORE_1e03', 'SCORE_5e03', 'SCORE_1e02', 'SCORE_5e02']

    df = pd.read_csv(data_path + "PRS_SCORE_Dataset.csv", low_memory=False, usecols=['eid', 'final_cluster']+NMR_shap_Cols+PRS_Cols)
    print(df.columns.tolist())

    for clt in [[0], [1], [2], [3], [1, 2, 3]]:
        corr(df, clt)

