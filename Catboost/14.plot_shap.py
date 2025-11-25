import shap
import matplotlib.pyplot as plt
from SUtils import *


def get_topN_reason(exp_Name, feature_importance, cols, model_name, top_num=20):
    #输出feature_importance值最高的N个标签
    feature_importance_dict = {}
    for i, f in zip(feature_importance, cols):
        feature_importance_dict[f] = i
    new_dict = dict(sorted(feature_importance_dict.items(), key=lambda e: e[1], reverse=True))
    f = open(img_path+"{}/{}_SHAP.csv".format(exp_Name, model_name), "w")
    print("key,value", file=f)
    for k, v in new_dict.items():
        if top_num > 0:
            print("{},{}".format(k, v), file=f)
            top_num -= 1
        else:
            break


def plot_shap(exp_Name, model_name, cols):
    # 绘制SHAP结果
    shap_values, X_test = [], []
    for fold in range(5):
        test = pd.read_csv(data_path + "cross_valid/{}/{}/X_test.csv".format(fold, exp_Name), usecols=cols)
        X_test.append(test)
        shap_result = np.load(result_path + "{}/{}/{}_test_shap.npy".format(fold, exp_Name, model_name))
        shap_values.append(shap_result)
    X_test = pd.concat(X_test, ignore_index=True, axis=0)
    shap_values = np.concatenate(shap_values, axis=0)
    print("exp: {}: X_test shape: {}, shap_values: {}".format(exp_Name, X_test.shape, shap_values.shape))

    plt.figure(figsize=(20, 8), dpi=720)
    plt.subplot(1, 2, 1)  # 1 row  2 col , the first fig
    # summarize the effects of all the features
    # shap.summary_plot(shap_values, X_test, show=False, max_display=20)
    shap.plots.violin(shap_values, features=X_test, feature_names=cols, plot_type="layered_violin", show=False, max_display=20)
    plt.xlabel('SHAP value')
    plt.subplot(1, 2, 2)  # 1 row  2 col , the second fig
    shap.summary_plot(shap_values, features=X_test, feature_names=cols, plot_type="bar", show=False, max_display=20)
    plt.xlabel('Mean(|SHAP|)')
    if not os.path.exists(img_path+"{}/".format(exp_Name)):
        os.makedirs(img_path+"{}/".format(exp_Name))
    plt.savefig(img_path+"{}/{}_SHAP.png".format(exp_Name, model_name), dpi=720)
    plt.clf()  # CLEAR FIGURE

    # Feature 重要性
    feature_importance = np.absolute(shap_values).mean(axis=0)
    get_topN_reason(exp_Name, feature_importance, cols, model_name, top_num=len(cols))


if __name__ == '__main__':
    m_names = ["catboost", "catboost_TRL", "catboost_TRL_TG", "catboost_TRL_C", "catboost_RC",
               "catboost_TRL_P", "catboost_TRL_PL", "catboost_TRL_FC", "catboost_TRL_CE",
               "catboost_TRL_L", "catboost_VLDL", "catboost_HDL", "catboost_LDL", "catboost_IDL",
               "catboost_C", "catboost_TG", "catboost_PL", "catboost_CE", "catboost_FC", "catboost_L",
               "catboost_P", "catboost_Size", "catboost_Total", "catboost_FA"]
    cols_l = [NMR_cols, TRL, TRL_TG, TRL_C, RC, TRL_P, TRL_PL, TRL_FC, TRL_CE, TRL_L, VLDL, HDL, LDL, IDL,
              C, TG, PL, CE, FC, L, P, Size, Total, FA]

    # m_names = ["catboost_FA"]
    # cols_l = [FA]

    for model_name, cols in zip(m_names, cols_l):
        print(model_name, "=" * 50)
        # plot_shap("cluster123_vs_0", model_name, cols)
        # plot_shap("cluster1_vs_other", model_name, cols)
        # plot_shap("cluster2_vs_other", model_name, cols)
        # plot_shap("cluster3_vs_other", model_name, cols)
        # plot_shap("NO_Chronic_cluster123_vs_0", model_name, cols)
        # plot_shap("NO_Chronic_cluster1_vs_other", model_name, cols)
        # plot_shap("NO_Chronic_cluster2_vs_other", model_name, cols)
        # plot_shap("NO_Chronic_cluster3_vs_other", model_name, cols)
        #
        # plot_shap("Male_cluster123_vs_0", model_name, cols)
        # plot_shap("Male_cluster1_vs_other", model_name, cols)
        # plot_shap("Male_cluster2_vs_other", model_name, cols)
        # plot_shap("Male_cluster3_vs_other", model_name, cols)
        # plot_shap("Female_cluster123_vs_0", model_name, cols)
        # plot_shap("Female_cluster1_vs_other", model_name, cols)
        # plot_shap("Female_cluster2_vs_other", model_name, cols)
        # plot_shap("Female_cluster3_vs_other", model_name, cols)
        # plot_shap("BMI_cluster123_vs_0", model_name, cols)
        # plot_shap("BMI_cluster1_vs_other", model_name, cols)
        # plot_shap("BMI_cluster2_vs_other", model_name, cols)
        # plot_shap("BMI_cluster3_vs_other", model_name, cols)

        plot_shap("immune_mediated_cluster123_vs_0", model_name, cols)
        plot_shap("immune_mediated_cluster1_vs_other", model_name, cols)
        plot_shap("immune_mediated_cluster2_vs_other", model_name, cols)
        plot_shap("immune_mediated_cluster3_vs_other", model_name, cols)
        plot_shap("Metabolic_cluster123_vs_0", model_name, cols)
        plot_shap("Metabolic_cluster1_vs_other", model_name, cols)
        plot_shap("Metabolic_cluster2_vs_other", model_name, cols)
        plot_shap("Metabolic_cluster3_vs_other", model_name, cols)
        plot_shap("3excluded_cluster123_vs_0", model_name, cols)
        plot_shap("3excluded_cluster1_vs_other", model_name, cols)
        plot_shap("3excluded_cluster2_vs_other", model_name, cols)
        plot_shap("3excluded_cluster3_vs_other", model_name, cols)

