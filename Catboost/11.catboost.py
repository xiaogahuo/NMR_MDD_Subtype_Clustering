from catboost import CatBoostClassifier
from SUtils import *
import csv
import shap

def run_cat(exp_Name, model_name, cols, n_estimators, class_weights):
    print("{}:{}".format(exp_Name, model_name))
    for fold in range(5):
        X_train = pd.read_csv(data_path + "cross_valid/{}/{}/X_train.csv".format(fold, exp_Name))[cols]
        X_test = pd.read_csv(data_path + "cross_valid/{}/{}/X_test.csv".format(fold, exp_Name))[cols]
        y_train = pd.read_csv(data_path + "cross_valid/{}/{}/y_train.csv".format(fold, exp_Name))['final_cluster']
        y_test = pd.read_csv(data_path + "cross_valid/{}/{}/y_test.csv".format(fold, exp_Name))['final_cluster']

        if not os.path.exists(result_path + "{}/{}/".format(fold, exp_Name)):
            os.makedirs(result_path + "{}/{}/".format(fold, exp_Name))
        # model = CatBoostClassifier(subsample=0.55, learning_rate=0.01, max_depth=5, n_estimators=param['n_estimators'], cat_features=param['cat_features'])
        # n_estimators = int(test.shape[0]*300/10000) # 按数据集大小等比例设置 n_estimators
        model = CatBoostClassifier(subsample=0.50, learning_rate=0.02, max_depth=5, eval_metric='F1', loss_function='Logloss', cat_features=None, n_estimators=n_estimators, class_weights=class_weights)
        model.fit(X_train, y_train, eval_set=[(X_test, y_test)], use_best_model=True, verbose=True)

        # Extract the loss values for training and validation datasets
        results = model.get_evals_result()
        train_loss = results['learn']['Logloss']
        test_loss = results['validation']['Logloss']
        pd.DataFrame({'loss': train_loss}).to_csv(result_path + "{}/{}/{}_train_loss.csv".format(fold, exp_Name, model_name), index=None)
        pd.DataFrame({'loss': test_loss}).to_csv(result_path + "{}/{}/{}_test_loss.csv".format(fold, exp_Name, model_name), index=None)

        predict_test_prob = model.predict_proba(X_test)  # 返回预测概率
        df_test = pd.DataFrame()
        df_test['label'] = y_test
        df_test['pred_prob'] = predict_test_prob[:, 1]
        df_test['pred_label'] = model.predict(X_test)
        df_test.to_csv(result_path + "{}/{}/{}_test_eval.csv".format(fold, exp_Name, model_name), index=None)

        predict_train_prob = model.predict_proba(X_train)  # 返回预测概率
        df_train = pd.DataFrame()
        df_train['label'] = y_train
        df_train['pred_prob'] = predict_train_prob[:, 1]
        df_train['pred_label'] = model.predict(X_train)
        df_train.to_csv(result_path + "{}/{}/{}_train_eval.csv".format(fold, exp_Name, model_name), index=None)

        # save SHAP_values
        explainer = shap.TreeExplainer(model)  # #这里的model在准备工作中已经完成建模，模型名称就是model
        shap_values = explainer.shap_values(X_test)  # 传入特征矩阵X，计算SHAP值
        print("shap_values shape: {}".format(shap_values.shape))
        np.save(result_path + "{}/{}/{}_test_shap.npy".format(fold, exp_Name, model_name), shap_values)


if __name__ == '__main__':
    m_names = ["catboost", "catboost_TRL", "catboost_TRL_TG", "catboost_TRL_C",
               "catboost_RC", "catboost_TRL_P", "catboost_TRL_PL", "catboost_TRL_FC",
               "catboost_TRL_CE", "catboost_TRL_L", "catboost_VLDL", "catboost_HDL", "catboost_LDL","catboost_IDL",
               "catboost_C", "catboost_TG", "catboost_PL", "catboost_CE", "catboost_FC", "catboost_L",
               "catboost_P", "catboost_Size", "catboost_Total", "catboost_FA"]

    cols_l = [NMR_cols, TRL, TRL_TG, TRL_C, RC, TRL_P, TRL_PL, TRL_FC, TRL_CE, TRL_L, VLDL, HDL, LDL, IDL,
              C, TG, PL, CE, FC, L, P, Size, Total, FA]

    for model_name, cols in zip(m_names, cols_l):
        print(model_name, "=" * 50)
        # run_cat("immune_mediated_cluster123_vs_0", model_name, cols, 700, {NC:1, MDD:20})
        # run_cat("immune_mediated_cluster1_vs_other", model_name, cols, 500, {NC:1, MDD:40})
        # run_cat("immune_mediated_cluster2_vs_other", model_name, cols, 500, {NC:1, MDD:40})
        # run_cat("immune_mediated_cluster3_vs_other", model_name, cols, 500, {NC:1, MDD:40})
        #
        # run_cat("Metabolic_cluster123_vs_0", model_name, cols, 700, {NC: 1, MDD: 20})
        # run_cat("Metabolic_cluster1_vs_other", model_name, cols, 500, {NC: 1, MDD: 40})
        # run_cat("Metabolic_cluster2_vs_other", model_name, cols, 500, {NC: 1, MDD: 40})
        # run_cat("Metabolic_cluster3_vs_other", model_name, cols, 500, {NC: 1, MDD: 40})

        run_cat("Drug_cluster123_vs_0", model_name, cols, 700, {NC: 1, MDD: 20})
        run_cat("Drug_cluster1_vs_other", model_name, cols, 500, {NC: 1, MDD: 50})
        run_cat("Drug_cluster2_vs_other", model_name, cols, 500, {NC: 1, MDD: 120})
        run_cat("Drug_cluster3_vs_other", model_name, cols, 500, {NC: 1, MDD: 60})
        #
        # run_cat("3excluded_cluster123_vs_0", model_name, cols, 700, {NC: 1, MDD: 20})
        # run_cat("3excluded_cluster1_vs_other", model_name, cols, 500, {NC: 1, MDD: 40})
        # run_cat("3excluded_cluster2_vs_other", model_name, cols, 500, {NC: 1, MDD: 200})
        # run_cat("3excluded_cluster3_vs_other", model_name, cols, 500, {NC: 1, MDD: 70})
        #
        # run_cat("cluster123_vs_0", model_name, cols, 700, {NC:1, MDD:20})
        # run_cat("cluster1_vs_other", model_name, cols, 500, {NC:1, MDD:40})
        # run_cat("cluster2_vs_other", model_name, cols, 500, {NC:1, MDD:40})
        # run_cat("cluster3_vs_other", model_name, cols, 500, {NC:1, MDD:40})
        #
        # run_cat("NO_Chronic_cluster123_vs_0", model_name, cols, 700, {NC: 1, MDD: 20})
        # run_cat("NO_Chronic_cluster1_vs_other", model_name, cols, 500, {NC: 1, MDD: 40})
        # run_cat("NO_Chronic_cluster2_vs_other", model_name, cols, 500, {NC: 1, MDD: 40})
        # run_cat("NO_Chronic_cluster3_vs_other", model_name, cols, 500, {NC: 1, MDD: 40})
        #
        # run_cat("Male_cluster123_vs_0", model_name, cols, 700, {NC: 1, MDD: 20})
        # run_cat("Male_cluster1_vs_other", model_name, cols, 500, {NC: 1, MDD: 40})
        # run_cat("Male_cluster2_vs_other", model_name, cols, 500, {NC: 1, MDD: 40})
        # run_cat("Male_cluster3_vs_other", model_name, cols, 500, {NC: 1, MDD: 40})
        #
        # run_cat("Female_cluster123_vs_0", model_name, cols, 700, {NC: 1, MDD: 20})
        # run_cat("Female_cluster1_vs_other", model_name, cols, 500, {NC: 1, MDD: 40})
        # run_cat("Female_cluster2_vs_other", model_name, cols, 500, {NC: 1, MDD: 40})
        # run_cat("Female_cluster3_vs_other", model_name, cols, 500, {NC: 1, MDD: 40})
        #
        # run_cat("BMI_cluster123_vs_0", model_name, cols, 700, {NC: 1, MDD: 20})
        # run_cat("BMI_cluster1_vs_other", model_name, cols, 500, {NC: 1, MDD: 40})
        # run_cat("BMI_cluster2_vs_other", model_name, cols, 500, {NC: 1, MDD: 40})
        # run_cat("BMI_cluster3_vs_other", model_name, cols, 500, {NC: 1, MDD: 40})

