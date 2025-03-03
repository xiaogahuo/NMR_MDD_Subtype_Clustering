
import pandas as pd
from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import label_binarize
import matplotlib.pyplot as plt
import seaborn as sns

from SUtils import *
import numpy as np

def get_score_label(fold, exp_Name, model_name, d_type):
    # get the raw scores and labels from the test csvfile for the ROC PR curves
    df = pd.read_csv(result_path + "{}/{}/{}_{}_eval.csv".format(fold, exp_Name, model_name, d_type))
    return df[['pred_label']].values, df['label'].values

def draw_confuse(exp_Name, model_name):
    lw = 2
    for fold in range(5):
        types = ['train', 'test']
        fig, ax = plt.subplots(2, 1, figsize=(8, 15), gridspec_kw={'height_ratios': [1, 1]})
        for k in range(2):
          pred_label, labels = get_score_label(fold, exp_Name, model_name, types[k])
          # pred_label=np.argmax(np.array(scores), axis=1)
          matrix = confusion_matrix(np.array(labels), np.array(pred_label))
          sns.heatmap(matrix, annot=True, fmt='d', cmap='Oranges',
                      xticklabels=['HC',"MDD"],
                      yticklabels=['HC',"MDD"],ax=ax[k])
          ax[k].set_xlabel('Predicted')
          ax[k].set_ylabel('Actual')
          ax[k].set_title('{}-Confusion Matrix'.format(types[k]))
        if not os.path.exists(img_path + "{}/".format(exp_Name)):
            os.makedirs(img_path + "{}/".format(exp_Name))
        plt.savefig(img_path + "{}/{}_fold{}_confusion_martix.png".format(exp_Name, model_name, fold), bbox_inches='tight')
        plt.clf()
        plt.close()


if __name__ == "__main__":
    m_names = ["catboost", "catboost_TRL", "catboost_TRL_TG", "catboost_TRL_C",
               "catboost_RC", "catboost_TRL_P", "catboost_TRL_PL", "catboost_TRL_FC",
               "catboost_TRL_CE", "catboost_TRL_L", "catboost_VLDL", "catboost_HDL", "catboost_LDL","catboost_IDL"]

    for model_name in m_names:
        # draw_confuse("cluster1_vs_0", model_name)
        # draw_confuse("cluster2_vs_0", model_name)
        # draw_confuse("cluster3_vs_0", model_name)
        draw_confuse("cluster123_vs_0", model_name)
        draw_confuse("cluster1_vs_other", model_name)
        draw_confuse("cluster2_vs_other", model_name)
        draw_confuse("cluster3_vs_other", model_name)
        draw_confuse("NO_Chronic_cluster123_vs_0", model_name)
        draw_confuse("NO_Chronic_cluster1_vs_other", model_name)
        draw_confuse("NO_Chronic_cluster2_vs_other", model_name)
        draw_confuse("NO_Chronic_cluster3_vs_other", model_name)
