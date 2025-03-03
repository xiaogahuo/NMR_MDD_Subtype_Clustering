import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn.manifold import TSNE
import umap
import pandas as pd
from SUtils import *

# algorithms = ['NMF', 'rbf_PCA', 'linear_PCA', 'poly_PCA', 'cosine_PCA']
algorithms = ['NMF']

dataset = "Diag_MDD"
df_MDD = pd.read_csv(data_path + "{}.csv".format(dataset))
print(df_MDD.head(20))
X = df_MDD[NMR_cols]
for algorithm in algorithms:
    if not os.path.exists(img_path + "{}/".format(algorithm)):
        os.makedirs(img_path + "{}/".format(algorithm))

    for n_components in range(2, 3):
        for cluster_num in range(3, 4):
            # 加载数据集
            rst_file = result_path + "{}/{}_components_cluster_rst_{}.csv".format(algorithm, n_components, cluster_num)
            img_file = img_path + "{}/{}_components_cluster_tsne_umap_{}.png".format(algorithm, n_components,
                                                                                     cluster_num)
            df_Y = pd.read_csv(rst_file)
            df_Y['final_cluster'] = df_Y['final_cluster'].replace({2: 3, 3: 2})
            y = df_Y['final_cluster']

            # UMAP降维
            umap_reducer = umap.UMAP(random_state=42)
            X_umap = umap_reducer.fit_transform(X)

            # t-SNE降维
            tsne_reducer = TSNE(random_state=42)
            X_tsne = tsne_reducer.fit_transform(X)

            colors = ["#ABDDA4", "#FDAE61", "#D7191C"]
            cmap = ListedColormap(colors)

            # 绘制UMAP结果
            plt.figure(figsize=(12, 6))
            plt.subplot(121)
            plt.scatter(X_umap[:, 0], X_umap[:, 1], c=y, cmap=cmap, s=4)
            # plt.colorbar()
            plt.title('UMAP')

            # 绘制t-SNE结果
            plt.subplot(122)
            plt.scatter(X_tsne[:, 0], X_tsne[:, 1], c=y, cmap=cmap, s=4)
            # plt.colorbar()
            plt.title('t-SNE')

            plt.savefig(img_file)
            plt.clf()
