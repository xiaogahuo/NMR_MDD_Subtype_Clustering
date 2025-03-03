import numpy as np
from SUtils import *
import seaborn as sns
import matplotlib.pyplot as plt

# algorithms = ['NMF', 'rbf_PCA', 'linear_PCA', 'poly_PCA', 'cosine_PCA']
algorithms = ['NMF']

for algorithm in algorithms:
    if not os.path.exists(img_path + "{}/".format(algorithm)):
        os.makedirs(img_path + "{}/".format(algorithm))
    for n_components in range(2, 9):
        for cluster_num in range(2, 9):
            # 按final_cluster对聚类结果排序，基于4个聚类算法的聚类结果生成加权similarity_matrix
            rst_file = result_path + "{}/{}_components_cluster_rst_{}.csv".format(algorithm, n_components, cluster_num)
            img_file = img_path + "{}/{}_components_cluster_rst_{}.png".format(algorithm, n_components, cluster_num)
            df = pd.read_csv(rst_file)
            df = df.sort_values(by='final_cluster')
            cluster_algorithms = df.columns.to_list()[1:5]
            n_algorithms = len(cluster_algorithms)
            m_samples = df.shape[0]
            print(cluster_algorithms, n_algorithms, m_samples)
            similarity_matrix = np.zeros((m_samples, m_samples))
            for col in cluster_algorithms:
                print(col)
                similarity_matrix += compute_similarity_matrix(df[col])
            similarity_matrix = similarity_matrix/n_algorithms

            sns.heatmap(similarity_matrix, annot=False, square=True)
            plt.savefig(img_file)
            plt.clf()
