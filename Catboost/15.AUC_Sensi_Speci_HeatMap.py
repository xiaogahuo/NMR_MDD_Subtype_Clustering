from SUtils import *
import matplotlib.pyplot as plt
import seaborn as sns


def plot_heatmap(file_name, sheet_name):
    df = pd.read_excel(result_path + "{}_AUC_Heatmap.xlsx".format(file_name), sheet_name=sheet_name, index_col=0)  # Estimate
    df_c = df.copy()
    for col in ['Subtype1', 'Subtype2', 'Subtype3', 'MDD']:
        df_c[col] = df_c[col]-df_c['MDD']
    for col in ['Subtype1_Sensi', 'Subtype2_Sensi', 'Subtype3_Sensi', 'MDD_Sensi']:
        df_c[col] = df_c[col]-df_c['MDD_Sensi']
    for col in ['Subtype1_Speci', 'Subtype2_Speci', 'Subtype3_Speci', 'MDD_Speci']:
        df_c[col] = df_c[col]-df_c['MDD_Speci']
    print(df_c.head())
    print("max: {}, min: {}".format(df_c.max().max(), df_c.min().min()))
    plt.figure(figsize=(12, 6))
    myPlot = sns.heatmap(df_c, annot=df, fmt='.3f', annot_kws={"size": 9, "color": "black"}, yticklabels=True, cmap="vlag", center=0, vmin=-0.45, vmax=0.80, linewidth=.5)
    myPlot.xaxis.tick_top()
    plt.savefig(img_path + "{}_{}_Heatmap.png".format(file_name, sheet_name))
    plt.show()
    plt.clf()


if __name__ == '__main__':
    # for file in ['UKB', 'DPUK']:
    for file in ['UKB']:
        for sheet in ['main', 'Chornic', 'Female', 'Male', 'BMI', 'Metabolic', 'immune_mediated', 'drug', '3_excluded']:
            plot_heatmap(file, sheet)