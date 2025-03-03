import gwaslab as gl
import pandas as pd

def plot():
    df=pd.read_csv("SAIGE_results.txt",sep=r'[ \t]',
                     engine="python",usecols=["CHR","POS","MarkerID","Allele1","Allele2","AF_Allele2",
                                              "BETA","SE","p.value"])
    df = df.rename(columns={
        'Allele1': 'nea',
        'Allele2': 'ea',
        'AF_Allele2': 'eaf',
        'MarkerID':'rsid'
    })

    df=gl.Sumstats(df,snpid="rsid",chrom="CHR",pos="POS",
             ea="ea",
             nea="nea",
             eaf="eaf",
             beta="BETA",
             se="SE",
             p="p.value")


    df.plot_mqq(skip=3,cut=20,anno=True,colors=["#ff0000","#fc4444",
                                           "#fc6404","#fcd444","#8cc43c",
                                           "#029658","#1abc9c","#5bc0de",
                                           "#6454ac","#fc8c84"],fontsize = 10,
              anno_fontsize = 10,title_fontsize = 13,marker_size=(5,25),
                  figargs={"figsize":(15,5),"dpi":300},save="./mqqplots.png")


if __name__=="__main__":
    plot()
