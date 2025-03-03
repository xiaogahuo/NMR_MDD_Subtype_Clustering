from collections import Counter

import pandas as pd
from plinkio import plinkfile
from SUtils import *

exp_name="GWAS_NMR_MDD_Subtype_Clustering_bgen"
result_name="final_result"

#生成每个染色体过滤SNP的txt文件
#174
def merge_txt(f_s):
    df = []
    for file_name in f_s:
        df_ch=pd.read_csv(bgen_path + "{}/{}/SAIGE_results.sig.csv".format(exp_name, file_name),usecols=["CHR","MarkerID","p.value"])  # sep='\t'
        df_ch=df_ch[df_ch["p.value"]<=5e-8]
        df.append(df_ch)
    df=pd.concat(df)
    df = df.drop_duplicates(subset=['CHR',"MarkerID"], keep='last').reset_index(drop=True)
    df.sort_values(by="p.value")
    df.to_csv(bgen_path+result_name+"/filter_rs_txt/all.csv",index=False)
    print("Total RS num:{}".format(len(df)))
    # 根据 CHR 列分组
    df_group = df.groupby('CHR')

    # 遍历每个分组
    for chr_num, group in df_group:
        # 获取该分组的 MarkerID 列，保存为 .txt 文件
        filename = bgen_path+result_name+"/filter_rs_txt/"+f'CHR_{chr_num}.txt'
        group['MarkerID'].to_csv(filename, index=False, header=False)
        print(f"数据已保存到 {filename}")

#读取基因型
def merge_genotypes_rs(chrs):
    # df_marker=pd.read_csv(bgen_path+result_name+"/filter_rs_txt/all.csv",usecols=["MarkerID"])
    # marker_ID=list(df_marker["MarkerID"])
    marker_ID=["rs964184"]
    eids=[]
    df=pd.DataFrame({})
    for chr in chrs:
        df_tmp=pd.DataFrame({})
        plink_file = plinkfile.open(bgen_path+result_name+"/filter_rs_bed/chr_{}".format(chr))
        if not plink_file.one_locus_per_row():
            print("This script requires that snps are rows and samples columns.")
            exit(1)

        #采样信息
        sample_list = plink_file.get_samples()
        locus_list = plink_file.get_loci()

        df_tmp["eid"]=[sample.iid for index,sample in zip(range(len(sample_list)),sample_list)]
        df_tmp["eid"]=df_tmp["eid"].astype(int)

        for locus, row in zip(locus_list, plink_file):
            #初始化基因每一列
            # 创建要添加的列数据
            id=str(locus.name)
            if id in marker_ID:
               column_data = [genotype for index, genotype in zip(range(len(row)), row)]
            # 使用pd.concat一次性添加所有列
               df_tmp = pd.concat([df_tmp] + [pd.DataFrame({"{}".format(id): column_data})], axis=1)
        if(len(df)==0):
            df=df_tmp
        else:
            df = reduce(lambda left, right: pd.merge(left, right, on=['eid'], how='inner'),[df, df_tmp])

            # print("Individual {0} has genotype {1} for snp {2}.".format(sample.iid, genotype, locus.name))
    df_social_MDD=pd.read_csv(path + "Social_MDD_format_new.csv", low_memory=False,usecols=['eid','MDD','Sex_cate','Age','survival_time'])
    print(len(df_social_MDD))
    #487409
    print("merge before sample num:{}".format(len(df)))
    df=reduce(lambda left, right: pd.merge(left, right, on=['eid'], how='inner'), [df, df_social_MDD])
    # 228154
    print("merge after sample num:{}".format(len(df)))
    #筛除正常死亡或失访
    df=df[df["survival_time"]!=-1]
    #217514

    print("filter after sample num:{}".format(len(df)))
    df.to_csv(bgen_path+result_name+"/filter_rs_csv/"+"sum_rs964184.csv")
#

if __name__=="__main__":
    # f_s = ["1_VS_3", "1_VS_2", "3_VS_2", "0_VS_123", "0_VS_1", "0_VS_2", "0_VS_3"]
    # merge_txt(f_s)
    chrs=[2,8,11,18]
    merge_genotypes_rs(chrs)

