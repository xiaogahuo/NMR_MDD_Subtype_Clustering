from SUtils import *


# 用于生成示例文件以便于查阅GWAS summary 数据确认A1 A2 beta 或者 or
def filter_gwas_summary_1K(file_name):
    print(file_name)
    df = pd.read_csv("/data/UKBioBank/LDScore_EUR/Downloads/{}".format(file_name), sep=r'\s+', engine='python')
    df_1000 = df.head(1000)
    df_1000.to_csv("/data/UKBioBank/LDScore_EUR/Format_GWAS_Summary_Data_1K/{}.txt".format(file_name), sep='\t', index=False)
    print("{} Done".format(file_name))


def format_gwas_summary(df, file_name, origin_cols, new_cols):
    print(file_name)
    # df_0 = pd.read_csv("/data/UKBioBank/LDScore_EUR/{}".format("w_hm3.snplist"), sep=r'\s+', engine='python')
    # df_f= pd.merge(df, df_0, left_on='variant_id', right_on='SNP')
    # print(df_f.head())
    df.rename(columns=dict(zip(origin_cols, new_cols)), inplace=True)
    # df[new_cols] = df[new_cols].applymap(lambda x: x.replace('"', '') if isinstance(x, str) else x)
    df[new_cols].to_csv("/data/UKBioBank/LDScore_EUR/Format_GWAS_Summary_Data/{}.txt".format(file_name), sep='\t', index=False)
    print(df.head())
    print("{} Done".format(file_name))


if __name__ == '__main__':
    # data = ['GCST90029070_buildGRCh37.tsv.gz', 'facialpain2_f6159_v3.bgenie.txt.gz',
    #             'neckshoulderpain2_f6159_v3.bgeniemodify.txt.gz', 'GCST90225528_buildGRCh37.tsv', 'IL5.data.gz',
    #             'SCGFb.data.gz', 'CTACK.data.gz', 'FGFBasic.data.gz', 'bNGF.data.gz','SDF1a.data.gz', 'HGF.data.gz', 'MCSF.data.gz',
    #         'GCST90042788_buildGRCh37.tsv.gz', 'IL10.data.gz', 'backpain2_f6159.2018-02-09.bgenie.gz',
    #         'headache2_f6159_WhiteBritishPCs_addPrunedRels_noPGC_noGenScot_v2.2017-08-30.bgenie.gz',
    #         'GCST90104898_buildGRCh37.tsv', 'PDGFbb.data.gz', 'abdominalpain2_f6159_v3.bgenie.txt.gz', 'VEGF.data.gz',
    #         'IL17.data.gz', 'GCST90042782_buildGRCh37.tsv.gz', 'MCP1.data.gz', 'painalloverbody_f6159_v3.bgenie.txt.gz',
    #         'MIP1a.data.gz', 'TNFb.data.gz', 'IL2.data.gz', 'IFNg.data.gz', 'IL13.data.gz', 'MIP1b.data.gz',
    #         'hippain2_f6159_v3.bgenie.txt.gz', 'Eotaxin.data.gz', 'IL2ra.data.gz', 'SCF.data.gz', 'MIG.data.gz',
    #         'TNFa.data.gz', 'IL7.data.gz', 'MIF.data.gz', 'IL4.data.gz', 'RANTES.data.gz', 'MCP3.data.gz',
    #         'IL8.data.gz', 'GROa.data.gz', 'kneepain2_f6159_v3_1812.bgenie.txt.gz', 'IL12p70.data.gz', 'IP10.data.gz',
    #         'IL9.data.gz', 'GCSF.data.gz', 'TRAIL.data.gz', 'IL1ra.data.gz', 'IL6.data.gz', 'IL1b.data.gz',
    #         'IL18.data.gz', 'IL16.data.gz']
    # print(len(data))
    # for file_name in data:
    #     filter_gwas_summary_1K(file_name)

    # cols = ['snpid', 'hg18chr', 'bp', 'a1', 'a2', 'beta', 'pval']
    #
    # file_name = 'headache2_f6159_WhiteBritishPCs_addPrunedRels_noPGC_noGenScot_v2.2017-08-30.bgenie.gz'
    # df = pd.read_csv("/data/UKBioBank/LDScore_EUR/Downloads/{}".format(file_name), sep=r'\s+', engine='python')
    # df['pval'] = 10 ** -df['headache2_f6159-log10p']
    # format_gwas_summary(df, file_name,['rsid', 'chr', 'pos', 'a_1', 'a_0', 'headache2_f6159_beta', 'pval'], cols) # beta or or
    #
    # for file_name in ['facialpain2_f6159_v3.bgenie.txt.gz', 'neckshoulderpain2_f6159_v3.bgeniemodify.txt.gz',
    #                   'backpain2_f6159.2018-02-09.bgenie.gz', 'abdominalpain2_f6159_v3.bgenie.txt.gz',
    #                   'painalloverbody_f6159_v3.bgenie.txt.gz', 'hippain2_f6159_v3.bgenie.txt.gz',
    #                   'kneepain2_f6159_v3_1812.bgenie.txt.gz']:
    #     df = pd.read_csv("/data/UKBioBank/LDScore_EUR/Downloads/{}".format(file_name), sep=r'\s+', engine='python')
    #     format_gwas_summary(df, file_name, ['rsid', 'chr', 'pos', 'a_1', 'a_0', 'beta', 'p'], cols)  # beta or or
    #
    # for file_name in ['GCST90029070_buildGRCh37.tsv.gz', 'GCST90225528_buildGRCh37.tsv', 'GCST90042788_buildGRCh37.tsv.gz',
    #                   'GCST90042782_buildGRCh37.tsv.gz', 'GCST90104898_buildGRCh37.tsv']:
    #     df = pd.read_csv("/data/UKBioBank/LDScore_EUR/Downloads/{}".format(file_name), sep=r'\s+', engine='python')
    #     format_gwas_summary(df, file_name, ['variant_id', 'chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'beta', 'p_value'], cols)  # beta or or

    #    if "data" in dat:
    #            format_gwas_summary(dat,
    #      ['MarkerName', 'Chromosome', 'Position', 'OtherAllele', 'EffectAllele', 'Effect', 'P.value'],
    #      ['snpid', 'hg18chr', 'bp', 'a1', 'a2', 'beta', 'pval']) # beta or or


    # df = pd.read_csv("/data/UKBioBank/LDScore_EUR/Downloads/sumstats_neuroticism_ctg_format.txt",
    #     sep=r'\s+', engine='python', usecols=['RSID', 'A1', 'A2', 'N', 'Z'])
    # df.rename(columns={'RSID': 'SNP'}, inplace=True)
    # print(df.head())
    # df.to_csv("/data/UKBioBank/LDScore_EUR/sumstats_ALL/neuroticism.sumstats.gz", sep='\t', index=False, compression='gzip')

    pass