# 排除 《20002self自我报告疾病排除.xlsx》 中指定疾病的受试者
# 按 '1286' 疾病代码定义自我报告 Depression    对应MergerAllData.py中 Dataset 4
from SUtils import *


def read_exclude_code():
    df = pd.DataFrame(pd.read_excel(root_path+'20002self自我报告疾病排除.xlsx', sheet_name=0), columns=['coding'])
    list = df['coding'].to_list()
    list = [str(elem) for elem in list]
    print('len: {}, list: {}'.format(len(list), list))
    return list


cols_names = ['20002-0.0', '20002-0.1', '20002-0.2', '20002-0.3', '20002-0.4', '20002-0.5', '20002-0.6', '20002-0.7', '20002-0.8', '20002-0.9', '20002-0.10', '20002-0.11', '20002-0.12', '20002-0.13', '20002-0.14', '20002-0.15', '20002-0.16', '20002-0.17', '20002-0.18', '20002-0.19', '20002-0.20', '20002-0.21', '20002-0.22', '20002-0.23', '20002-0.24', '20002-0.25', '20002-0.26', '20002-0.27', '20002-0.28', '20002-0.29', '20002-0.30', '20002-0.31', '20002-0.32', '20002-0.33', '20002-1.0', '20002-1.1', '20002-1.2', '20002-1.3', '20002-1.4', '20002-1.5', '20002-1.6', '20002-1.7', '20002-1.8', '20002-1.9', '20002-1.10', '20002-1.11', '20002-1.12', '20002-1.13', '20002-1.14', '20002-1.15', '20002-1.16', '20002-1.17', '20002-1.18', '20002-1.19', '20002-1.20', '20002-1.21', '20002-1.22', '20002-1.23', '20002-1.24', '20002-1.25', '20002-1.26', '20002-1.27', '20002-1.28', '20002-1.29', '20002-1.30', '20002-1.31', '20002-1.32', '20002-1.33', '20002-2.0', '20002-2.1', '20002-2.2', '20002-2.3', '20002-2.4', '20002-2.5', '20002-2.6', '20002-2.7', '20002-2.8', '20002-2.9', '20002-2.10', '20002-2.11', '20002-2.12', '20002-2.13', '20002-2.14', '20002-2.15', '20002-2.16', '20002-2.17', '20002-2.18', '20002-2.19', '20002-2.20', '20002-2.21', '20002-2.22', '20002-2.23', '20002-2.24', '20002-2.25', '20002-2.26', '20002-2.27', '20002-2.28', '20002-2.29', '20002-2.30', '20002-2.31', '20002-2.32', '20002-2.33', '20002-3.0', '20002-3.1', '20002-3.2', '20002-3.3', '20002-3.4', '20002-3.5', '20002-3.6', '20002-3.7', '20002-3.8', '20002-3.9', '20002-3.10', '20002-3.11', '20002-3.12', '20002-3.13', '20002-3.14', '20002-3.15', '20002-3.16', '20002-3.17', '20002-3.18', '20002-3.19', '20002-3.20', '20002-3.21', '20002-3.22', '20002-3.23', '20002-3.24', '20002-3.25', '20002-3.26', '20002-3.27', '20002-3.28', '20002-3.29', '20002-3.30', '20002-3.31', '20002-3.32', '20002-3.33']
def single_row(x):
    result_list = [str(int(elem)) if not np.isnan(elem) else None for elem in x[cols_names]]
    result_list = list(filter(None, result_list))  # 去除列表空值
    t_1243 = 1 if '1243' in result_list else 0
    t_1286 = 1 if '1286' in result_list else 0
    t_1287 = 1 if '1287' in result_list else 0
    t_1288 = 1 if '1288' in result_list else 0
    t_1289 = 1 if '1289' in result_list else 0
    t_1290 = 1 if '1290' in result_list else 0
    t_1291 = 1 if '1291' in result_list else 0
    t_1461 = 1 if '1461' in result_list else 0
    t_1462 = 1 if '1462' in result_list else 0
    t_1463 = 1 if '1463' in result_list else 0
    t_1482 = 1 if '1482' in result_list else 0
    exclude_diagnosis = list(set(exclude_code) & set(result_list))  # 求交集
    return t_1243, t_1286, t_1287, t_1288, t_1289, t_1290, t_1291, t_1461, t_1462, t_1463, t_1482, exclude_diagnosis, 1 if len(exclude_diagnosis) > 0 else 0


def self_diagnosis(exclude_code):
    # df = pd.read_csv(root_path + "Cate_3001.csv")
    df = pd.read_csv(root_path + "Cate_3001.csv", usecols=['eid']+cols_names)
    # df_columns = df.columns.tolist()
    # print(df.columns.tolist())
    row = df.apply(lambda x: single_row(x), axis=1)
    t_1243, t_1286, t_1287, t_1288, t_1289, t_1290, t_1291, t_1461, t_1462, t_1463, t_1482, exclude_diagnosis, exclude_tab = zip(*row)
    df['t_1243'] = pd.Series(t_1243)
    df['t_1286'] = pd.Series(t_1286)
    df['t_1287'] = pd.Series(t_1287)
    df['t_1288'] = pd.Series(t_1288)
    df['t_1289'] = pd.Series(t_1289)
    df['t_1290'] = pd.Series(t_1290)
    df['t_1291'] = pd.Series(t_1291)
    df['t_1461'] = pd.Series(t_1461)
    df['t_1462'] = pd.Series(t_1462)
    df['t_1463'] = pd.Series(t_1463)
    df['t_1482'] = pd.Series(t_1482)
    df['exclude_diagnosis'] = pd.Series(exclude_diagnosis)
    df['exclude_tab'] = pd.Series(exclude_tab)
    if not os.path.exists(data_path +"Self_report_diagnosis/"):
        os.makedirs(data_path +"Self_report_diagnosis/")
    df.to_csv(data_path + "Self_report_diagnosis/" + "Self_report_diagnosis.csv", columns=['eid', 't_1243', 't_1286', 't_1287', 't_1288', 't_1289', 't_1290', 't_1291', 't_1461', 't_1462', 't_1463', 't_1482', 'exclude_diagnosis', 'exclude_tab'], index=None)


if __name__ == '__main__':
    exclude_code = read_exclude_code()
    self_diagnosis(exclude_code)
