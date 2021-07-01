import pandas as pd
from functools import reduce
import os
import pickle

path = f'E:/YandexDisk/Work/dnamvae/data/annotation'

fn_pkl = f"{path}/ann.pkl"
if os.path.isfile(fn_pkl):

    f = open(fn_pkl, 'rb')
    annotations = pickle.load(f)
    f.close()

else:

    df_27k = pd.read_excel(f'{path}/ann27k.xlsx')
    df_450k = pd.read_excel(f'{path}/ann450k.xlsx')
    df_850k = pd.read_excel(f'{path}/ann850k.xlsx')

    annotations = [df_27k, df_450k, df_850k]

    f = open(fn_pkl, 'wb')
    pickle.dump(annotations, f, pickle.HIGHEST_PROTOCOL)
    f.close()


df_27k = annotations[0]
df_450k = annotations[1]
df_850k = annotations[2]

add_450k_cpgs = list(set(df_450k['CpG'].to_list()) - set(df_850k['CpG'].to_list()))
df_450k_add = df_450k[df_450k['CpG'].isin(add_450k_cpgs)]
df_450k_add.sort_values(['chr', 'Position'], ascending=[True, True])

ololo = set.intersection(set(df_450k_add['CpG'].to_list()), set(df_850k['CpG'].to_list()))

df_res = pd.concat([df_850k, df_450k_add], ignore_index=True)

add_27k_cpgs = list(set(df_27k['CpG'].to_list()) - set(df_res['CpG'].to_list()))
df_27k_add = df_27k[df_27k['CpG'].isin(add_27k_cpgs)]
df_27k_add.sort_values(['chr', 'Position'], ascending=[True, True])

df_res = pd.concat([df_res, df_27k_add], ignore_index=True)


df_res.to_excel(f'{path}/ann_full.xlsx', index=False)

f = open(f'{path}/cpgs.pkl', 'wb')
pickle.dump(df_res['CpG'].to_list(), f, pickle.HIGHEST_PROTOCOL)
f.close()

