import pandas as pd
from functools import reduce
import os
import pickle


save_path = f'E:/YandexDisk/Work/dnamvae/data/datasets/unn'

load_fn = "E:/YandexDisk/Work/pydnameth/unn_epic/observables_part(v1).xlsx"
df = pd.read_excel(load_fn)
df['Sample_Name'] = 'X' + df['Sample_Name']
pheno = df.set_index('Sample_Name')

load_fn = "E:/YandexDisk/Work/pydnameth/unn_epic/betas_part(v1)_config(0.01_0.10_0.10)_norm(fun).txt"
df = pd.read_csv(load_fn, delimiter = "\t", index_col='IlmnID')
beta = df.T

d = {'beta': beta, 'pheno':pheno}

f = open(f'{save_path}/data.pkl', 'wb')
pickle.dump(d, f, pickle.HIGHEST_PROTOCOL)
f.close()