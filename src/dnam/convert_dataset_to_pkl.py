import pandas as pd
from functools import reduce
import os
import pickle


save_path = f'E:/YandexDisk/Work/dnamvae/data/datasets/GSE87571'

load_fn = "E:/YandexDisk/Work/pydnameth/script_datasets/GPL13534/filtered/blood(whole)/GSE87571/observables_part(wo_missedFeatures).xlsx"
df = pd.read_excel(load_fn)
#df['Sample_Name'] = 'X' + df['Sample_Name']
pheno = df.set_index('geo_accession')

load_fn = "E:/YandexDisk/Work/pydnameth/script_datasets/GPL13534/filtered/blood(whole)/GSE87571/betas_part(wo_missedFeatures)_config(0.01_0.10_0.10)_norm(fun).txt"
df = pd.read_csv(load_fn, delimiter = "\t", index_col='IlmnID')
beta = df.T

pheno_ids = pheno.index.values.tolist()
beta_ids = beta.index.values.tolist()

print(f"is equal ids: {pheno_ids == beta_ids}")


d = {'beta': beta, 'pheno':pheno}

f = open(f'{save_path}/data.pkl', 'wb')
pickle.dump(d, f, pickle.HIGHEST_PROTOCOL)
f.close()