import pandas as pd
from functools import reduce
import os
import pickle


dataset = "GSE55763"

save_path = f'E:/YandexDisk/Work/dnamvae/data/datasets/{dataset}'

load_fn = f"E:/YandexDisk/Work/pydnameth/{dataset}/observables.xlsx"
df = pd.read_excel(load_fn)
#df['id']= df['id'].astype(str)
#df['id'] = 'X' + df['id']
pheno = df.set_index('geo_accession')

load_fn = f"E:/YandexDisk/Work/pydnameth/{dataset}/betas.txt"
df = pd.read_csv(load_fn, delimiter = "\t", index_col='ID_REF')
beta = df.T

pheno_ids = pheno.index.values.tolist()
#tmps = pheno['id'].tolist()
beta_ids = beta.index.values.tolist()

print(f"is equal ids: {pheno_ids == beta_ids}")

#beta.index = pheno_ids

d = {'beta': beta, 'pheno': pheno}

f = open(f'{save_path}/data.pkl', 'wb')
pickle.dump(d, f, pickle.HIGHEST_PROTOCOL)
f.close()
