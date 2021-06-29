import pandas as pd
from functools import reduce
import os
import pickle


dataset = "GSE84727"
path = f"E:/YandexDisk/Work/pydnameth/datasets/450K/{dataset}"

fn = f"{path}/observables.xlsx"
df = pd.read_excel(fn)
#df['id']= df['id'].astype(str)
#df['id'] = 'X' + df['id']
pheno = df.set_index('geo_accession')

fn = f"{path}/raw/GSE84727_normalisedBetas.csv"
df = pd.read_csv(fn, delimiter=",")
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
