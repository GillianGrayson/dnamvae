import pandas as pd
import os
import pickle
import numpy as np


th = 0.000
lme = "more"
m = "vt_score"

path = f"E:/YandexDisk/Work/dnamvae/data/datasets"
outcome = "age"
datasets = ["GSE40279", "GSE87571", "EPIC", "GSE55763"]

save_path = f"{path}/combo/{'_'.join(datasets)}/{m}_{lme}_{str(th)}"

target_cpgs = ['cg16867657']

f = open(f'{save_path}/data_nn.pkl', 'rb')
data = pickle.load(f)
f.close()

d = {outcome: data['pheno'][outcome].to_numpy()}
for cpg in target_cpgs:
    d[cpg] = data['beta'][cpg].to_numpy()

df = pd.DataFrame(d)
df.to_excel(f"{save_path}/check.xlsx", index=False)