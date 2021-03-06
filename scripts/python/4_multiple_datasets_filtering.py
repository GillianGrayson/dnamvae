import pandas as pd
import os
import pickle
import numpy as np
from sklearn.preprocessing import MinMaxScaler, StandardScaler

m = "vt_score"
lme = "more"
th = 0.005
scaler_type = 'none'

path = f"E:/YandexDisk/Work/dnamvae/data/datasets"
outcome = "age"
datasets = ["GSE40279", "GSE87571", "EPIC", "GSE55763"]

df = pd.read_excel(f"{path}/combo/{'_'.join(datasets)}.xlsx")

metrics = df[m].to_numpy()
if lme == "more":
    idx = np.where(metrics > th)[0]
    mask = metrics > th

cpgs = df["CpG"].to_numpy()
cpgs = cpgs[idx]

df_save = df.loc[idx, :]

save_path = f"{path}/combo/{'_'.join(datasets)}/{m}_{lme}_{str(th)}_{scaler_type}"
if not os.path.exists(save_path):
    os.makedirs(save_path)

df_save.to_excel(f"{save_path}/features.xlsx", index=False)

for d_id, d in enumerate(datasets):
    f = open(f"{path}/{d}/data_nn.pkl", 'rb')
    data = pickle.load(f)
    f.close()
    if d_id == 0:
        beta = data["beta"]
        pheno = data["pheno"].filter([outcome], axis=1)
    else:
        beta = beta.append(data["beta"])
        pheno = pheno.append(data["pheno"].filter([outcome], axis=1))

beta = beta[cpgs]

if scaler_type == 'minmax':
    scaler = MinMaxScaler()
    beta = pd.DataFrame(scaler.fit_transform(beta.values), columns=beta.columns, index=beta.index)
elif scaler_type == 'standard':
    scaler = StandardScaler()
    beta = pd.DataFrame(scaler.fit_transform(beta.values), columns=beta.columns, index=beta.index)

data = {"beta": beta, "pheno": pheno}

f = open(f'{save_path}/data_nn.pkl', 'wb')
pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)
f.close()