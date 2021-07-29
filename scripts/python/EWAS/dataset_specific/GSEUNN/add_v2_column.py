import pandas as pd
import numpy as np

dataset = "GSEUNN"
platform = "GPL21145"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

pheno_v1 = pd.read_excel(f"{path}/{platform}/{dataset}/raw/result/part(v1)_config(0.01_0.10_0.10)/pheno.xlsx", index_col="Sample_Name")
pheno_v2 = pd.read_excel(f"{path}/{platform}/{dataset}/raw/result/part(v2)_config(0.01_0.10_0.10)/pheno.xlsx", index_col="Sample_Name")

cols_to_use = pheno_v2.columns.difference(pheno_v1.columns)
pheno = pd.merge(pheno_v1, pheno_v2[cols_to_use], left_index=True, right_index=True, how='outer', indicator='is_v2')
pheno['is_v2'] = np.where(pheno['is_v2'] == 'both', True, False)

pheno.to_excel(f"{path}/{platform}/{dataset}/pheno.xlsx", index=True)
