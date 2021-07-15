import numpy as np
import pandas as pd
import os
from src.dnam.routines.datasets_features import *


dataset = "GSE144858"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

age_col = get_column_name(dataset, 'Age')
sex_col = get_column_name(dataset, 'Sex')
sex_dict = get_sex_dict(dataset)

save_path = f"{path}/{platform}/{dataset}/calculator"
if not os.path.exists(save_path):
    os.makedirs(save_path)

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno.pkl")
pheno = pheno[[age_col, sex_col]]
pheno[sex_col] = pheno[sex_col].map({sex_dict["F"]: 1, sex_dict["M"]: 0})
pheno.rename(columns={age_col: 'Age', sex_col: 'Female'}, inplace=True)
pheno["Tissue"] = "Blood WB"
pheno.to_csv(f"{save_path}/pheno.csv", na_rep="NA")

with open(f"{path}/calculator/cpgs_horvath_calculator.txt") as f:
    cpgs_h = f.read().splitlines()
betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")
cpgs_na = list(set(cpgs_h) - set(betas.columns.values))
betas = betas[betas.columns.intersection(cpgs_h)]
betas[cpgs_na] = np.nan
betas = betas.T
betas.to_csv(f"{save_path}/betas.csv", na_rep="NA")
