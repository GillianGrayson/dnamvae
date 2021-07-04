import numpy as np
import pandas as pd
import os
from src.dnam.routines.datasets_features import *


dataset = "GSE168739"
platform = "GPL21145"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

status_pair = get_status_pair(dataset)
age_pair = get_age_pair(dataset)
sex_pair = get_sex_pair(dataset)
status_vals_pairs = get_status_vals_pairs(dataset)
sex_vals_pairs = get_sex_vals_pairs(dataset)

save_path = f"{path}/{platform}/{dataset}/calculator"
if not os.path.exists(save_path):
    os.makedirs(save_path)

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno.pkl")
pheno = pheno[[age_pair[0], sex_pair[0]]]
pheno[sex_pair[0]] = pheno[sex_pair[0]].map({sex_vals_pairs[0][0]: 1, sex_vals_pairs[1][0]: 0})
pheno.rename(columns={age_pair[0]: 'Age', sex_pair[0]: 'Female'}, inplace=True)
pheno["Tissue"] = "Blood WB"
pheno.to_csv(f"{save_path}/pheno.csv", na_rep="NaN")

with open(f"{path}/calculator/cpgs_horvath_calculator.txt") as f:
    cpgs_h = f.read().splitlines()
betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")
cpgs_na = list(set(cpgs_h) - set(betas.columns.values))
betas = betas[betas.columns.intersection(cpgs_h)]
betas[cpgs_na] = np.nan
betas = betas.T
betas.to_csv(f"{save_path}/betas.csv", na_rep="NaN")
