import numpy as np
import pandas as pd
import os


dataset = "GSE42861"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

save_path = f"{path}/{platform}/{dataset}/calculator"
if not os.path.exists(save_path):
    os.makedirs(save_path)

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno.pkl")
pheno = pheno[["age", "gender"]]
pheno.replace({"gender": {"f": 1, "m": 0}}, inplace=True)
pheno.rename(columns={'age': 'Age', 'gender': 'Female'}, inplace=True)
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
