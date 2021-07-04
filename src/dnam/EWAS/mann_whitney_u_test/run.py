import pandas as pd
from src.dnam.routines.manifest import get_manifest
from src.dnam.EWAS.mann_whitney_u_test.routines.process import perform_mann_whitney_u_test
from src.dnam.EWAS.mann_whitney_u_test.routines.plot import plot_mann_whitney_u_test
from src.dnam.routines.datasets_features import *


dataset = "GSE152027"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

status_pair = get_status_pair(dataset)
age_pair = get_age_pair(dataset)
sex_pair = get_sex_pair(dataset)
status_vals_pairs = get_status_vals_pairs(dataset)
sex_vals_pairs = get_sex_vals_pairs(dataset)

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno.pkl")
betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")

df_global = pd.merge(pheno, betas, left_index=True, right_index=True)
df_1 = df_global.loc[df_global[status_pair[0]] == status_vals_pairs[0][0], :]
df_2 = df_global.loc[df_global[status_pair[0]] == status_vals_pairs[1][0], :]

cpgs = betas.columns.values

manifest = get_manifest(platform)

result = perform_mann_whitney_u_test(df_1, df_2, cpgs, manifest, f"{path}/{platform}/{dataset}")
plot_mann_whitney_u_test(df_1, df_2, result, f"{path}/{platform}/{dataset}", 10, [status_vals_pairs[0][1], status_vals_pairs[1][1]])
