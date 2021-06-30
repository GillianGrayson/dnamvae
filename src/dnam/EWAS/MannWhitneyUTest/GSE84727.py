import pandas as pd
from src.dnam.EWAS.routines.manifest import get_manifest
from src.dnam.EWAS.MannWhitneyUTest.routines.process import perform_mann_whitney_u_test


dataset = "GSE84727"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno.pkl")
betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")

df_global = pd.merge(pheno, betas, left_index=True, right_index=True)
df_1 = df_global.loc[df_global['disease_status'] == 1, :]
df_2 = df_global.loc[df_global['disease_status'] == 2, :]

cpgs = betas.columns.values

manifest = get_manifest(platform)

perform_mann_whitney_u_test(df_1, df_2, cpgs, manifest, f"{path}/{platform}/{dataset}")
