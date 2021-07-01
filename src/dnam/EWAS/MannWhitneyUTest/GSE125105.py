import pandas as pd
from src.dnam.routines.manifest import get_manifest
from src.dnam.EWAS.MannWhitneyUTest.routines.process import perform_mann_whitney_u_test
from src.dnam.EWAS.MannWhitneyUTest.routines.plot import plot_mann_whitney_u_test


dataset = "GSE125105"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno.pkl")
betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")

df_global = pd.merge(pheno, betas, left_index=True, right_index=True)
df_1 = df_global.loc[df_global['diagnosis'] == "control", :]
df_2 = df_global.loc[df_global['diagnosis'] == "case", :]

cpgs = betas.columns.values

manifest = get_manifest(platform)

result = perform_mann_whitney_u_test(df_1, df_2, cpgs, manifest, f"{path}/{platform}/{dataset}")
plot_mann_whitney_u_test(df_1, df_2, result, f"{path}/{platform}/{dataset}", 10, ["Status: Control", "Status: Depression"])
