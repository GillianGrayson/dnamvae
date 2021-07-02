import pandas as pd
from src.dnam.routines.manifest import get_manifest
from src.dnam.EWAS.regression.routines.process import perform_regression
from src.dnam.EWAS.regression.routines.plot import plot_regression_scatter
import os

dataset = "GSE80417"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

is_recalc = False

formula = "age * C(disease_status)"
terms = ["age:C(disease_status)[T.2]"]
aim = "age_status"

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno.pkl")
betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")

df = pd.merge(pheno, betas, left_index=True, right_index=True)
df.rename(columns={'disease status': 'disease_status'}, inplace=True)
df = df[df['age'].notnull()]

cpgs = betas.columns.values

manifest = get_manifest(platform)

if is_recalc or not os.path.isfile(f"{path}/{platform}/{dataset}/EWAS/regression/{aim}/table.xlsx"):
    result = perform_regression(df, cpgs, manifest, formula, terms, f"{path}/{platform}/{dataset}/EWAS/regression/{aim}")
else:
    result = pd.read_excel(f"{path}/{platform}/{dataset}/EWAS/regression/{aim}/table.xlsx", index_col="CpG")

plot_regression_scatter(df, ("age", "Age"), "subject", {"Status: Control": 1, "Status: Schizophrenia": 2}, result, 10, f"{path}/{platform}/{dataset}/EWAS/regression/{aim}")
