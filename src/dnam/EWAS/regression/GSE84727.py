import pandas as pd
from src.dnam.routines.manifest import get_manifest
from src.dnam.EWAS.regression.routines.process import perform_regression
from src.dnam.EWAS.regression.routines.plot import plot_regression_scatter

dataset = "GSE84727"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

formula = "age * C(disease_status)"
terms = ["age:C(disease_status)[T.2]"]
aim = "age_status"

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno.pkl")
betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")

df = pd.merge(pheno, betas, left_index=True, right_index=True)
df = df[df['age'].notnull()]

cpgs = betas.columns.values

manifest = get_manifest(platform)

#result = perform_regression(df, cpgs, manifest, formula, terms, f"{path}/{platform}/{dataset}/EWAS/regression/{aim}")
result = pd.read_excel(f"{path}/{platform}/{dataset}/EWAS/regression/{aim}/test.xlsx", index_col="CpG")
plot_regression_scatter(df, ("age", "Age"), "disease_status", {"Status: Control": 1, "Status: Schizophrenia": 2}, result, 10, f"{path}/{platform}/{dataset}/EWAS/regression/{aim}")
