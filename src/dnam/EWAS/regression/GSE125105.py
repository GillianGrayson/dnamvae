import pandas as pd
from src.dnam.routines.manifest import get_manifest
from src.dnam.EWAS.regression.routines.process import perform_regression
from src.dnam.EWAS.regression.routines.plot import plot_regression_scatter
import os


dataset = "GSE125105"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

is_recalc = False

formula = "age * C(diagnosis)"
terms = ["age:C(diagnosis)[T.control]"]
aim = "age_status"

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno.pkl")
betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")

df = pd.merge(pheno, betas, left_index=True, right_index=True)
df = df[df['age'].notnull()]

cpgs = betas.columns.values

manifest = get_manifest(platform)

if is_recalc or not os.path.isfile(f"{path}/{platform}/{dataset}/EWAS/regression/{aim}/table.xlsx"):
    result = perform_regression(df, cpgs, manifest, formula, terms, f"{path}/{platform}/{dataset}/EWAS/regression/{aim}")
else:
    result = pd.read_excel(f"{path}/{platform}/{dataset}/EWAS/regression/{aim}/table.xlsx", index_col="CpG")

plot_regression_scatter(df, ("age", "Age"), "diagnosis", {"Status: Control": "control", "Status: Depression": "case"}, result, 10, f"{path}/{platform}/{dataset}/EWAS/regression/{aim}")
