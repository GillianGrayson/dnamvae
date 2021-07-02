import pandas as pd
from src.dnam.routines.manifest import get_manifest
from src.dnam.EWAS.regression.routines.process import perform_regression
from src.dnam.EWAS.regression.routines.plot import plot_regression_scatter
import os


dataset = "GSE152027"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

is_recalc = False

formula = "ageatbloodcollection * C(status)"
terms = ["ageatbloodcollection:C(status)[T.FEP]"]
aim = "age_status"

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno.pkl")
betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")

df = pd.merge(pheno, betas, left_index=True, right_index=True)
df = df[df['ageatbloodcollection'].notnull()]
df = df.loc[df['status'].isin(["CON", "FEP"]), :]

cpgs = betas.columns.values

manifest = get_manifest(platform)

if is_recalc or not os.path.isfile(f"{path}/{platform}/{dataset}/EWAS/regression/{aim}/table.xlsx"):
    result = perform_regression(df, cpgs, manifest, formula, terms, f"{path}/{platform}/{dataset}/EWAS/regression/{aim}")
else:
    result = pd.read_excel(f"{path}/{platform}/{dataset}/EWAS/regression/{aim}/table.xlsx", index_col="CpG")

plot_regression_scatter(df, ("ageatbloodcollection", "Age"), "status", {"Status: Control": "CON", "Status: First Epizode": "FEP"}, result, 10, f"{path}/{platform}/{dataset}/EWAS/regression/{aim}")
