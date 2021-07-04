import pandas as pd
from src.dnam.routines.manifest import get_manifest
from src.dnam.EWAS.regression.routines.process import perform_regression
from src.dnam.EWAS.regression.routines.plot import plot_regression_scatter
from src.dnam.routines.datasets_features import *
import os

dataset = "GSE125105"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

is_recalc = True

status_pair = tuple([x.replace(' ','_') for x in get_status_pair(dataset)])
age_pair = tuple([x.replace(' ','_') for x in get_age_pair(dataset)])
sex_pair = tuple([x.replace(' ','_') for x in get_sex_pair(dataset)])
status_vals_pairs = get_status_vals_pairs(dataset)
sex_vals_pairs = get_sex_vals_pairs(dataset)

formula = f"{age_pair[0]} * C({status_pair[0]})"
status_vals = sorted([x for (x,y) in status_vals_pairs])
terms = [f"{age_pair[0]}:C({status_pair[0]})[T.{status_vals[-1]}]"]
aim = "age_status"

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno.pkl")
pheno.columns = pheno.columns.str.replace(' ','_')
betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")

df = pd.merge(pheno, betas, left_index=True, right_index=True)
df = df[df[age_pair[0]].notnull()]
df = df.loc[df[status_pair[0]].isin(status_vals), :]

cpgs = betas.columns.values

manifest = get_manifest(platform)

if is_recalc or not os.path.isfile(f"{path}/{platform}/{dataset}/EWAS/regression/{aim}/table.xlsx"):
    result = perform_regression(df, cpgs, manifest, formula, terms, f"{path}/{platform}/{dataset}/EWAS/regression/{aim}")
else:
    result = pd.read_excel(f"{path}/{platform}/{dataset}/EWAS/regression/{aim}/table.xlsx", index_col="CpG")

plot_regression_scatter(df, age_pair, status_pair[0], status_vals_pairs, result, 10, f"{path}/{platform}/{dataset}/EWAS/regression/{aim}")
