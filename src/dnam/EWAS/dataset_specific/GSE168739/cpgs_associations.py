import pandas as pd
from src.dnam.routines.manifest import get_manifest
from src.dnam.EWAS.regression.routines.process import perform_regression
from src.dnam.EWAS.regression.routines.plot import plot_regression_scatter
from src.dnam.routines.plot.routines import save_figure
from src.dnam.routines.plot.scatter import add_scatter_trace
from src.dnam.routines.plot.layout import add_layout
from src.dnam.routines.datasets_features import *
import os
import plotly.graph_objects as go


dataset = "GSE168739"
platform = "GPL21145"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

is_recalc = True

age_pair = tuple([x.replace(' ','_') for x in get_age_pair(dataset)])
sex_pair = tuple([x.replace(' ','_') for x in get_sex_pair(dataset)])
sex_vals_pairs = get_sex_vals_pairs(dataset)

formula = f"{age_pair[0]}"
terms = [f"{age_pair[0]}"]
aim = "age"

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
pheno.columns = pheno.columns.str.replace(' ','_')
betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")

supp = pd.read_excel(f"{path}/{platform}/{dataset}/paper/suppl/mmc4.xls", skiprows=1, index_col="CpG ID")
cpgs = list(set.intersection(set(supp.index.values), set(betas.columns.values)))
betas = betas[cpgs]

df = pd.merge(pheno, betas, left_index=True, right_index=True)
df = df[df[age_pair[0]].notnull()]

manifest = get_manifest(platform)

if is_recalc or not os.path.isfile(f"{path}/{platform}/{dataset}/EWAS/regression/{aim}/table.xlsx"):
    result = perform_regression(df, cpgs, manifest, formula, terms, f"{path}/{platform}/{dataset}/EWAS/regression/{aim}")
else:
    result = pd.read_excel(f"{path}/{platform}/{dataset}/EWAS/regression/{aim}/table.xlsx", index_col="CpG")

save_path = f"{path}/{platform}/{dataset}/EWAS/regression/{aim}/figs"
if not os.path.exists(save_path):
    os.makedirs(save_path)
for cpg_id, (cpg, row) in enumerate(result.iterrows()):
    fig = go.Figure()
    add_scatter_trace(fig, df[age_pair[0]].values, df[cpg].values, "")
    add_layout(fig, age_pair[1], 'Methylation Level', f"{cpg} ({row['Gene']})")
    fig.update_layout({'colorway': ['blue']})
    save_figure(fig, f"{save_path}/{cpg_id}_{cpg}")

