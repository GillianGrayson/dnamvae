import pandas as pd
from src.dnam.routines.manifest import get_manifest
from src.dnam.routines.datasets_features import *
import os
import plotly.graph_objects as go
from src.dnam.routines.plot.save import save_figure
from src.dnam.routines.plot.scatter import add_scatter_trace
from src.dnam.routines.plot.layout import add_layout
import statsmodels.formula.api as smf


dataset = "GSE147221"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

path_save = f"{path}/{platform}/{dataset}/EWAS/cpgs_plot"
path_save = f"{path_save}/figs"
if not os.path.exists(path_save):
    os.makedirs(path_save)

status_pair = tuple([x.replace(' ','_') for x in get_status_pair(dataset)])
age_pair = tuple([x.replace(' ','_') for x in get_age_pair(dataset)])
sex_pair = tuple([x.replace(' ','_') for x in get_sex_pair(dataset)])
status_vals_pairs = get_status_vals_pairs(dataset)
status_vals = sorted([x for (x,y) in status_vals_pairs])
sex_vals_pairs = get_sex_vals_pairs(dataset)

x_feat = age_pair[0]
x_name = age_pair[1]

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
pheno.columns = pheno.columns.str.replace(' ','_')
betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")
df = pd.merge(pheno, betas, left_index=True, right_index=True)
df = df.loc[df[status_pair[0]].isin(status_vals), :]
df_1 = df.loc[(df[status_pair[0]] == status_vals_pairs[0][0]) & (df[x_feat].notnull()), :]
df_2 = df.loc[(df[status_pair[0]] == status_vals_pairs[1][0]) & (df[x_feat].notnull()), :]

with open(f"cpgs_to_plot.txt") as f:
    cpgs = f.read().splitlines()

manifest = get_manifest(platform)

for cpg_id, cpg in enumerate(cpgs):
    reg_1 = smf.ols(formula=f"{cpg} ~ {x_feat}", data=df_1).fit()
    reg_2 = smf.ols(formula=f"{cpg} ~ {x_feat}", data=df_2).fit()
    fig = go.Figure()
    add_scatter_trace(fig, df_1[x_feat].values, df_1[cpg].values, status_vals_pairs[0][1])
    add_scatter_trace(fig, df_1[x_feat].values, reg_1.fittedvalues.values, "", "lines")
    add_scatter_trace(fig, df_2[x_feat].values, df_2[cpg].values, status_vals_pairs[1][1])
    add_scatter_trace(fig, df_2[x_feat].values, reg_2.fittedvalues.values, "", "lines")
    add_layout(fig, x_name, 'Methylation Level', f"{cpg} ({manifest.loc[cpg, 'Gene']})")
    fig.update_layout({'colorway': ['blue', 'blue', 'red', 'red']})
    save_figure(fig, f"{path_save}/{cpg_id}_{cpg}")
