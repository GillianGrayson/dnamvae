import pandas as pd
from src.dnam.routines.manifest import get_manifest
import plotly.graph_objects as go
from src.dnam.routines.datasets_features import *
import numpy as np
import os
from tqdm import tqdm
from scipy.stats import mannwhitneyu
from src.dnam.EWAS.routines.correction import correct_pvalues
from src.dnam.routines.plot.save import save_figure
from src.dnam.routines.plot.layout import add_layout
from src.dnam.routines.plot.box import add_box_trace


dataset = "GSE125105"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

is_rerun = True
num_cpgs_to_plot = 10

status_pair = get_status_pair(dataset)
age_pair = get_age_pair(dataset)
sex_pair = get_sex_pair(dataset)
status_vals_pairs = get_status_vals_pairs(dataset)
sex_vals_pairs = get_sex_vals_pairs(dataset)

path_save = f"{path}/{platform}/{dataset}/EWAS/mann_whitney_u_test"
if not os.path.exists(f"{path_save}/figs"):
    os.makedirs(f"{path_save}/figs")

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")

df_global = pd.merge(pheno, betas, left_index=True, right_index=True)
df_1 = df_global.loc[df_global[status_pair[0]] == status_vals_pairs[0][0], :]
df_2 = df_global.loc[df_global[status_pair[0]] == status_vals_pairs[1][0], :]

cpgs = betas.columns.values

manifest = get_manifest(platform)

if is_rerun:
    result = {'CpG': cpgs}
    result['Gene'] = np.zeros(len(cpgs), dtype=object)
    metrics = ['statistic', 'pval']
    for m in metrics:
        result[m] = np.zeros(len(cpgs))

    for cpg_id, cpg in tqdm(enumerate(cpgs), desc='Mann-Whitney U test', total=len(cpgs)):
        result['Gene'][cpg_id] = manifest.loc[cpg, 'Gene']
        data_1 = df_1[cpg].values
        data_2 = df_2[cpg].values
        statistic, pvalue = mannwhitneyu(data_1, data_2)
        result['statistic'][cpg_id] = statistic
        result['pval'][cpg_id] = pvalue

    result = correct_pvalues(result, ['pval'])
    result = pd.DataFrame(result)
    result.set_index("CpG", inplace=True)
    result.sort_values(['pval'], ascending=[True], inplace=True)
    result.to_excel(f"{path_save}/table.xlsx", index=True)
else:
    result = pd.read_excel(f"{path_save}/table.xlsx", index_col="CpG")

result = result.head(num_cpgs_to_plot)
for cpg_id, (cpg, row) in enumerate(result.iterrows()):

    fig = go.Figure()
    add_box_trace(fig, df_1[cpg].values, status_vals_pairs[0][1])
    add_box_trace(fig, df_2[cpg].values, status_vals_pairs[1][1])
    add_layout(fig, '', "Methylation Level", f"{cpg} ({manifest.loc[cpg, 'Gene']}): {row['pval']:0.4e}")
    fig.update_layout({'colorway': ['blue', "red"]})
    save_figure(fig, f"{path_save}/figs/{cpg_id}_{cpg}")
