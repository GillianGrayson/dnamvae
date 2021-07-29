import pandas as pd
from scripts.python.routines.manifest import get_manifest
import plotly.graph_objects as go
import numpy as np
import os
from tqdm import tqdm
from scipy.stats import mannwhitneyu
from scripts.python.EWAS.routines.correction import correct_pvalues
from scripts.python.routines.plot.save import save_figure
from scripts.python.routines.plot.layout import add_layout
from scripts.python.routines.plot.box import add_box_trace
from scripts.python.routines.filter.pheno import filter_pheno


dataset = "GSE53740"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

is_rerun = True
num_cpgs_to_plot = 10

status_col = get_column_name(dataset, 'Status').replace(' ','_')
age_col = get_column_name(dataset, 'Age').replace(' ','_')
sex_col = get_column_name(dataset, 'Sex').replace(' ','_')
status_dict = get_status_dict(dataset)
status_names_dict = get_status_names_dict(dataset)
sex_dict = get_sex_dict(dataset)

path_save = f"{path}/{platform}/{dataset}/EWAS/mann_whitney_u_test"
if not os.path.exists(f"{path_save}/figs"):
    os.makedirs(f"{path_save}/figs")

continuous_vars = {'Age': age_col}
categorical_vars = {status_col: status_dict, sex_col: sex_dict}
pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
pheno = filter_pheno(dataset, pheno, continuous_vars, categorical_vars)

betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")

df = pd.merge(pheno, betas, left_index=True, right_index=True)
df_1 = df.loc[(df[status_col] == status_dict['Control']), :]
df_2 = df.loc[(df[status_col] == status_dict['Case']), :]

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
    add_box_trace(fig, df_1[cpg].values, status_names_dict['Control'])
    add_box_trace(fig, df_2[cpg].values, status_names_dict['Case'])
    add_layout(fig, '', "Methylation Level", f"{cpg} ({manifest.loc[cpg, 'Gene']}): {row['pval']:0.4e}")
    fig.update_layout({'colorway': ['blue', "red"]})
    save_figure(fig, f"{path_save}/figs/{cpg_id}_{cpg}")
