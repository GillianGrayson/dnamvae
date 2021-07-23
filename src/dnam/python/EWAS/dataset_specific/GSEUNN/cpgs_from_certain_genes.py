import pandas as pd
from src.dnam.python.routines.manifest import get_manifest
import plotly.graph_objects as go
import numpy as np
import os
from tqdm import tqdm
from scipy.stats import mannwhitneyu
from src.dnam.python.EWAS.routines.correction import correct_pvalues
from src.dnam.python.routines.plot.save import save_figure
from src.dnam.python.routines.plot.layout import add_layout
from src.dnam.python.routines.plot.box import add_box_trace
from src.dnam.python.routines.datasets_features import *
from src.dnam.python.routines.filter.pheno import filter_pheno
from scipy.stats import spearmanr
import seaborn as sns
import matplotlib.pyplot as plt


dataset = "GSEUNN"
platform = "GPL21145"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

features_type = ['immuno', 'cytokines']
features = []
for ft in features_type:
    with open(f"{path}/{platform}/{dataset}/features/{ft}.txt") as f:
        f = f.read().splitlines()
        features.extend(f)
features = ['Age', 'DNAmAge', 'DNAmAgeHannum', 'DNAmPhenoAge', 'DNAmGrimAge', 'PhenoAge', 'ImmunoAge'] + features

manifest = get_manifest(platform)
genes = ['GDF15', 'FGF21', 'CXCL9']
manifest_trgt = manifest.loc[manifest['Gene'].isin(genes), :]

status_col = get_column_name(dataset, 'Status').replace(' ','_')
age_col = get_column_name(dataset, 'Age').replace(' ','_')
sex_col = get_column_name(dataset, 'Sex').replace(' ','_')
status_dict = get_status_dict(dataset)
status_names_dict = get_status_names_dict(dataset)
sex_dict = get_sex_dict(dataset)

path_save = f"{path}/{platform}/{dataset}/special/cpgs_from_certain_genes"
if not os.path.exists(f"{path_save}/figs"):
    os.makedirs(f"{path_save}/figs")

continuous_vars = {'Age': age_col}
categorical_vars = {status_col: status_dict, sex_col: sex_dict}
pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno.pkl")
pheno = filter_pheno(dataset, pheno, continuous_vars, categorical_vars)
betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")

df = pd.merge(pheno, betas, left_index=True, right_index=True)
df_ctrl = df.loc[(df[status_col] == status_dict['Control']), :]
df_case = df.loc[(df[status_col] == status_dict['Case']), :]

cpgs = list(set(manifest_trgt.index.values).intersection(set(betas.columns.values)))
manifest_trgt = manifest_trgt.loc[manifest_trgt.index.isin(cpgs), :]
cpgs = list(manifest_trgt.index.values)
cpgs_to_show = [f"{cpg}({manifest_trgt.loc[cpg, 'Gene']})" for cpg in cpgs]

corr_mtx = pd.DataFrame(data=np.zeros(shape=(len(features), len(cpgs))), index=features, columns=cpgs_to_show)
for f in features:
    for cpg_id, cpg in enumerate(cpgs):
        corr, _ = spearmanr(df[f], df[cpg])
        corr_mtx.loc[f, cpgs_to_show[cpg_id]] = corr

fig, ax = plt.subplots(figsize=(13,15))
sns_plot = sns.heatmap(corr_mtx, annot=True, vmin=-1, vmax=1, center=0, linewidths=0.5, ax=ax, annot_kws={'fontsize': 5}, cbar_kws={"shrink": 1.0, 'orientation': 'horizontal'})
sns_plot.figure.savefig(f"{path_save}/figs/corr_mtx.png")
sns_plot.figure.savefig(f"{path_save}/figs/corr_mtx.pdf")
