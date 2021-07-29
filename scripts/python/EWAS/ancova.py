import pandas as pd
from src.dnam.python.routines.manifest import get_manifest
from tqdm import tqdm
from src.dnam.python.EWAS.routines.correction import correct_pvalues
import plotly.graph_objects as go
from src.dnam.python.routines.plot.save import save_figure
from src.dnam.python.routines.plot.scatter import add_scatter_trace
from src.dnam.python.routines.plot.layout import add_layout
import os
import numpy as np
from pingouin import ancova
from src.dnam.python.routines.datasets_features import *
from src.dnam.python.routines.filter.pheno import filter_pheno

platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"
datasets = ["GSE53740"]

is_rerun = True
num_cpgs_to_plot = 10

for dataset in datasets:
    print(dataset)

    status_col = get_column_name(dataset, 'Status').replace(' ', '_')
    age_col = get_column_name(dataset, 'Age').replace(' ', '_')
    sex_col = get_column_name(dataset, 'Sex').replace(' ', '_')
    status_dict = get_status_dict(dataset)
    status_vals = sorted(list(status_dict.values()))
    status_names_dict = get_status_names_dict(dataset)
    sex_dict = get_sex_dict(dataset)

    terms = [status_col, age_col]
    aim = f"Age_Status"

    path_save = f"{path}/{platform}/{dataset}/EWAS/ancova/{aim}"
    if not os.path.exists(f"{path_save}/figs"):
        os.makedirs(f"{path_save}/figs")

    continuous_vars = {'Age': age_col}
    categorical_vars = {status_col: status_dict, sex_col: sex_dict}
    pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
    pheno = filter_pheno(dataset, pheno, continuous_vars, categorical_vars)
    betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")
    df = pd.merge(pheno, betas, left_index=True, right_index=True)

    cpgs = betas.columns.values

    manifest = get_manifest(platform)

    if is_rerun:
        result = {'CpG': cpgs}
        result['Gene'] = np.zeros(len(cpgs), dtype=object)
        for t in terms:
            result[f"{t}_pval"] = np.zeros(len(cpgs))

        for cpg_id, cpg in tqdm(enumerate(cpgs), desc='from_formula', total=len(cpgs)):
            result['Gene'][cpg_id] = manifest.loc[cpg, 'Gene']
            res = ancova(data=df, dv=cpg, covar=age_col, between=status_col)
            for t in terms:
                result[f"{t}_pval"][cpg_id] = res.loc[res['Source'] == t, 'p-unc'].values[0]

        result = correct_pvalues(result, [f"{t}_pval" for t in terms])
        result = pd.DataFrame(result)
        result.set_index("CpG", inplace=True)
        result.sort_values([f"{t}_pval" for t in terms], ascending=[True] * len(terms), inplace=True)
        result.to_excel(f"{path_save}/table.xlsx", index=True)
    else:
        result = pd.read_excel(f"{path_save}/table.xlsx", index_col="CpG")

    result = result.head(num_cpgs_to_plot)
    for cpg_id, (cpg, row) in enumerate(result.iterrows()):
        fig = go.Figure()
        add_scatter_trace(fig,  df.loc[df[status_col] == status_dict['Control'], age_col].values, df.loc[df[status_col] == status_dict['Control'], cpg].values, status_names_dict['Control'])
        add_scatter_trace(fig, df.loc[df[status_col] == status_dict['Case'], age_col].values, df.loc[df[status_col] == status_dict['Case'], cpg].values, status_names_dict['Case'])
        add_layout(fig, "Age", 'Methylation Level', f"{cpg} ({manifest.loc[cpg, 'Gene']})")
        fig.update_layout({'colorway': ['blue', "red"]})
        save_figure(fig, f"{path_save}/figs/{cpg_id}_{cpg}")
