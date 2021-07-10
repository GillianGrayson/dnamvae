import pandas as pd
from src.dnam.routines.manifest import get_manifest
from src.dnam.routines.datasets_features import *
from tqdm import tqdm
from src.dnam.EWAS.routines.correction import correct_pvalues
import statsmodels.formula.api as smf
import plotly.graph_objects as go
from src.dnam.routines.plot.save import save_figure
from src.dnam.routines.plot.scatter import add_scatter_trace
from src.dnam.routines.plot.layout import add_layout
import os
import numpy as np


platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"
datasets = ["GSE80417", "GSE42861", "GSE84727", "GSE125105", "GSE147221"]

is_rerun = True
num_cpgs_to_plot = 10

for dataset in datasets:
    print(dataset)

    status_pair = tuple([x.replace(' ','_') for x in get_status_pair(dataset)])
    age_pair = tuple([x.replace(' ','_') for x in get_age_pair(dataset)])
    sex_pair = tuple([x.replace(' ','_') for x in get_sex_pair(dataset)])
    status_vals_pairs = get_status_vals_pairs(dataset)
    sex_vals_pairs = get_sex_vals_pairs(dataset)

    cont_feat = "DNAmPhenoAgeAcc"
    cont_show = "DNAmPhenoAgeAcc"

    formula = f"{cont_feat} * C({status_pair[0]})"
    status_vals = sorted([x for (x,y) in status_vals_pairs])
    terms = [f"{cont_feat}:C({status_pair[0]})[T.{status_vals[-1]}]", f"{cont_feat}", f"C({status_pair[0]})[T.{status_vals[-1]}]"]
    aim = f"{cont_feat}_status"

    path_save = f"{path}/{platform}/{dataset}/EWAS/from_formula/{aim}"
    if not os.path.exists(f"{path_save}/figs"):
        os.makedirs(f"{path_save}/figs")

    pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
    pheno.columns = pheno.columns.str.replace(' ', '_')
    betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")

    df = pd.merge(pheno, betas, left_index=True, right_index=True)
    df = df[df[cont_feat].notnull()]
    df = df.loc[df[status_pair[0]].isin(status_vals), :]

    cpgs = betas.columns.values

    manifest = get_manifest(platform)

    if is_rerun:
        result = {'CpG': cpgs}
        result['Gene'] = np.zeros(len(cpgs), dtype=object)
        metrics = ['R2', 'R2_adj']
        for m in metrics:
            result[m] = np.zeros(len(cpgs))
        for t in terms:
            result[f"{t}_pvalue"] = np.zeros(len(cpgs))

        for cpg_id, cpg in tqdm(enumerate(cpgs), desc='from_formula', total=len(cpgs)):
            result['Gene'][cpg_id] = manifest.loc[cpg, 'Gene']
            reg = smf.ols(formula=f"{cpg} ~ {formula}", data=df).fit()
            pvalues = dict(reg.pvalues)
            result['R2'][cpg_id] = reg.rsquared
            result['R2_adj'][cpg_id] = reg.rsquared_adj
            for t in terms:
                result[f"{t}_pvalue"][cpg_id] = pvalues[t]

        result = correct_pvalues(result, [f"{t}_pvalue" for t in terms])
        result = pd.DataFrame(result)
        result.set_index("CpG", inplace=True)
        result.sort_values([f"{t}_pvalue" for t in terms], ascending=[True] * len(terms), inplace=True)
        result.to_excel(f"{path_save}/table.xlsx", index=True)
    else:
        result = pd.read_excel(f"{path_save}/table.xlsx", index_col="CpG")

    result = result.head(num_cpgs_to_plot)
    for cpg_id, (cpg, row) in enumerate(result.iterrows()):
        fig = go.Figure()
        for (real, show) in status_vals_pairs:
            df_curr = df.loc[df[status_pair[0]] == real, :]
            reg = smf.ols(formula=f"{cpg} ~ {cont_feat}", data=df_curr).fit()
            add_scatter_trace(fig, df_curr[cont_feat].values, df_curr[cpg].values, show)
            add_scatter_trace(fig, df_curr[cont_feat].values, reg.fittedvalues.values, "", "lines")
        add_layout(fig, cont_show, 'Methylation Level', f"{cpg} ({manifest.loc[cpg, 'Gene']})")
        fig.update_layout({'colorway': ['blue', 'blue', "red", "red"]})
        save_figure(fig, f"{path_save}/figs/{cpg_id}_{cpg}")
