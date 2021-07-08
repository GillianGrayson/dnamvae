import pandas as pd
from src.dnam.routines.manifest import get_manifest
from src.dnam.routines.datasets_features import *
import os
import numpy as np
from tqdm import tqdm
import plotly.graph_objects as go
import statsmodels.formula.api as smf
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from src.dnam.EWAS.routines.correction import correct_pvalues
from src.dnam.routines.plot.save import save_figure
from src.dnam.routines.plot.scatter import add_scatter_trace
from src.dnam.routines.plot.layout import add_layout

platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"
datasets = ["GSE42861", "GSE80417", "GSE84727", "GSE125105", "GSE147221"]

is_rerun = True
num_cpgs_to_plot = 10

feats = {
    "DNAmPhenoAgeAcc": "DNAmPhenoAgeAcc",
    "DNAmGrimAgeAcc": "DNAmGrimAgeAcc"
}

for dataset in datasets:

    path_save = f"{path}/{platform}/{dataset}/EWAS/cpg_vs_continuous/control"

    status_pair = tuple([x.replace(' ','_') for x in get_status_pair(dataset)])
    age_pair = tuple([x.replace(' ','_') for x in get_age_pair(dataset)])
    sex_pair = tuple([x.replace(' ','_') for x in get_sex_pair(dataset)])
    status_vals_pairs = get_status_vals_pairs(dataset)
    status_vals = sorted([x for (x,y) in status_vals_pairs])
    sex_vals_pairs = get_sex_vals_pairs(dataset)

    pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
    pheno.columns = pheno.columns.str.replace(' ','_')
    betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")
    df = pd.merge(pheno, betas, left_index=True, right_index=True)
    df = df.loc[df[status_pair[0]].isin(status_vals), :]

    df_1 = df.loc[df[status_pair[0]] == status_vals_pairs[0][0], :]
    df_2 = df.loc[df[status_pair[0]] == status_vals_pairs[1][0], :]

    cpgs = betas.columns.values

    manifest = get_manifest(platform)

    for k, v in feats.items():

        df_1_curr = df_1[df_1[k].notnull()]
        df_2_curr = df_2[df_2[k].notnull()]

        path_curr = f"{path_save}/{v}/figs"
        if not os.path.exists(path_curr):
            os.makedirs(path_curr)

        if is_rerun:
            result = {'CpG': cpgs}
            result['Gene'] = np.zeros(len(cpgs), dtype=object)
            metrics = ['R2', 'R2_adj', f"{v}_pval", 'pearson_r', 'pearson_pval', 'spearman_r', 'spearman_pval']
            for m in metrics:
                result[m] = np.zeros(len(cpgs))

            for cpg_id, cpg in tqdm(enumerate(cpgs), desc='Regression', total=len(cpgs)):
                result['Gene'][cpg_id] = manifest.loc[cpg, 'Gene']
                reg = smf.ols(formula=f"{cpg} ~ {k}", data=df_1_curr).fit()
                pvalues = dict(reg.pvalues)
                result['R2'][cpg_id] = reg.rsquared
                result['R2_adj'][cpg_id] = reg.rsquared_adj
                result[f"{v}_pval"][cpg_id] = pvalues[k]
                pearson_r, pearson_pval = pearsonr(df_1_curr[cpg].values, df_1_curr[k].values)
                result['pearson_r'][cpg_id] = pearson_r
                result['pearson_pval'][cpg_id] = pearson_pval
                spearman_r, spearman_pval = spearmanr(df_1_curr[cpg].values, df_1_curr[k].values)
                result['spearman_r'][cpg_id] = spearman_r
                result['spearman_pval'][cpg_id] = spearman_pval

            result = correct_pvalues(result, [f"{v}_pval", 'pearson_pval', 'spearman_pval'])
            result = pd.DataFrame(result)
            result.set_index("CpG", inplace=True)
            result.sort_values([f"{v}_pval"], ascending=[True], inplace=True)
            result.to_excel(f"{path_save}/{v}/table.xlsx", index=True)
        else:
            result = pd.read_excel(f"{path_save}/{v}/table.xlsx", index_col="CpG")

        result = result.head(num_cpgs_to_plot)

        for cpg_id, (cpg, row) in enumerate(result.iterrows()):
            reg = smf.ols(formula=f"{cpg} ~ {k}", data=df_1_curr).fit()
            fig = go.Figure()
            add_scatter_trace(fig, df_1_curr[k].values, df_1_curr[cpg].values, status_vals_pairs[0][1])
            add_scatter_trace(fig, df_1_curr[k].values, reg.fittedvalues.values, "", "lines")
            add_scatter_trace(fig, df_2_curr[k].values, df_2_curr[cpg].values, status_vals_pairs[1][1])
            add_layout(fig, f"{v}", 'Methylation Level', f"{cpg} ({row['Gene']})")
            fig.update_layout({'colorway': ['blue', 'blue', "red"]})
            save_figure(fig, f"{path_curr}/{cpg_id}_{cpg}")
