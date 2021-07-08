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
from pingouin import ancova


platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"
datasets = ["GSE42861", "GSE80417", "GSE84727", "GSE125105", "GSE147221"]

is_rerun = True
num_cpgs_to_plot = 10

for dataset in datasets:
    print(dataset)

    status_pair = tuple([x.replace(' ','_') for x in get_status_pair(dataset)])
    age_pair = tuple([x.replace(' ','_') for x in get_age_pair(dataset)])
    sex_pair = tuple([x.replace(' ','_') for x in get_sex_pair(dataset)])
    status_vals_pairs = get_status_vals_pairs(dataset)
    sex_vals_pairs = get_sex_vals_pairs(dataset)

    cont_feat = age_pair[0]
    cont_show = age_pair[1]
    cat_feat = status_pair[0]
    cat_show = status_pair[1]

    status_vals = sorted([x for (x,y) in status_vals_pairs])
    terms = [cat_feat, cont_feat]
    aim = f"{cont_show}_{cat_show}"

    path_save = f"{path}/{platform}/{dataset}/EWAS/ancova/{aim}"
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
        for t in terms:
            result[f"{t}_pval"] = np.zeros(len(cpgs))

        for cpg_id, cpg in tqdm(enumerate(cpgs), desc='from_formula', total=len(cpgs)):
            result['Gene'][cpg_id] = manifest.loc[cpg, 'Gene']
            res = ancova(data=df, dv=cpg, covar=cont_feat, between=cat_feat)
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
        for (real, show) in status_vals_pairs:
            add_scatter_trace(fig,  df.loc[df[status_pair[0]] == real, cont_feat].values, df.loc[df[status_pair[0]] == real, cpg].values, show)
        add_layout(fig, cont_show, 'Methylation Level', f"{cpg} ({manifest.loc[cpg, 'Gene']})")
        fig.update_layout({'colorway': ['blue', "red"]})
        save_figure(fig, f"{path_save}/figs/{cpg_id}_{cpg}")
