import pandas as pd
from scripts.python.routines.manifest import get_manifest
from scripts.python.routines.plot.save import save_figure
from scripts.python.routines.plot.scatter import add_scatter_trace
from scripts.python.routines.plot.layout import add_layout
import os
import plotly.graph_objects as go
import statsmodels.formula.api as smf
import numpy as np
from tqdm import tqdm
from scripts.python.EWAS.routines.correction import correct_pvalues
from scipy.stats import pearsonr
from scipy.stats import spearmanr


dataset = "GSE168739"
platform = "GPL21145"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

is_recalc = True

age_pair = tuple([x.replace(' ','_') for x in get_age_pair(dataset)])
sex_pair = tuple([x.replace(' ','_') for x in get_sex_pair(dataset)])
sex_vals_pairs = get_sex_vals_pairs(dataset)

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
pheno.columns = pheno.columns.str.replace(' ','_')
betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")

supp = pd.read_excel(f"{path}/{platform}/{dataset}/paper/suppl/mmc4.xls", skiprows=1, index_col="CpG ID")
cpgs = list(set.intersection(set(supp.index.values), set(betas.columns.values)))
betas = betas[cpgs]

df = pd.merge(pheno, betas, left_index=True, right_index=True)
df = df[df[age_pair[0]].notnull()]

manifest = get_manifest(platform)


ages_types = ["DNAmAge", "DNAmAgeHannum", "DNAmPhenoAge", "DNAmGrimAge"]
for at in ages_types:
    df[f"{at}Acc"] = df[at] - df[age_pair[0]]

metrics = {
    age_pair[0]: age_pair[1],
    "DNAmAgeAcc": "DNAmAgeAcc",
    "DNAmAgeHannumAcc": "DNAmAgeHannumAcc",
    "DNAmPhenoAgeAcc": "DNAmPhenoAgeAcc",
    "DNAmGrimAgeAcc": "DNAmGrimAgeAcc"
}

for k, v in metrics.items():
    formula = k
    terms = [k]
    aim = v

    save_path = f"{path}/{platform}/{dataset}/EWAS/regression/{aim}/figs"
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    result = {'CpG': cpgs}
    result['Gene'] = np.zeros(len(cpgs), dtype=object)
    metrics = ['R2', 'R2_adj', 'pearson_r', 'pearson_pval', 'spearman_r', 'spearman_pval']
    for m in metrics:
        result[m] = np.zeros(len(cpgs))
    for t in terms:
        result[f"{t}_pvalue"] = np.zeros(len(cpgs))

    for cpg_id, cpg in tqdm(enumerate(cpgs), desc='Regression', total=len(cpgs)):
        result['Gene'][cpg_id] = manifest.loc[cpg, 'Gene']
        reg = smf.ols(formula=f"{cpg} ~ {formula}", data=df).fit()
        pvalues = dict(reg.pvalues)
        result['R2'][cpg_id] = reg.rsquared
        result['R2_adj'][cpg_id] = reg.rsquared_adj
        for t in terms:
            result[f"{t}_pvalue"][cpg_id] = pvalues[t]
        pearson_r, pearson_pval = pearsonr(df[cpg].values, df[k].values)
        result['pearson_r'][cpg_id] = pearson_r
        result['pearson_pval'][cpg_id] = pearson_pval
        spearman_r, spearman_pval = spearmanr(df[cpg].values, df[k].values)
        result['spearman_r'][cpg_id] = spearman_r
        result['spearman_pval'][cpg_id] = spearman_pval

    result = correct_pvalues(result, [f"{t}_pvalue" for t in terms] + ['pearson_pval', 'spearman_pval'])
    result = pd.DataFrame(result)
    result.set_index("CpG", inplace=True)
    result.sort_values([f"{t}_pvalue" for t in terms], ascending=[True] * len(terms), inplace=True)
    result.to_excel(f"{path}/{platform}/{dataset}/EWAS/regression/{aim}/table.xlsx", index=True)

    for cpg_id, (cpg, row) in enumerate(result.iterrows()):
        reg = smf.ols(formula=f"{cpg} ~ {formula}", data=df).fit()
        fig = go.Figure()
        add_scatter_trace(fig, df[k].values, df[cpg].values, "")
        add_scatter_trace(fig, df[k].values, reg.fittedvalues.values, "", "lines")
        add_layout(fig, f"{v}", 'Methylation Level', f"{cpg} ({row['Gene']})")
        fig.update_layout({'colorway': ['blue', 'blue']})
        save_figure(fig, f"{save_path}/{cpg_id}_{cpg}")

