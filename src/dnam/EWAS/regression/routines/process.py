import numpy as np
import pandas as pd
from tqdm import tqdm
from src.dnam.EWAS.routines.correction import correct_pvalues
import statsmodels.formula.api as smf
import os


def perform_regression(df: pd.DataFrame, cpgs: list, manifest: pd.DataFrame, formula: str, terms: list, path: str):

    result = {'CpG': cpgs}
    result['Gene'] = np.asarray(['non-genic'] * len(cpgs))
    metrics = ['R2', 'R2_adj']
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

    result = correct_pvalues(result, [f"{t}_pvalue" for t in terms])

    if not os.path.exists(path):
        os.makedirs(path)

    result = pd.DataFrame(result)
    result.set_index("CpG", inplace=True)
    result.sort_values([f"{t}_pvalue" for t in terms], ascending=[True] * len(terms), inplace=True)
    result.to_excel(f"{path}/table.xlsx", index=True)

    return result
