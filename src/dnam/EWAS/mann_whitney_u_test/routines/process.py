import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import mannwhitneyu
from src.dnam.EWAS.routines.correction import correct_pvalues
import os


def perform_mann_whitney_u_test(df_1: pd.DataFrame, df_2: pd.DataFrame, cpgs, manifest: pd.DataFrame, path: str):

    result = {'CpG': cpgs}
    result['Gene'] = np.zeros(len(cpgs), dtype=object)
    metrics = ['statistic', 'pvalue']
    for m in metrics:
        result[m] = np.zeros(len(cpgs))

    for cpg_id, cpg in tqdm(enumerate(cpgs), desc='Mann-Whitney U test', total=len(cpgs)):
        result['Gene'][cpg_id] = manifest.loc[cpg, 'Gene']
        data_1 = df_1[cpg].values
        data_2 = df_2[cpg].values
        statistic, pvalue = mannwhitneyu(data_1, data_2)
        result['statistic'][cpg_id] = statistic
        result['pvalue'][cpg_id] = pvalue

    result = correct_pvalues(result, ['pvalue'])

    save_path = f"{path}/EWAS/mann_whitney_u_test"
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    result = pd.DataFrame(result)
    result.set_index("CpG", inplace=True)
    result.sort_values(['pvalue'], ascending=[True], inplace=True)
    result.to_excel(f"{save_path}/table.xlsx", index=True)

    return result
