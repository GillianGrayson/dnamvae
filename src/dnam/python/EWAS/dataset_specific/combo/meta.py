import pandas as pd
from src.dnam.python.routines.manifest import get_manifest
from src.dnam.python.EWAS.routines.correction import correct_pvalues
import os
from scipy.stats import norm
import numpy as np
from src.dnam.python.routines.datasets_features import *


platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"
datasets = ["GSE147221", "GSE84727", "GSE125105", "GSE111629", "GSE128235", "GSE72774", "GSE53740", "GSE144858"]

dnam_acc_type = 'DNAmGrimAgeAcc'

target = f"Status_{dnam_acc_type}"
path_save = f"{path}/{platform}/combo/EWAS/meta/{target}"
if not os.path.exists(f"{path_save}"):
    os.makedirs(f"{path_save}")

tables_aim = f"Age_Status_{dnam_acc_type}"

pval_suff = ''

manifest = get_manifest(platform)
tables = manifest[['Gene']]

sizes = {}
signs = {}
cols = {'meta':
    [
        f'Age_pvalue{pval_suff}',
        f'Status_pvalue{pval_suff}',
        f"{dnam_acc_type}_pvalue{pval_suff}"
    ]
}
for dataset in datasets:
    print(dataset)
    status_col = get_column_name(dataset, 'Status').replace(' ', '_')
    age_col = get_column_name(dataset, 'Age').replace(' ', '_')
    sex_col = get_column_name(dataset, 'Sex').replace(' ', '_')
    status_dict = get_status_dict(dataset)
    status_vals = sorted(list(status_dict.values()))
    case_name = get_status_names_dict(dataset)
    sex_dict = get_status_dict(dataset)

    cols[dataset] = [
        f"{age_col}_pvalue{pval_suff}",
        f"C({status_col})[T.{status_vals[-1]}]_pvalue{pval_suff}",
        f"{dnam_acc_type}_pvalue{pval_suff}"
    ]

    continuous_vars = {'Age': age_col, dnam_acc_type: dnam_acc_type}
    categorical_vars = {status_col: status_dict}
    pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
    pheno.columns = pheno.columns.str.replace(' ', '_')
    for name, feat in continuous_vars.items():
        pheno = pheno[pheno[feat].notnull()]
    for feat, groups in categorical_vars.items():
        pheno = pheno.loc[pheno[feat].isin(list(groups.values())), :]
    sizes[dataset] = pheno.shape[0]

    signs[dataset] = 1

    path_load = f"{path}/{platform}/{dataset}/EWAS/from_formula/{tables_aim}"
    tbl = pd.read_excel(f"{path_load}/table.xlsx", index_col="CpG")
    tbl = tbl[cols[dataset]]
    new_cols = [x + f"_{dataset}" for x in cols[dataset]]
    cols[dataset] = new_cols
    tbl = tbl.add_suffix(f"_{dataset}")

    tables = tables.merge(tbl, how='inner', left_index=True, right_index=True)

nums = dict((col, np.zeros(tables.shape[0])) for col in cols['meta'])
dens = dict((col, 0) for col in cols['meta'])
for col_id, col in enumerate(cols['meta']):
    for dataset in datasets:
        if signs[dataset] < 0:
            zi = -norm.ppf(tables[cols[dataset][col_id]].values * 0.5)
        else:
            zi = norm.ppf(tables[cols[dataset][col_id]].values * 0.5)
        wi = np.sqrt(sizes[dataset])
        nums[col] += zi * wi
        dens[col] += wi * wi
    z = nums[col] / np.sqrt(dens[col])
    pvals = 2.0 * norm.cdf(-np.abs(z))
    tables[col] = pvals
result = tables[['Gene'] + cols['meta']]
result = correct_pvalues(result, cols['meta'])
result.sort_values(cols['meta'], ascending=[True] * len(cols['meta']), inplace=True)
result.to_excel(f"{path_save}/meta.xlsx", index=True)

# for dataset in datasets:
#     tables = tables.loc[(tables[cols[dataset][4]] < pval_thld) & (tables[cols[dataset][5]] <= pval_thld)]
# print(f"Number of CpGs: {tables.shape[0]}")
# tables.to_excel(f"{path_save}/table.xlsx", index=True)



