import pandas as pd
from src.dnam.python.routines.manifest import get_manifest
from src.dnam.python.EWAS.routines.correction import correct_pvalues
import os
from scipy.stats import norm
import numpy as np
from src.dnam.python.routines.datasets_features import *
from src.dnam.python.routines.filter.pheno import filter_pheno
import upsetplot as upset
from matplotlib import pyplot


platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"
#datasets = ["GSE84727", "GSE147221", "GSE125105", "GSE111629", "GSE128235", "GSE72774", "GSE53740", "GSE144858"]
datasets = ["GSE84727", "GSE147221", "GSE125105", "GSE111629", "GSE128235", "GSE72774"]

dnam_acc_type = 'DNAmGrimAgeAcc'

target = f"Age_Status"
path_save = f"{path}/{platform}/combo/EWAS/meta/{target}"
if not os.path.exists(f"{path_save}"):
    os.makedirs(f"{path_save}")

tables_aim = f"Age_Status_{dnam_acc_type}"

pval_suff = '_fdr_bh'
pval_thld = 0.05

manifest = get_manifest(platform)
tables_meta = manifest[['Gene']]
tables_single = manifest[['Gene']]

sizes = {}
signs = {}
meta_cols = {'meta':
    [
        f'Age_pvalue',
        f'Status_pvalue',
    ]
}
single_cols = {}
for dataset in datasets:
    print(dataset)
    status_col = get_column_name(dataset, 'Status').replace(' ', '_')
    age_col = get_column_name(dataset, 'Age').replace(' ', '_')
    sex_col = get_column_name(dataset, 'Sex').replace(' ', '_')
    status_dict = get_status_dict(dataset)
    status_vals = sorted(list(status_dict.values()))
    status_names_dict = get_status_names_dict(dataset)
    sex_dict = get_sex_dict(dataset)

    meta_cols[dataset] = [
        f"{age_col}_pvalue",
        f"C({status_col})[T.{status_vals[-1]}]_pvalue",
    ]

    single_cols[dataset] = [
        f"{age_col}_pvalue{pval_suff}",
        f"C({status_col})[T.{status_vals[-1]}]_pvalue{pval_suff}",
    ]

    #single_cols[dataset] = [single_cols[dataset][-1]]

    continuous_vars = {'Age': age_col, dnam_acc_type: dnam_acc_type}
    categorical_vars = {status_col: status_dict, sex_col: sex_dict}
    pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
    pheno = filter_pheno(pheno, continuous_vars, categorical_vars)
    sizes[dataset] = pheno.shape[0]
    print(f'Number of subjects: {pheno.shape[0]}')
    print(f"Number of Control: {pheno.loc[pheno[status_col] == status_dict['Control'], :].shape[0]}")
    print(f"Number of Case: {pheno.loc[pheno[status_col] == status_dict['Case'], :].shape[0]}")

    signs[dataset] = 1

    path_load = f"{path}/{platform}/{dataset}/EWAS/from_formula/{tables_aim}"
    tbl = pd.read_excel(f"{path_load}/table.xlsx", index_col="CpG")

    tbl_meta = tbl[meta_cols[dataset]]
    new_cols = [x + f"_{dataset}" for x in meta_cols[dataset]]
    meta_cols[dataset] = new_cols
    tbl_meta = tbl_meta.add_suffix(f"_{dataset}")
    tables_meta = tables_meta.merge(tbl_meta, how='inner', left_index=True, right_index=True)

    tbl_single = tbl[single_cols[dataset]]
    new_cols = [x + f"_{dataset}" for x in single_cols[dataset]]
    single_cols[dataset] = new_cols
    tbl_single = tbl_single.add_suffix(f"_{dataset}")
    tables_single = tables_single.merge(tbl_single, how='inner', left_index=True, right_index=True)

print(f'Total number of subjects: {sum(sizes.values())}')

upset_df = pd.DataFrame(index=tables_single.index)
for dataset in datasets:
    upset_df[dataset] = (tables_single[single_cols[dataset][0]] < pval_thld)
    for col in  single_cols[dataset][1::]:
        upset_df[dataset] =  upset_df[dataset] & (tables_single[col] < pval_thld)
upset_df = upset_df.set_index(datasets)
plt = upset.UpSet(upset_df, subset_size='count', show_counts=True).plot()
pyplot.savefig(f"{path_save}/single.png", bbox_inches='tight')
pyplot.savefig(f"{path_save}/single.pdf", bbox_inches='tight')

for dataset in datasets:
    for col in single_cols[dataset]:
        tables_single = tables_single.loc[(tables_single[col] < pval_thld), :]
        print(f"Number of CpGs: {tables_single.shape[0]}")
tables_single.to_excel(f"{path_save}/single.xlsx", index=True)

nums = dict((col, np.zeros(tables_meta.shape[0])) for col in meta_cols['meta'])
dens = dict((col, 0) for col in meta_cols['meta'])
for col_id, col in enumerate(meta_cols['meta']):
    for dataset in datasets:
        if signs[dataset] < 0:
            zi = -norm.ppf(tables_meta[meta_cols[dataset][col_id]].values * 0.5)
        else:
            zi = norm.ppf(tables_meta[meta_cols[dataset][col_id]].values * 0.5)
        wi = np.sqrt(sizes[dataset])
        nums[col] += zi * wi
        dens[col] += wi * wi
    z = nums[col] / np.sqrt(dens[col])
    pvals = 2.0 * norm.cdf(-np.abs(z))
    tables_meta[col] = pvals
result = tables_meta[['Gene'] + meta_cols['meta']]
result = correct_pvalues(result, meta_cols['meta'])
result.sort_values(meta_cols['meta'], ascending=[True] * len(meta_cols['meta']), inplace=True)
result.to_excel(f"{path_save}/meta.xlsx", index=True)










