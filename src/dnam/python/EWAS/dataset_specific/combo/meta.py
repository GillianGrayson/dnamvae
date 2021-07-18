import pandas as pd
from src.dnam.python.routines.manifest import get_manifest
from tqdm import tqdm
from src.dnam.python.EWAS.routines.correction import correct_pvalues
import statsmodels.formula.api as smf
import plotly.graph_objects as go
from src.dnam.python.routines.plot.save import save_figure
from src.dnam.python.routines.plot.scatter import add_scatter_trace
from src.dnam.python.routines.plot.layout import add_layout
import os
import numpy as np
from src.dnam.python.routines.datasets_features import *


platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"
datasets = ["GSE147221", "GSE84727",]

dnam_acc_type = 'DNAmGrimAgeAcc'

target = f"Status_{dnam_acc_type}"
path_save = f"{path}/{platform}/combo/EWAS/tables_filter/{target}"
if not os.path.exists(f"{path_save}"):
    os.makedirs(f"{path_save}")

tables_aim = f"Age_Status_{dnam_acc_type}"

pval_suff = '_fdr_bh'
pval_thld = 1e-4

for d_id, dataset in enumerate(datasets):
    print(dataset)

    status_col = get_column_name(dataset, 'Status').replace(' ', '_')
    age_col = get_column_name(dataset, 'Age').replace(' ', '_')
    sex_col = get_column_name(dataset, 'Sex').replace(' ', '_')
    status_dict = get_status_dict(dataset)
    status_vals = sorted(list(status_dict.values()))
    case_name = get_status_case_name(dataset)
    sex_dict = get_status_dict(dataset)

    cols = [f"{age_col}_pvalue{pval_suff}", f"C({status_col})[T.{status_vals[-1]}]_pvalue{pval_suff}", f"{dnam_acc_type}_pvalue{pval_suff}"]

    path_load = f"{path}/{platform}/{dataset}/EWAS/from_formula/{tables_aim}"
    tbl = pd.read_excel(f"{path_load}/table.xlsx", index_col="CpG")
    tbl = tbl[cols]
    tbl = tbl.add_suffix(f"_{dataset}")

    if d_id == 0:
        tables = tbl
    else:
        tables = tables.merge(tbl, how='inner', left_index=True, right_index=True)

for dataset in datasets:

    status_col = get_column_name(dataset, 'Status').replace(' ', '_')
    age_col = get_column_name(dataset, 'Age').replace(' ', '_')
    sex_col = get_column_name(dataset, 'Sex').replace(' ', '_')
    status_dict = get_status_dict(dataset)
    status_vals = sorted(list(status_dict.values()))
    case_name = get_status_case_name(dataset)
    sex_dict = get_status_dict(dataset)

    tables = tables.loc[(tables[f"C({status_col})[T.{status_vals[-1]}]_pvalue{pval_suff}"] < pval_thld) & (tables[f"{dnam_acc_type}_pvalue{pval_suff}"] <= pval_thld)]

print(f"Number of CpGs: {tables.shape[0]}")
tables.to_excel(f"{path_save}/table.xlsx", index=True)



