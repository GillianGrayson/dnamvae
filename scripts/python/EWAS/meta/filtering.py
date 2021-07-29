import pandas as pd
import numpy as np


platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

target = f"Age_Status"
percent_limit = 1
pval_suff = '_fdr_bh'
metrics = ['Age', 'Status']

path_load = f"{path}/{platform}/combo/EWAS/meta/{target}"
tbl = pd.read_excel(f"{path_load}/meta.xlsx", index_col="CpG")

intersection_tbl = tbl
for m in metrics:
    thld = np.percentile(tbl[f"{m}_pvalue{pval_suff}"].values, percent_limit)
    intersection_tbl = intersection_tbl.loc[tbl[f"{m}_pvalue{pval_suff}"] < thld, :]

intersection_tbl.to_excel(f"{path_load}/intersection_tbl.xlsx", index=True)
