import pandas as pd
from src.dnam.python.routines.manifest import get_manifest
from src.dnam.python.routines.plot.save import save_figure
from src.dnam.python.routines.plot.scatter import add_scatter_trace
from src.dnam.python.routines.plot.layout import add_layout
import os
import plotly.graph_objects as go
import statsmodels.formula.api as smf
import numpy as np
from tqdm import tqdm
from src.dnam.python.EWAS.routines.correction import correct_pvalues
from scipy.stats import pearsonr
from scipy.stats import spearmanr


dataset = "GSEUNN"
platform = "GPL21145"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

pheno_v1 = pd.read_excel(f"{path}/{platform}/{dataset}/raw/result/part(v1)_config(0.01_0.10_0.10)/pheno.xlsx", index_col="Sample_Name")
pheno_v2 = pd.read_excel(f"{path}/{platform}/{dataset}/raw/result/part(v2)_config(0.01_0.10_0.10)/pheno.xlsx", index_col="Sample_Name")

cols_to_use = pheno_v2.columns.difference(pheno_v1.columns)
pheno = pd.merge(pheno_v1, pheno_v2[cols_to_use], left_index=True, right_index=True, how='outer', indicator='is_v2')
pheno['is_v2'] = np.where(pheno['is_v2'] == 'both', True, False)

pheno.to_excel(f"{path}/{platform}/{dataset}/pheno.xlsx", index=True)
