import pandas as pd
import plotly.graph_objects as go
from src.dnam.routines.plot.save import save_figure
from src.dnam.routines.plot.histogram import add_histogram_trace
from src.dnam.routines.plot.layout import add_layout
from src.dnam.routines.datasets_features import *
import os


dataset = "GSE80417"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

is_update = True

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
pheno.columns = pheno.columns.str.replace(' ','_')

save_path = f"{path}/{platform}/{dataset}/pheno/distribution"
fig_path = f"{save_path}/figs"
if not os.path.exists(fig_path):
    os.makedirs(fig_path)

status_pair = tuple([x.replace(' ','_') for x in get_status_pair(dataset)])
age_pair = tuple([x.replace(' ','_') for x in get_age_pair(dataset)])
sex_pair = tuple([x.replace(' ','_') for x in get_sex_pair(dataset)])
status_vals_pairs = get_status_vals_pairs(dataset)
sex_vals_pairs = get_sex_vals_pairs(dataset)

df_1 = pheno.loc[(pheno[status_pair[0]] == status_vals_pairs[0][0]) & (pheno[age_pair[0]].notnull()), :]
df_2 = pheno.loc[(pheno[status_pair[0]] == status_vals_pairs[1][0]) & (pheno[age_pair[0]].notnull()), :]

fig = go.Figure()
add_histogram_trace(fig, df_1[age_pair[0]].values, f"{status_vals_pairs[0][1]} ({df_1.shape[0]})")
add_histogram_trace(fig, df_2[age_pair[0]].values, f"{status_vals_pairs[1][1]} ({df_2.shape[0]})")
add_layout(fig, age_pair[1], "Count", "")
fig.update_layout(colorway = ['blue', 'red'], barmode = 'overlay')
save_figure(fig, f"{fig_path}/histogram_{age_pair[1]}_{status_pair[1]}")
