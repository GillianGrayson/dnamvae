import pandas as pd
import plotly.graph_objects as go
from src.dnam.routines.plot.save import save_figure
from src.dnam.routines.plot.histogram import add_histogram_trace
from src.dnam.routines.plot.layout import add_layout
from src.dnam.routines.datasets_features import *
import os


dataset = "GSE144858"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

is_update = True

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
pheno.columns = pheno.columns.str.replace(' ','_')

save_path = f"{path}/{platform}/{dataset}/pheno/distribution"
fig_path = f"{save_path}/figs"
if not os.path.exists(fig_path):
    os.makedirs(fig_path)

status_col = get_column_name(dataset, 'Status').replace(' ','_')
age_col = get_column_name(dataset, 'Age').replace(' ','_')
sex_col = get_column_name(dataset, 'Sex').replace(' ','_')
status_dict = get_status_dict(dataset)
case_name = get_status_case_name(dataset)
sex_dict = get_status_dict(dataset)

df_1 = pheno.loc[(pheno[status_col] == status_dict['Control']) & (pheno[age_col].notnull()), :]
df_2 = pheno.loc[(pheno[status_col] == status_dict['Case']) & (pheno[age_col].notnull()), :]

fig = go.Figure()
add_histogram_trace(fig, df_1[age_col].values, f"Control ({df_1.shape[0]})")
add_histogram_trace(fig, df_2[age_col].values, f"{case_name} ({df_2.shape[0]})")
add_layout(fig, "Age", "Count", "")
fig.update_layout(colorway = ['blue', 'red'], barmode = 'overlay')
save_figure(fig, f"{fig_path}/histogram_Age_Status")
