import pandas as pd
import plotly.graph_objects as go
from scripts.python.routines.plot.save import save_figure
from scripts.python.routines.plot.histogram import add_histogram_trace
from scripts.python.routines.plot.layout import add_layout
import os
from scripts.python.routines.filter.pheno import filter_pheno


dataset = "GSE53740"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

is_update = True

save_path = f"{path}/{platform}/{dataset}/pheno/distribution"
fig_path = f"{save_path}/figs"
if not os.path.exists(fig_path):
    os.makedirs(fig_path)

status_col = get_column_name(dataset, 'Status').replace(' ','_')
age_col = get_column_name(dataset, 'Age').replace(' ','_')
sex_col = get_column_name(dataset, 'Sex').replace(' ','_')
status_dict = get_status_dict(dataset)
get_status_names = get_status_names_dict(dataset)
sex_dict = get_sex_dict(dataset)

continuous_vars = {'Age': age_col}
categorical_vars = {status_col: status_dict, sex_col: sex_dict}
pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
pheno = filter_pheno(dataset, pheno, continuous_vars, categorical_vars)

df_1 = pheno.loc[(pheno[status_col] == status_dict['Control']), :]
df_2 = pheno.loc[(pheno[status_col] == status_dict['Case']), :]

fig = go.Figure()
add_histogram_trace(fig, df_1[age_col].values, f"{get_status_names['Control']} ({df_1.shape[0]})")
add_histogram_trace(fig, df_2[age_col].values, f"{get_status_names['Case']} ({df_2.shape[0]})")
add_layout(fig, "Age", "Count", "")
fig.update_layout(colorway = ['blue', 'red'], barmode = 'overlay')
save_figure(fig, f"{fig_path}/histogram_Age_Status")
