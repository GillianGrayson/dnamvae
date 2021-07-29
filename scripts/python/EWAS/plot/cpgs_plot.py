import pandas as pd
from scripts.python.routines.manifest import get_manifest
import os
import plotly.graph_objects as go
from scripts.python.routines.plot.save import save_figure
from scripts.python.routines.plot.scatter import add_scatter_trace
from scripts.python.routines.plot.layout import add_layout
import statsmodels.formula.api as smf
from scripts.python.routines.filter.pheno import filter_pheno


dataset = "GSE147221"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

path_save = f"{path}/{platform}/{dataset}/EWAS/cpgs_plot"
path_save = f"{path_save}"
if not os.path.exists(path_save):
    os.makedirs(path_save)

status_col = get_column_name(dataset, 'Status').replace(' ', '_')
age_col = get_column_name(dataset, 'Age').replace(' ', '_')
sex_col = get_column_name(dataset, 'Sex').replace(' ', '_')
status_dict = get_status_dict(dataset)
status_vals = sorted(list(status_dict.values()))
get_status_names = get_status_names_dict(dataset)
sex_dict = get_sex_dict(dataset)

continuous_vars = {'Age': age_col, 'DNAmGrimAgeAcc': 'DNAmGrimAgeAcc'}
categorical_vars = {status_col: status_dict, sex_col: sex_dict}
pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
pheno = filter_pheno(dataset, dataset, pheno, continuous_vars, categorical_vars)
betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")
df = pd.merge(pheno, betas, left_index=True, right_index=True)

with open(f"cpgs_to_plot.txt") as f:
    cpgs = f.read().splitlines()

manifest = get_manifest(platform)

for cpg_id, cpg in enumerate(cpgs):
    for name_cont, feat_cont in continuous_vars.items():
        fig = go.Figure()
        for feat, groups in categorical_vars.items():
            for group_show, group_val in groups.items():
                df_curr = df.loc[df[feat] == group_val, :]
                reg = smf.ols(formula=f"{cpg} ~ {feat_cont}", data=df_curr).fit()
                add_scatter_trace(fig, df_curr[feat_cont].values, df_curr[cpg].values, group_show)
                add_scatter_trace(fig, df_curr[feat_cont].values, reg.fittedvalues.values, "", "lines")
            add_layout(fig, name_cont, 'Methylation Level', f"{cpg} ({manifest.loc[cpg, 'Gene']})")

            fig.update_layout({'colorway': ['blue', 'blue', "red", "red"]})

            if not os.path.exists(f"{path_save}/figs/{name_cont}"):
                os.makedirs(f"{path_save}/figs/{name_cont}")
            save_figure(fig, f"{path_save}/figs/{name_cont}/{cpg_id}_{cpg}")
