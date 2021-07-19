import pandas as pd
from src.dnam.python.routines.manifest import get_manifest
import statsmodels.formula.api as smf
import plotly.graph_objects as go
from src.dnam.python.routines.plot.save import save_figure
from src.dnam.python.routines.plot.scatter import add_scatter_trace
from src.dnam.python.routines.plot.layout import add_layout
from src.dnam.python.routines.plot.violin import add_violin_trace
import os
from src.dnam.python.routines.datasets_features import *
from src.dnam.python.routines.filter.pheno import filter_pheno


platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"
datasets = ["GSE147221", "GSE84727", "GSE125105", "GSE111629", "GSE128235", "GSE72774", "GSE53740", "GSE144858"]

dnam_acc_type = 'DNAmGrimAgeAcc'

target = f"Status_{dnam_acc_type}"
path_save = f"{path}/{platform}/combo/EWAS/meta/{target}"
if not os.path.exists(f"{path_save}"):
    os.makedirs(f"{path_save}")

num_cpgs_to_plot = 10

manifest = get_manifest(platform)

with open(f"cpgs.txt") as f:
    cpgs = f.read().splitlines()

for dataset in datasets:
    print(dataset)

    status_col = get_column_name(dataset, 'Status').replace(' ', '_')
    age_col = get_column_name(dataset, 'Age').replace(' ', '_')
    sex_col = get_column_name(dataset, 'Sex').replace(' ', '_')
    status_dict = get_status_dict(dataset)
    status_vals = sorted(list(status_dict.values()))
    status_names_dict = get_status_names_dict(dataset)
    sex_dict = get_status_dict(dataset)

    continuous_vars = {'Age': age_col, dnam_acc_type: dnam_acc_type}
    categorical_vars = {status_col: status_dict, sex_col: sex_dict}
    pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
    pheno = filter_pheno(pheno, continuous_vars, categorical_vars)
    betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")
    df = pd.merge(pheno, betas, left_index=True, right_index=True)

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
                if not os.path.exists(f"{path_save}/figs/{dataset}/{name_cont}"):
                    os.makedirs(f"{path_save}/figs/{dataset}/{name_cont}")
            save_figure(fig, f"{path_save}/figs/{dataset}/{name_cont}/{cpg_id}_{cpg}")

        fig = go.Figure()
        for k, v in status_dict.items():
            add_violin_trace(fig, df.loc[df[status_col] == v, cpg].values, status_names_dict[k])
        add_layout(fig, '', "Methylation Level", f"{cpg} ({manifest.loc[cpg, 'Gene']})")
        fig.update_layout({'colorway': ['blue', "red"]})
        if not os.path.exists(f"{path_save}/figs/{dataset}/status"):
            os.makedirs(f"{path_save}/figs/{dataset}/status")
        save_figure(fig, f"{path_save}/figs/{dataset}/status/{cpg_id}_{cpg}")
