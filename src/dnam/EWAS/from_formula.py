import pandas as pd
from src.dnam.routines.manifest import get_manifest
from src.dnam.routines.datasets_features import *
from tqdm import tqdm
from src.dnam.EWAS.routines.correction import correct_pvalues
import statsmodels.formula.api as smf
import plotly.graph_objects as go
from src.dnam.routines.plot.save import save_figure
from src.dnam.routines.plot.scatter import add_scatter_trace
from src.dnam.routines.plot.layout import add_layout
import os
import numpy as np


platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"
datasets = ["GSE42861", "GSE53740", "GSE72774", "GSE80417", "GSE84727", "GSE87648", "GSE111629", "GSE125105", "GSE128235", "GSE144858", "GSE147221"]

is_rerun = True
num_cpgs_to_plot = 10

for dataset in datasets:
    print(dataset)

    status_col = get_column_name(dataset, 'Status').replace(' ', '_')
    age_col = get_column_name(dataset, 'Age').replace(' ', '_')
    sex_col = get_column_name(dataset, 'Sex').replace(' ', '_')
    status_dict = get_status_dict(dataset)
    status_vals = sorted(list(status_dict.values()))
    case_name = get_status_case_name(dataset)
    sex_dict = get_status_dict(dataset)

    dnam_acc_type = 'DNAmGrimAgeAcc'

    continuous_vars = {'Age': age_col, dnam_acc_type: dnam_acc_type}
    categorical_vars = {status_col: status_dict}

    formula = f"{age_col} + C({status_col}) + {dnam_acc_type}"
    terms = [f"{age_col}", f"C({status_col})[T.{status_vals[-1]}]", f"{dnam_acc_type}"]
    aim = f"Age_Status_{dnam_acc_type}"

    path_save = f"{path}/{platform}/{dataset}/EWAS/from_formula/{aim}"
    if not os.path.exists(f"{path_save}/figs"):
        os.makedirs(f"{path_save}/figs")

    pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
    pheno.columns = pheno.columns.str.replace(' ', '_')
    betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")

    df = pd.merge(pheno, betas, left_index=True, right_index=True)
    for name, feat in continuous_vars.items():
        df = df[df[feat].notnull()]
    for feat, groups in categorical_vars.items():
        df = df.loc[df[feat].isin(list(groups.values())), :]

    cpgs = betas.columns.values

    manifest = get_manifest(platform)

    if is_rerun:
        result = {'CpG': cpgs}
        result['Gene'] = np.zeros(len(cpgs), dtype=object)
        metrics = ['R2', 'R2_adj']
        for m in metrics:
            result[m] = np.zeros(len(cpgs))
        for t in terms:
            result[f"{t}_pvalue"] = np.zeros(len(cpgs))

        for cpg_id, cpg in tqdm(enumerate(cpgs), desc='from_formula', total=len(cpgs)):
            result['Gene'][cpg_id] = manifest.loc[cpg, 'Gene']
            reg = smf.ols(formula=f"{cpg} ~ {formula}", data=df).fit()
            pvalues = dict(reg.pvalues)
            result['R2'][cpg_id] = reg.rsquared
            result['R2_adj'][cpg_id] = reg.rsquared_adj
            for t in terms:
                result[f"{t}_pvalue"][cpg_id] = pvalues[t]

        result = correct_pvalues(result, [f"{t}_pvalue" for t in terms])
        result = pd.DataFrame(result)
        result.set_index("CpG", inplace=True)
        result.sort_values([f"{t}_pvalue" for t in terms], ascending=[True] * len(terms), inplace=True)
        result.to_excel(f"{path_save}/table.xlsx", index=True)
    else:
        result = pd.read_excel(f"{path_save}/table.xlsx", index_col="CpG")

    result = result.head(num_cpgs_to_plot)
    for cpg_id, (cpg, row) in enumerate(result.iterrows()):
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
