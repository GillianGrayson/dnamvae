import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import statsmodels.formula.api as smf
import plotly.graph_objects as go
from src.dnam.routines.plot.save import save_figure
from src.dnam.routines.plot.scatter import add_scatter_trace
from src.dnam.routines.plot.box import add_box_trace
from src.dnam.routines.plot.layout import add_layout
from src.dnam.routines.datasets_features import *
import os


dataset = "GSE84727"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

is_update = True

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
pheno.columns = pheno.columns.str.replace(' ','_')

save_path = f"{path}/{platform}/{dataset}/pheno/acceleration"
fig_path = f"{save_path}/figs"
if not os.path.exists(fig_path):
    os.makedirs(fig_path)

status_pair = tuple([x.replace(' ','_') for x in get_status_pair(dataset)])
age_pair = tuple([x.replace(' ','_') for x in get_age_pair(dataset)])
sex_pair = tuple([x.replace(' ','_') for x in get_sex_pair(dataset)])
status_vals_pairs = get_status_vals_pairs(dataset)
sex_vals_pairs = get_sex_vals_pairs(dataset)

x_feat = age_pair[0]
x_name = age_pair[1]

y_feats = ["DNAmAge", "DNAmAgeHannum", "DNAmPhenoAge", "DNAmGrimAge"]
y_names = ["DNAmAge", "DNAmAgeHannum", "DNAmPhenoAge", "DNAmGrimAge"]
res_names = ["DNAmAgeAcc", "DNAmAgeHannumAcc", "DNAmPhenoAgeAcc", "DNAmGrimAgeAcc"]

df_1 = pheno.loc[(pheno[status_pair[0]] == status_vals_pairs[0][0]) & (pheno[x_feat].notnull()), :]
df_2 = pheno.loc[(pheno[status_pair[0]] == status_vals_pairs[1][0]) & (pheno[x_feat].notnull()), :]

metrics = ['R2', 'R2_adj', 'MW_statistic', 'MW_pvalue']
res_dict = {"metric": y_names}
for m in metrics:
    res_dict[m] = np.zeros(len(y_feats))

for y_id, y in enumerate(y_feats):
    formula = f"{y} ~ {x_feat}"
    reg = smf.ols(formula=formula, data=df_1).fit()
    res_dict['R2'][y_id] = reg.rsquared
    res_dict['R2_adj'][y_id] = reg.rsquared_adj

    if is_update:
        pheno[res_names[y_id]] = pheno[y] - reg.predict(pheno)
        df_1 = pheno.loc[(pheno[status_pair[0]] == status_vals_pairs[0][0]) & (pheno[x_feat].notnull()), :]
        df_2 = pheno.loc[(pheno[status_pair[0]] == status_vals_pairs[1][0]) & (pheno[x_feat].notnull()), :]

    scatter = go.Figure()
    add_scatter_trace(scatter, df_1[x_feat].values, df_1[y].values, status_vals_pairs[0][1])
    add_scatter_trace(scatter, df_1[x_feat].values, reg.fittedvalues.values, "", "lines")
    add_scatter_trace(scatter, df_2[x_feat].values, df_2[y].values, status_vals_pairs[1][1])
    add_layout(scatter, x_name, y_names[y_id], "")
    scatter.update_layout({'colorway': ['blue', 'blue', 'red']})
    save_figure(scatter, f"{fig_path}/scatter_{x_name}_{y}")

    statistic, pvalue = mannwhitneyu(df_1[res_names[y_id]].values, df_2[res_names[y_id]].values)
    res_dict['MW_statistic'][y_id] = statistic
    res_dict['MW_pvalue'][y_id] = pvalue

    box = go.Figure()
    add_box_trace(box, df_1[res_names[y_id]].values, status_vals_pairs[0][1])
    add_box_trace(box, df_2[res_names[y_id]].values, status_vals_pairs[1][1])
    add_layout(box, "", res_names[y_id], f"{res_names[y_id]}: {pvalue:0.4e}")
    box.update_layout({'colorway': ['blue', 'red']})
    save_figure(box, f"{fig_path}/box_{res_names[y_id]}")

res_df = pd.DataFrame(res_dict)
res_df.set_index("metric", inplace=True)
res_df.to_excel(f"{save_path}/table.xlsx", index=True)

pheno.to_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
pheno.to_excel(f"{path}/{platform}/{dataset}/pheno_xtd.xlsx", index=True)

