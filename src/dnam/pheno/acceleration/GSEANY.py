import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import statsmodels.formula.api as smf
import plotly.graph_objects as go
from src.dnam.routines.plot.routines import save_figure
from src.dnam.routines.plot.scatter import add_scatter_trace
from src.dnam.routines.plot.box import add_box_trace
from src.dnam.routines.plot.layout import add_layout
import os


dataset = "GSE42861"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")

save_path = f"{path}/{platform}/{dataset}/pheno/acceleration"
fig_path = f"{save_path}/figs"
if not os.path.exists(fig_path):
    os.makedirs(fig_path)

is_update = True

x_feat = "age"
x_name = "Age"

y_feats = ["DNAmAge", "DNAmAgeHannum", "DNAmPhenoAge", "DNAmGrimAge"]
y_names = ["DNAmAge", "DNAmAgeHannum", "DNAmPhenoAge", "DNAmGrimAge"]

res_names = ["DNAmAgeAcc", "DNAmAgeHannumAcc", "DNAmPhenoAgeAcc", "DNAmGrimAgeAcc"]

status_cat = ["Normal", "rheumatoid arthritis"]
status_val = ["Status: Normal", "Status: Rheumatoid Arthritis"]
pheno = pheno[pheno[x_feat].notnull()]
df_1 = pheno.loc[pheno['disease state'] == status_cat[0], :]
df_2 = pheno.loc[pheno['disease state'] == status_cat[1], :]

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
        pheno[res_names[y_id]] = reg.predict(pheno) - pheno[y]
        df_1 = pheno.loc[pheno['disease state'] == status_cat[0], :]
        df_2 = pheno.loc[pheno['disease state'] == status_cat[1], :]

    scatter = go.Figure()
    add_scatter_trace(scatter, df_1[x_feat].values, df_1[y].values, status_val[0])
    add_scatter_trace(scatter, df_1[x_feat].values, reg.fittedvalues.values, "", "lines")
    add_scatter_trace(scatter, df_2[x_feat].values, df_2[y].values, status_val[1])
    add_layout(scatter, x_name, y_names[y_id], "")
    scatter.update_layout({'colorway': ['blue', 'blue', 'red']})
    save_figure(scatter, f"{fig_path}/scatter_{x_name}_{y}")

    statistic, pvalue = mannwhitneyu(df_1[res_names[y_id]].values, df_2[res_names[y_id]].values)
    res_dict['MW_statistic'][y_id] = statistic
    res_dict['MW_pvalue'][y_id] = pvalue

    box = go.Figure()
    add_box_trace(box, df_1[res_names[y_id]].values, status_val[0])
    add_box_trace(box, df_2[res_names[y_id]].values, status_val[1])
    add_layout(box, "", y_names[y_id], f"{res_names[y_id]}: {pvalue:0.4e}")
    box.update_layout({'colorway': ['blue', 'red']})
    save_figure(box, f"{fig_path}/box_{res_names[y_id]}")

res_df = pd.DataFrame(res_dict)
res_df.set_index("metric", inplace=True)
res_df.to_excel(f"{save_path}/table.xlsx", index=True)

