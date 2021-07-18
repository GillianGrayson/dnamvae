import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import statsmodels.formula.api as smf
import plotly.graph_objects as go
from src.dnam.python.routines.plot.save import save_figure
from src.dnam.python.routines.plot.scatter import add_scatter_trace
from src.dnam.python.routines.plot.box import add_box_trace
from src.dnam.python.routines.plot.layout import add_layout
import os
from src.dnam.python.routines.datasets_features import *


dataset = "GSE156994"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

is_update = True

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
pheno.columns = pheno.columns.str.replace(' ', '_')

save_path = f"{path}/{platform}/{dataset}/pheno/acceleration"
fig_path = f"{save_path}/figs"
if not os.path.exists(fig_path):
    os.makedirs(fig_path)

status_col = get_column_name(dataset, 'Status').replace(' ','_')
age_col = get_column_name(dataset, 'Age').replace(' ','_')
sex_col = get_column_name(dataset, 'Sex').replace(' ','_')
status_dict = get_status_dict(dataset)
case_name = get_status_case_name(dataset)
sex_dict = get_status_dict(dataset)

y_feats = ["DNAmAge", "DNAmAgeHannum", "DNAmPhenoAge", "DNAmGrimAge"]
res_names = ["DNAmAgeAcc", "DNAmAgeHannumAcc", "DNAmPhenoAgeAcc", "DNAmGrimAgeAcc"]

df_1 = pheno.loc[(pheno[status_col] == status_dict['Control']) & (pheno[age_col].notnull()), :]
df_2 = pheno.loc[(pheno[status_col] == status_dict['Case']) & (pheno[age_col].notnull()), :]

metrics = ['R2', 'R2_adj', 'MW_statistic', 'MW_pvalue']
res_dict = {"metric": y_feats}
for m in metrics:
    res_dict[m] = np.zeros(len(y_feats))

for y_id, y in enumerate(y_feats):
    formula = f"{y} ~ {age_col}"
    reg = smf.ols(formula=formula, data=df_1).fit()
    res_dict['R2'][y_id] = reg.rsquared
    res_dict['R2_adj'][y_id] = reg.rsquared_adj

    if is_update:
        pheno[res_names[y_id]] = pheno[y] - reg.predict(pheno)
        df_1 = pheno.loc[(pheno[status_col] == status_dict['Control']) & (pheno[age_col].notnull()), :]
        df_2 = pheno.loc[(pheno[status_col] == status_dict['Case']) & (pheno[age_col].notnull()), :]

    scatter = go.Figure()
    add_scatter_trace(scatter, df_1[age_col].values, df_1[y].values, 'Control')
    add_scatter_trace(scatter, df_1[age_col].values, reg.fittedvalues.values, "", "lines")
    add_scatter_trace(scatter, df_2[age_col].values, df_2[y].values, case_name)
    add_layout(scatter, "Age", y, "")
    scatter.update_layout({'colorway': ['blue', 'blue', 'red']})
    save_figure(scatter, f"{fig_path}/scatter_Age_{y}")

    statistic, pvalue = mannwhitneyu(df_1[res_names[y_id]].values, df_2[res_names[y_id]].values)
    res_dict['MW_statistic'][y_id] = statistic
    res_dict['MW_pvalue'][y_id] = pvalue

    box = go.Figure()
    add_box_trace(box, df_1[res_names[y_id]].values, 'Control')
    add_box_trace(box, df_2[res_names[y_id]].values, case_name)
    add_layout(box, "", res_names[y_id], f"{res_names[y_id]}: {pvalue:0.4e}")
    box.update_layout({'colorway': ['blue', 'red']})
    save_figure(box, f"{fig_path}/box_{res_names[y_id]}")

res_df = pd.DataFrame(res_dict)
res_df.set_index("metric", inplace=True)
res_df.to_excel(f"{save_path}/table.xlsx", index=True)

pheno.to_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
pheno.to_excel(f"{path}/{platform}/{dataset}/pheno_xtd.xlsx", index=True)

