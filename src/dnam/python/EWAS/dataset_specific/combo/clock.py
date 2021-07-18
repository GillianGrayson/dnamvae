import pandas as pd
from src.dnam.python.routines.manifest import get_manifest
import numpy as np
from src.dnam.python.routines.datasets_features import *
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import RepeatedKFold
from sklearn.metrics import mean_squared_error, mean_absolute_error
import os
from scipy.stats import mannwhitneyu
import statsmodels.formula.api as smf
import plotly.graph_objects as go
from src.dnam.python.routines.plot.save import save_figure
from src.dnam.python.routines.plot.scatter import add_scatter_trace
from src.dnam.python.routines.plot.violin import add_violin_trace
from src.dnam.python.routines.plot.layout import add_layout
import os


platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"
datasets = ["GSE147221", "GSE84727", "GSE125105", "GSE111629", "GSE128235", "GSE72774", "GSE53740", "GSE144858"]

dnam_acc_type = 'DNAmGrimAgeAcc'

target = f"Status_{dnam_acc_type}"
path_save = f"{path}/{platform}/combo/EWAS/meta/{target}"
if not os.path.exists(f"{path_save}/clock"):
    os.makedirs(f"{path_save}/clock")

manifest = get_manifest(platform)

pheno_all = pd.DataFrame(columns=['Age', 'Status'])
pheno_all.index.name = 'subject_id'
for d_id, dataset in enumerate(datasets):
    print(dataset)
    status_col = get_column_name(dataset, 'Status').replace(' ', '_')
    age_col = get_column_name(dataset, 'Age').replace(' ', '_')
    sex_col = get_column_name(dataset, 'Sex').replace(' ', '_')
    status_dict = get_status_dict(dataset)
    status_vals = sorted(list(status_dict.values()))
    status_names_dict = get_status_names_dict(dataset)
    sex_dict = get_sex_dict(dataset)

    continuous_vars = {'Age': age_col}
    categorical_vars = {status_col: status_dict}

    pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
    pheno.columns = pheno.columns.str.replace(' ', '_')
    betas = pd.read_pickle(f"{path}/{platform}/{dataset}/betas.pkl")

    df = pd.merge(pheno, betas, left_index=True, right_index=True)
    for name, feat in continuous_vars.items():
        df = df[df[feat].notnull()]
    for feat, groups in categorical_vars.items():
        df = df.loc[df[feat].isin(list(groups.values())), :]

    pheno = df[[age_col, status_col]]
    status_dict_inverse = dict((v, k) for k, v in status_dict.items())
    pheno[status_col] = pheno[status_col].map(status_dict_inverse)
    pheno.rename(columns={age_col: 'Age', status_col: 'Status'}, inplace=True)
    pheno_all = pheno_all.append(pheno, verify_integrity=True)

    cpgs = betas.columns.values
    betas = df[cpgs].T
    if d_id == 0:
        betas_all = betas
    else:
        betas_all = betas_all.merge(betas, how='inner', left_index=True, right_index=True)

betas_all = betas_all.T
betas_all.index.name = "subject_id"
df_all = pd.merge(pheno_all, betas_all, left_index=True, right_index=True)

with open(f"cpgs.txt") as f:
    cpgs_target = f.read().splitlines()

X_target = df_all.loc[df_all['Status'] == 'Control', cpgs_target].to_numpy()
y_target = df_all.loc[df_all['Status'] == 'Control', 'Age'].to_numpy()

X_all = df_all.loc[:, cpgs_target].to_numpy()
y_all = df_all.loc[:, 'Age'].to_numpy()

cv = RepeatedKFold(n_splits=5, n_repeats=5, random_state=1337)
model = ElasticNetCV(n_alphas=100, cv=cv, n_jobs=2)
model.fit(X_target, y_target)

model_dict = {'feature': ['Intercept'], 'coef': [model.intercept_]}
num_features = 0
for cpg_id, cpg in enumerate(cpgs_target):
    coef = model.coef_[cpg_id]
    if abs(coef) > 0:
        model_dict['feature'].append(cpg)
        model_dict['coef'].append(coef)
        num_features += 1
model_df = pd.DataFrame(model_dict)
if not os.path.exists(f"{path_save}/clock/{num_features}"):
    os.makedirs(f"{path_save}/clock/{num_features}")
model_df.to_excel(f"{path_save}/clock/{num_features}/clock.xlsx", index=False)

metrics_dict = {'alpha': model.alpha_, 'l1_ratio': model.l1_ratio_, 'num_features': num_features}
y_target_pred = model.predict(X_target)
metrics_dict['R2_Control'] = model.score(X_target, y_target)
metrics_dict['RMSE_Control'] = np.sqrt(mean_squared_error(y_target_pred, y_target))
metrics_dict['MAE_Control'] = mean_absolute_error(y_target_pred, y_target)
y_all_pred = model.predict(X_all)
metrics_dict['R2_All'] = model.score(X_all, y_all)
metrics_dict['RMSE_All'] = np.sqrt(mean_squared_error(y_all_pred, y_all))
metrics_dict['MAE_All'] = mean_absolute_error(y_all_pred, y_all)
metrics_df = pd.DataFrame(metrics_dict, index=[0])
metrics_df.to_excel(f"{path_save}/clock/{num_features}/metrics.xlsx", index=False)

pheno_all[f'AgeEST'] = y_all_pred

formula = f"AgeEST ~ Age"
reg = smf.ols(formula=formula, data=pheno_all.loc[pheno_all['Status'] == 'Control', :]).fit()
res_dict = {'R2': reg.rsquared, 'R2_adj': reg.rsquared_adj}
pheno_all['Acceleration'] = pheno_all[f'AgeEST'] - reg.predict(pheno_all)
pheno_all.to_excel(f"{path_save}/clock/{num_features}/pheno.xlsx", index=True)

scatter = go.Figure()
add_scatter_trace(scatter, pheno_all.loc[pheno_all['Status'] == 'Control', 'Age'].values, pheno_all.loc[pheno_all['Status'] == 'Control', 'AgeEST'].values, 'Control')
add_scatter_trace(scatter, pheno_all.loc[pheno_all['Status'] == 'Control', 'Age'].values, reg.fittedvalues.values, "", "lines")
add_scatter_trace(scatter, pheno_all.loc[pheno_all['Status'] == 'Case', 'Age'].values, pheno_all.loc[pheno_all['Status'] == 'Case', 'AgeEST'].values, 'Case')
add_layout(scatter, "Age", "AgeEST", "")
scatter.update_layout({'colorway': ['blue', 'blue', 'red']})
save_figure(scatter, f"{path_save}/clock/{num_features}/scatter_Age_AgeEST")

statistic, pvalue = mannwhitneyu(pheno_all.loc[pheno_all['Status'] == 'Control', 'Acceleration'].values, pheno_all.loc[pheno_all['Status'] == 'Case', 'Acceleration'].values)
res_dict['MW_statistic'] = statistic
res_dict['MW_pvalue'] = pvalue
res_df = pd.DataFrame(res_dict, index=[0])
res_df.to_excel(f"{path_save}/clock/{num_features}/res.xlsx", index=False)

box = go.Figure()
add_violin_trace(box, pheno_all.loc[pheno_all['Status'] == 'Control', 'Acceleration'].values, 'Control')
add_violin_trace(box, pheno_all.loc[pheno_all['Status'] == 'Case', 'Acceleration'].values, 'Case')
add_layout(box, "", 'Acceleration', f"")
box.update_layout({'colorway': ['blue', 'red']})
save_figure(box, f"{path_save}/clock/{num_features}/box_Acceleration")
