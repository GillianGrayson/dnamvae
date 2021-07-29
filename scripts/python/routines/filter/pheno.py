import pandas as pd


def filter_dataset(dataset: str, df: pd.DataFrame):
    if dataset == "GSEUNN":
        df = df.loc[df['is_v2'] == True, :]
    else:
        pass
    return df


def filter_pheno(dataset, pheno, continuous_vars, categorical_vars):
    pheno = filter_dataset(dataset, pheno)
    pheno.columns = pheno.columns.str.replace(' ', '_')
    for name, feat in continuous_vars.items():
        pheno = pheno[pheno[feat].notnull()]
    for feat, groups in categorical_vars.items():
        pheno = pheno.loc[pheno[feat].isin(list(groups.values())), :]
    return pheno