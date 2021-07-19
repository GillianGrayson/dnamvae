def filter_pheno(pheno, continuous_vars, categorical_vars):
    pheno.columns = pheno.columns.str.replace(' ', '_')
    for name, feat in continuous_vars.items():
        pheno = pheno[pheno[feat].notnull()]
    for feat, groups in categorical_vars.items():
        pheno = pheno.loc[pheno[feat].isin(list(groups.values())), :]
    return pheno