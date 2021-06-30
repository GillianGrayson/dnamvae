import pandas as pd


def get_pheno_betas_with_common_subjects(pheno: pd.DataFrame, betas: pd.DataFrame):

    pheno_subj_ids = list(pheno.index.values)
    betas_subj_ids = list(betas.index.values)
    if set(pheno_subj_ids) == set(betas_subj_ids):
        print(f"In pheno and betas subjects are the same")
    else:
        print(f"Warning! In pheno and betas subjects are not the same")
    if not pheno_subj_ids == betas_subj_ids:
        print(f"Warning! In pheno and betas subjects have different order")
        intersection = list(set(pheno_subj_ids).intersection(betas_subj_ids))
        pheno = pheno.loc[intersection, :]
        betas = betas.loc[intersection, :]
    else:
        print(f"In pheno and betas subjects have the same order")

    return pheno, betas
