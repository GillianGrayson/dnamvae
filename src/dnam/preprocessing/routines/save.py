import pandas as pd


def save_pheno_betas_to_pkl(pheno: pd.DataFrame, betas: pd.DataFrame, path: str):
    pheno.to_pickle(f"{path}/pheno.pkl")
    betas.to_pickle(f"{path}/betas.pkl")
