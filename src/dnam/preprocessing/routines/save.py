import pandas as pd
import pickle


def save_pheno_betas_to_pkl(pheno: pd.DataFrame, betas: pd.DataFrame, path: str):
    d = {'betas': betas, 'pheno': pheno}
    f = open(f'{path}/data.pkl', 'wb')
    pickle.dump(d, f, pickle.HIGHEST_PROTOCOL)
    f.close()
