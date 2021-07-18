import pandas as pd
from tqdm import tqdm

def get_forbidden_cpgs(path, types):
    forbidden_cpgs = set()
    for type in types:
        df = pd.read_excel(f"{path}/{type}.xlsx", header=None)
        tmp = df[0].values.tolist()
        forbidden_cpgs |= set(tmp)
    return list(forbidden_cpgs)


def betas_pvals_filter(betas: pd.DataFrame, pvals: pd.DataFrame, det_pval=0.01, det_cpg_cutoff=0.1):
    num_subjects = betas.shape[0]
    num_cpgs = betas.shape[1]
    passed_cols = []
    for col_id, col in tqdm(enumerate(betas.columns)):
        num_failed = (pvals[col].values > det_pval).sum()
        if num_failed < det_cpg_cutoff * num_subjects:
            passed_cols.append(col)
    print(f"Removing {num_cpgs - len(passed_cols)} failed CpGs with detection p-value above {det_pval}")
    betas = betas.loc[:, passed_cols]
    return betas


def manifest_filter(betas: pd.DataFrame, manifest: pd.DataFrame):
    betas_cpgs = betas.columns.values
    manifest_cpgs = manifest.index.values
    common_cpgs = list(set(betas_cpgs).intersection(set(manifest_cpgs)))
    betas = betas.loc[:, common_cpgs]
    return betas
