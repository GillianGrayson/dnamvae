import pandas as pd
import pickle
from src.dnam.preprocessing.routines.filter import get_forbidden_cpgs


dataset = "GSE84727"
array_type = "450K"
path = f"E:/YandexDisk/Work/pydnameth/datasets"
forbidden_types = ["NoCG", "SNP", "MultiHit", "XY"]

fn = f"{path}/{array_type}/{dataset}/pheno.xlsx"
df = pd.read_excel(fn)
pheno = df.set_index('sentrixids')
pheno.index.name = "subject_id"

fn = f"{path}/{array_type}/{dataset}/raw/GSE84727_normalisedBetas.csv"
df = pd.read_csv(fn, delimiter=",")
df.rename(columns={df.columns[0]: 'CpG'}, inplace=True)
df.set_index('CpG', inplace=True)
betas = df.T
betas.index.name = "subject_id"
forbidden_cpgs = get_forbidden_cpgs(f"{path}/{array_type}/forbidden_cpgs", forbidden_types)
betas = betas.loc[:, ~betas.columns.isin(forbidden_cpgs)]

pheno_subj_ids = pheno.index.values.tolist()
betas_subj_ids = betas.index.values.tolist()
is_equal_ids = pheno_subj_ids == betas_subj_ids
print(f"Is equal ids: {is_equal_ids}")
if not is_equal_ids:
    intersection = list(set(pheno_subj_ids).intersection(betas_subj_ids))
    pheno = pheno.loc[intersection, :]
    betas = betas.loc[intersection, :]

d = {'betas': betas, 'pheno': pheno}
f = open(f'{path}/{array_type}/{dataset}/data.pkl', 'wb')
pickle.dump(d, f, pickle.HIGHEST_PROTOCOL)
f.close()
