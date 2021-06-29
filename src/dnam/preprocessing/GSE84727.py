import pandas as pd
import pickle


dataset = "GSE84727"
path = f"E:/YandexDisk/Work/pydnameth/datasets/450K/{dataset}"

fn = f"{path}/observables.xlsx"
df = pd.read_excel(fn)
pheno = df.set_index('sentrixids')
pheno.index.name = "subject_id"

fn = f"{path}/raw/GSE84727_normalisedBetas.csv"
df = pd.read_csv(fn, delimiter=",")
df.rename(columns={df.columns[0]: 'CpG'}, inplace=True)
df.set_index('CpG', inplace=True)
betas = df.T
betas.index.name = "subject_id"

pheno_subj_ids = pheno.index.values.tolist()
betas_subj_ids = betas.index.values.tolist()
is_equal_ids = pheno_subj_ids == betas_subj_ids
print(f"Is equal ids: {is_equal_ids}")
if not is_equal_ids:
    intersection = list(set(pheno_subj_ids).intersection(betas_subj_ids))
    pheno = pheno.loc[intersection, :]
    betas = betas.loc[intersection, :]

d = {'beta': betas, 'pheno': pheno}

f = open(f'{path}/data.pkl', 'wb')
pickle.dump(d, f, pickle.HIGHEST_PROTOCOL)
f.close()
