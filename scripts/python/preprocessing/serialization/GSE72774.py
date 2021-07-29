import pandas as pd
from src.dnam.python.preprocessing.serialization.routines.filter import get_forbidden_cpgs, manifest_filter
from src.dnam.python.preprocessing.serialization.routines.pheno_betas_checking import get_pheno_betas_with_common_subjects
from src.dnam.python.preprocessing.serialization.routines.save import save_pheno_betas_to_pkl
from src.dnam.python.routines.manifest import get_manifest


dataset = "GSE72774"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"
forbidden_types = ["NoCG", "SNP", "MultiHit", "XY"]

manifest = get_manifest(platform)

fn = f"{path}/{platform}/{dataset}/pheno.xlsx"
df = pd.read_excel(fn)
df['subject_id'] = df['title']
pheno = df.set_index('subject_id')
pheno.index.name = "subject_id"

fn = f"{path}/{platform}/{dataset}/raw/GSE72774_datBetaNormalized.csv"
df = pd.read_csv(fn, delimiter=",")
df.rename(columns={df.columns[0]: 'CpG'}, inplace=True)
df.set_index('CpG', inplace=True)
betas = df.T
betas.index.name = "subject_id"
betas.dropna(axis='columns', how='any', inplace=True)
betas = manifest_filter(betas, manifest)
forbidden_cpgs = get_forbidden_cpgs(f"{path}/{platform}/manifest/forbidden_cpgs", forbidden_types)
betas = betas.loc[:, ~betas.columns.isin(forbidden_cpgs)]

pheno, betas = get_pheno_betas_with_common_subjects(pheno, betas)
save_pheno_betas_to_pkl(pheno, betas, f"{path}/{platform}/{dataset}")
