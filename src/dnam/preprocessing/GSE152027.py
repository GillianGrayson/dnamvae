import pandas as pd
from src.dnam.preprocessing.routines.filter import get_forbidden_cpgs, betas_pvals_filter
from src.dnam.preprocessing.routines.pheno_betas_checking import get_pheno_betas_with_common_subjects
from src.dnam.preprocessing.routines.save import save_pheno_betas_to_pkl


dataset = "GSE152027"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"
forbidden_types = ["NoCG", "SNP", "MultiHit", "XY"]

fn = f"{path}/{platform}/{dataset}/pheno.xlsx"
df = pd.read_excel(fn)
df[['subject_id', 'Sample']] = df['title'].str.split(' ',expand=True)
pheno = df.set_index('subject_id')
pheno.index.name = "subject_id"

fn = f"{path}/{platform}/{dataset}/raw/GSE152027_IOP_processed_signals.csv"
df = pd.read_csv(fn, delimiter=",")
df.rename(columns={df.columns[0]: 'CpG'}, inplace=True)
df.set_index('CpG', inplace=True)
betas = df.iloc[:, 0::2]
pvals = df.iloc[:, 1::2]
betas = betas.T
pvals = pvals.T
pvals.index = betas.index.values.tolist()
betas.index.name = "subject_id"
pvals.index.name = "subject_id"
betas = betas_pvals_filter(betas, pvals, 0.01, 0.1)
forbidden_cpgs = get_forbidden_cpgs(f"{path}/{platform}/manifest/forbidden_cpgs", forbidden_types)
betas = betas.loc[:, ~betas.columns.isin(forbidden_cpgs)]

pheno, betas = get_pheno_betas_with_common_subjects(pheno, betas)
save_pheno_betas_to_pkl(pheno, betas, f"{path}/{platform}/{dataset}")
