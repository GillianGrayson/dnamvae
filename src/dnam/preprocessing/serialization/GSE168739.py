import pandas as pd
from src.dnam.preprocessing.serialization.routines.filter import get_forbidden_cpgs, manifest_filter
from src.dnam.preprocessing.serialization.routines.pheno_betas_checking import get_pheno_betas_with_common_subjects
from src.dnam.preprocessing.serialization.routines.save import save_pheno_betas_to_pkl
from src.dnam.routines.manifest import get_manifest


dataset = "GSE168739"
platform = "GPL21145"
path = f"E:/YandexDisk/Work/pydnameth/datasets"
forbidden_types = ["NoCG", "SNP", "MultiHit", "XY"]

manifest = get_manifest(platform)

fn = f"{path}/{platform}/{dataset}/pheno.xlsx"
df = pd.read_excel(fn)
pheno = df.set_index('Sample_Name')
pheno.index.name = "subject_id"

fn = f"{path}/{platform}/{dataset}/raw/result/part(all)_config(0.01_0.10_0.10)/beta_funnorm_filtered.txt"
df = pd.read_csv(fn, delimiter="\t", index_col='CpG')
betas = df.T
betas.index.name = "subject_id"
betas = manifest_filter(betas, manifest)
forbidden_cpgs = get_forbidden_cpgs(f"{path}/{platform}/manifest/forbidden_cpgs", forbidden_types)
betas = betas.loc[:, ~betas.columns.isin(forbidden_cpgs)]

pheno, betas = get_pheno_betas_with_common_subjects(pheno, betas)
save_pheno_betas_to_pkl(pheno, betas, f"{path}/{platform}/{dataset}")
