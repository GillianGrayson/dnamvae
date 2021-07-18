import pandas as pd
from src.dnam.python.preprocessing.serialization.routines.filter import get_forbidden_cpgs, manifest_filter, betas_pvals_filter
from src.dnam.python.preprocessing.serialization.routines.pheno_betas_checking import get_pheno_betas_with_common_subjects
from src.dnam.python.preprocessing.serialization.routines.save import save_pheno_betas_to_pkl
from src.dnam.python.preprocessing.serialization.routines.download import download_betas_and_pvals_from_gsms
from src.dnam.python.routines.manifest import get_manifest


dataset = "GSE144858"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"
forbidden_types = ["NoCG", "SNP", "MultiHit", "XY"]

manifest = get_manifest(platform)

fn = f"{path}/{platform}/{dataset}/pheno.xlsx"
df = pd.read_excel(fn)
pheno = df.set_index('geo_accession')
pheno.index.name = "subject_id"

betas, pvals = download_betas_and_pvals_from_gsms(pheno.index.values, f"{path}/{platform}/{dataset}/raw")
betas = betas_pvals_filter(betas, pvals, 0.01, 0.1)
betas = manifest_filter(betas, manifest)
forbidden_cpgs = get_forbidden_cpgs(f"{path}/{platform}/manifest/forbidden_cpgs", forbidden_types)
betas = betas.loc[:, ~betas.columns.isin(forbidden_cpgs)]

pheno, betas = get_pheno_betas_with_common_subjects(pheno, betas)
save_pheno_betas_to_pkl(pheno, betas, f"{path}/{platform}/{dataset}")
