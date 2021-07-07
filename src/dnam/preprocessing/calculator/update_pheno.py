import pandas as pd


dataset = "GSE80417"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

pheno = pd.read_pickle(f"{path}/{platform}/{dataset}/pheno.pkl")

calcs = pd.read_csv(f"{path}/{platform}/{dataset}/calculator/betas.output.csv", delimiter=",", index_col='subject_id')
calcs.drop(columns=["SampleID", "Age", "Female", "Tissue"])
features = ["DNAmAge", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "DNAmAgeHannum", "DNAmPhenoAge", "DNAmGDF15", "DNAmGrimAge", "IEAA", "EEAA",	"IEAA.Hannum"]
calcs = calcs[features]

df = pd.merge(pheno, calcs, left_index=True, right_index=True)
df.to_excel(f"{path}/{platform}/{dataset}/pheno_xtd.xlsx")
df.to_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")
