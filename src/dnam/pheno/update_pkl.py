import pandas as pd


dataset = "GSE168739"
platform = "GPL21145"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

pheno = pd.read_excel(f"{path}/{platform}/{dataset}/pheno_xtd.xlsx", index_col="subject_id")
pheno.to_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")