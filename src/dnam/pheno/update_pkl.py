import pandas as pd


dataset = "GSE147221"
platform = "GPL13534"
path = f"E:/YandexDisk/Work/pydnameth/datasets"

pheno = pd.read_excel(f"{path}/{platform}/{dataset}/pheno_xtd.xlsx", index_col="subject_id")
pheno.to_pickle(f"{path}/{platform}/{dataset}/pheno_xtd.pkl")