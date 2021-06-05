import pandas as pd
import os
import pickle
from sklearn.feature_selection import VarianceThreshold
import numpy as np
from skfeature.function.similarity_based import lap_score
from skfeature.function.similarity_based import SPEC
from skfeature.function.sparse_learning_based import UDFS
from skfeature.utility.sparse_learning import feature_ranking


th = 0.0

path = f"E:/YandexDisk/Work/dnamvae/data/datasets"
outcome = "age"
datasets = ["GSE40279", "GSE87571", "EPIC", "GSE55763"]

for d_id, d in enumerate(datasets):
    f = open(f"{path}/{d}/data_nn.pkl", 'rb')
    data = pickle.load(f)
    f.close()
    if d_id == 0:
        beta = data["beta"]
        pheno = data["pheno"].filter([outcome], axis=1)
    else:
        beta = beta.append(data["beta"])
        pheno = pheno.append(data["pheno"].filter([outcome], axis=1))

save_path = f"{path}/combo"
if not os.path.exists(save_path):
    os.makedirs(save_path)

vt = VarianceThreshold(th)
vt.fit(beta)
vt_metrics = vt.variances_
vt_idx = np.where(vt_metrics > th)[0]
vt_mask = vt_metrics > th
vt_bool = vt.get_support()

cpgs = beta.columns[vt_bool]
d = {'CpG': cpgs.to_list()}
d["vt_score"] = vt.variances_[vt_bool]

df = pd.DataFrame(d)
df.to_excel(f"{save_path}/{'_'.join(datasets)}.xlsx", index=False)

beta = beta[cpgs]

lap_metrics = lap_score.lap_score(beta.to_numpy())
lap_idx = lap_score.feature_ranking(lap_metrics)
d["lap_score"] = lap_metrics

df = pd.DataFrame(d)
df.to_excel(f"{save_path}/{'_'.join(datasets)}.xlsx", index=False)

spec_metrics = SPEC.spec(beta.to_numpy())
spec_idx = SPEC.feature_ranking(spec_metrics)
d["spec_metrics"] = spec_metrics

df = pd.DataFrame(d)
df.to_excel(f"{save_path}/{'_'.join(datasets)}.xlsx", index=False)

udfs_metrics = UDFS.udfs(beta.to_numpy())
udfs_idx = feature_ranking(udfs_metrics)
d["udfs_metrics"] = udfs_metrics

df = pd.DataFrame(d)
df.to_excel(f"{save_path}/{'_'.join(datasets)}.xlsx", index=False)




