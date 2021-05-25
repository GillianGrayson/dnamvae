import numpy as np
import pandas as pd
from functools import reduce
import os
import pickle
import tqdm


path = f'E:/YandexDisk/Work/dnamvae/data'

f = open(f'{path}/annotation/cpgs.pkl', 'rb')
cpgs = pickle.load(f)
f.close()

f = open(f'{path}/datasets/GSE87571/data.pkl', 'rb')
raw_data = pickle.load(f)
f.close()

raw_data['beta'] = raw_data['beta'].astype(np.float32)

num_subjects = raw_data['beta'].shape[0]
num_features = len(cpgs)

data = np.zeros((num_subjects, num_features), dtype=np.float32)

X = pd.DataFrame(data=data, index=raw_data['beta'].index, columns=cpgs)
for cpg in tqdm.tqdm(list(raw_data['beta'].columns.values), mininterval=1.0, desc='DNAmDataset creating'):
    X[cpg] = raw_data['beta'][cpg]

raw_data['beta'] = X

f = open(f'{path}/datasets/GSE87571/data_nn.pkl', 'wb')
pickle.dump(raw_data, f, pickle.HIGHEST_PROTOCOL)
f.close()

