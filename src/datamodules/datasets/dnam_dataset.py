import torch
import numpy as np
from torch.utils.data import Dataset
import pandas as pd
import tqdm


class DNAmDataset(Dataset):

    def __init__(
            self,
            raw_data: dict,
            cpgs: list,
            outcome: str = 'Age'
    ):
        self.raw_data = raw_data
        self.cpgs = cpgs
        self.outcome = outcome
        self.num_subjects = self.raw_data['beta'].shape[0]
        self.num_features = len(self.cpgs)

        data = np.zeros((self.num_subjects, self.num_features), dtype=np.float32)

        self.X = pd.DataFrame(data=data, index=self.raw_data['beta'].index, columns=self.cpgs)
        for cpg in tqdm.tqdm(list(self.raw_data['beta'].columns.values), mininterval=1.0, desc='DNAmDataset creating'):
            self.X[cpg] = self.raw_data['beta'][cpg]

        self.y = self.raw_data['pheno'].loc[:, self.outcome]

    def __getitem__(self, idx: int):
        x = self.X.iloc[idx].to_numpy()
        y = self.y[idx]

        return (x, y)

    def __len__(self):
        return self.num_subjects
