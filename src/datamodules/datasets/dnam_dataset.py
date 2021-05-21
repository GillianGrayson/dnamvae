import torch
import numpy as np
from torch.utils.data import Dataset
import pandas as pd
import tqdm


class DNAmDataset(Dataset):

    def __init__(
            self,
            data: dict,
            outcome: str = 'Age'
    ):
        self.data = data
        self.outcome = outcome
        self.num_subjects = self.data['beta'].shape[0]
        self.num_features = self.data['beta'].shape[1]

        self.X = data['beta']
        self.y = self.data['pheno'].loc[:, self.outcome]

    def __getitem__(self, idx: int):
        x = self.X.iloc[idx].to_numpy()
        y = self.y[idx]

        return (x, y)

    def __len__(self):
        return self.num_subjects
