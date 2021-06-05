import torch
import numpy as np
from torch.utils.data import Dataset
import pandas as pd
import tqdm
from pandas.api.types import is_string_dtype
from pandas.api.types import is_numeric_dtype


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
        if is_numeric_dtype(self.y):
            self.y = self.y.astype('float32')

    def __getitem__(self, idx: int):
        x = self.X.iloc[idx].to_numpy()
        y = self.y[idx]

        return (x, y)

    def __len__(self):
        return self.num_subjects
