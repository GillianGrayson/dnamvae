"""
    Cross validation for Pytorch Lightning Data Modules
"""

import os
from abc import abstractmethod, ABC
from typing import Tuple

import pytorch_lightning as pl
from sklearn.model_selection import KFold
from torch.utils.data import DataLoader, ConcatDataset, Subset


class CVDataModule(ABC):

    def __init__(self,
                 data_module: pl.LightningDataModule,
                 n_splits: int = 10,
                 shuffle: bool = True):
        self.data_module = data_module
        self._n_splits = n_splits
        self._shuffle = shuffle

    @abstractmethod
    def split(self):
        pass


class KFoldCVDataModule(CVDataModule):
    """
        K-fold cross-validation data module

    Args:
        data_module: data module containing data to be split
        n_splits: number of k-fold iterations/data splits
    """

    def __init__(self,
                 data_module: pl.LightningDataModule,
                 n_splits: int = 10):
        super().__init__(data_module, n_splits)
        self._k_fold = KFold(n_splits=self._n_splits, shuffle=self._shuffle)

    def get_data(self):
        """
            Extract and concatenate training and validation datasets from data module.
        """
        self.data_module.setup()
        train_ds = self.data_module.train_dataloader().dataset
        val_ds = self.data_module.val_dataloader().dataset
        return ConcatDataset([train_ds, val_ds])

    def get_test_dataloader(self):
        return self.data_module.test_dataloader()

    def split(self) -> Tuple[DataLoader, DataLoader]:
        """
            Split data into k-folds and yield each pair
        """
        # 0. Get data to split
        data = self.get_data()

        # 1. Iterate through splits
        for train_idx, val_idx in self._k_fold.split(range(len(data))):

            train_dl = DataLoader(Subset(data, train_idx),
                                  batch_size=self.data_module.batch_size,
                                  num_workers=self.data_module.num_workers,
                                  pin_memory=self.data_module.pin_memory,
                                  shuffle=self._shuffle)
            val_dl = DataLoader(Subset(data, val_idx),
                                batch_size=self.data_module.batch_size,
                                num_workers=self.data_module.num_workers,
                                pin_memory=self.data_module.pin_memory,
                                shuffle=self._shuffle)

            yield train_dl, val_dl