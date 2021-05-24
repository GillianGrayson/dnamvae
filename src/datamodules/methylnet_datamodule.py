from typing import Optional, Tuple
from .datasets.dnam_dataset import DNAmDataset
from pytorch_lightning import LightningDataModule
from torch.utils.data import DataLoader, Dataset
import numpy as np
import pickle


class MethylNetDataModule(LightningDataModule):

    def __init__(
            self,
            data_path: str = "E:/YandexDisk/Work/dnamvae/data/datasets/methylnet",
            outcome: str = 'Age',
            batch_size: int = 64,
            num_workers: int = 0,
            pin_memory: bool = False,
            **kwargs,
    ):
        super().__init__()

        self.data_path = data_path
        self.outcome = outcome
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.pin_memory = pin_memory

        self.data_train: Optional[Dataset] = None
        self.data_val: Optional[Dataset] = None
        self.data_test: Optional[Dataset] = None

    def prepare_data(self):
        """Download data if needed. This method is called only from a single GPU.
        Do not use it to assign state (self.x = y)."""
        pass

    def setup(self, stage: Optional[str] = None):
        """Load data. Set variables: self.data_train, self.data_val, self.data_test."""
        f = open(f"{self.data_path}/train_methyl_array.pkl", 'rb')
        tmp = pickle.load(f)
        tmp['beta'] =  tmp['beta'].astype(np.float32)
        self.train_methyl_array = DNAmDataset(tmp, self.outcome)
        f.close()

        f = open(f"{self.data_path}/val_methyl_array.pkl", 'rb')
        tmp = pickle.load(f)
        tmp['beta'] = tmp['beta'].astype(np.float32)
        self.val_methyl_array = DNAmDataset(tmp, self.outcome)
        f.close()

        f = open(f"{self.data_path}/test_methyl_array.pkl", 'rb')
        tmp = pickle.load(f)
        tmp['beta'] = tmp['beta'].astype(np.float32)
        self.test_methyl_array = DNAmDataset(tmp, self.outcome)
        f.close()

    def train_dataloader(self):
        return DataLoader(
            dataset=self.train_methyl_array,
            batch_size=self.batch_size,
            num_workers=self.num_workers,
            pin_memory=self.pin_memory,
            shuffle=True,
        )

    def val_dataloader(self):
        return DataLoader(
            dataset=self.val_methyl_array,
            batch_size=self.batch_size,
            num_workers=self.num_workers,
            pin_memory=self.pin_memory,
            shuffle=False,
        )

    def test_dataloader(self):
        return DataLoader(
            dataset=self.test_methyl_array,
            batch_size=self.batch_size,
            num_workers=self.num_workers,
            pin_memory=self.pin_memory,
            shuffle=False,
        )
