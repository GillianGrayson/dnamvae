from typing import Optional, Tuple
from .datasets.dnam_dataset import DNAmDataset, DNAmDatasetIndexRetrieve
from pytorch_lightning import LightningDataModule
from torch.utils.data import DataLoader, Dataset, random_split
import pickle
import numpy as np
from src.utils import utils

log = utils.get_logger(__name__)

class DNAmDataModule(LightningDataModule):

    def __init__(
            self,
            data_fn: str = "E:/YandexDisk/Work/dnamvae/data/datasets/unn/data_nn.pkl",
            outcome: str = 'Age',
            train_val_test_split: Tuple[float, float, float] = (0.7, 0.2, 0.1),
            batch_size: int = 64,
            num_workers: int = 0,
            pin_memory: bool = False,
            **kwargs,
    ):
        super().__init__()

        self.data_fn = data_fn
        self.outcome = outcome
        self.train_val_test_split = train_val_test_split
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
        f = open(self.data_fn, 'rb')
        self.data = pickle.load(f)
        f.close()

        # self.dims is returned when you call datamodule.size()
        self.dims = (1, self.data['beta'].shape[1])

        dataset = DNAmDataset(self.data, self.outcome)

        if self.train_val_test_split[2] == 0:
            total_count = self.data['beta'].shape[0]
            train_count = int(np.floor(total_count * self.train_val_test_split[0]))
            valid_count = total_count - train_count
            test_count = 0
        else:
            total_count = self.data['beta'].shape[0]
            train_count = int(np.floor(total_count * self.train_val_test_split[0]))
            valid_count = int(np.floor(total_count * self.train_val_test_split[1]))
            test_count = total_count - train_count - valid_count

        log.info(f"total_count: {total_count}")
        log.info(f"train_count: {train_count}")
        log.info(f"valid_count: {valid_count}")
        log.info(f"test_count: {test_count}")

        self.data_train, self.data_val, self.data_test = random_split(
            dataset, [train_count, valid_count, test_count]
        )

    def train_dataloader(self):
        return DataLoader(
            dataset=self.data_train,
            batch_size=self.batch_size,
            num_workers=self.num_workers,
            pin_memory=self.pin_memory,
            shuffle=True,
        )

    def val_dataloader(self):
        return DataLoader(
            dataset=self.data_val,
            batch_size=self.batch_size,
            num_workers=self.num_workers,
            pin_memory=self.pin_memory,
            shuffle=False,
        )

    def test_dataloader(self):
        return DataLoader(
            dataset=self.data_test,
            batch_size=self.batch_size,
            num_workers=self.num_workers,
            pin_memory=self.pin_memory,
            shuffle=False,
        )
