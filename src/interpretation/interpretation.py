from src.models.fcmlp_model import FCMLPModel
import pickle
from src.datamodules.datasets.dnam_dataset import DNAmDataset
from torch.utils.data import DataLoader
from tqdm import tqdm
import pandas as pd
import dotenv
from pytorch_lightning import LightningModule, LightningDataModule
from torch.utils.data import DataLoader, ConcatDataset, Subset
from pytorch_lightning import seed_everything
import hydra
from src.utils import utils
from omegaconf import DictConfig
import shap

log = utils.get_logger(__name__)

dotenv.load_dotenv(override=True)

@hydra.main(config_path="../../configs/", config_name="config.yaml")
def main(config: DictConfig):

    if "seed" in config:
        seed_everything(config.seed)

    model_path = "E:/YandexDisk/Work/dnamvae/models/fcmlp_model/logs/runs/2021-06-09/01-45-07"
    ckpt_path = f"{model_path}/checkpoints/epoch=223_fold_1.ckpt"

    model_type = "FCMLPModel"
    if model_type == "FCMLPModel":
        model = FCMLPModel.load_from_checkpoint(checkpoint_path=ckpt_path)
    else:
        raise ValueError("Unsupported model type!")

    # switch to evaluation mode
    model.eval()
    model.freeze()

    # Init Lightning datamodule
    log.info(f"Instantiating datamodule <{config.datamodule._target_}>")
    datamodule: LightningDataModule = hydra.utils.instantiate(config.datamodule)
    datamodule.setup()

    train_dataloader = datamodule.train_dataloader()
    val_dataloader = datamodule.val_dataloader()
    test_dataloader = datamodule.test_dataloader()

    batch = next(iter(train_dataloader))
    background, _ = batch

    e = shap.DeepExplainer(model, background)
    shap_values = e.shap_values(background)

    shap.summary_plot(shap_values, background)


if __name__ == "__main__":
    main()
