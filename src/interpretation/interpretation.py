import numpy as np
from src.models.fcmlp_model import FCMLPModel
import os
import matplotlib.pyplot as plt
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

    num_top_features = 100

    if "seed" in config:
        seed_everything(config.seed)

    save_dir = os.path.dirname(config.datamodule.data_fn)

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

    train_dataloader = datamodule.train_dataloader().dataset
    val_dataloader = datamodule.val_dataloader().dataset
    test_dataloader = datamodule.test_dataloader().dataset
    dataset = ConcatDataset([train_dataloader, val_dataloader, test_dataloader])
    dataloader = DataLoader(
        dataset,
        batch_size=config.datamodule.batch_size,
        num_workers=config.datamodule.num_workers,
        pin_memory=config.datamodule.pin_memory,
        shuffle=True
    )

    batch_id = 0
    d = {}
    for background, _, indexes in tqdm(dataloader):

        outs = model(background).flatten()

        if batch_id == 0:
            e = shap.DeepExplainer(model, background)

        shap_values = e.shap_values(background)

        if batch_id == 0:
            shap_abs = np.absolute(shap_values)
            shap_mean_abs = np.mean(shap_abs, axis=0)
            order = np.argsort(shap_mean_abs)[::-1]
            features = datamodule.data["beta"].columns.values
            features_best = features[order[0:num_top_features]]

        subject_indices = indexes.flatten().cpu().detach().numpy()
        subjects = datamodule.data["beta"].index.values[subject_indices]
        outcomes = datamodule.data["pheno"].loc[subjects, config.datamodule.outcome].to_numpy()

        betas = background.cpu().detach().numpy()
        preds = outs.cpu().detach().numpy()

        if batch_id == 0:
            d['subject'] = subjects
            d['outcome'] = outcomes
            d['preds'] = preds

            for f_id in range(0, num_top_features):
                feat = features_best[f_id]
                curr_beta = betas[:, order[f_id]]
                curr_shap = shap_values[:, order[f_id]]
                d[f"{feat}_beta"] = curr_beta
                d[f"{feat}_shap"] = curr_shap
        else:
            d['subject'] = np.append(d['subject'], subjects)
            d['outcome'] = np.append(d['outcome'], outcomes)
            d['preds'] = np.append(d['preds'], preds)

            for f_id in range(0, num_top_features):
                feat = features_best[f_id]
                curr_beta = betas[:, order[f_id]]
                curr_shap = shap_values[:, order[f_id]]
                d[f"{feat}_beta"] = np.append(d[f"{feat}_beta"], curr_beta)
                d[f"{feat}_shap"] = np.append(d[f"{feat}_shap"], curr_shap)

        batch_id += 1

    df_features = pd.DataFrame(d)
    df_features.to_excel(f"{save_dir}/interpretation/shap_values_{config.datamodule.batch_size}_{num_top_features}.xlsx", index=False)


if __name__ == "__main__":
    main()
