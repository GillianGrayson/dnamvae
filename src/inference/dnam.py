from src.models.fcmlp_model import FCMLPModel
import pickle
from src.datamodules.datasets.dnam_dataset import DNAmDataset
from torch.utils.data import DataLoader
from tqdm import tqdm
import pandas as pd


def inference():

    model_type = "FCMLPModel"

    model_path = "E:/YandexDisk/Work/dnamvae/models/fcmlp_model/logs/runs/2021-06-09/01-45-07"
    ckpt_path = f"{model_path}/checkpoints/epoch=223_fold_1.ckpt"

    data_path = "E:/YandexDisk/Work/dnamvae/data/datasets/combo/GSE40279_GSE87571_EPIC_GSE55763/vt_score_more_0.005_none/data_nn.pkl"
    outcome = "age"

    if model_type == "FCMLPModel":
        model = FCMLPModel.load_from_checkpoint(checkpoint_path=ckpt_path)
    else:
        raise ValueError("Unsupported model type!")

    # print model hyperparameters
    print(model.hparams)

    # switch to evaluation mode
    model.eval()
    model.freeze()

    # load data
    f = open(data_path, 'rb')
    data = pickle.load(f)
    f.close()

    dataset = DNAmDataset(data, outcome)
    dataloader = DataLoader(
        dataset=dataset,
        batch_size=5,
        num_workers=0,
        pin_memory=False,
        shuffle=False,
    )

    y_real = []
    y_pred = []

    for x, y in tqdm(dataloader):
        y_p = model(x)
        y_real.extend(y.tolist())
        y_pred.extend(y_p.flatten().tolist())

    d = {
        f"{outcome}_real": y_real,
        f"{outcome}_pred": y_pred
    }
    df = pd.DataFrame(d)
    df.to_excel(f'{model_path}/inference.xlsx', index=False)


if __name__ == "__main__":
    inference()
