from src.models.fcmlp_model import FCMLPModel
from src.models.fcvae_model_v2 import FCVAEModelV2
from src.models.extractor_mlp_model import ExtractorFCMLPModel
import pickle
from src.datamodules.datasets.dnam_dataset import DNAmDataset
from torch.utils.data import DataLoader


def save_onnx():

    model_type = "ExtractorFCMLPModel"


    model_path = "E:/YandexDisk/Work/dnamvae/models/fcvae_mlp_model/logs/runs/2021-06-11/00-34-12"
    ckpt_name = "epoch=77_fold_2"
    ckpt_path = f"{model_path}/checkpoints/{ckpt_name}.ckpt"

    data_path = "E:/YandexDisk/Work/dnamvae/data/datasets/combo/GSE40279_GSE87571_EPIC_GSE55763/vt_score_more_0.005_none/data_nn.pkl"
    outcome = "age"

    if model_type == "FCMLPModel":
        model = FCMLPModel.load_from_checkpoint(checkpoint_path=ckpt_path)
    elif model_type == "FCVAEModelV2":
        model = FCVAEModelV2.load_from_checkpoint(checkpoint_path=ckpt_path)
    elif model_type == "ExtractorFCMLPModel":
        model = ExtractorFCMLPModel.load_from_checkpoint(checkpoint_path=ckpt_path)
    else:
        raise ValueError("Unsupported model type!")

    # print model hyperparameters
    print(model.hparams)

    print(model)

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
        batch_size=1,
        num_workers=0,
        pin_memory=False,
        shuffle=False,
    )

    x, y = next(iter(dataloader))

    model.to_onnx(f"{model_path}/checkpoints/{ckpt_name}.onnx", x, export_params=True)

if __name__ == "__main__":
    save_onnx()
