from typing import Any, List
from pytorch_lightning import LightningModule
from src.models.fcvae_model_v1 import FCVAEModelV1
from src.models.fcvae_model_v2 import FCVAEModelV2
from src.models.fcae_model import FCAEModel
from torch import nn
import torch
from torchmetrics.classification.accuracy import Accuracy


class ExtractorFCMLPModel(LightningModule):
    """
    A LightningModule organizes your PyTorch code into 5 sections:
        - Computations (init).
        - Train loop (training_step)
        - Validation loop (validation_step)
        - Test loop (test_step)
        - Optimizers (configure_optimizers)

    Read the docs:
        https://pytorch-lightning.readthedocs.io/en/latest/common/lightning_module.html
    """

    def __init__(
            self,
            extractor_path: str = "",
            task: str = "regression",
            n_output: int = 1,
            topology: List[int] = None,
            dropout: float = 0.1,
            num_unfreeze_epochs = 10,
            loss_type: str = "MSE",
            extractor_type: str = "FCVAEModelV2",
            lr: float = 0.001,
            weight_decay: float = 0.0005,
            **kwargs
    ):
        super().__init__()
        self.save_hyperparameters()

        self.extractor_type = extractor_type
        if self.extractor_type == "FCVAEModelV1":
            self.feature_extractor = FCVAEModelV1.load_from_checkpoint(extractor_path)
        elif self.extractor_type == "FCVAEModelV2":
            self.feature_extractor = FCVAEModelV2.load_from_checkpoint(extractor_path)
        elif self.extractor_type == "FCAEModel":
            self.feature_extractor = FCAEModel.load_from_checkpoint(extractor_path)
        else:
            raise ValueError("Unsupported extractor_type")

        self.feature_extractor.freeze()

        self.task = task
        self.n_output = n_output
        self.topology = [self.feature_extractor.model.n_latent] + list(topology)

        self.num_unfreeze_epochs = num_unfreeze_epochs

        self.mlp_layers = []
        for i in range(len(self.topology) - 1):
            layer = nn.Linear(self.topology[i], self.topology[i + 1])
            self.mlp_layers.append(nn.Sequential(layer, nn.ReLU(), nn.BatchNorm1d(self.topology[i + 1]), nn.Dropout(dropout)))
        self.mlp_layers.append(nn.Linear(self.topology[-1], self.n_output))

        if task == "classification":
            self.loss_fn = torch.nn.CrossEntropyLoss(reduction='mean')
            if n_output < 2:
                raise ValueError(f"Classification with {n_output} classes")
        elif task == "regression":
            if self.hparams.loss_type == "MSE":
                self.loss_fn = torch.nn.MSELoss(reduction='mean')
            elif self.hparams.loss_type == "L1Loss":
                self.loss_fn = torch.nn.L1Loss(reduction='mean')
            else:
                raise ValueError("Unsupported loss_type")

        self.mlp = nn.Sequential(*self.mlp_layers)

        self.accuracy = Accuracy()

    def on_epoch_end(self):
        # a hook is cleaner (but a callback is much better)
        if self.trainer.current_epoch == self.num_unfreeze_epochs:
            self.feature_extractor.unfreeze()

    def forward(self, x: torch.Tensor):
        z = self.feature_extractor.get_latent(x)
        return self.mlp(z)

    def get_probabilities(self, x: torch.Tensor):
        x = self.feature_extractor.get_latent(x)
        x = self.mlp(x)
        return torch.softmax(x, dim=1)

    def step(self, batch: Any):
        x, y = batch
        out = self.forward(x)
        batch_size = x.size(0)
        y = y.view(batch_size, -1)
        loss = self.loss_fn(out, y)

        logs = {"loss": loss}
        if self.task == "classification":
            out_tag = torch.argmax(out, dim=1)
            acc = self.accuracy(out_tag, y)
            logs["acc"] = acc

        return loss, logs

    def training_step(self, batch: Any, batch_idx: int):
        loss, logs = self.step(batch)
        d = {f"train/{k}": v for k, v in logs.items()}
        self.log_dict(d, on_step=False, on_epoch=True, logger=True)
        return logs

    def training_epoch_end(self, outputs: List[Any]):
        pass

    def validation_step(self, batch: Any, batch_idx: int):
        loss, logs = self.step(batch)
        d = {f"val/{k}": v for k, v in logs.items()}
        self.log_dict(d, on_step=False, on_epoch=True, logger=True)
        return logs

    def validation_epoch_end(self, outputs: List[Any]):
        pass

    def test_step(self, batch: Any, batch_idx: int):
        loss, logs = self.step(batch)
        d = {f"test_{k}": v for k, v in logs.items()}
        self.log_dict(d, on_step=False, on_epoch=True, logger=True)
        return logs

    def test_epoch_end(self, outputs: List[Any]):
        pass

    def configure_optimizers(self):
        """Choose what optimizers and learning-rate schedulers to use in your optimization.
        Normally you'd need one. But in the case of GANs or similar you might have multiple.

        See examples here:
            https://pytorch-lightning.readthedocs.io/en/latest/common/lightning_module.html#configure-optimizers
        """
        return torch.optim.Adam(
            params=self.parameters(), lr=self.hparams.lr, weight_decay=self.hparams.weight_decay
        )
