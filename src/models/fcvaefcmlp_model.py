from typing import Any, List
from pytorch_lightning import LightningModule
from torch.nn import functional as F
from src.models.fcvae_model_v1 import FCVAEModelV1
from torch import nn
import torch
from torchmetrics.classification.accuracy import Accuracy


class FCVAEFCMLPModel(LightningModule):
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
            fcvae_path: str = "",
            task: str = "regression",
            n_output: int = 1,
            topology: List[int] = None,
            dropout: float = 0.1,
            is_softmax: bool = False,
            kl_coeff: int = 1.0,
            lr: float = 0.001,
            weight_decay: float = 0.0005,
            **kwargs
    ):
        super().__init__()
        self.save_hyperparameters()

        self.feature_extractor = FCVAEModelV1.load_from_checkpoint(fcvae_path)
        self.feature_extractor.freeze()

        self.task = task
        self.n_output = n_output
        self.topology = [self.feature_extractor.model.n_latent] + list(topology)

        self.mlp_layers = []
        for i in range(len(self.topology) - 1):
            layer = nn.Linear(self.topology[i], self.topology[i + 1])
            self.mlp_layers.append(nn.Sequential(layer, nn.ReLU(), nn.Dropout(dropout)))
        self.mlp_layers.append(nn.Linear(self.topology[-1], self.n_output))

        if task == "classification":
            if n_output == 1:
                self.mlp_layers.append(nn.Sigmoid())
                self.loss_fn = torch.nn.BCEWithLogitsLoss(reduction='mean')
            elif n_output > 2:
                self.loss_fn = torch.nn.CrossEntropyLoss(reduction='mean')
            else:
                raise ValueError(f"Classification with {n_output} classes")
        elif task == "regression":
            self.loss_fn = torch.nn.MSELoss(reduction='mean')

        self.mlp = nn.Sequential(*self.mlp_layers)

        self.train_accuracy = Accuracy()
        self.val_accuracy = Accuracy()
        self.test_accuracy = Accuracy()

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
        loss = self.loss_fn(out, y)

        if self.task == "classification":
            if self.n_output == 1:
                out_tag = torch.round(torch.sigmoid(out))
            elif self.n_output > 2:
                out_tag = torch.argmax(out, dim=1)

        x_hat = self.model.decoder(z)

        recon_loss = F.mse_loss(x_hat, x, reduction='mean')

        log_qz = q.log_prob(z)
        log_pz = p.log_prob(z)

        kl = log_qz - log_pz
        kl_mean = kl.mean()
        kl_final = kl_mean * self.hparams.kl_coeff

        loss = kl_final + recon_loss

        logs = {
            "loss": loss,
            "recon_loss": recon_loss,
            "kl_loss": kl_final
        }
        return loss, logs

    def training_step(self, batch: Any, batch_idx: int):
        loss, logs = self.step(batch)
        d = {f"train_{k}": v for k, v in logs.items()}
        self.log_dict(d)

        # we can return here dict with any tensors
        # and then read it in some callback or in training_epoch_end() below
        # remember to always return loss from training_step, or else backpropagation will fail!
        return logs

    def training_epoch_end(self, outputs: List[Any]):
        # `outputs` is a list of dicts returned from `training_step()`
        pass

    def validation_step(self, batch: Any, batch_idx: int):
        loss, logs = self.step(batch)
        d = {f"val_{k}": v for k, v in logs.items()}
        self.log_dict(d)
        return logs

    def validation_epoch_end(self, outputs: List[Any]):
        pass

    def test_step(self, batch: Any, batch_idx: int):
        loss, logs = self.step(batch)
        d = {f"test_{k}": v for k, v in logs.items()}
        self.log_dict(d)
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
