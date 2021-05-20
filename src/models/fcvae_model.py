from typing import Any, List
import torch
from pytorch_lightning import LightningModule
from src.models.modules.fcvae_net import FCVAENet
from torch.nn import functional as F


class FCVAEModel(LightningModule):
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
        n_input: int = 784,
        n_latent: int = 256,
        topology: int = 256,
        kl_coeff: int = 256,
        lr: float = 0.001,
        weight_decay: float = 0.0005,
        **kwargs
    ):
        super().__init__()

        # this line ensures params passed to LightningModule will be saved to ckpt
        # it also allows to access params with 'self.hparams' attribute
        self.save_hyperparameters()

        self.model = FCVAENet(hparams=self.hparams)

    def forward(self, x: torch.Tensor):
        return self.model(x)

    def step(self, batch: Any):
        x, y = batch

        mu, log_var = self.model.encode(x)
        p, q, z =  self.model.sample(mu, log_var)
        x_hat =  self.model.decoder(z)

        recon_loss = F.mse_loss(x_hat, x, reduction='mean')

        log_qz = q.log_prob(z)
        log_pz = p.log_prob(z)

        kl = log_qz - log_pz
        kl_mean = kl.mean()
        kl_final = kl_mean * self.kl_coeff

        kl_vanilla_1 = (-0.5 * (1 + log_var - mu ** 2 - torch.exp(log_var)).sum(dim=1)).mean(dim=0)
        kl_vanilla_2 = torch.mean(-0.5 * torch.sum(1 + log_var - mu ** 2 - log_var.exp(), dim=1), dim=0)
        z_vanilla = self.model.alt_reparametrize(mu, log_var)
        x_hat_vanilla = self.model.decoder(z_vanilla)

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
