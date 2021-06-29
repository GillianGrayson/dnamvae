from typing import Any, List
import torch
from pytorch_lightning import LightningModule
from src.models.modules.fcvae_net import FCVAENet


class FCVAEModelV2(LightningModule):
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
            loss_type: str = "MSE",
            kl_coeff: int = 256,
            lr: float = 0.001,
            weight_decay: float = 0.0005,
            **kwargs
    ):
        super().__init__()

        # this line ensures params passed to LightningModule will be saved to ckpt
        # it also allows to access params with 'self.hparams' attribute
        self.save_hyperparameters()

        if self.hparams.loss_type == "MSE":
            self.loss_fn = torch.nn.MSELoss(reduction='mean')
        elif self.hparams.loss_type == "BCE":
            self.loss_fn = torch.nn.BCELoss(reduction='mean')
        elif self.hparams.loss_type == "L1Loss":
            self.loss_fn = torch.nn.L1Loss(reduction='mean')
        else:
            raise ValueError("Unsupported loss_type")

        self.model = FCVAENet(hparams=self.hparams)

    def forward(self, x: torch.Tensor):
        return self.model.forward_v2(x)

    def get_latent(self, x: torch.Tensor):
        mu, log_var = self.model.encode(x)
        z = self.model.v2_reparametrize(mu, log_var)
        return z

    def step(self, batch: Any):
        x, y, ind = batch

        mu, log_var = self.model.encode(x)

        kl_final = (-0.5 * (1 + log_var - mu ** 2 - torch.exp(log_var)).sum(dim=1)).mean(dim=0)

        z = self.model.v2_reparametrize(mu, log_var)
        x_hat = self.model.decoder(z)

        recon_loss = self.loss_fn(x_hat, x)

        loss = kl_final + recon_loss

        logs = {
            "loss": loss,
            "recon_loss": recon_loss,
            "kl_loss": kl_final
        }
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
        d = {f"test/{k}": v for k, v in logs.items()}
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
