from torch import nn
import torch
from torch.nn import functional as F


class FCVAE(nn.Module):

    def __init__(self, hparams: dict):
        super().__init__()

        self.n_input = hparams['n_input']
        self.n_latent = hparams['n_latent']
        self.topology = hparams['topology']
        self.kl_coeff = hparams['kl_coeff']

        self.encoder_topology = [self.n_inputn_input] + self.topology
        self.decoder_topology = [ self.n_latent] + self.topology[::-1]

        self.encoder_layers = []
        for i in range(len(self.encoder_topology) - 1):
            layer = nn.Linear(self.encoder_topology[i], self.encoder_topology[i + 1])
            nn.init.xavier_uniform(layer.weight)
            self.encoder_layers.append(nn.Sequential(layer, nn.BatchNorm1d(self.encoder_topology[i + 1]), nn.ReLU()))
        self.encoder = nn.Sequential(*self.encoder_layers)

        self.hidden_mu = nn.Linear(self.encoder_topology[-1], self.n_latent)
        self.hidden_log_var = nn.Linear(self.encoder_topology[-1], self.n_latent)

        self.decoder_layers = []
        for i in range(len(self.decoder_topology) - 1):
            layer = nn.Linear(self.decoder_topology[i], self.decoder_topology[i + 1])
            nn.init.xavier_uniform(layer.weight)
            self.decoder_layers.append(nn.Sequential(layer, nn.BatchNorm1d(self.decoder_topology[i + 1]), nn.ReLU()))
        self.decoder_layers = nn.Sequential(*self.decoder_layers)
        self.output_layer = nn.Sequential(nn.Linear(self.decoder_topology[-1], self.n_input), nn.Sigmoid())
        self.decoder = nn.Sequential(*[self.decoder_layers, self.output_layer])

    def encode(self,x):
        hidden = self.encoder(x)
        mu = self.hidden_mu(hidden)
        log_var = self.hidden_log_var(hidden)
        return mu,log_var

    def forward(self, x):
        mu, log_var = self.encode(x)
        p, q, z = self.sample(mu, log_var)
        return self.decoder(z)

    def _run_step(self, x):
        mu, log_var = self.encode(x)
        p, q, z = self.sample(mu, log_var)
        return z, self.decoder(z), p, q

    def sample(self, mu, log_var):
        std = torch.exp(0.5 * log_var)
        p = torch.distributions.Normal(torch.zeros_like(mu), torch.ones_like(std))
        q = torch.distributions.Normal(mu, std)
        z = q.rsample()
        return p, q, z

    def step(self, batch):
        x, y = batch

        mu, log_var = self.encode(x)
        p, q, z = self.sample(mu, log_var)
        x_hat = self.decoder(z)

        recon_loss = F.mse_loss(x_hat, x, reduction='mean')

        log_qz = q.log_prob(z)
        log_pz = p.log_prob(z)

        kl = log_qz - log_pz
        kl_mean = kl.mean()
        kl_final = kl_mean * self.kl_coeff

        kl_vanilla_1 = (-0.5 * (1 + log_var - mu ** 2 - torch.exp(log_var)).sum(dim=1)).mean(dim=0)
        kl_vanilla_2 = torch.mean(-0.5 * torch.sum(1 + log_var - mu ** 2 - log_var.exp(), dim=1), dim=0)
        z_vanilla = self.alt_reparametrize(mu, log_var)
        x_hat_vanilla = self.decoder(z_vanilla)

        loss = kl_final + recon_loss

        logs = {
            "recon_loss": recon_loss,
            "kl_loss": kl_final,
            "loss": loss,
        }
        return loss, logs

    def training_step(self, batch, batch_idx):
        loss, logs = self.step(batch)
        self.log_dict({f"train_{k}": v for k, v in logs.items()})
        return loss

    def validation_step(self, batch, batch_idx):
        loss, logs = self.step(batch)
        self.log_dict({f"val_{k}": v for k, v in logs.items()})
        return loss

    def alt_reparametrize(self, mu, log_var):
        # std = torch.exp(0.5 * log_var)
        # eps = torch.randn(size=(mu.size(0), mu.size(1)))
        # eps = eps.type_as(mu)

        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std)

        return mu + std * eps
