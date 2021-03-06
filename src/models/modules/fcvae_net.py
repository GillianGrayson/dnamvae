from torch import nn
import torch


class FCVAENet(nn.Module):

    def __init__(self, hparams: dict):
        super().__init__()

        self.n_input = hparams['n_input']
        self.n_latent = hparams['n_latent']
        self.topology = list(hparams['topology'])
        self.kl_coeff = hparams['kl_coeff']
        self.loss_type = hparams['loss_type']

        self.encoder_topology = [self.n_input] + self.topology
        self.decoder_topology = [ self.n_latent] + self.topology[::-1]

        self.encoder_layers = []
        for i in range(len(self.encoder_topology) - 1):
            layer = nn.Linear(self.encoder_topology[i], self.encoder_topology[i + 1])
            self.encoder_layers.append(nn.Sequential(layer, nn.BatchNorm1d(self.encoder_topology[i + 1]), nn.ReLU()))
        self.encoder = nn.Sequential(*self.encoder_layers)

        self.hidden_mu = nn.Linear(self.encoder_topology[-1], self.n_latent)
        self.hidden_log_var = nn.Linear(self.encoder_topology[-1], self.n_latent)

        self.decoder_layers = []
        for i in range(len(self.decoder_topology) - 1):
            layer = nn.Linear(self.decoder_topology[i], self.decoder_topology[i + 1])
            self.decoder_layers.append(nn.Sequential(layer, nn.BatchNorm1d(self.decoder_topology[i + 1]), nn.ReLU()))
        self.decoder_layers = nn.Sequential(*self.decoder_layers)
        if self.loss_type == "BCE":
            self.output_layer = nn.Sequential(nn.Linear(self.decoder_topology[-1], self.n_input), nn.Sigmoid())
        elif self.loss_type == "MSE" or self.loss_type == "L1Loss":
            self.output_layer = nn.Linear(self.decoder_topology[-1], self.n_input)
        else:
            raise ValueError("Unsupported loss_type")
        self.decoder = nn.Sequential(*[self.decoder_layers, self.output_layer])

    def encode(self,x):
        hidden = self.encoder(x)
        mu = self.hidden_mu(hidden)
        log_var = self.hidden_log_var(hidden)
        return mu,log_var

    def forward_v1(self, x):
        mu, log_var = self.encode(x)
        p, q, z = self.v1_reparametrize(mu, log_var)
        return self.decoder(z)

    def forward_v2(self, x):
        mu, log_var = self.encode(x)
        z = self.v2_reparametrize(mu, log_var)
        return self.decoder(z)

    def v1_reparametrize(self, mu, log_var):
        std = torch.exp(0.5 * log_var)
        p = torch.distributions.Normal(torch.zeros_like(mu), torch.ones_like(std))
        q = torch.distributions.Normal(mu, std)
        z = q.rsample()
        return p, q, z

    def v2_reparametrize(self, mu, log_var):
        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std)
        return mu + std * eps
