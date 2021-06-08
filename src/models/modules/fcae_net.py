from torch import nn
import torch


class FCAENet(nn.Module):

    def __init__(self, hparams: dict):
        super().__init__()

        self.n_input = hparams['n_input']
        self.n_latent = hparams['n_latent']
        self.topology = list(hparams['topology'])
        self.loss_type = hparams['loss_type']

        self.encoder_topology = [self.n_input] + self.topology
        self.decoder_topology = [ self.n_latent] + self.topology[::-1]

        self.encoder_layers = []
        for i in range(len(self.encoder_topology) - 1):
            layer = nn.Linear(self.encoder_topology[i], self.encoder_topology[i + 1])
            self.encoder_layers.append(nn.Sequential(layer, nn.BatchNorm1d(self.encoder_topology[i + 1]), nn.ReLU()))
        self.encoder_layers.append(nn.Linear(self.encoder_topology[-1], self.n_latent))
        self.encoder = nn.Sequential(*self.encoder_layers)

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
        x = self.encoder(x)
        return x

    def forward(self, x):
        x = self.encode(x)
        x = self.decoder(x)
        return x
