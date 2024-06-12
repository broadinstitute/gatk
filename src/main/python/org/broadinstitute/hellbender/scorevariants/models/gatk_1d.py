#/usr/bin/python3

import torch
from torch import nn
import numpy as np

class GATK_CNN_1D(nn.Module):
    def __init__(self):
        """Architecture of PyTorch 1D CNN model
        """
        super(GATK_CNN_1D, self).__init__()
        self.conv = nn.Sequential(nn.Conv1d(in_channels=4, out_channels=256, kernel_size=12),
                                  nn.ReLU(),
                                  nn.Conv1d(in_channels=256, out_channels=256, kernel_size=12),
                                  nn.ReLU(),
                                  nn.Dropout(0.1),
                                  nn.Conv1d(in_channels=256, out_channels=128, kernel_size=12),
                                  nn.ReLU(),
                                  nn.Dropout(0.1))
        self.normalize = nn.BatchNorm1d(7, eps=0.001, momentum=0.99)
        self.dense_1 = nn.Sequential(nn.Linear(7, 40),
                                     nn.ReLU())
        self.dense_2 = nn.Sequential(nn.Linear(12200, 40),
                                     nn.ReLU(),
                                     nn.Dropout(0.2))
        self.dense_3 = nn.Linear(47, 4)
        self.loss_fn = nn.CrossEntropyLoss()

    def forward(self, batch):
        """Forwrd path
        Args:
            batch: input minibatch
        Returns
            output of the network after the last fully-connected layer
        """
        conv = self.conv(torch.transpose(batch['reference'], 1, 2))
        conv = torch.flatten(torch.transpose(conv, 1, 2), 1)
        normalized_annotations = self.normalize(batch['best_practices'])
        dense_1 = self.dense_1(normalized_annotations)
        cat_1 = torch.cat((conv, dense_1), -1)
        dense_2 = self.dense_2(cat_1)
        cat_2 = torch.cat((dense_2, normalized_annotations), -1)
        dense_3 = self.dense_3(cat_2)
        return dense_3

    def port_params(self, model):
        """Port parameters (weights, biases, etc.) from the Keras 1D model to PyTorch 1D model
        Args:
            model: Keras model
        Returns
            PyTorch model
        """
        weights = model.get_weights()
        self.conv[0].weight.data=torch.from_numpy(np.transpose(weights[0]))
        self.conv[0].bias.data=torch.from_numpy(np.transpose(weights[1]))
        self.conv[2].weight.data=torch.from_numpy(np.transpose(weights[2]))
        self.conv[2].bias.data=torch.from_numpy(np.transpose(weights[3]))
        self.conv[5].weight.data=torch.from_numpy(np.transpose(weights[4]))
        self.conv[5].bias.data=torch.from_numpy(np.transpose(weights[5]))
        self.normalize.weight.data=torch.from_numpy(weights[6])
        self.normalize.bias.data=torch.from_numpy(weights[7])
        self.normalize.running_mean.data = torch.from_numpy(weights[8])
        self.normalize.running_var.data = torch.from_numpy(weights[9])
        self.dense_1[0].weight.data=torch.from_numpy(np.transpose(weights[10]))
        self.dense_1[0].bias.data=torch.from_numpy(np.transpose(weights[11]))
        self.dense_2[0].weight.data=torch.from_numpy(np.transpose(weights[12]))
        self.dense_2[0].bias.data=torch.from_numpy(np.transpose(weights[13]))
        self.dense_3.weight.data=torch.from_numpy(np.transpose(weights[14]))
        self.dense_3.bias.data=torch.from_numpy(np.transpose(weights[15]))
        return self
