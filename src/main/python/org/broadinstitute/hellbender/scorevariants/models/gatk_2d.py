#/usr/bin/python3

import torch
from torch import nn
import numpy as np

class GATK_CNN_2D(nn.Module):
    def __init__(self):
        """Architecture of PyTorch 2D CNN model
        """
        super(GATK_CNN_2D, self).__init__()
        self.conv = nn.Sequential(nn.Conv2d(in_channels=15, out_channels=64, kernel_size=(25,1)),
                                  nn.ReLU(),
                                  nn.Dropout(0.1),
                                  nn.Conv2d(in_channels=64, out_channels=48, kernel_size=(1,25)),
                                  nn.ReLU(),
                                  nn.Dropout2d(0.1),
                                  nn.Conv2d(in_channels=48, out_channels=32, kernel_size=(25,1)),
                                  nn.ReLU(),
                                  nn.Dropout2d(0.1),
                                  nn.MaxPool2d(kernel_size=(3,1), stride=(3,1)),
                                  nn.Conv2d(in_channels=32, out_channels=24, kernel_size=(1,25)),
                                  nn.ReLU(),
                                  nn.Dropout2d(0.1),
                                  nn.MaxPool2d(kernel_size=(3,1), stride=(3,1)))
        self.dense_1 = nn.Sequential(nn.BatchNorm1d(7, eps=0.001, momentum=0.99),
                                     nn.Linear(7, 64),
                                     nn.ReLU())
        self.dense_2 = nn.Sequential(nn.Linear(15424, 24),
                                     nn.ReLU(),
                                     nn.Dropout(0.3),
                                     nn.Linear(24, 4))
        self.loss_fn = nn.CrossEntropyLoss()

    def forward(self, batch):
        """Forwrd path
        Args:
            batch: input minibatch
        Returns
            output of the network after the last fully-connected layer
        """
        conv = self.conv(batch['read_tensor'].permute((0, 3, 1, 2)))
        conv = torch.flatten(conv.permute((0, 2, 3, 1)), 1)
        dense_1 = self.dense_1(batch['best_practices'])
        cat_1 = torch.cat((conv, dense_1), -1)
        dense_2 = self.dense_2(cat_1)
        return dense_2

    def port_params(self, model):
        """Port parameters (weights, biases, etc.) from the Keras 2D model to PyTorch 2D model
        Args:
            model: Keras model
        Returns
            PyTorch model
        """
        weights = model.get_weights()
        self.conv[0].weight.data=torch.from_numpy(np.transpose(weights[0],(3,2,0,1)))
        self.conv[0].bias.data=torch.from_numpy(np.transpose(weights[1]))
        self.conv[3].weight.data=torch.from_numpy(np.transpose(weights[2],(3,2,0,1)))
        self.conv[3].bias.data=torch.from_numpy(np.transpose(weights[3]))
        self.conv[6].weight.data=torch.from_numpy(np.transpose(weights[4],(3,2,0,1)))
        self.conv[6].bias.data=torch.from_numpy(np.transpose(weights[5]))
        self.conv[10].weight.data=torch.from_numpy(np.transpose(weights[6],(3,2,0,1)))
        self.conv[10].bias.data=torch.from_numpy(np.transpose(weights[7]))
        self.dense_1[0].weight.data=torch.from_numpy(weights[8])
        self.dense_1[0].bias.data=torch.from_numpy(weights[9])
        self.dense_1[0].running_mean.data=torch.from_numpy(weights[10])
        self.dense_1[0].running_var.data=torch.from_numpy(weights[11])
        self.dense_1[1].weight.data=torch.from_numpy(np.transpose(weights[12]))
        self.dense_1[1].bias.data=torch.from_numpy(np.transpose(weights[13]))
        self.dense_2[0].weight.data=torch.from_numpy(np.transpose(weights[14]))
        self.dense_2[0].bias.data=torch.from_numpy(np.transpose(weights[15]))
        self.dense_2[3].weight.data=torch.from_numpy(np.transpose(weights[16]))
        self.dense_2[3].bias.data=torch.from_numpy(np.transpose(weights[17]))
        return self
