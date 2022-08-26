import torch.nn.functional as F
import pytorch_lightning as pl
from scorevariants.models.common import predictions_to_score
import numpy as np

class LightningWrapper(pl.LightningModule):
    def __init__(self, model,
                 tmp_file=None,
                 **additional_hparams,
                ):
        """PyTorch Lighting wrapper for training and evaluation
        
        Args:
            model: nn.Module, pytorch model to train/evaluate
            additional_hparams: saved in hparams
        
        """
        super().__init__()
        self.model = model
        
        # save all hyperparameters
        # in addition to provenance, this helps
        # with model checkpointing by saving the model
        self.save_hyperparameters()

        if tmp_file:
            self.tmp_file = open(tmp_file, 'w')

    def forward(self, batch):
        outputs = self.model(batch)
        return outputs

    def on_test_end(self):
        self.tmp_file.close()
        
    def test_step(self, batch, batch_idx):
        predictions = self(batch)
        detached_predictions = predictions.detach()
        prob_predictions = F.softmax(detached_predictions)
        scores = predictions_to_score(prob_predictions.detach(), batch['type'])
        if self.tmp_file:
            for i in range(len(scores)):
                self.tmp_file.write('%s\t%d\t%s\t[%s]\t%.3f\n' % (batch['chrom'][i],
                                                                      batch['pos'][i],
                                                                      batch['ref'][i],
                                                                      batch['alt'][i],
                                                                      scores[i]))
