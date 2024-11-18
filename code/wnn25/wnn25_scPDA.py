import random
import torch
import numpy as np
import pandas as pd
from scPDA import model
import h5py

torch.manual_seed(1)
random.seed(1)
np.random.seed(1)

wnn25_counts = pd.read_csv("../../results/wnn25/wnn25_raw.csv", index_col=0).transpose()
wnn25_counts_tensor = torch.tensor(wnn25_counts.to_numpy())
wnn25_mu1 = pd.read_csv("../../results/wnn25/wnn25_GMM_mu1.csv", index_col=0).transpose()
wnn25_mu1_tensor = torch.tensor(wnn25_mu1.to_numpy())

scPDA = model(raw_counts=wnn25_counts_tensor, bg_mean=wnn25_mu1_tensor)
scPDA.train()
scPDA.inference()

with h5py.File('../../results/wnn25/wnn25_scPDA.h5', 'w') as f:
    f.create_dataset('mu1', data=scPDA.mu1.numpy())
    f.create_dataset('mu2', data=(scPDA.mu1*scPDA.alpha).numpy())
    f.create_dataset('theta1', data=scPDA.theta1.numpy())
    f.create_dataset('theta2', data=scPDA.theta2.numpy())
    f.create_dataset('pi', data=scPDA.pi.numpy())