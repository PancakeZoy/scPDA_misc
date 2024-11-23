import os
import pandas as pd
from scar import model
import h5py
import torch
import random
import numpy as np

torch.manual_seed(1)
random.seed(1)
np.random.seed(1)

raw_ADT = pd.read_csv("../../results/tea/tea_raw.csv", index_col=0).transpose()
ambient_profile = pd.read_csv("../../results/tea/tea_EmptyProfile.csv", index_col=0)

ADT_scar_amb = model(raw_count = raw_ADT,
                     ambient_profile = ambient_profile,  
                     feature_type = 'ADT',
                     count_model = 'binomial',
                     device = 'cpu',
                     )
ADT_scar_amb.train(epochs=300,
                   batch_size=128,
                   verbose=True
                   )
ADT_scar_amb.inference()

denoised_ADT = pd.DataFrame(ADT_scar_amb.native_counts, index=raw_ADT.index, columns=raw_ADT.columns)
noise_ratio = ADT_scar_amb.noise_ratio[:,0]

with h5py.File('../../results/tea/tea_scAR.h5', 'w') as f:
    f.create_dataset('count', data=denoised_ADT)
    f.create_dataset('NoiseRatio', data=noise_ratio)