# CFS: Correlation-based feature selector
#
import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob

DATAPATH = './data/'
PARAMS_FILE = './params_std'
ele = 'Si'
nF  = 26
CORR = 0.98
N_iteration = 10
sampling_size = 500

pt_list = np.array(glob.glob(f'{DATAPATH}/*'))
with open('./scale_factor','rb') as f:
    scale = torch.load(f)

CSF = []
for N in range(N_iteration):
    random_idx = np.random.randint(low=1,high=len(pt_list),size=sampling_size)
    random_idx.sort()
    sample = pt_list[random_idx]
    feature_list = np.empty((0,nF))
    for idx, item in enumerate(sample):
        with open(item,'rb') as f:
            tmp = torch.load(f)
        feature_list = np.vstack((feature_list,tmp['x'][ele]))
    scaled = (feature_list - scale[ele][0,:]) / scale[ele][1,:]
    
    df = pd.DataFrame(scaled)
    corr = df.corr()
    columns = np.full((corr.shape[0],), True, dtype=bool)
    for i in range(corr.shape[0]):
        for j in range(i+1, corr.shape[0]):
            if corr.iloc[i,j] >= CORR:
                if columns[j]:
                    columns[j] = False
    selected_cols = df.columns[columns]
    select_df = df[selected_cols]
    print(len(selected_cols), selected_cols)
    CSF += selected_cols.to_list()
CSF = np.array(list(set(CSF)))

params = open(PARAMS_FILE,'r').readlines()
opt_params = open(f'params_opt_{len(CSF)}','w')
for i in range(len(CSF)):
    opt_params.write(f'{params[CSF[i]]}')
opt_params.close()
