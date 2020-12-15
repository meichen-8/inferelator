# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 10:12:03 2020

@author: fang

"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, average_precision_score

#%%
def aupr_(outputdir, figfile):
    AUPRs = []
    i = -1
    fig,ax = plt.subplots(10, 1, figsize = (6,60))
    random_seed_folders = os.listdir(outputdir)
    
    for rs in random_seed_folders:
        if re.match("random_seed",rs):
            i = i + 1
            edges = pd.read_csv(os.path.join(outputdir, rs, 'network.tsv'), sep='\t', header=0)
            y = edges.loc[:,['gold_standard','combined_confidences']]
            y.replace('', np.nan, inplace=True)
            y.dropna(inplace=True)
            y = y.values
            precision, recall, thresholds = precision_recall_curve(y[:,0], y[:,1])    
            ap = average_precision_score(y[:,0], y[:,1])
            AUPRs.append(ap)
            ax[i].plot(recall, precision, '.')
            ax[i].set_ylabel('precision')
            ax[i].set_xlabel('recall')
            ax[i].set_title('AUPR={0:0.2f}'.format(ap))
        
    plt.savefig(figfile)
    plt.close()
    return AUPRs

#%%
AUPRs = []
outputdirs = ['./output/figure_5d_stl_2', './output_gene/figure_5d_stl_1']
figfiles = ['./figure/ori-5d-stl-2.png', './figure/gene-5d-stl-1.png']

for i in range(len(figfiles)):
    auprs = aupr_(outputdirs[i], figfiles[i])
    AUPRs.append(auprs)

#%% aupr_tcc(outputdir, figfile)
outputdir = './output_tcc/figure_5d_stl_3'
figfile = './figure/tcc-5d-stl-3.png'
gs = pd.read_csv('./data/gold_standard_tcc.tsv', sep = '\t', header = 0, index_col = 0)
gs_col = list(gs.columns)
gs_idx = list(gs.index)
auprs = []
i = -1
fig,ax = plt.subplots(10, 1, figsize = (6,60))
random_seed_folders = os.listdir(outputdir)

for rs in random_seed_folders:
    if re.match("random_seed",rs):
        i = i + 1
        edges = pd.read_csv(os.path.join(outputdir, rs, 'network.tsv'), sep='\t', header=0)
        target = edges.loc[:,'target']
        regulator = edges.loc[:,'regulator']
        score = edges.loc[:,'combined_confidences']
        new_target = []
        new_regulator = []
        in_gs =[]
        new_score = []
        for k in range(edges.shape[0]):
            tfs = regulator[k].split(',')
            genes = target[k].split(',')
            for tf in tfs:
                if tf in gs_col:
                    for gene in genes:
                        if gene in gs_idx:
                            new_target.append(gene)
                            new_regulator.append(tf)
                            in_gs.append(abs(gs.loc[gene,tf]))
                            new_score.append(k)
                            
        new_edges = {'target': new_target, 'regulator':new_regulator, 'gold_standard': in_gs, 'combined_rank': new_score}
        new_edges = pd.DataFrame(data = new_edges)
        new_edges = new_edges.groupby(['target','regulator']).mean()
        new_edges.sort_values(by='combined_rank', axis=0, ascending=True, inplace=True)
        new_edges['combined_rank'] = (np.max(new_edges['combined_rank'])-new_edges['combined_rank'])/np.max(new_edges['combined_rank'])
        
        precision, recall, thresholds = precision_recall_curve(new_edges['gold_standard'], new_edges['combined_rank'])    
        ap = average_precision_score(new_edges['gold_standard'], new_edges['combined_rank'])
        auprs.append(ap)
        ax[i].plot(recall, precision, '.')
        ax[i].set_ylabel('precision')
        ax[i].set_xlabel('recall')
        ax[i].set_title('AUPR={0:0.2f}'.format(ap))
    
plt.savefig(figfile)
plt.close()

AUPRs.append(auprs)

#%%
plt.close('all')
fig, ax = plt.subplots(1,1, figsize = (6,5))
ax.plot(AUPRs, 'k.', alpha = 0.2)
ax.scatter(['paper', 'kl-gene', 'kl-Tcc-3'], np.mean(AUPRs,axis = 1), c='r')
plt.savefig("./figure/AUPRs.pdf")
plt.close()