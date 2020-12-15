#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 16:35:47 2020

@author: fang

Generate gold standard and prior network for TCCs


"""
import anndata as ad
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt


def dict_g2ec(all_tccs, t2g):
    g2tcc = {}   
    for tccs in all_tccs:
        tlist = tccs.split(',')
        for t in tlist:
            if t in t2g:
                gene = t2g[t]
            else:
                gene = t.split('_')[0]
            if gene in g2tcc:
                temp = g2tcc[gene]
                g2tcc[gene] = np.append(temp, tccs)
            else:
                g2tcc[gene] = [tccs]
    return g2tcc

def list_g2ec(gene_filename, tcc_filename, g2tcc):
    genes = np.loadtxt(gene_filename, dtype = str, delimiter = '\t')
    idx = genes[0]
    tccs = []
    for gene in idx:
        if gene in g2tcc:
            tccs = np.append(tccs, g2tcc[gene])
    np.savetxt(tcc_filename, tccs, fmt = '%s' , delimiter = '\t')
    return tccs


#%%

if __name__ == '__main__':
    
    tccs = ad.read_h5ad("./data/expression_data/all_tcc.h5ad")
    all_tccs = list(tccs.var_names)
#%%   
    t2g = {}
    g2t = {}
    iso = 0
    with open("./data/t2g_101.txt") as f:
        for line in f:
            l = line.split()
            t2g[l[0]] = l[1]
            if l[1] in g2t:
                g2t[l[1]] = np.append(g2t[l[1]], l[0])
                iso = iso + 1
            else:
                g2t[l[1]] = l[0]
#%%                
                
    g2ec = dict_g2ec(all_tccs, t2g)

    # prior
    gene_filename = './data/YEASTRACT_20190713_BOTH.tsv'
    tcc_filename = './data/YEASTRACT_20190713_BOTH_tcc-t.tsv' 
    
    gene_network = pd.read_csv(gene_filename, delimiter = '\t', header = 0, index_col = 0)
    col = list(gene_network.columns)
    idx = list(gene_network.index)
    tcc_idx = []
    tcc_col = []
    idx_map = []
    #col_map = []
    # build index and map of index from genes_network to tccs_network
    
    j = 0
    for i,gene in enumerate(col):
        tcc_col.append(gene+'_mRNA')
        #if gene in g2ec:
            #tccs = g2ec[gene]
            #tcc_col = np.append(tcc_col,tccs)
            #col_map.append([j,j+len(tccs)])
            #j = j + len(tccs)
        #else:
            #col_map.append([j,j])
            #print(gene)
    
    j = 0
    for i,gene in enumerate(idx):
        if gene in g2ec:
            tccs = g2ec[gene]
            tcc_idx = np.append(tcc_idx,tccs)
            idx_map.append([j,j+len(tccs)])
            j = j + len(tccs)
        else:
            idx_map.append([j,j])
            #print(gene)
            
    
    # generate tcc_network
    tcc_network = np.zeros([len(tcc_idx),len(tcc_col)])
    value = 1
    index = np.argwhere(gene_network.values==value)
    
    for k in index:
        i,j = k
        if idx_map[i][0] != idx_map[i][1]:
            tcc_network[ idx_map[i][0]:idx_map[i][1], j ] = value/(idx_map[i][1]-idx_map[i][0])
        
    tcc_network = pd.DataFrame(data = tcc_network, index = tcc_idx, columns = tcc_col)
    tcc_network_ = tcc_network.groupby(by = str).sum()
    #tcc_network_ = tcc_network_.groupby(by = str, axis=1).sum()
    tcc_network_.to_csv(tcc_filename, sep = '\t')
    
    tf_tcc = list(tcc_network.columns)
    np.savetxt('./data/tf_names_yeastract_tcc.tsv', tf_tcc, fmt = '%s' , delimiter = '\t')
   
    #tf_ = list(gene_network_.columns)
    #np.savetxt('./data/tf_names_yeastract.tsv', tf_, fmt = '%s' , delimiter = '\t')

#%%   
    gs = pd.read_csv('./data/gold_standard.tsv', delimiter = '\t', header = 0, index_col = 0)
    new_index = []
    index_mask = []
    for g in gs.index:
        if g in g2t:
            index_mask.append(True)
            new_index.append(g2t[g])
        else:
            index_mask.append(False)
    new_col = []
    col_mask = []
    for g in gs.columns:
        if g in g2t:
            col_mask.append(True)
            new_col.append(g2t[g])
        else:
            col_mask.append(False)        
            
    new_gs = gs.values[index_mask,:]
    new_gs = new_gs[:,col_mask]
    new_gs = pd.DataFrame(data = new_gs , index = new_index, columns = new_col)
    new_gs.to_csv('./data/gold_standard_tcc.tsv', sep = '\t')
    

    
#%%
    priors = tcc_network_
    gs_label = list(new_gs.index)+list(new_gs.columns)
    droplist = []
    for names in list(priors.index):
        name = names.split(',')
        flag = 0
        for t in name:
            if t in gs_label:
                flag = 1
                break
        if flag == 1:
            droplist.append(names)
                            
    split_axis = 0
    new_priors = priors.drop(droplist, axis=split_axis, errors='ignore')
    