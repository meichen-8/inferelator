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


def dict_g2tcc(all_tccs, t2g):
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
                g2tcc[gene] = tccs
    return g2tcc

def network_g2tcc(gene_filename, tcc_filename, g2tcc):
    gene_network = pd.read_csv(gene_filename, delimiter = '\t', header = 0, index_col = 0)
    col = list(gene_network.columns)
    idx = list(gene_network.index)
    tcc_idx = []
    tcc_col = []
    idx_map = []
    col_map = []
    # build index and map of index from genes_network to tccs_network
    j = 0
    for i,gene in enumerate(col):
        if gene in g2tcc:
            tccs = g2tcc[gene]
            tcc_col = np.append(tcc_col,tccs)
            col_map.append([j,j+len(tccs)])
            j = j + len(tccs)
        else:
            col_map.append([j,j])
            #print(gene)
    
    j = 0
    for i,gene in enumerate(idx):
        if gene in g2tcc:
            tccs = g2tcc[gene]
            tcc_idx = np.append(tcc_idx,tccs)
            idx_map.append([j,j+len(tccs)])
            j = j + len(tccs)
        else:
            idx_map.append([j,j])
            #print(gene)
            
    
    # generate tcc_network
    tcc_network = np.zeros([len(tcc_idx),len(tcc_col)])
    for value in [1,-1]:
        index = np.argwhere(gene_network.values==value)
        
        for k in index:
            i,j = k
            if idx_map[i][0] != idx_map[i][1] and col_map[j][0] != col_map[j][1]:
                tcc_network[ idx_map[i][0]:idx_map[i][1], col_map[j][0]:col_map[j][1] ] = value
            
    tcc_network = pd.DataFrame(data = tcc_network, index = tcc_idx, columns = tcc_col)
    tcc_network.to_csv(tcc_filename, sep = '\t')
    
    return tcc_network, gene_network



if __name__ == '__main__':
    
    tccs = ad.read_h5ad("./data/expression_data/all_tcc.h5ad")
    all_tccs = list(tccs.var_names)
    
    t2g = {}
    with open("./data/t2g_101.txt") as f:
        for line in f:
            l = line.split()
            t2g[l[0]] = l[1]
    g2tcc = dict_g2tcc(all_tccs, t2g)
    
    # gold_standard
    gene_filename = './data/gold_standard.tsv'
    tcc_filename = './data/gold_standard_tcc.tsv' 
    tcc_network, gene_network = network_g2tcc(gene_filename, tcc_filename, g2tcc)
    tcc_tf = list(tcc_network.columns)
    
    # prior
    gene_filename = './data/YEASTRACT_20190713_BOTH.tsv'
    tcc_filename = './data/YEASTRACT_20190713_BOTH_tcc.tsv' 
    tcc_network_, gene_network_ = network_g2tcc(gene_filename, tcc_filename, g2tcc)
    tcc_tf = np.intersect1d(tcc_tf, list(tcc_network_.columns))
    
    np.savetxt('./data/tf_names_tcc.tsv', tcc_tf, fmt = '%s', delimiter='\t')
    


    
    
    