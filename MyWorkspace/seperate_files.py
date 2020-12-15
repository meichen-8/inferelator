#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 29 15:36:13 2020

@author: fang
"""
import numpy as np
import pandas as pd

#%% generate expression data files
# growth conditions in GEO and folder naming
growth_conditions=["DIAUXY", "YPETOH", "MMD", "MMETOH",  "NLIMNH4", "NLIMUREA", "NLIMGLN", \
                   "NLIMPRO", "YPD", "RAPA", "CSTARVE"]

# condition names in above data matrix
conditions=["YPDDiauxic", "YPEtOH", "MinimalGlucose", "MinimalEtOH", "AmmoniumSulfate", "Urea", "Glutamine", \
            "Proline", "YPD", "YPDRapa", "CStarve"]

c2gc=dict(zip(conditions, growth_conditions))

ori = pd.read_csv('./data/103118_SS_data.tsv', delimiter = '\t', header = 0, index_col = 0)

for c in conditions:
    data = ori[ori['Condition']==c]
    data.to_csv('./data/103118_SS_'+c2gc[c]+'.tsv', sep = '\t')
    
    
#%% TF

tf = pd.read_csv("./data/tf_names.tsv", sep = '\t', header = None)
tfr = pd.read_csv("./data/tf_names_restrict.tsv", sep = '\t', header = None)


    
gold_standard = pd.read_csv('./data/gold_standard.tsv', delimiter = '\t', header = 0, index_col = 0)
prior = pd.read_csv('./data/YEASTRACT_20190713_BOTH.tsv', delimiter = '\t', header = 0, index_col = 0)
prior_tf = list(prior.columns)

ctf = np.intersect1d(prior_tf, tfr[0])

