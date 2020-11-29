#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 28 23:45:01 2020

@author: fang

process inferelator results

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

f = pd.read_csv("./output/figure_5d_stl_2/crossvalidation_performance.tsv", sep = '\t', header = 0)

figure5table = pd.read_csv( "../../../../data/source-code1/data/STable5.tsv", sep = '\t', header = 0)
figure5d = figure5table[ figure5table['Figure'] == 'Figure_5d']
figure5d_stl = figure5d[ figure5d['Tasking']=='STL' ]


tf = pd.read_csv("./data/tf_names.tsv", sep = '\t', header = 0)
tfr = pd.read_csv("./data/tf_names_restrict.tsv", sep = '\t', header = 0)
