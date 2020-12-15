#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 00:00:43 2020

@author: fang
"""
import numpy as np
import pandas as pd

my = np.loadtxt('./data/tf_names_yeastract.tsv', dtype = str, delimiter = '\t')
cj = np.loadtxt('./data/tf_names_yeastract_cj.txt', dtype = str, delimiter = '\t')

common = np.intersect1d(my, cj)

old = pd.read_csv('./data/YEASTRACT_Both_20181118.tsv', delimiter = '\t', header = 0, index_col = 0)
new = pd.read_csv('./data/YEASTRACT_20190713_BOTH.tsv', delimiter = '\t', header = 0, index_col = 0)
