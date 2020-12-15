#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 11:09:18 2020

@author: fang
"""

# Load modules
from inferelator import inferelator_workflow, inferelator_verbose_level, MPControl, CrossValidationManager
import pandas as pd

# Set verbosity level to "Talky"
inferelator_verbose_level(1)

# Set the location of the input data and the desired location of the output files
DATA_DIR = './data/'
OUTPUT_DIR = './output_gene/'

growth_conditions=["DIAUXY", "YPETOH", "MMD", "MMETOH",  "NLIMNH4", "NLIMUREA", "NLIMGLN", "NLIMPRO", "YPD", "RAPA", "CSTARVE"]


EXPRESSION_FILE_NAME = 'expression_data/all_gene.h5ad'
#GENE_METADATA_FILE_NAME = 'orfs.tsv'
GOLD_STANDARD_FILE_NAME = 'gold_standard.tsv'
METADATA_COLUMNS = ['TF', 'strain', 'date', 'restriction', 'mechanism', 'time']

YEASTRACT_PRIOR = "YEASTRACT_20190713_BOTH.tsv"

TF_NAMES = "tf_names_yeastract.tsv"


n_cores_local = 16
local_engine = True

def set_up_workflow(wkf):
    wkf.set_file_paths(input_dir = DATA_DIR,
                       output_dir = OUTPUT_DIR,
                       gold_standard_file = GOLD_STANDARD_FILE_NAME,
                       priors_file = YEASTRACT_PRIOR,
                       tf_names_file=TF_NAMES)
    
    wkf.set_expression_file(h5ad = EXPRESSION_FILE_NAME)
    
    wkf.set_file_properties(extract_metadata_from_expression_matrix=True,
                            expression_matrix_metadata=METADATA_COLUMNS,
                            expression_matrix_columns_are_genes=True)
    
    wkf.set_crossvalidation_parameters(split_gold_standard_for_crossvalidation=True,cv_split_ratio=0.5, cv_split_axis=0)                               

    wkf.set_run_parameters(num_bootstraps=5)
    wkf.set_count_minimum(0.05)
    wkf.add_preprocess_step("log2")
    return wkf


def set_up_cv(wkf):
    cv_wrap = CrossValidationManager(wkf)
    cv_wrap.add_gridsearch_parameter('random_seed', list(range(52, 62)))
    return cv_wrap

def set_up_task(wkf, name):
    task = wkf.create_task(task_name = name,
                           input_dir = DATA_DIR,
                           tf_names_file = TF_NAMES,
                           priors_file = YEASTRACT_PRIOR,
                           workflow_type = "single-cell")
                           #meta_data_file = GENE_METADATA_FILE_NAME,
    
    task.set_expression_file(h5ad = 'expression_data/'+name+'_gene.h5ad')
    
    task.set_file_properties(extract_metadata_from_expression_matrix=True,
                            expression_matrix_metadata=METADATA_COLUMNS,
                            expression_matrix_columns_are_genes=True)
                            #gene_list_index="SystematicName")
    return wkf

# %%
if __name__ == '__main__' and local_engine:
    MPControl.set_multiprocess_engine("multiprocessing")
    MPControl.client.processes = n_cores_local
    MPControl.connect()
# %%
    # Figure 5D: STL
    wkfname = 'figure_5d_stl_20201214'
    worker = set_up_workflow(inferelator_workflow(regression="bbsr", workflow="single-cell"))
    worker.append_to_path('output_dir', wkfname)
   
    cv_wrap = set_up_cv(worker)
    cv_wrap.run()
    
    del worker
    del cv_wrap
    
    info = {'expression_fiie': EXPRESSION_FILE_NAME, 'gold_standard_file': GOLD_STANDARD_FILE_NAME ,'priors': YEASTRACT_PRIOR, 'tf_names_file': TF_NAMES, 'cv_split_axis': 0, 'cv_split_ratio': 0.5, 'num_bootstraps': 5}
    info = pd.DataFrame.from_dict(info, orient='index', columns = ['value'])

    info.to_csv(OUTPUT_DIR + wkfname + '/information.csv')
  
#%%
    # Figure 5D: MTL_bbsr
    wkfname = 'figure_5d_mtl_bbsr_20201214'
    worker = inferelator_workflow(regression="bbsr", workflow="multitask")
    worker.set_file_paths(input_dir = DATA_DIR,
                       output_dir = OUTPUT_DIR, 
                       gold_standard_file = GOLD_STANDARD_FILE_NAME)
    
    for gc in growth_conditions:
        worker = set_up_task(worker, gc)
        
    worker.set_crossvalidation_parameters(split_gold_standard_for_crossvalidation=True, cv_split_ratio=0.5, cv_split_axis=0)
    worker.set_run_parameters(num_bootstraps=5)
    worker.set_count_minimum(0.05)
    worker.add_preprocess_step("log2")
    worker.append_to_path("output_dir", wkfname)
    
    cv_wrap = set_up_cv(worker)
    cv_wrap.run()
    
    del worker
    del cv_wrap
    
    info = {'expression_fiie': EXPRESSION_FILE_NAME, 'gold_standard_file': GOLD_STANDARD_FILE_NAME ,'priors': YEASTRACT_PRIOR, 'tf_names_file': TF_NAMES, 'cv_split_axis': 0, 'cv_split_ratio': 0.5, 'num_bootstraps': 5}
    info = pd.DataFrame.from_dict(info, orient='index', columns = ['value'])
    info.to_csv(OUTPUT_DIR + wkfname + '/information.csv')
    
#%% Figure 5D: MTL_bbsr
    wkfname = 'figure_5d_mtl_amusr_20201214'
    worker = inferelator_workflow(regression="amusr", workflow="multitask")
    worker.set_file_paths(input_dir = DATA_DIR,
                       output_dir = OUTPUT_DIR, 
                       gold_standard_file = GOLD_STANDARD_FILE_NAME)
    
    for gc in growth_conditions:
        worker = set_up_task(worker, gc)
        
    worker.set_crossvalidation_parameters(split_gold_standard_for_crossvalidation=True, cv_split_ratio=0.5, cv_split_axis=0)
    worker.set_run_parameters(num_bootstraps=5)
    worker.set_count_minimum(0.05)
    worker.add_preprocess_step("log2")
    worker.append_to_path("output_dir", wkfname)
    
    cv_wrap = set_up_cv(worker)
    cv_wrap.run()
    
    del worker
    del cv_wrap
    
    info = {'expression_fiie': EXPRESSION_FILE_NAME, 'gold_standard_file': GOLD_STANDARD_FILE_NAME ,'priors': YEASTRACT_PRIOR, 'tf_names_file': TF_NAMES, 'cv_split_axis': 0, 'cv_split_ratio': 0.5, 'num_bootstraps': 5}
    info = pd.DataFrame.from_dict(info, orient='index', columns = ['value'])
    info.to_csv(OUTPUT_DIR + wkfname + '/information.csv')