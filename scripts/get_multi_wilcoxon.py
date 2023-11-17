'''
Runs TISSUE multiple imputation Mann-Whitney/Wilcoxon on specified dataset and testing setup (conditions)

Example: python get_multi_wilcoxon.py SCPI_k4_k1 Dataset15 100 Region none none knn_spage_tangram --non-symmetric
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq
from scipy.stats import pearsonr, spearmanr
from sklearn.decomposition import PCA, NMF
import anndata as ad
import warnings
import time
import pickle
import os

from TISSUE.tissue.main import *
from TISSUE.tissue.downstream import *
from TISSUE.tissue.experiments import *
from TISSUE.tissue.utils import *

import logging
logging.getLogger("imported_module").setLevel(logging.WARNING)

import argparse

# set up arguments
parser = argparse.ArgumentParser()
parser.add_argument("savedir", help="path to directory with intermediate results")
parser.add_argument("dataset", help="name of dataset folder in DataUpload/")
parser.add_argument("n_imputations", help="number of imputations", type=int)
parser.add_argument("condition", help="condition")
parser.add_argument("group1", help="group1")
parser.add_argument("group2", help="group2")
parser.add_argument("prediction_models", help="prediction model strings separated by '_'", type=str)
parser.add_argument('--symmetric', action='store_true')
parser.add_argument('--non-symmetric', dest='symmetric', action='store_false')
parser.set_defaults(symmetric=True)
args = parser.parse_args()

# load parameters from arguments
savedir = args.savedir
dataset_name = args.dataset
n_imputations = args.n_imputations
condition = args.condition
group1 = args.group1
group2 = args.group2
symmetric = args.symmetric
methods = list(args.prediction_models.split("_"))

if group1.lower() == "none":
    group1 = None
if group2.lower() == "none":
    group2 = None

# run MI testing
adata = group_multiple_imputation_testing_from_intermediate(dataset_name, methods, symmetric, condition, n_imputations=n_imputations,
                                                        group1=group1, group2=group2, savedir=savedir, test="wilcoxon_greater")
                
# save results in anndata
preprocess_data(adata, standardize=False, normalize=False) # to keep consistent with predictions
# if error loading (i.e. metadata too large), then large_save instead
try:
    adata.write(savedir+"/"+dataset_name+"_"+"_".join(methods)+"_MI_WILCOXON.h5ad")
    adata2 = sc.read_h5ad(savedir+"/"+dataset_name+"_"+"_".join(methods)+"_MI_WILCOXON.h5ad")
except:
    large_save(adata, savedir+"/"+dataset_name+"_"+"_".join(methods)+"_MI_WILCOXON")
    os.remove(savedir+"/"+dataset_name+"_"+"_".join(methods)+"_MI_WILCOXON.h5ad")