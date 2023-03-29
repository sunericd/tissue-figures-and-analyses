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
parser = argparse.ArgumentParser()
parser.add_argument("savedir", help="path to directory with intermediate results")
parser.add_argument("dataset", help="name of dataset folder in DataUpload/")
parser.add_argument('--symmetric', action='store_true')
parser.add_argument('--non-symmetric', dest='symmetric', action='store_false')
parser.set_defaults(symmetric=True)
args = parser.parse_args()

savedir = args.savedir
dataset_name = args.dataset
symmetric = args.symmetric
methods = ["knn", "spage", "tangram"]
alpha_levels = np.linspace(0.01, 0.99, 1000)

# run calibration
res_dict, adata = group_conformalize_from_intermediate(dataset_name, methods, symmetric, alpha_levels,
                                                       save_alpha=[0.05,0.1,0.2,0.5,0.33], savedir=savedir)  

# pickle res_dict
with open(savedir+"/"+dataset_name+"_conformal_dict.pkl", "wb") as f:
    pickle.dump(res_dict, f)
                
# save results in anndata
preprocess_data(adata, standardize=False, normalize=False) # to keep consistent with predictions
adata.write(savedir+"/"+dataset_name+"_"+"_".join(methods)+"_SCPI.h5ad")

