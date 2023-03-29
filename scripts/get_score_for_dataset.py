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
parser.add_argument("dataset", help="name of dataset folder in DataUpload/")
parser.add_argument("k_gene", help="number of gene groups", type=int)
parser.add_argument("k_cell", help="number of cell groups", type=int)
parser.add_argument('--symmetric', action='store_true')
parser.add_argument('--non-symmetric', dest='symmetric', action='store_false')
parser.set_defaults(symmetric=True)
args = parser.parse_args()

dataset_name = args.dataset
k_gene = args.k_gene
k_cell = args.k_cell
symmetric = args.symmetric
methods = ["knn", "spage", "tangram"]

savedir = "SCPI_k"+str(k_gene)+"_k"+str(k_cell)

# read in spatial and RNAseq datasets
if os.path.isfile("DataUpload/"+dataset_name+"/Metadata.txt"):
    adata, RNAseq_adata = load_paired_datasets("DataUpload/"+dataset_name+"/Spatial_count.txt",
                                                "DataUpload/"+dataset_name+"/Locations.txt",
                                                "DataUpload/"+dataset_name+"/scRNA_count.txt",
                                                spatial_metadata = "DataUpload/"+dataset_name+"/Metadata.txt")
else:
    adata, RNAseq_adata = load_paired_datasets("DataUpload/"+dataset_name+"/Spatial_count.txt",
                                                "DataUpload/"+dataset_name+"/Locations.txt",
                                                "DataUpload/"+dataset_name+"/scRNA_count.txt")
adata.var_names = [x.lower() for x in adata.var_names]
RNAseq_adata.var_names = [x.lower() for x in RNAseq_adata.var_names]

# preprocess RNAseq data
preprocess_data(RNAseq_adata, standardize=False, normalize=True)

# subset spatial data into shared genes
gene_names = np.intersect1d(adata.var_names, RNAseq_adata.var_names)
adata = adata[:, gene_names]

# filler/token gene to use for prediction task
target_genes = np.array([np.setdiff1d(RNAseq_adata.var_names, gene_names)[0]])

# build spatial graph
build_spatial_graph(adata, method="fixed_radius", n_neighbors=15)
#calc_adjacency_weights(adata, method="cosine")


# iterate through different methods and: (1) impute gene expression, (2) compute calibration scores, (3) run DGEA with MI
for method in methods:
    
    # (1) impute gene expression
    if method == "spage":
        if len(adata.var_names) < 40:
            n_pv = 20
        else:
            n_pv = round(np.min([len(adata.var_names), len(adata.obs_names)])/2)
        
        predict_gene_expression (adata, RNAseq_adata, target_genes,
                                 method=method, n_folds=10, n_pv=n_pv)
    elif method == "knn":
        predict_gene_expression (adata, RNAseq_adata, target_genes,
                                 method=method, n_folds=10, n_neighbors=10)
    else:
        predict_gene_expression (adata, RNAseq_adata, target_genes,
                                 method=method, n_folds=10)
    
    # (2) compute calibration scores
    predicted = method+"_predicted_expression"
    calib_genes = gene_names
    conformalize_spatial_uncertainty(adata, predicted, calib_genes, weight="exp_cos", mean_normalized=False, add_one=True,
                                     grouping_method="kmeans_gene_cell", k=k_gene, k2=k_cell, n_pc=15)

# Save results
adata.write(savedir+"/"+dataset_name+"_"+"_".join(methods)+"_NOCV.h5ad")