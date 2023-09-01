#!/bin/bash
#SBATCH -n 1 # Number of cores requested 
#SBATCH -J conformal9 # job name
#SBATCH -N 1 # Ensure that all cores are on one machine 
#SBATCH -t 3600 # Runtime in minutes 
#SBATCH -p batch # Partition to submit to 
#SBATCH --account=abrunet1
#SBATCH --mem=128G # Memory per cpu in MB (see also --mem-per-cpu) 
#SBATCH -o conformal9_%A.out # Standard out goes to this file 
#SBATCH -e conformal9_%A.err # Standard err goes to this filehostname 

source activate geneimputation
python spatial_conformal_uncertainty.py Dataset9 10 10 4 1 knn_spage_tangram --save_intermediate --non-symmetric
python get_calibration.py SCPI_k4_k1 Dataset9 knn_spage_tangram --non-symmetric
source deactivate