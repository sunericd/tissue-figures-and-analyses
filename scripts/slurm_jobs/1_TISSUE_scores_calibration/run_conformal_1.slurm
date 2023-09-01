#!/bin/bash
#SBATCH -n 1 # Number of cores requested 
#SBATCH -J conformal1 # job name
#SBATCH -N 1 # Ensure that all cores are on one machine 
#SBATCH -t 3600 # Runtime in minutes 
#SBATCH -p batch # Partition to submit to 
#SBATCH --account=abrunet1
#SBATCH --mem=256G # Memory per cpu in MB (see also --mem-per-cpu) 
#SBATCH -o conformal1_%A.out # Standard out goes to this file 
#SBATCH -e conformal1_%A.err # Standard err goes to this filehostname 

source activate geneimputation
python spatial_conformal_uncertainty.py Dataset1 10 10 4 1 knn_spage_tangram --save_intermediate --non-symmetric
python get_calibration.py SCPI_k4_k1 Dataset1 knn_spage_tangram --non-symmetric --num_alphas 100
conda deactivate