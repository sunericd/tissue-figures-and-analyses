#!/bin/bash
#SBATCH -n 1 # Number of cores requested 
#SBATCH -J conformalAxolotl # job name
#SBATCH -N 1 # Ensure that all cores are on one machine 
#SBATCH -t 10080 # Runtime in minutes 
#SBATCH -p batch # Partition to submit to 
#SBATCH --account=abrunet1
#SBATCH --mem=64G # Memory per cpu in MB (see also --mem-per-cpu) 
#SBATCH -o conformalAxolotl_%A.out # Standard out goes to this file 
#SBATCH -e conformalAxolotl_%A.err # Standard err goes to this filehostname 

source activate geneimputation
python spatial_conformal_uncertainty.py Axolotl 10 10 4 1 knn_spage_tangram --save_intermediate --non-symmetric
python get_calibration.py SCPI_k4_k1 Axolotl knn_spage_tangram --non-symmetric
conda deactivate