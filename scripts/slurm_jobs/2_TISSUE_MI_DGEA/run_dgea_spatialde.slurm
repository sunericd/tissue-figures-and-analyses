#!/bin/bash
#SBATCH -n 1 # Number of cores requested 
#SBATCH -J dgea_spatialde # job name
#SBATCH -N 1 # Ensure that all cores are on one machine 
#SBATCH -t 10000 # Runtime in minutes 
#SBATCH -p batch # Partition to submit to 
#SBATCH --account=abrunet1
#SBATCH --mem=64G # Memory per cpu in MB (see also --mem-per-cpu) 
#SBATCH -o dgea_spatialde_%A.out # Standard out goes to this file 
#SBATCH -e dgea_spatialde_%A.err # Standard err goes to this filehostname 

source activate geneimputation
python get_multi_spatialde.py SCPI_k4_k1 Dataset11 10 none none none knn_spage_tangram --non-symmetric
python get_multi_spatialde.py SCPI_k4_k1 Dataset17 10 none none none knn_spage_tangram --non-symmetric
source deactivate