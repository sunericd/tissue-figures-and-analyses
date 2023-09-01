#!/bin/bash
#SBATCH -n 1 # Number of cores requested 
#SBATCH -J score9 # job name
#SBATCH -N 1 # Ensure that all cores are on one machine 
#SBATCH -t 3600 # Runtime in minutes 
#SBATCH -p batch # Partition to submit to 
#SBATCH --account=abrunet1
#SBATCH --mem=64G # Memory per cpu in MB (see also --mem-per-cpu) 
#SBATCH -o score9_%A.out # Standard out goes to this file 
#SBATCH -e score9_%A.err # Standard err goes to this filehostname 

source activate geneimputation
python get_score_for_dataset.py Dataset9 4 1 knn_spage_tangram --non-symmetric
source deactivate