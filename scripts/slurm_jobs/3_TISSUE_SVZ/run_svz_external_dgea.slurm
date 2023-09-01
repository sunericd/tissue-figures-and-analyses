#!/bin/bash
#SBATCH -n 1 # Number of cores requested 
#SBATCH -J svz_dgea # job name
#SBATCH -N 1 # Ensure that all cores are on one machine 
#SBATCH -t 1000 # Runtime in minutes 
#SBATCH -p batch # Partition to submit to 
#SBATCH --account=abrunet1
#SBATCH --mem=64G # Memory per cpu in MB (see also --mem-per-cpu) 
#SBATCH -o svz_dgea_%A.out # Standard out goes to this file 
#SBATCH -e svz_dgea_%A.err # Standard err goes to this filehostname 

source activate geneimputation
python get_external_multi_ttest.py SVZsubambiguous lineage_marker_genes.txt 100 celltype none none 4 1 knn_spage_tangram --non-symmetric
source deactivate