#!/bin/bash
#SBATCH -n 1 # Number of cores requested 
#SBATCH -J dgea_wil # job name
#SBATCH -N 1 # Ensure that all cores are on one machine 
#SBATCH -t 10000 # Runtime in minutes 
#SBATCH -p batch # Partition to submit to 
#SBATCH --account=abrunet1
#SBATCH --mem=256G # Memory per cpu in MB (see also --mem-per-cpu) 
#SBATCH -o dgea_wil_%A.out # Standard out goes to this file 
#SBATCH -e dgea_wil_%A.err # Standard err goes to this filehostname 

source activate geneimputation
python get_multi_wilcoxon.py SCPI_k4_k1 Dataset14 100 region none none knn_spage_tangram --non-symmetric
python get_multi_wilcoxon.py SCPI_k4_k1 Dataset13 100 celltype_lowres none none knn_spage_tangram --non-symmetric
python get_multi_wilcoxon.py SCPI_k4_k1 Dataset15 100 Region none none knn_spage_tangram --non-symmetric
python get_multi_wilcoxon.py SCPI_k4_k1 Axolotl 100 Annotation none none knn_spage_tangram --non-symmetric
source deactivate