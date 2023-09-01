#!/bin/bash
#SBATCH -J conformal11 # job name
#SBATCH -G 1
#SBATCH --cpus-per-gpu=4
#SBATCH -t 2880 # Runtime in minutes 
#SBATCH -p owners # Partition to submit to 
#SBATCH --account=jamesz
#SBATCH --mem=128G # Memory per cpu in MB (see also --mem-per-cpu) 
#SBATCH -o conformal11_gimvi_%A.out # Standard out goes to this file 
#SBATCH -e conformal11_gimvi_%A.err # Standard err goes to this filehostname 

source activate geneimputation
python spatial_conformal_uncertainty.py Dataset11 10 10 4 1 gimvi --save_intermediate --non-symmetric --no-preprocess_RNA
python get_calibration.py SCPI_k4_k1 Dataset11 gimvi --non-symmetric
conda deactivate