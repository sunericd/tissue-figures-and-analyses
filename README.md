# TISSUE Figures and Analyses
This repository contains all Jupyter notebooks and Python scripts for generating data and figures associated with the TISSUE manuscript. If you use this code or find it useful, we would appreciate it if you cite the relevant publication:

Sun ED, Ma R, Navarro Negredo P, Brunet A, Zou J. TISSUE: uncertainty-calibrated prediction of single-cell spatial transcriptomics improves downstream analyses.
Preprint at https://doi.org/10.1101/2023.04.25.538326 (2023).


To set up a conda environment for these analyses, we recommend installing all dependencies in a new conda environment and then setting that new environment as a jupyter kernel for use in the notebooks. Refer to ```environment.yml``` for the main dependencies that we used in these notebooks. In some select notebooks, a different conda environment was used (e.g. ```environment_dynamicviz.yml```, ```environment_merfish.yml```) and this will be noted in the header of the notebooks.

Jupyter notebooks containing code for making figures and running simulations for all results can be found in ```notebooks```. They are generally listed in order of chronology within the TISSUE manuscript with ```00_...``` corresponding to data preprocessing, ```01-05_...``` corresponding to panels for main figures, and ```Supplemental_...``` corresponding to supplementary analysis and figures. In detail, these notebooks are:
- ```00_format_srtsim_data.ipynb``` - formats SRTsim simulated data into TISSUE-readable format for experiments
- ```00_generate_srtsim_data.R``` - R script for generating simulated spatial transcriptomics data using SRTsim
- ```00_metadata_retrieval.ipynb``` - mapping of existing metadata to spatial transcriptomics datasets (or construction of new metadata)
- ```01_prediction_model_evaluation_and_TISSUE.ipynb``` - analyses corresponding to Figure 1 (performance of prediction methods, cell-centric variability)
- ```02_TISSUE_prediction_intervals.ipynb``` - analyses corresponding to Figure 2 (TISSUE calibration performance)
- ```03_TISSUE_multiple_imputation_for_differential_gene_expression.ipynb``` - analyses corresponding to Figure 3 (TISSUE multiple imputation and hypothesis testing)
- ```04_TISSUE_filtering_for_downstream_real_data.ipynb``` - analyses corresponding to Figure 4 for real data (TISSUE cell filtering for downstream tasks)
- ```04_TISSUE_filtering_for_downstream_simulated_data.ipynb``` - analyses corresponding to Figure 4 for simulated data (TISSUE cell filtering for downstream tasks)
- ```04_dynamicviz_pca_plots.ipynb``` - generates DynamicViz plots using multiple simulated datasets
- ```05_SVZ_MERFISH_case_study.ipynb``` - analyses corresponding to Figure 5 for our MERFISH adult mouse SVZ dataset
- ```05_SVZ_MERFISH_initial_processing.ipynb``` - initial processing of our MERFISH adult mouse SVZ dataset (pre-processing, clustering, annotation)
- ```Supplemental_TISSUE_WPCA_clustering_visualization.ipynb``` - Weighted principal component analysis using TISSUE prediction intervals
- ```Supplemental_gimvi_analysis.ipynb``` - gimVI prediction method results
- ```Supplemental_multiple_imputation_SpatialDE.ipynb``` - SpatialDE framework with TISSUE multiple imputation
- ```Supplemental_multiple_imputation_wilcoxon.ipynb``` - Mann-Whitney/Wilcoxon framework with TISSUE multiple imputation
- ```Supplemental_number_of_neighbors_sensitivity_analysis.ipynb``` - Sensitivity of TISSUE performance to the number of neighbors used in cell-centric variability
- ```Supplemental_predicted_expression_and_prediction_error_correlation.ipynb``` - Correlation between pairwise similarities of predicted expression and predicted error for motivating stratified grouping approach
- ```Supplemental_replicates_analysis.ipynb``` - reproducibility analysis using technical replicate of mouse gastrulation dataset
- ```Supplemental_runtime_analysis.ipynb``` - computational runtimes for TISSUE and prediction
- ```Supplemental_spatial_data_visualization.ipynb``` - visualization of datasets along spatial coordinates
- ```Supplemental_spatial_proximity_assumption_validation.ipynb``` - validation of TISSUE assumption on similarity in expression between neighbors
- ```Supplemental_sprod_denoising.ipynb``` - experiments involving Sprod de-noised data and using the Sprod cell similarity graph within TISSUE

Python scripts for running all TISSUE analyses (i.e. generating imputations and uncertainties through cross-validation) can be found in ```scripts```. These are general pipelines for generating results using the TISSUE package and their particular uses/relations to the figures and analyses are outlined in the Jupyter notebooks. In detail, these scripts are:
- ```get_calibration.py``` - generates ```_SCPI.h5ad``` AnnData objects containing the TISSUE prediction intervals and cross-validated gene expression predictions, and ```.pkl``` pickle files containing calibration curve values
- ```get_external_multi_ttest.py``` - generates ```_MI_EXTERNAL_TTEST.h5ad``` AnnData objects containing the TISSUE multiple imputation t-test statistics for differential gene expression on genes outside the spatial transcriptomics gene panel
- ```get_multi_spatialde.py``` - generates ```_MI_SPATIALDE.h5ad``` AnnData objects containing the TISSUE multiple imputation SpatialDE statistics for spatially variable gene expression on genes within the spatial transcriptomics gene panel
- ```get_multi_ttest.py``` - generates ```_MI_TTEST.h5ad``` AnnData objects containing the TISSUE multiple imputation t-test statistics for differential gene expression on genes within the spatial transcriptomics gene panel
- ```get_multi_wilcoxon.py``` - generates ```_MI_WILCOXON.h5ad``` AnnData objects containing the TISSUE multiple imputation Mann-Whitney/Wilcoxon statistics for differential gene expression on genes within the spatial transcriptomics gene panel
- ```get_score_for_dataset.py``` - generates ```_NOCV.h5ad``` AnnData objects containing predicted gene expression, cell-centric variability, stratified groups, etc for a single pass over the dataset
- ```spatial_conformal_uncertainty.py``` - generates intermediate gene expression predictions in a cross-validated manner; NOTE: this script should be run before any others
- ```spatial_conformal_uncertainty_sprod.py``` - variation using Sprod graph

Finally, this repository also contains several slurm job scripts that were used to run the associated Python scripts with different sets of parameters. These can be found in ```scripts/slurm_jobs``` with labeled subdirectories, and their associated use is noted in the corresponding notebooks. These job scripts will provide the parameter/settings that were used to run the analyses in the manuscript and also provide an (over-)estimate of the computational resources required for each job (e.g. memory, time, number of CPU/GPUs).