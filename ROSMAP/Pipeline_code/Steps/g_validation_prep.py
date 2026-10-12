#Load packages
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import zscore
import matplotlib
matplotlib.use('Agg')
import sys
from Scripts.metacell_proc_utils import calculate_clusters

def val_set_prep(val_data, gene_clus_list, cell_type, main_dir):
    """
    ** STEP 7 **
    Prepares the validation set by applying gene cluster data to filter and analyze gene expression.

    This function: 
    1. Gets validation data and gene clust list passed from previous functions
    2. According to predefined gene clusters, calculating average cluster expressions, merging with phenotype data,
    filtering outliers based on z-scores, and saving the processed data.

    Parameters:
        val_data (DataFrame): The DataFrame containing validation set data with gene expressions and phenotypes.
        gene_clus_list (dict): A dictionary where keys are subtypes and values are DataFrames with gene cluster information.
        cell_type (str): The major cell type name, used in naming output files.
        output_dir (str): The directory path where the output files will be saved.

    Returns:
        val_set_dict (dict): A dictionary where keys are subtypes and values are DataFrames of the processed validation data for each subtype.

    Each subtype's DataFrame includes gene expression data filtered by cluster presence and phenotypic information.
    Outliers based on z-scores are removed from this analysis. The function also prints a confirmation message for each processed subtype.
    """
    
    # Make output dir
    output_dir = os.path.join(main_dir, 'Gene_clusters')
    os.makedirs(output_dir, exist_ok=True)
    
    val_set_dict = {}
    # Group data by cell subtype
    for subtype, group in val_data.groupby('cell_subtype'): 
        if subtype in gene_clus_list:
            clusters_df = gene_clus_list[subtype]

            # Filter the data by the genes listed in the clusters
            filt_df_genes = group[clusters_df['Gene'].values] if 'Gene' in clusters_df.columns else group

            # Calculate average clusters expression per subtype
            results = calculate_clusters(filt_df_genes, clusters_df, cluster_col='Cluster_k')
            subtype_df = pd.DataFrame(results)
            
            # Assume phenotype data follows after gene columns      
            index_of_projid = group.columns.get_loc('sparsity')
            pheno_cols = group.columns[index_of_projid:]  
            filt_df_pheno = group[pheno_cols]

            # Merge the cluster expression data and pheno data
            comb_pheno_data = subtype_df.merge(filt_df_pheno, left_index=True, right_index=True)

            #Calculate zscores and filter outliers
            cluster_columns = [col for col in comb_pheno_data.columns if 'Cluster_' in col]
            z_scores = comb_pheno_data[cluster_columns].apply(zscore)
            outliers = (z_scores.abs() > 3).any(axis=1)
            filtered_data = comb_pheno_data[~outliers]
            # Save output
            filtered_output_file = os.path.join(output_dir, f'{cell_type}_{subtype}_validation_gene_clust_metacell_zscorefilt.csv')
            filtered_data.to_csv(filtered_output_file, index=True)
            
            val_set_dict[subtype]= filtered_data

            print(f"Processed and saved data for {subtype}")
    
    return val_set_dict
