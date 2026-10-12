#Load packages
import os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import zscore

import sys
def calculate_clusters(df_reduced, gene_cluster_df, cluster_col='Cluster_k'):
    """
    Calculates the average gene cluster expression in each cell.

    Args:
        df_reduced (pandas.DataFrame):
            DataFrame containing the gene expression matrix with cells as rows and genes as columns.
        gene_cluster_df (pandas.DataFrame): 
            DataFrame containing gene information with a column for gene names and a column for cluster labels.
        cluster_col (str):
            The name of the column in `gene_cluster_df` that holds the cluster labels.

    Returns:
        results_dict (dict):
            A dictionary where keys are cluster labels (e.g., 'Cluster_0', 'Cluster_1', etc.) and values are pandas Series
            representing the mean expression values of genes in each cluster for each cell.
    """
    results_dict = {}
    
   
    # loop thru all the different clusters in k17 clustering
    for cluster_label in gene_cluster_df[cluster_col].unique():
        # get the gene lists if each cluster, and then use those list of genes to calc the mean of 
        #genes exp of the whole module per donor
        genes = gene_cluster_df[gene_cluster_df[cluster_col] == cluster_label]['Gene'].values
        mean_values = df_reduced[genes].mean(axis=1)
        
        # get the results dictionary
        results_dict[f'Cluster_{cluster_label}'] = mean_values
        
    return results_dict

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
    
    #MAKE OUTPUT DIR
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
            #pheno_cols = group.columns[len(filt_df_genes.columns):]
            filt_df_pheno = group[pheno_cols]

            # Merge the cluster expression data and pheno data
            comb_pheno_data = subtype_df.merge(filt_df_pheno, left_index=True, right_index=True)
            
            #Calculate zscores and filter outliers
            cluster_columns = [col for col in comb_pheno_data.columns if 'Cluster_' in col]
            z_scores = comb_pheno_data[cluster_columns].apply(zscore)
            outliers = (z_scores.abs() > 3).any(axis=1)
            filtered_data = comb_pheno_data[~outliers]
            # Save output
            filtered_output_file = os.path.join(output_dir, f'{cell_type}_{subtype}_test_gene_clust_metacell_zscorefilt.csv')
            filtered_data.to_csv(filtered_output_file, index=True)
            
            val_set_dict[subtype]= filtered_data

            print(f"Processed and saved data for {subtype}")
    
    return val_set_dict


# USER VARIABLES
main_dir = '/Users/ga0028vi/Documents/MCcharactertization/ROSMAP/Top_sets/Set1'
cell_types = ['Astrocytes', 'Excitatory', 'Inhibitory', 'Immune', 'Oligodendrocytes', 'OPCs', 'Vasculature']

for cell_type in cell_types:
    test_data_path = os.path.join(main_dir, 'Test_set', f'{cell_type}_metacells_test_subset_sparcfilt.csv')
    gene_clustering_dir = os.path.join(main_dir, 'Gene_clusters')  # Where your *_gene_clustering_list.csv files are stored

    # LOAD TEST DATA
    val_data = pd.read_csv(test_data_path)  # Reused as val_data in function, even though it’s the test set

    # LOAD GENE CLUSTER FILES INTO A DICTIONARY
    gene_clus_list = {}
    pattern = os.path.join(gene_clustering_dir, f'{cell_type}_*_gene_clustering_list.csv')
    for file in glob.glob(pattern):
        base = os.path.basename(file)
        subtype = base.replace(f'{cell_type}_', '').replace('_gene_clustering_list.csv', '')
        df = pd.read_csv(file)
        gene_clus_list[subtype] = df

    #PROCESS TEST SET
    test_set_dict = val_set_prep(val_data=val_data, 
                                gene_clus_list=gene_clus_list, 
                                cell_type=cell_type, 
                                main_dir=main_dir)



