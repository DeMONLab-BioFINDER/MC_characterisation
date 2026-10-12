#Load packages
import os
import pandas as pd
import numpy as np
from sklearn.cluster import AgglomerativeClustering
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import zscore
import sys
import matplotlib
matplotlib.use('Agg')

from Scripts.metacell_proc_utils import get_cmats_df, calculate_clusters

#Analysing each cell type file (iterate per each major cell type)
def cluster_genes_main(train_metacells, main_dir, cell_type, subtype_k):
    """
    ** STEP 6 **
    Processes gene clusters for each cell subtype based on optimal clustering resolution and performs post-processing 
    including z-score filtering. Additionally, it generates heatmaps correlating phenotype data with gene clusters.

    Parameters:
        train_metacells (DataFrame): The DataFrame containing metacell data for training.
        main_dir (str): The main directory path for output.
        cell_type (str): The cell type name used in naming output files.
        subtype_k (dict): A dictionary containing the optimal number of clusters for each subtype.

    Returns:
        Tuple: A tuple containing two dictionaries:
        - subtype_clusters (dict): A dictionary mapping each cell subtype to its DataFrame after z-score filtering.
        - subtype_geneclus_list (dict): A dictionary mapping each cell subtype to its DataFrame of genes and their cluster labels.

    Outputs:
       Several CSV files for each cell subtype with gene clusters and phenotype data, filtered by z-score.
       Correlation heatmaps between gene clusters and phenotype data are saved as PNG files.
    """

    # Make output dir
    output_dir = os.path.join(main_dir, 'Gene_clusters')
    os.makedirs(output_dir, exist_ok=True)

    print(f"Starting cluster_genes_main data for {cell_type}")


    # Get cmats
    last_col = train_metacells.columns.get_loc('sparsity') #Get column which splits genes vs pheno
    # Group data by 'cell_subtype' and store each group as a separate DataFrame including phenotype data
    grouped_data = train_metacells.groupby('cell_subtype')
    subtype_dfs = {name: group.copy() for name, group in grouped_data}

    filtered_subtype_dfs = {}
    for subtype, df in subtype_dfs.items():
        # Select the gene expression columns
        gene_data = df.iloc[:, :last_col]
        # Filter to keep only columns with at least one non-zero value
        expressed_genes = gene_data.loc[:, (gene_data != 0).any(axis=0)]
        # Combine the filtered gene expression data with phenotype data
        phenotype_data = df.iloc[:, last_col:]
        filtered_subtype_dfs[subtype] = pd.concat([expressed_genes, phenotype_data], axis=1)

    # Generate correlation matrices for the new groups using function
    correlation_matrices = {}
    for subtype, subtype_df in filtered_subtype_dfs.items():
        cmat = get_cmats_df(subtype_df, 'sparsity')
        correlation_matrices[subtype] = cmat

    #Loop through each cell subtype
    subtype_geneclus_list = {}
    subtype_clusters = {}

    #var subtype contains both cell type and subtype in it
    for subtype, cmat in correlation_matrices.items():
        print(f"Starting cluster_genes_main data for {cell_type}_{subtype}")

        #check if the subtype is in the k value list, if so proceed with processing
        if subtype in subtype_k:
            k = subtype_k[subtype][0] #get the k value for current cell subtype
            print(f"Processing {cell_type}_{subtype} with k={k}")

            cluster = AgglomerativeClustering(n_clusters=k, metric='euclidean', linkage='ward')
            cluster_labels = cluster.fit_predict(cmat)
            border = filtered_subtype_dfs[subtype].columns.get_loc('sparsity') #Get column which splits genes vs pheno
            gene_list = list(filtered_subtype_dfs[subtype].columns[:border]) 

            #Get gene, cluster df for the subtype
            clusters_df = pd.DataFrame({'Gene': gene_list, #this needs to be a data table with specific query of columns up till variable 'last_column'
                                       'Cluster_k': cluster_labels})
            clusters_df.to_csv(os.path.join(output_dir, f"{cell_type}_{subtype}_gene_clustering_list.csv"), index=False)
            subtype_geneclus_list[subtype] = clusters_df #store in dict
            
            #filter the data to keep only gene data_frame (for specific cell subtype)
            filt_df_genes = filtered_subtype_dfs[subtype].iloc[:, :border]
            
            # Calculate average clusters expression per subtype and store the results in a DataFrame
            results = calculate_clusters(filt_df_genes, clusters_df, cluster_col='Cluster_k')
            subtype_df = pd.DataFrame(results)
            
            #Filter the data to only keep pheno data (for specific cell subtype)
            filt_df_pheno = filtered_subtype_dfs[subtype].iloc[:, border:]
            
            #Merge the cluster expression data and pheno data for each cell subtype as a DataFrame
            comb_pheno_data = subtype_df.merge(filt_df_pheno, left_index = True, right_index=True)
            
            # Calculate z-scores for clustering columns and filter out outliers
            clusters = [col for col in comb_pheno_data.columns if 'Cluster_' in col]
            z_scores = comb_pheno_data[clusters].apply(zscore)
            outliers = (z_scores.abs() > 3).any(axis=1)
            filtered_data = comb_pheno_data[~outliers]
            
            filtered_output_file = os.path.join(output_dir, f'{cell_type}_{subtype}_train_gene_clust_metacell_zscorefilt.csv')
            filtered_data.to_csv(filtered_output_file, index=True)
            subtype_clusters[subtype] = filtered_data #store in dict

            #get heatmaps between clusters and pheno for zscore filtered data
            clusters = [col for col in filtered_data.columns if 'Cluster_' in col]
            pheno_cols = ['E2', 'E4', 'new_cogdx', 'age_death', 'msex', 'ad_reagan', 
                          'braaksc', 'ceradsc', 'gpath', 'niareagansc', 'amyloid', 'tangles', 'arteriol_scler', 'caa_4gp', 
                          'cvda_4gp2', 'ci_num2_gct', 'ci_num2_mct', 'sparsity']       

            # Calculate correlation matrix
            gene_pheno_corr = filtered_data[clusters + pheno_cols].corr()
            fig, ax = plt.subplots(figsize=(10, 8))
            heatmap = sns.heatmap(gene_pheno_corr.loc[pheno_cols, clusters], annot=False, center=0, vmin=-0.5, vmax=0.5, cmap='RdBu_r', ax=ax)
            heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=45, horizontalalignment='right', fontsize=8)
            heatmap.set_yticklabels(heatmap.get_yticklabels(), fontsize=10)
            #save fig
            heatmap_filename = os.path.join(output_dir, f"{cell_type}_{subtype}_clusters_pheno_corr_zscorefilt.png")
            fig.savefig(heatmap_filename, bbox_inches='tight')
            plt.close(fig)
            print(f"Finished and saved cluster_genes_main data for {subtype}")

    print(f"Finished cluster_genes_main data for {cell_type}")

    return subtype_clusters, subtype_geneclus_list
