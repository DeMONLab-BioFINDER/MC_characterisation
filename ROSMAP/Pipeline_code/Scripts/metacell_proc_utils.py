"""
metacell_process_utils.py

Utility functions for preprocessing, filtering, and analyzing metacell data and gene clusters.

Functions:
- group_data(data, grouping_col, boundary_col): Groups data by a specified column and calculates the boundary column index.
- get_cmats(grouped_data, last_col): Computes correlation matrices for grouped data.
- get_cmats_df(data, border_col): Calculates correlation matrices for a single subtype.
- calculate_silhouette_scores(k_max, cmat): Computes silhouette scores for agglomerative clustering.
- cluster_genes(k, cmat, gene_names): Clusters genes into groups using agglomerative clustering.
- calculate_clusters(df_reduced, gene_cluster_df, cluster_col='Cluster_k'): Calculates the average gene cluster expression for each cell.
- plot_clustermaps(correlation_matrices, cell_type, cmap='RdBu_r', center=0, figsize=(10, 10), save_path=None): Generates clustermaps for correlation matrices.
"""

import numpy as np
import os
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score
import seaborn as sns
import matplotlib.pyplot as plt


def clean_cell_names(adata):
    #Rename all the cell types to remove spaces in between and remove the dashes too
    adata.obs.columns = adata.obs.columns.str.replace(' ', '_')
    adata.obs['cell_type_high_resolution'] = adata.obs['cell_type_high_resolution'].str.replace(' ', '_')
    adata.obs['cell_type_high_resolution'] = adata.obs['cell_type_high_resolution'].str.replace('-', '')
    adata.obs['cell_type_high_resolution'] = adata.obs['cell_type_high_resolution'].str.replace(r'\(', '', regex=True)
    adata.obs['cell_type_high_resolution'] = adata.obs['cell_type_high_resolution'].str.replace(r'\)', '', regex=True)
    adata.obs['cell_type_high_resolution'] = adata.obs['cell_type_high_resolution'].str.replace('/', '')

    #create a new adata.obs column for new subtype groupings
    adata.obs['cell_subtypes'] = adata.obs['cell_type_high_resolution']
    
    return adata

def group_data(data, grouping_col, boundary_col):
    """
    Groups a DataFrame by a specified column and identifies the boundary column index.

    Args:
        data (pandas.DataFrame): The input DataFrame to be grouped.
        grouping_col (str): The column name to group the data by (e.g., 'cell_subtype').
        boundary_col (str): The column indicating the boundary between phenotype and gene data (e.g., 'projid').

    Returns:
        int: The index of the boundary column.
        pandas.core.groupby.DataFrameGroupBy: Grouped data object for further analysis.
    """
    grouped_data = data.groupby(grouping_col)
    last_col = data.columns.get_loc(boundary_col)
    return last_col, grouped_data


def get_cmats(grouped_data, last_col):
    """
    Calculates correlation matrices for all groups in the grouped data.

    Args:
        grouped_data (pandas.core.groupby.DataFrameGroupBy): Grouped data object created using `groupby`.
        last_col (int): The index of the boundary column, indicating where the gene data ends.

    Returns:
        dict: A dictionary where keys are group labels (e.g., cell subtypes) and values are
            correlation matrices (numpy.ndarray) for the gene expression data within each group.
    """
    correlation_matrices = {}

    for cell_group, group in grouped_data:
        # Extract the gene expression data for the current group
        df_gene_data = group.iloc[:, :last_col]
    
        # Compute the correlation matrix
        cmat = np.corrcoef(df_gene_data.T)  # Transpose to correlate columns not rows
    
        # Store the correlation matrix with the corresponding cell_subtype key
        correlation_matrices[cell_group] = cmat

    return correlation_matrices
    #return cmat #or correlation_matrices dict?
 
   
def get_cmats_df(data, border_col):
    """
    Calculates a correlation matrix for a single group, using only gene expression data.

    Args:
        data (pandas.DataFrame): The input DataFrame.
        border_col (str): The column name indicating the boundary between phenotype and gene data
            (e.g., 'projid').

    Returns:
        numpy.ndarray: A correlation matrix computed for the gene expression data.
    """
    df_gene_data = data.loc[:, data.columns[:data.columns.get_loc(border_col)]]
    cmat = np.corrcoef(df_gene_data.T)
    return cmat


def calculate_silhouette_scores(k_max, cmat):
    """
    Calculates silhouette scores for agglomerative clustering for a range of cluster counts.

    Args:
        k_max (int): The maximum number of clusters to evaluate (starts from 2 to `k_max`).
        cmat (numpy.ndarray): The input data, typically a correlation matrix or feature matrix.

    Returns:
        tuple:
            - list: Silhouette scores for each value of k (number of clusters).
            - range: The range of cluster counts evaluated (2 to `k_max`).
    """
    silhouette_scores = []
    k_values = range(2, k_max)
     
    for k in k_values:
        clusterer = AgglomerativeClustering(n_clusters=k)
        cluster_labels = clusterer.fit_predict(cmat)
  
        # Calculate Sil score
        sil_score = silhouette_score(cmat, cluster_labels)
        silhouette_scores.append(sil_score)
        
        # Log the silhouette score
        print(f'For n_clusters = {k}, the silhouette score is {sil_score:.3f}')
    

    
    return silhouette_scores, k_values

def cluster_genes(k, cmat, gene_names):
    """
    Clusters genes into groups using agglomerative clustering.

    Args:
        k (int): The number of clusters to form.
        cmat (numpy.ndarray): A square matrix representing the correlation coefficients between genes.
        gene_names (list of str): A list containing the names of the genes.

    Returns:
        pandas.DataFrame: A DataFrame with two columns:
            - 'Gene': Gene names.
            - 'Cluster_k': Cluster labels corresponding to each gene.
    """
    cluster = AgglomerativeClustering(n_clusters=k, metric='euclidean', 
                              linkage='ward')
    cluster_labels = cluster.fit_predict(cmat)
  
    gene_cluster_df = pd.DataFrame({
    'Gene': gene_names, #this needs to be a data table with specific query of columns up till variable 'last_column'
    'Cluster_k': cluster_labels})
    
    return gene_cluster_df


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

def plot_clustermaps(correlation_matrices, cell_type, cmap='RdBu_r', center=0, figsize=(10, 10), save_path=None):
    """
    Plots clustermaps for gene correlation matrices for each cell subtype.

    Args:
        correlation_matrices (dict): A dictionary of correlation matrices for each cell subtype.
        cell_type (str): The type of cell being analyzed.
        cmap (str, optional): The colormap for the plot. Defaults to 'RdBu_r'.
        center (int, optional): The center value for the colormap. Defaults to 0.
        figsize (tuple, optional): The figure size for the plot. Defaults to (10, 10).
        save_path (str, optional): Path to save the plots. If None, plots are not saved.

    Returns:
        None
    """
    for cell_subtype, cmat in correlation_matrices.items():
        # Create a clustermap for each correlation matrix
        cell_type=cell_type
        sns.clustermap(cmat, cmap=cmap, center=center, figsize=figsize)  #can tack on variable at the front if I decide to mod the fig
        plt.title(f"Clustermap for {cell_subtype}")
        if save_path:
            plt.savefig(os.path.join(save_path, f"{cell_type}_{cell_subtype}_gene_clustermap.png"))
