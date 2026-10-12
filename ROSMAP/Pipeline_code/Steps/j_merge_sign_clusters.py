import pandas as pd
import os
import glob
import numpy as np
from sklearn.impute import KNNImputer
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import zscore
import matplotlib
matplotlib.use('Agg')

def extract_significant_clusters(main_dir):
    """
    Extracts significant gene clusters across all cell types from summary files.

    Args:
        summary_dir (str): Directory containing summary CSV files for all major cell types.

    Returns:
        tuple: Two lists:
            - Clusters that passed the validation p-value threshold.
            - Clusters that passed the validation adjusted p-value threshold.
    """
    
    summary_dir = os.path.join(main_dir, 'Summary_files_nosex')
    os.makedirs(summary_dir, exist_ok=True)    

    val_pval_clusters = []
    val_adj_pval_clusters = []

    # Find all summary files
    summary_files = glob.glob(os.path.join(summary_dir, '*_cluster_summary_nosex.csv'))

    for file_path in summary_files:
        # Load the summary file
        df = pd.read_csv(file_path, index_col=0)
        # Skip empty files (no rows or columns)
        if df.empty or df.shape[1] == 0:
            print(f"⚠️ Skipping empty summary file: {file_path}")
            continue

        # Extract clusters passing the validation p-value threshold
        val_pval_clusters.extend(df.index[df['Pass_val_pval'] == 1].tolist())

        # Extract clusters passing the validation adjusted p-value threshold
        val_adj_pval_clusters.extend(df.index[df['Pass_val_adj_pval'] == 1].tolist())


    val_pval_df = pd.DataFrame(val_pval_clusters, columns=["Cluster"])
    val_adj_pval_df = pd.DataFrame(val_adj_pval_clusters, columns=["Cluster"])
    val_pval_df.to_csv(os.path.join(summary_dir, 'all_clusters_pass_pval_nosex.csv'), index=False)
    val_adj_pval_df.to_csv(os.path.join(summary_dir, 'all_clusters_pass_adj_pval_nosex.csv'),index=False)

    return val_pval_clusters, val_adj_pval_clusters




def aggregate_significant_clusters(main_dir, significant_clusters):
    """
    Aggregates significant clusters across training and validation files into a single CSV file.

    Args:
        main_dir (str): Main directory path where the final aggregated file will be saved.
        cluster_folder (str): Directory containing the training and validation cluster files.
        significant_clusters (list): List of significant clusters to aggregate.

    Returns:
        None
    """
    cluster_folder = os.path.join(main_dir, 'Gene_clusters')
    
    final_data = []

    for celltype_cluster in significant_clusters:
        subtype = '_'.join(celltype_cluster.split('_')[:2])
        cluster = '_'.join(celltype_cluster.split('_')[2:])

        # Define paths for train and validation files
        train_file_path = os.path.join(cluster_folder, f'{subtype}_train_gene_clust_metacell_zscorefilt.csv')
        validation_file_path = os.path.join(cluster_folder, f'{subtype}_validation_gene_clust_metacell_zscorefilt.csv')

        # Ensure both files exist
        if not os.path.exists(train_file_path):
            print(f"Train file missing for {subtype}: {train_file_path}")
            continue

        if not os.path.exists(validation_file_path):
            print(f"Validation file missing for {subtype}: {validation_file_path}")
            continue

        # Read files
        train_file = pd.read_csv(train_file_path, index_col=0)
        validation_file = pd.read_csv(validation_file_path, index_col=0)

        # Combine train and validation data
        merged_data = pd.concat([train_file, validation_file], axis=0)

        # Filter to keep only projid and the specific cluster column
        if cluster not in merged_data.columns:
            print(f"Cluster column {cluster} not found in files for {subtype}.")
            continue

        cluster_data = merged_data[['projid', cluster]]

        # Group by projid and calculate the mean for the cluster
        avg_cluster_data = cluster_data.groupby('projid').mean().reset_index()
        avg_cluster_data = avg_cluster_data.rename(columns={cluster: celltype_cluster})  # Rename column

        # Append to the final_data list with projid as the index
        final_data.append(avg_cluster_data.set_index('projid'))

    if final_data:
            final_data_df = pd.concat(final_data, axis=1)
            final_output_path = os.path.join(main_dir, "Clusters_train_PLS_bypval_nosexlmer.csv")
            final_data_df.to_csv(final_output_path, index=True)
            print(f"Aggregated clusters saved to {final_output_path}")
            return final_data_df
    else:
        print("No significant clusters were aggregated.")
        return None

def clean_impute_minmax_gene_data(final_df, filter_threshold, knn_neighbors=5):
    """
    Cleans, filters, and imputes gene expression data.

    Args:
        final_df (pd.DataFrame): Aggregated gene cluster data.
        pheno_file (str): Path to the phenotype file.
        filter_threshold (float): Number (i.e. )
        knn_neighbors (int): Number of neighbors for KNN imputation (default 5).

    Returns:
        pd.DataFrame: A DataFrame with imputed gene data.
    """
    # Filter out genes with too many NaNs
    valid_gene_threshold = int(final_df.shape[1] * filter_threshold)  # donor row filter
    genes_filtered = final_df.dropna(thresh=valid_gene_threshold, axis=0)
    
    valid_col_thresh = int(genes_filtered.shape[0] * filter_threshold)  # cluster column filter
    genes_filtered = genes_filtered.dropna(thresh=valid_col_thresh, axis=1)

    # Impute missing values using KNN
    imputer = KNNImputer(n_neighbors=knn_neighbors)
    imputed_array = imputer.fit_transform(genes_filtered)
    imputed_genes = pd.DataFrame(imputed_array, index=genes_filtered.index, columns=genes_filtered.columns)
    
    scaler = MinMaxScaler()
    imputed_minmax_genes = pd.DataFrame(scaler.fit_transform(imputed_genes), columns=imputed_genes.columns, index=imputed_genes.index)

    return imputed_minmax_genes

def impute_minmax_pheno(pheno_file, imputed_genes, pheno_cols_drop=None, knn_neighbors=5):
    """
    Processes and imputes phenotype data, then scales it.

    Args:
        pheno_file (str): Path to the phenotype file.
        imputed_genes (pd.DataFrame): DataFrame of imputed gene data to align indices.
        pheno_cols_drop (list): List of columns to drop from the phenotype data (default None).
        knn_neighbors (int): Number of neighbors for KNN imputation (default 5).

    Returns:
        pd.DataFrame: A DataFrame with scaled and imputed phenotype data.
    """
    
    pheno = pd.read_excel(pheno_file, index_col='projid')

    # Drop unnecessary columns
    if pheno_cols_drop:
        pheno.drop(columns=pheno_cols_drop, inplace=True, errors='ignore')

    # Align indices with gene data
    pheno = pheno[pheno.index.isin(imputed_genes.index)]

    # Impute missing values using KNN
    imputer = KNNImputer(n_neighbors=knn_neighbors)
    imputed_pharray = imputer.fit_transform(pheno)
    imputed_pheno = pd.DataFrame(imputed_pharray, index=pheno.index, columns=pheno.columns)
    
    
    # Log-transform continuous variables
    log_cols = ['tangles', 'amyloid', 'gpath']
    for col in log_cols:
        imputed_pheno[col] = np.log1p(imputed_pheno[col])

    # Scale the phenotype data
    scaler = MinMaxScaler()
    scaled_pheno_array = scaler.fit_transform(imputed_pheno)
    scaled_pheno = pd.DataFrame(scaled_pheno_array, index=imputed_pheno.index, columns=imputed_pheno.columns)

    return scaled_pheno