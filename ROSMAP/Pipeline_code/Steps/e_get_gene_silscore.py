import os
import contextlib
import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
from scipy.signal import find_peaks
from Scripts.metacell_proc_utils import plot_clustermaps, calculate_silhouette_scores, get_cmats_df

def calculate_sil_scores(train_metacells, main_dir, cell_type, k_max = 26, min_peak= 10, max_peak=20):
    """
    **STEP 5**
    
    Analyzes silhouette scores for various cell subtypes based on correlation matrices,
    identifying optimal clustering configurations and producing visualizations.

    This function:
    1. Groups input metacell data by 'cell_subtype' and filters out non-expressed genes.
    2. Calculates correlation matrices for each subtype.
    3. Computes silhouette scores and identifies peak scores within specified cluster ranges.
    4. Visualizes silhouette scores and peaks for each subtype and saves these plots.

    Parameters:
        train_metacells (DataFrame): DataFrame containing metacell data for a specific cell type.
        main_dir (str): Main directory path where the silhouette score results directory will be created.
        cell_type (str): Type of cell to process which is used for logging and file naming.
        k_max (int): Max number of k used when calculating sil scores (go from 2 to k_max).
        min_peak (int): Optimal peak cluster selection floor.
        max_peak (int): Optimal peak cluster selection ceiling.

    Returns:
    - subtype_peak_scores (dict): A dictionary mapping each cell subtype to its best cluster size (k)
      and the corresponding silhouette score.

    Outputs:
    - Silhouette score plots saved to the output directory for each cell subtype.
    - Log files capturing the processing details for each subtype.
    """

    
    # Make output directory
    sil_dir = os.path.join(main_dir, 'Gene_sil_scores')
    os.makedirs(sil_dir, exist_ok=True)
    
    print(f"Starting calculate_sil_scores for {cell_type}")
    
    #Get column which splits genes vs pheno
    last_col = train_metacells.columns.get_loc('sparsity')
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
    
    # Generate correlation matrices for the new groups
    correlation_matrices = {}
    for subtype, subtype_df in filtered_subtype_dfs.items():
        cmat = get_cmats_df(subtype_df, 'sparsity')
        correlation_matrices[subtype] = cmat

    #plot and save clustermaps of each cell subtype
    plot_clustermaps(correlation_matrices=correlation_matrices, cell_type=cell_type, save_path=sil_dir)
    
    subtype_peak_scores = {}
    log_path = os.path.join(sil_dir, f'{cell_type}_clustering_sil_scores_log.txt')

    
    # Loop through each group in 'group_col' for silhouette score calculation
    for subtype, cmat in correlation_matrices.items():
        with open(log_path, 'a') as f:
            with contextlib.redirect_stdout(f):
                print(f"Processing {cell_type}_{subtype}")
                sil_scores, k_values = calculate_silhouette_scores(k_max, cmat)  # Calculate silhouette scores
        peaks, _ = find_peaks(sil_scores) # find peaks
    
        # Filter peaks to find the highest within the defined k range
        filtered_peaks = [p for p in peaks if min_peak <= k_values[p] <= max_peak]
        if filtered_peaks:
            best_peak_idx = max(filtered_peaks, key=lambda x: sil_scores[x])
            best_peak_k = k_values[best_peak_idx]
            best_peak_score = sil_scores[best_peak_idx]
    
            # Save the best peak information in the dictionary
            subtype_peak_scores[f"{subtype}"] = (best_peak_k, best_peak_score)
        
    
            # Plot silhouette scores with marked peaks
            plt.figure(figsize=(10, 5))
            k_values = range(2, 2 + len(sil_scores))
            plt.scatter(k_values, sil_scores, color='blue')
            plt.plot(k_values, sil_scores, label='Silhouette Scores')
            plt.title(f"Silhouette Scores for {cell_type}_{subtype}")
            plt.xlabel('Number of Clusters k')
            plt.ylabel('Sil score')
            for peak in peaks:
                color = 'green' if peak == best_peak_idx else 'red'
                marker = 'P' if peak == best_peak_idx else 'o'
                plt.scatter(k_values[peak], sil_scores[peak], color=color, s=100, marker=marker)  # Mark the peak
                plt.text(k_values[peak], sil_scores[peak], f'k={k_values[peak]}', fontsize=12, ha='center')
            plt.grid(True)
            plt.savefig(os.path.join(sil_dir, f'{cell_type}_{subtype}_sil_scores.png'))  # Save the plot
        
            # Log peak silhouette scores
            print(f"Peaks for {cell_type}_{subtype}:")
            for peak in peaks:
                print(f'Peak at k={k_values[peak]} with silhouette score={sil_scores[peak]:.3f}')
        
            print(f"Completed calculate_sil_scores for {cell_type}_{subtype}")  # Print to terminal for completion
        
        
        else:
            # Log subtypes with no peaks and skip further processing
            no_peaks_log_path = os.path.join(sil_dir, f'{cell_type}_{subtype}_no_peaks_log.txt')  # Log for subtypes with no peaks

            with open(no_peaks_log_path, 'a') as error_file:
                error_file.write(f"No silhouette score peaks found for {cell_type}_{subtype}\n")
            print(f"Skipped {cell_type}_{subtype} due to no peaks.")
            continue  # Skip to next subtype
    
    print(f'Finished calculate_sil_scores for {cell_type}')
    return subtype_peak_scores
