import pandas as pd
import os
import glob
import numpy as np
from statsmodels.stats.multitest import multipletests

def summarize_lmer_files(main_dir, cell_type, sig_thresh = 0.05):
    """
    ** STEP 9 **
    Processes LMER results for each cell subtype within a major cell type, 
    applies FDR correction, and summarizes the results in a CSV file.

    Args:
        main_dir (str): Main directory path where LMER results are stored.
        cell_type (str): The major cell type being analyzed.
        sig_thresh (float, optional): Significance threshold for FDR correction. Default is 0.05.

    Returns:
        pd.DataFrame: A summary DataFrame containing cluster analysis results for all subtypes.
    """
    def create_base_key(filename): #gets type_subtype
        """
        Extracts the subtype base key from a filename.
        
        Args:
            filename (str): The filename from which to extract the subtype.
        
        Returns:
            str: A subtype identifier extracted from the filename.
        """
        return '_'.join(filename.split('_')[:2])

    def read_lmer_files(file_pattern, path):
        """
        Reads LMER result files into a dictionary.

        Args:
            file_pattern (str): The pattern to match files in the directory.
            path (str): The directory path to search for files.
        
        Returns:
            dict: A dictionary of DataFrames containing LMER results, keyed by subtype.
        """
        data_type = "training" if "train" in file_pattern else "validation"

        
        file_paths = glob.glob(os.path.join(path, file_pattern))
        dict = {}
        for file_path in file_paths:
            subtype = create_base_key(os.path.basename(file_path))
            df = pd.read_csv(file_path, index_col=0)
            if df.isnull().all().all():  # Check if all values in the DataFrame are NaN
                print(f"{subtype} ({data_type} set) have no lmer results (all values are NaN).")
                df = df.fillna(1) #replace all NaNs with one only in the sets that are completely empty
            
            elif df.isnull().any().any():
                print(f"{subtype} ({data_type} set) have some NaNs.")
            
            dict[subtype] = df
        return dict

    def add_summary_col(dictionary, threshold, summary_key, data_summary):
        """
        Adds a summary column to the data summary dictionary based on the specified threshold.

        Args:
            dictionary (dict): Dictionary of DataFrames from LMER results.
            threshold (float): The significance threshold.
            summary_key (str): Key for the new column in the data summary.
            data_summary (dict): Dictionary to store the summary results.
        """
        for subtype, df in dictionary.items():
            df = (df <= threshold).astype(int)  # Apply threshold check
            for cluster in df.index:
                key = f"{subtype}_{cluster}"
                if key not in data_summary:
                    data_summary[key] = {}
                data_summary[key][summary_key] = df.loc[cluster].max()
    
    
    train_path = os.path.join(main_dir, 'LMER_res_nosex')
    val_path = os.path.join(main_dir, 'LMER_res_nosex')
    train_pval_dict = read_lmer_files(f'{cell_type}_*_train_lmer_pvals_nosex.csv', train_path)
    train_coef_dict = read_lmer_files(f'{cell_type}_*_train_lmer_coef_nosex.csv', train_path)
    val_pval_dict = read_lmer_files(f'{cell_type}_*_validation_lmer_pvals_nosex.csv', val_path)
    val_coef_dict = read_lmer_files(f'{cell_type}_*_validation_lmer_coef_nosex.csv', val_path)


    # Calculate adj pval for training set
    train_adj_pval_dict = {}
    
    for subtype, pval_df in train_pval_dict.items():
        # Flatten p-values and adjust using FDR
        pvals_flat = pval_df.values.flatten()
        adj_pvals = multipletests(pvals_flat, alpha=0.05, method='fdr_bh')[1]
    
        # Reshape adjusted p-values back into DataFrame
        adjp_df = pd.DataFrame(adj_pvals.reshape(pval_df.shape), index=pval_df.index, columns=pval_df.columns)
        train_adj_pval_dict[subtype] = adjp_df
        # Save adjusted p-values to a new CSV file (replace 'pvals.csv' with 'adj_pvals.csv' in the path)
        adjp_filename = f"{subtype}_train_adj_pvals_nosex.csv"
        adjp_df.to_csv(os.path.join(train_path, adjp_filename))


    # Filter validation p-values based on directionalities and on significant training adj pvals
    val_adj_pval_dict = {}
    for subtype, val_pval_df in val_pval_dict.items():
        train_adj_pval_df = train_adj_pval_dict[subtype]
        train_coef_df = train_coef_dict[subtype]
        val_coef_df = val_coef_dict[subtype]
        

        if not (train_adj_pval_df <= sig_thresh).any().any():
            # If no significant p-values in the training set, print message
            print(f"{subtype}: No significant adjusted p-values in the training set.")
            # Set all relationships in validation to 1
            val_pval_df[:] = 1  # Assign 1 to all validation p values
            fdr_val_df = val_pval_df.copy()
            val_adj_pval_dict[subtype] = fdr_val_df #Keep the 1's in validation adj p val
            
            # Save adjusted validation p-values and filtered p-values
            fdr_val_df.to_csv(os.path.join(val_path, f"{subtype}_validation_adj_pvals_nosex.csv"))
            val_pval_df.to_csv(os.path.join(val_path, f'{subtype}_validation_lmer_pvals_nosex.csv'))
            continue  # Skip to the next subtype, no need for further processing
        
        
        # Create a directionality matrix: 1 if direction matches, 0 otherwise      
        direction_match = (np.sign(train_coef_df) == np.sign(val_coef_df)).astype(int)
        # Combine the conditions: train adj pval < 0.05 and directionality matches (i.e., directionality_matrix == 1)
        sig_and_matching = (train_adj_pval_df <= sig_thresh) & (direction_match == 1)
        
        #Set validation pvalues to 1 if they don't meet the train q-val and direction
        val_pval_df.values[np.where(~sig_and_matching)] = 1
        fdr_val_df = val_pval_df.copy()
        
        if sig_and_matching.any().any():
            fdr = multipletests(val_pval_df.values[np.where(sig_and_matching)], method='fdr_bh')[1]
            fdr_val_df.values[np.where(sig_and_matching)] = fdr
            
        val_adj_pval_dict[subtype] = fdr_val_df
        
        # Save adjusted validation p-values and filtered p-values
        adjp_filename = f"{subtype}_validation_adj_pvals_nosex.csv"
        fdr_val_df.to_csv(os.path.join(val_path, adjp_filename))

        val_pval_filename = f"{subtype}_validation_lmer_pvals_nosex.csv"
        val_pval_df.to_csv(os.path.join(val_path, val_pval_filename))

    
    # Process p-value and adjusted p-value files
    data_summary = {}
    add_summary_col(train_pval_dict, sig_thresh, 'Pass_train_pval', data_summary)
    add_summary_col(train_adj_pval_dict, sig_thresh, 'Pass_train_adj_pval', data_summary)
    add_summary_col(val_pval_dict, sig_thresh, 'Pass_val_pval', data_summary)
    add_summary_col(val_adj_pval_dict, sig_thresh, 'Pass_val_adj_pval', data_summary)
   
    summary_dir = os.path.join(main_dir, 'Summary_files_nosex')
    os.makedirs(summary_dir, exist_ok=True)
    
    summary_df = pd.DataFrame.from_dict(data_summary, orient='index')
    summary_df.to_csv(os.path.join(summary_dir, f'{cell_type}_cluster_summary_nosex.csv'), index=True)

    return summary_df