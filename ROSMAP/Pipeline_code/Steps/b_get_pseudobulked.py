import os
import scanpy as sc
import pandas as pd
import scipy as sp
import numpy as np
from Scripts.metacell_proc_utils import clean_cell_names

def get_pseudobulk(main_dir, adata=None,  h5ad_file=None, input_dir=None, base_file_name=None):
    """
    Pseudobulk = mean expression per (projid, cell_subtype).
    Adds sparsity computed on gene columns only.

    Output CSV:
      Raw_metacell_files/{base_file_name}_pseudobulk_expr.csv

    Index:
      projid_cellsubtype (e.g. 12345_AstrocytesAstD)

    Columns:
      [genes..., projid, cell_subtype, sparsity]
    """
    #############################################
    ### Set up directories and data
    #############################################
    # Adata can be either passed down form previous func as a variable or read-in from provided file
    # Ensure either adata or h5ad_path is provided
    if adata is None and h5ad_file is None:
        raise ValueError("Either `adata` or `h5ad_path` must be provided.")
    
    # If adata is None, read from file and derive base_file_name
    if adata is None:
        h5ad_path = os.path.join(input_dir, h5ad_file)
        adata = sc.read_h5ad(h5ad_path)
        if base_file_name is None:
            base_file_name = os.path.splitext(h5ad_file)[0].split("_")[0]  # e.g., 'Immune' or 'Excitatory'
    
    print(f"Running b_get_pseudobulk for {base_file_name}.")

    # Ensure the directory structure
    metacells_dir = os.path.join(main_dir, 'Raw_metacell_files')
    os.makedirs(metacells_dir, exist_ok=True)    

    #############################################
    ### Clean-up cell names BUT DONT GROUP
    #############################################
    adata = clean_cell_names(adata)
    adata.obs['cell_subtypes'] = adata.obs['cell_subtypes'].str.replace('_', '', regex=True)


    #############################################
    ### Pseudobulk the data
    #############################################
    pseudobulk_means_list = []
    groups = adata.obs.groupby(["projid", "cell_subtypes"]).groups  # (donor, subtype) -> row indices

    for (donor, subtype), idx in groups.items():
        subset_name = f"{donor}_{subtype}"

        # subset AnnData
        sub = adata[idx].copy()

        # build expression df (cells x genes)
        X = sub.X.toarray() if sp.sparse.issparse(sub.X) else np.asarray(sub.X)
        expr_df = pd.DataFrame(X, index=sub.obs_names, columns=sub.var_names)

        # one label for ALL cells in this donor-by-subtype subset
        labels = pd.Series(subset_name, index=sub.obs_names, dtype="object")

        # mean expression for this one pseudobulk
        mean_expr = expr_df.groupby(labels).mean()
        pseudobulk_means_list.append(mean_expr)

        del sub

    ##############################################
    ### Concatenate all pseudobulks, calculate sparsity
    ##############################################
    pseudobulk_df = pd.concat(pseudobulk_means_list, axis=0)
    pseudobulk_df.index.name = "metacell"

    # sparsity across gene columns
    pseudobulk_df["sparsity"] = (pseudobulk_df == 0).sum(axis=1) / pseudobulk_df.shape[1]

    ##############################################
    ### Save
    ##############################################
    out_csv = os.path.join(metacells_dir, f"{base_file_name}_pseudobulk_expr.csv")
    pseudobulk_df.to_csv(out_csv)

    print(f"Finished running b_get_pseudobulk for {base_file_name}, saved output.")

    if h5ad_file is not None:
        return pseudobulk_df, base_file_name
    return pseudobulk_df
