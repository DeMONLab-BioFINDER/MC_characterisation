import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import scipy as sp
import numpy as np
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA

sys.path.insert(0, '/proj/nobackup/sens2023026/wharf/gabriv/gabriv-sens2023026/MCcharacterization/Scripts/Pipeline_code/') 

from Scripts.metacell_proc_utils import clean_cell_names
import Scripts.cell_mappings as subtype_mapping

def get_metacells(main_dir, subtype_map, adata=None, h5ad_file=None, input_dir=None, k=25, base_file_name=None):
    """
    Generate metacells from preprocessed AnnData objects or `.h5ad` files.

    Args:
        adata (AnnData): (optional) A preprocessed AnnData object. If provided, the function uses this object directly.
        h5ad_path (str): (optional) Path to a preprocessed `.h5ad` file. Required if `adata` is not provided.
        k (int): Number of neighbors for the KNN algorithm to define metacells.
        main_dir (str): (optional) Root directory for outputs (e.g., Metacells_raw/).

    Returns:
        metacells_df (pandas.DataFrame): A DataFrame containing average expression values for each metacell.
    """
    
    #############################################
    ### Set up directories and data
    #############################################
    # adata can be either passed down form previous func as a variable or read-in from provided file
    # Ensure either adata or h5ad_path is provided
    if adata is None and h5ad_file is None:
        raise ValueError("Either `adata` or `h5ad_path` must be provided.")
    
    # If adata is None, read from file and derive base_file_name
    if adata is None:
        h5ad_path = os.path.join(input_dir, h5ad_file)
        adata = sc.read_h5ad(h5ad_path)
        if base_file_name is None:
            base_file_name = os.path.splitext(h5ad_file)[0].split("_")[0]  # e.g., 'Immune' or 'Excitatory'
    
    print(f"Running b_get_metacells for {base_file_name}.")

    # Ensure the directory structure
    metacells_dir = os.path.join(main_dir, 'Raw_metacell_files')
    os.makedirs(metacells_dir, exist_ok=True)    

    #############################################
    ### Clean-up cell names and group if necessary
    #############################################
    adata = clean_cell_names(adata)
    adata.obs['cell_subtypes'] = adata.obs['cell_subtypes'].replace(subtype_map)
    adata.obs['cell_subtypes'] = adata.obs['cell_subtypes'].str.replace('_', '', regex=True)    
    
    min_c = k * 1.5
    
    # collect mean-expression-per-metacell for each donor-by-subtype subset
    metacell_means_list = []

    
    #############################################
    ### Subset data by cell SUBTYPE and PROJID
    #############################################
    print(f'Start metacel clustering for {base_file_name} cells.')
    groups = adata.obs.groupby(["projid", "cell_subtypes"]).groups  # dict: (donor, subtype) -> row indices

    for (donor, subtype), idx in groups.items():
        subset_name = f"{donor}_{subtype}"

        # subset AnnData
        sub = adata[idx].copy()
        n_cells = sub.n_obs

        # metacell label per cell (index = cell barcodes)
        metacell_labels = pd.Series(index=sub.obs_names, dtype="object")

        # Assign metacell labels
        if n_cells < k:
            metacell_labels.loc[:] = f"{subset_name}_jnk_cell_small"

        elif n_cells < min_c:
            metacell_labels.loc[:] = f"{subset_name}_jnk_cell_big"

        else:
            # build expression matrix (dense) with cells as rows
            X = sub.X.toarray() if sp.sparse.issparse(sub.X) else np.asarray(sub.X)
            current_x = pd.DataFrame(X, index=sub.obs_names, columns=sub.var_names)

            knn = NearestNeighbors(n_neighbors=k).fit(current_x)

            ind = 0
            # while enough cells remain to form another "normal" metacell
            while len(current_x) >= min_c:
                
                # PCA cell choosing
                # PCA only to choose the seed (PC1 extreme)
                pca = PCA(n_components=1, random_state=123)
                pc1 = pca.fit_transform(current_x.values).ravel()
                seed_pos = int(np.argmax(np.abs(pc1)))

                # take neighbors of the FIRST cell in current_x
                _, nei_matrix = knn.kneighbors(current_x.iloc[[seed_pos]])
                chosen_cells = current_x.index[nei_matrix[0]]

                metacell_labels.loc[chosen_cells] = f"{subset_name}_{ind}"
                current_x.drop(index=chosen_cells, inplace=True)

                ind += 1
                
                if len(current_x) >= min_c:
                    knn.fit(current_x)

            # leftovers -> jnk buckets
            if len(current_x) < k:
                metacell_labels.loc[current_x.index] = f"{subset_name}_jnk_cell_leftovers_small"
            elif len(current_x) < min_c:
                metacell_labels.loc[current_x.index] = f"{subset_name}_jnk_cell_leftovers_big"

        
        # Compute mean expression per metacell for this subset
        X = sub.X.toarray() if sp.sparse.issparse(sub.X) else np.asarray(sub.X)
        expr_df = pd.DataFrame(X, index=sub.obs_names, columns=sub.var_names)

        mean_expr = expr_df.groupby(metacell_labels).mean()
        metacell_means_list.append(mean_expr)

        # free memory
        del sub

    ##############################################
    ### Concatenate all metacells, calculate sparsity and add plot
    ##############################################
    all_metacells_df = pd.concat(metacell_means_list, axis=0) #concat
    all_metacells_df.index.name = "metacell"

    all_metacells_df["sparsity"] = (all_metacells_df == 0).sum(axis=1) / all_metacells_df.shape[1] #calc sparsity
    
    # Sparsity line plot
    sparsity = all_metacells_df["sparsity"]
    categories = {
    "jnk_cell_small": "jnk_cell_small",
    "jnk_cell_big": "jnk_cell_big",
    "jnk_cell_leftovers_small": "jnk_cell_leftovers_small",
    "jnk_cell_leftovers_big": "jnk_cell_leftovers_big",
    "normal_metacells": r"\d+$",  # index ends with digits
    }
    plt.figure(figsize=(12, 8))
    for label, pattern in categories.items():
        # regex only for normal metacells
        use_regex = pattern.endswith("$")
        filtered = sparsity[sparsity.index.str.contains(pattern, regex=use_regex, na=False)]
        if filtered.empty:
            continue

        filtered = filtered.sort_values(ascending=False)
        plt.plot(
            filtered.values,
            marker="o",
            linestyle="-",
            label=f"{label} (n={len(filtered)})"
        )

    plt.title("Sparsity Across Metacell Categories")
    plt.xlabel("Metacells (sorted within category)")
    plt.ylabel("Sparsity (fraction zeros)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(
    os.path.join(metacells_dir, f"{base_file_name}_meta_sparsity_lines.png"),
    dpi=150
    )
    plt.close()

    ##################################
    ### Save
    ##################################
    out_csv = os.path.join(metacells_dir, f"{base_file_name}_metacell_expr.csv")
    all_metacells_df.to_csv(out_csv)
    
    print(f"Finished running b_get_metacells for {base_file_name}, saved output.")
    
    if h5ad_file is not None:
        return all_metacells_df, base_file_name
    return all_metacells_df
