import scanpy as sc
import pandas as pd
import numpy as np
import scipy as sp
import anndata as ad
import os
import matplotlib.pyplot as plt
import gc

plt.switch_backend('Agg')  # Non-interactive backend safe for parallel processing

def preprocess_sc(
    data_file, 
    input_dir, 
    output_dir, 
    hvg_min_mean, 
    hvg_max_mean, 
    hvg_min_disp, 
    normz_target=1e4, 
    save_h5ad=False, 
    save_figures=False, 
    return_output_dir=False
):
    """
    Preprocess a single-cell RNA-seq dataset stored in an h5ad file.

    Args:
        data_file (str): Name of the input h5ad file containing single-cell RNA-seq data.
        input_dir (str): Directory containing raw AnnData files to be processed.
        output_dir (str): Directory where the pipeline results (all) will be stored.
        hvg_min_mean (float): Minimum mean expression threshold for identifying 
                               highly variable genes.
        hvg_max_mean (float): Maximum mean expression threshold for identifying 
                               highly variable genes.
        hvg_min_disp (float): Minimum dispersion threshold for identifying 
                               highly variable genes.
        normz_target (float, optional): Total count target for library size 
                                        normalization. Defaults to 10,000.
        save_h5ad (bool, optional): Whether to save the processed h5ad file. Defaults to True.
        save_figures (bool, optional): Whether to save QC and HVG figures. Defaults to True.
        return_output_dir (bool, optional): Whether to return the output directory path. Defaults to False.

    Returns:
        adata (AnnData): Preprocessed AnnData object for downstream processing.
        base_file_name (str): The base name of the input file (without extensions).
        output_dir (str, optional): The output directory path (only if return_output_dir is True).

    Generates:
        - A preprocessed .h5ad file with filtered and normalized single-cell RNA-seq data 
          (if save_h5ad is True).
        - QC figures as PNG files (if save_figures is True), including:
            * Violin plots for QC metrics (n_genes_by_counts, total_counts, pct_counts_mt).
            * Scatter plot of total counts vs. number of genes by counts.
            * Bar plot of top 20 highest-expressed genes.
            * Highly variable genes plot.

    Notes:
        - The function creates necessary subdirectories within the specified output_dir.
        - Filtering steps include removal of ribosomal and mitochondrial genes as well as MALAT1.
    """
    
    #########################################
    ### Set up directories
    ########################################

    #Define directories used for this function
    preproc_output_dir = os.path.join(output_dir, 'preproc_h5ad') #processed h5ad files to be used for downstream
    data_path = os.path.join(input_dir, data_file)
    base_file_name = os.path.splitext(data_file)[0].split("_")[0]  # e.g., 'Immune' or 'Excitatory'
    
    print(base_file_name)
    os.makedirs(preproc_output_dir, exist_ok=True)    #Create dirs if necessary
    
    # Make fig folder and save QC figures if save_figures = True
    if save_figures:
        qc_dir = os.path.join(preproc_output_dir, 'QC_figs_sc')
        os.makedirs(qc_dir, exist_ok=True)


    #########################################
    ### Read-in and prep the data
    ########################################
    adata = sc.read(data_path)
    
    if base_file_name != "Excitatory":
        adata.X = sp.sparse.csr_matrix(adata.layers['counts']) #Excitatories dont have counts layer
    if "counts" in adata.layers:
        del adata.layers["counts"]
        gc.collect()  # Free memory now
        
    # Convert to sparse if dense
    if not sp.sparse.issparse(adata.X):
        adata.X = sp.sparse.csr_matrix(adata.X)
        
    # Convert only stored values in sparse format
    if sp.sparse.issparse(adata.X):
        adata.X.data = adata.X.data.astype(np.int32)  # Convert only nonzero values
    else:
        adata.X = adata.X.astype(np.int32)  # Dense case  
        
    #########################################
    ### Basic QC metrics
    ########################################
    
    #FILTERING
    sc.pp.filter_cells(adata, min_genes=100) #filter out cells with no genes
    sc.pp.filter_genes(adata, min_cells=3) #filter out genes that were not detected in these cell types

    if save_figures:
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        # ribosomal genes
        adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, log1p=False)

        #QC violin plots
        violin_filename = f'{base_file_name}_QCviolin.png'
        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True
        )
        plt.savefig(os.path.join(qc_dir, violin_filename))
        plt.close()

        #QC scatterplot (with raw counts)
        scatter_filename = f'{base_file_name}_QCscatter.png'
        sc.pl.scatter(adata, "total_counts", "n_genes_by_counts",
                    color="pct_counts_mt")
        plt.savefig(os.path.join(qc_dir, scatter_filename))
        plt.close()

        #Top highest expression genes n=20 plot
        high_exp_genes_filename = f'{base_file_name}_high_exp_genes.png'
        sc.pl.highest_expr_genes(adata, n_top=20)
        plt.savefig(os.path.join(qc_dir, high_exp_genes_filename))
        plt.close()


    #########################################
    ### Filter mito and ribo proteins
    ########################################
    #Find all ribosomal genes (both RPS and RPL) IN THE ADATA.var (not .raw)
    
    ribomito_mask = (
        adata.var_names.str.startswith(('RPS','RPL','MT-'))  # ribo/mito 
        | (adata.var_names == 'MALAT1')
    )
    
    #Count the values in ribomito_genes_mask
    num_true = np.sum(ribomito_mask)
    num_false = ribomito_mask.size - num_true
    print(f"Number of Ribosomal and Mitochondrial genes: {num_true}")
    print(f"Number of other genes: {num_false}")
    
    #Filter out the ribosomal genes
    adata = adata[:, ~ribomito_mask]
    del ribomito_mask  # Remove the mask
    gc.collect()


    #########################################
    ### Normalization and HVGs selection
    ########################################
    #NORMALIZATION
    sc.pp.normalize_total(adata, target_sum = normz_target)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean = hvg_min_mean, max_mean = hvg_max_mean, min_disp = hvg_min_disp)

    # If save_figures = True, save HVG plot
    if save_figures:
        HVGplot_filename = f'{base_file_name}_HVG.png'
        sc.pl.highly_variable_genes(adata)
        plt.savefig(os.path.join(qc_dir, HVGplot_filename))
        plt.close()
        print('adata_HVG_count', adata.var['highly_variable'].value_counts())

    # Filter out genes so only HVGs are left in the dataset
    adata._inplace_subset_var(adata.var['highly_variable'])
    
    #SAVE preprocessed single-cell files
    if save_h5ad:
        output_h5ad_file = os.path.join(preproc_output_dir, f'{base_file_name}_preprocessed.h5ad')
        adata.write_h5ad(output_h5ad_file)
    
    if return_output_dir:
        return adata, base_file_name, output_dir
    return adata, base_file_name
