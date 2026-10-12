"""
This script coordinates:
  1) The single-cell pipeline for multiple .h5ad files in parallel
  2) The PLS pipeline, after single-cell pipelines finish
"""

import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
import matplotlib.pyplot as plt
import pandas as pd
import multiprocessing
import tracemalloc

# Add the pipeline folder to Python path
sys.path.insert(0, '/Users/ga0028vi/Documents/MCcharactertization/ROSMAP/Pipeline_param_code/') 

from Mains.main_param_preproc import run_celltype_pipeline
from Mains.main_param_pls import run_pls_pipe

def master_pipeline(
    main_dir,
    raw_h5ad_dir,
    pheno_path,
    projid_split_dir,
    r_script_path,
    celltype_files,
    pheno_cols_drop,
    # The variable parameters (from CSV):
    sparc_filt_thresh,
    k,                # for metacell clustering
    min_peak,
    max_peak,
    hvg_min_mean,
    hvg_max_mean,
    hvg_min_disp,
    # mostly-constant with defaults:
    k_max=26,         # how many gene clusters max for silhouette score
    normz_target=1e4,
    # PLS parameters:
    gene_filter_threshold=0.7,
    n_bsiterations=1000,
    gene_knn_neighbours=5,
    n_test_pls_comp=10,
    max_workers_inner=4

):
    """
    Main pipeline for one parameter set:
    1) Runs single-cell pipelines in parallel across multiple .h5ad files.
    2) Runs a single PLS pipeline after all single-cell pipelines finish.
    """
    plt.switch_backend('Agg')  # Non-interactive backend


    print(f"\n=== Starting single-cell pipelines in parallel. main_dir={main_dir} ===")
    
    effective_workers = min(len(celltype_files), max_workers_inner)
    with ProcessPoolExecutor(max_workers=effective_workers) as executor:
        future_to_celltype = {}
        for ctf in celltype_files:
            fut = executor.submit(
                run_celltype_pipeline,
                cell_type_file=ctf,
                main_dir=main_dir,
                raw_h5ad_dir=raw_h5ad_dir,
                pheno_path=pheno_path,
                projid_split_dir=projid_split_dir,
                r_script_path=r_script_path,
                # The parameters to pass along:
                sparc_filt_thresh=sparc_filt_thresh,
                hvg_min_mean=hvg_min_mean,
                hvg_max_mean=hvg_max_mean,
                hvg_min_disp=hvg_min_disp,
                normz_target=normz_target,
                k=k,
                k_max=k_max,
                min_peak=min_peak,
                max_peak=max_peak
            )
            future_to_celltype[fut] = ctf
            
        #Collecting results
        for fut in as_completed(future_to_celltype):
            celltype = future_to_celltype[fut]
            try:
                result = fut.result()
                print(f"[{celltype}] pipeline finished. Summary returned: {result}")
            except Exception as exc:
                print(f"[{celltype}] pipeline failed with exception: {exc}")
                


    # After all single-cell pipelines finish, run PLS once
    print("\nAll single-cell pipelines completed. Now running PLS analysis.")
    plt.switch_backend('Agg')

    run_pls_pipe(
        main_dir=main_dir,
        pheno_path=pheno_path,
        pheno_cols_drop=pheno_cols_drop,
        gene_filter_threshold=gene_filter_threshold,
        n_bsiterations=n_bsiterations,
        gene_knn_neighbours=gene_knn_neighbours,
        n_test_pls_comp=n_test_pls_comp
    )
    print("PLS pipeline finished. Entire process complete.")


def run_full_pipeline_for_paramrow(row, 
                                   base_main_dir,
                                   raw_h5ad_dir, 
                                   pheno_path,
                                   projid_split_dir, 
                                   r_script_path,
                                   celltype_files, 
                                   pheno_cols_drop,
                                   max_workers_inner=4
                                   ):
    """
    A helper function that:
      1) Reads the relevant columns from the CSV row
      2) Creates a unique main_dir for this parameter set
      3) Calls master_pipeline
    """
    
    # Create a subfolder for each set, so each run has its own outputs
    set_id = row["SetID"]
    set_main_dir = os.path.join(base_main_dir, set_id)
    os.makedirs(set_main_dir, exist_ok=True)

    # Extract variable parameters from the CSV row
    # (cast them to int or float)
    min_disp_val = float(row["min_disp"])
    k_raw = row["k"]
    # conditional parse
    k_val = "pseudobulk" if str(k_raw).lower() == "pseudobulk" else int(k_raw)
    
    sparc_filt_thresh_val = float(row["sparsity_filt"])
    min_disp_val = float(row["min_disp"])
    min_mean_val = float(row["min_mean"])
    max_mean_val = float(row["max_mean"])
    min_peak_val = int(row["min_peak"])
    max_peak_val = int(row["max_peak"])

    # call pipeline
    master_pipeline(
        main_dir=set_main_dir,
        raw_h5ad_dir=raw_h5ad_dir,
        pheno_path=pheno_path,
        projid_split_dir=projid_split_dir,
        r_script_path=r_script_path,
        celltype_files=celltype_files,
        pheno_cols_drop=pheno_cols_drop,
        # The CSV-provided variables:
        sparc_filt_thresh=sparc_filt_thresh_val,
        k=k_val,
        min_peak=min_peak_val,
        max_peak=max_peak_val,
        hvg_min_mean=min_mean_val,
        hvg_max_mean=max_mean_val,
        hvg_min_disp=min_disp_val,
        # pass concurrency setting for inner tasks:
        max_workers_inner=max_workers_inner
    )

    return f"Completed pipeline for param set {set_id}."


def main():
    """
    Main entry point:
      1) Define constants.
      2) Read param_config.csv.
      3) For each row, run the pipeline with that row's parameters (serially or in parallel).
    """

    # Define paths and lists
    base_main_dir = "/Users/ga0028vi/Documents/MCcharactertization/ROSMAP/test/"
    raw_h5ad_dir = os.path.join(base_main_dir, "raw_h5ad")
    pheno_path = os.path.join(base_main_dir, "other_data", "pheno_filt_touse.xlsx")
    projid_split_dir = os.path.join(base_main_dir, "other_data")
    r_script_path = "/Users/ga0028vi/Documents/MCcharactertization/ROSMAP/Pipeline_code/Steps/h_train_val_lmer_log_summary.R"

    celltype_files = ["Vasculature_cells.h5ad", "Astrocytes.h5ad", 'Inhibitory_neurons.h5ad', ]
    pheno_cols_drop = ['ad_reagan','educ','nft','ci_num2_gct']

    # Read CSV of parameter sets
    param_csv = "/Users/ga0028vi/Documents/MCcharactertization/ROSMAP/test/other_data/parameters_combos.csv"
    df_params = pd.read_csv(param_csv)  

    # figure out how many total CPU cores available
    num_cores = multiprocessing.cpu_count()
    print(f"Detected {num_cores} CPU cores on this machine.")

    # define how many param sets can run in parallel
    max_workers_outer = min(1, num_cores // 2)

    # how many .h5ad tasks in parallel per param set
    # e.g. if we do 2 param sets in parallel, each set might do 4 .h5ad tasks
    max_workers_inner = min(1, num_cores // 2)


    # Parallel param sets
    with ProcessPoolExecutor(max_workers=max_workers_outer) as outer_executor:
        future_to_setid = {}
        for idx, row in df_params.iterrows():
            sid = row["SetID"]
            fut2 = outer_executor.submit(
                run_full_pipeline_for_paramrow,
                row,
                base_main_dir,
                raw_h5ad_dir,
                pheno_path,
                projid_split_dir,
                r_script_path,
                celltype_files,
                pheno_cols_drop,
                max_workers_inner  # pass concurrency setting for the .h5ad tasks
            )
            future_to_setid[fut2] = sid

        for fut2 in as_completed(future_to_setid):
            sid = future_to_setid[fut2]
            try:
                res = fut2.result()
                print(f"Param set {sid} completed. {res}")
            except Exception as e:
                print(f"Param set {sid} failed with error: {e}")

if __name__ == "__main__":
    main()