'''
wrapper_single.py
'''

import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
import matplotlib.pyplot as plt
import pandas as pd
import multiprocessing
import tracemalloc
import argparse


# Add the pipeline folder to Python path
sys.path.insert(0, '/proj/nobackup/sens2023026/wharf/gabriv/gabriv-sens2023026/MCcharacterization/Scripts/Pipeline_code/') 

# Import pipeline functions
from Mains.main_param_preproc import run_celltype_pipeline
from Mains.main_param_pls import run_pls_pipe


def master_pipeline(
    main_dir,
    raw_h5ad_dir,
    pheno_path,
    projid_split_dir,
    r_script_lmer_path,
    r_script_lm_path,
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
    max_workers_inner=4,
    pls_results_folder='PLS_results',
    fixed_n_components=None
    ):
    """
    Main pipeline for one parameter set:
      1) Run per-celltype processing in parallel
      2) Run the PLS step across the full dataset
    """
    plt.switch_backend('Agg')  # Non-interactive backend


    print(f"\n===Starting single-cell pipelines in parallel. main_dir={main_dir} ===")
    
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
                r_script_lmer_path=r_script_lmer_path,
                r_script_lm_path=r_script_lm_path,
                # parameters to pass along:
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
    print("\n=== All single-cell pipelines completed. Now running PLS analysis ===")
    plt.switch_backend('Agg')

    run_pls_pipe(
        main_dir=main_dir,
        pheno_path=pheno_path,
        pheno_cols_drop=pheno_cols_drop,
        gene_filter_threshold=gene_filter_threshold,
        n_bsiterations=n_bsiterations,
        gene_knn_neighbours=gene_knn_neighbours,
        n_test_pls_comp=n_test_pls_comp,
        pls_results_folder=pls_results_folder,
        fixed_n_components=fixed_n_components
    )
    print("PLS pipeline finished. Entire process complete.")


def run_full_pipeline_for_paramrow(row, 
                                   base_main_dir,
                                   raw_h5ad_dir, 
                                   pheno_path,
                                   projid_split_dir, 
                                   r_script_lmer_path,
                                   r_script_lm_path,
                                   celltype_files, 
                                   pheno_cols_drop,
                                   pls_results_folder,
                                   max_workers_inner=4,
                                   fixed_n_components=None 
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
    # read raw value from CSV
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
        r_script_lmer_path=r_script_lmer_path,
        r_script_lm_path=r_script_lm_path,
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
        max_workers_inner=max_workers_inner,
        pls_results_folder=pls_results_folder,
        fixed_n_components=fixed_n_components
    )

    return f"Completed pipeline for param set {set_id}."


def main():
    # Parse command-line argument
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=int, default=0,
                        help="Which row of the CSV to use (0-based).")
    args = parser.parse_args()

    # Hard-coded paths
    base_main_dir = "/proj/nobackup/sens2023026/wharf/gabriv/gabriv-sens2023026/MCcharacterization/data/Results/"
    raw_h5ad_dir = '/proj/nobackup/sens2023026/wharf/gabriv/gabriv-sens2023026/MCcharacterization/data/core_data/'
    pheno_path = os.path.join(base_main_dir, "other_data", "pheno_filt_touse.xlsx")
    projid_split_dir = os.path.join(base_main_dir, "other_data")
    r_script_lmer_path = "/proj/nobackup/sens2023026/wharf/gabriv/gabriv-sens2023026/MCcharacterization/Scripts/Pipeline_code/Steps/h_train_val_lmer_log_summary.R"
    r_script_lm_path = "/proj/nobackup/sens2023026/wharf/gabriv/gabriv-sens2023026/MCcharacterization/Scripts/Pipeline_code/Steps/h_train_val_lm_log_summary.R"

    celltype_files = ['Immune_cells.h5ad']
    pheno_cols_drop = ['ad_reagan', 'ci_num2_gct', 'niareagansc', 'braaksc', 'ceradsc', 'nft' ]
    pls_results_folder = 'PLS_results_logt_filtpheno' 
    fixed_n_components = None

    param_csv = "/proj/nobackup/sens2023026/wharf/gabriv/gabriv-sens2023026/MCcharacterization/data/Results/other_data/parameters_combos_expanded.csv"
    
    # Load parameter combinations from CSV 
    df_params = pd.read_csv(param_csv)

    # Validate that 'args.index' is within range
    if args.index < 0 or args.index >= len(df_params):
        print(f"Index {args.index} is out of range (0 to {len(df_params)-1}).")
        sys.exit(1)

    # Select the correct row for this run
    row = df_params.iloc[args.index]
    sid = row["SetID"]
    print(f"Selected param set: {sid} (index={args.index})")

    # Run pipeline for the selected parameter row
    num_cores = multiprocessing.cpu_count()
    max_workers_inner = min(15, num_cores)

    result_msg = run_full_pipeline_for_paramrow(
        row,
        base_main_dir=base_main_dir,
        raw_h5ad_dir=raw_h5ad_dir,
        pheno_path=pheno_path,
        projid_split_dir=projid_split_dir,
        r_script_lmer_path=r_script_lmer_path,
        r_script_lm_path=r_script_lm_path,
        celltype_files=celltype_files,
        pheno_cols_drop=pheno_cols_drop,
        max_workers_inner=max_workers_inner,
        pls_results_folder=pls_results_folder,
        fixed_n_components=fixed_n_components
        )
    print(result_msg)


if __name__ == "__main__":
    main()
