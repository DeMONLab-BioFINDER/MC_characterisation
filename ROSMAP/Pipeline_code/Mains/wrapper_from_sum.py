import os
import sys
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import traceback

# Add the pipeline folder to Python path
sys.path.insert(0, '/proj/nobackup/sens2023026/wharf/gabriv/gabriv-sens2023026/MCcharacterization/Scripts/Pipeline_code/') 

from Steps.i_summarize_files import summarize_lmer_files
from Mains.main_param_pls import run_pls_pipe

plt.switch_backend("Agg")


def master_pipeline_summary_pls(
    main_dir,
    pheno_path,
    pheno_cols_drop,
    major_celltypes,
    # PLS parameters:
    gene_filter_threshold=0.7,
    n_bsiterations=1000,
    gene_knn_neighbours=5,
    n_test_pls_comp=10,
    pls_results_folder="PLS_results",
    fixed_n_components=None,
):
    """
    Truncated pipeline for one parameter set:
      1) Re-run summarization (Summary_files/*_cluster_summary.csv)
      2) Run PLS (reads Summary_files)
    """

    print(f"\n Re-running summarization for majors. main_dir={main_dir}")
    for major in major_celltypes:
        print(f"  -> summarizing {major}")
        summarize_lmer_files(main_dir=main_dir, cell_type=major)

    print("\nSummaries done. Now running PLS analysis")
    run_pls_pipe(
        main_dir=main_dir,
        pheno_path=pheno_path,
        pheno_cols_drop=pheno_cols_drop,
        gene_filter_threshold=gene_filter_threshold,
        n_bsiterations=n_bsiterations,
        gene_knn_neighbours=gene_knn_neighbours,
        n_test_pls_comp=n_test_pls_comp,
        pls_results_folder=pls_results_folder,
        fixed_n_components=fixed_n_components,
    )

    print("=== Summary->PLS finished for this set. ===")


def run_pipeline_for_paramrow(row,
                              base_main_dir,
                              pheno_path,
                              pheno_cols_drop,
                              pls_results_folder,
                              major_celltypes,
                              fixed_n_components=None):
    """
    Same as original helper:
      1) create set_main_dir from SetID
      2) run the truncated pipeline for that set
    """
    set_id = row["SetID"]
    set_main_dir = os.path.join(base_main_dir, set_id)
    os.makedirs(set_main_dir, exist_ok=True)

    master_pipeline_summary_pls(
        main_dir=set_main_dir,
        pheno_path=pheno_path,
        pheno_cols_drop=pheno_cols_drop,
        major_celltypes=major_celltypes,
        pls_results_folder=pls_results_folder,
        fixed_n_components=fixed_n_components,
    )

    return f"Completed rerun (summary->PLS) for param set {set_id}."


def main():
    # Parse command-line argument
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=int, default=0,
                        help="Which row of the CSV to use (0-based).")
    args = parser.parse_args()

    # Hard-coded paths
    base_main_dir = "/proj/nobackup/sens2023026/wharf/gabriv/gabriv-sens2023026/MCcharacterization/data/Results/"
    pheno_path = os.path.join(base_main_dir, "other_data", "pheno_filt_touse.xlsx")
    param_csv = "/proj/nobackup/sens2023026/wharf/gabriv/gabriv-sens2023026/MCcharacterization/data/Results/other_data/parameters_combos_expanded.csv"

    pheno_cols_drop = ['ad_reagan', 'ci_num2_gct', 'niareagansc', 'braaksc', 'ceradsc', 'nft']
    pls_results_folder = 'PLS_results_logt_filtpheno' 
    fixed_n_components = None

    # Major cell types:
    major_celltypes = ["Inhibitory", "Excitatory", "Astrocytes", "OPCs",
                       "Oligodendrocytes", "Immune", "Vasculature"]

    # Load parameter combinations from CSV
    df_params = pd.read_csv(param_csv)

    # Validate index
    if args.index < 0 or args.index >= len(df_params):
        print(f"Index {args.index} is out of range (0 to {len(df_params)-1}).")
        sys.exit(1)

    # Select row
    row = df_params.iloc[args.index]
    sid = row["SetID"]
    print(f"Selected param set: {sid} (index={args.index})")

    # Run truncated pipeline
    result_msg = run_pipeline_for_paramrow(
        row,
        base_main_dir=base_main_dir,
        pheno_path=pheno_path,
        pheno_cols_drop=pheno_cols_drop,
        pls_results_folder=pls_results_folder,
        major_celltypes=major_celltypes,
        fixed_n_components=fixed_n_components,
    )
    print(result_msg)



if __name__ == "__main__":
    main()

