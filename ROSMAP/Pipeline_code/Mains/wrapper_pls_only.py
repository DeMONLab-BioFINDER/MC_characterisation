import os
import sys
import argparse
import pandas as pd

# Add the pipeline folder to Python path
sys.path.insert(0, '/proj/nobackup/sens2023026/wharf/gabriv/gabriv-sens2023026/MCcharacterization/Scripts/Pipeline_code/') 

# Import the PLS function
from Mains.main_param_pls import run_pls_pipe


def run_pls_for_paramrow(row, 
                         base_main_dir, 
                         pheno_path, 
                         pheno_cols_drop,
                         pls_results_folder,
                         gene_filter_threshold=0.7, 
                         n_bsiterations=1000, 
                         gene_knn_neighbours=5, 
                         n_test_pls_comp=10,
                         fixed_n_components=None,
                         projid_prune_list=None):
    
    """
    Run the PLS pipeline for a single parameter row.
    """
    
    set_id = row["SetID"]
    set_main_dir = os.path.join(base_main_dir, set_id)
    print(f"\n=== Running PLS only for {set_id} ===")

    run_pls_pipe(
        main_dir=set_main_dir,
        pheno_path=pheno_path,
        pheno_cols_drop=pheno_cols_drop,
        pls_results_folder=pls_results_folder,
        gene_filter_threshold=gene_filter_threshold,
        n_bsiterations=n_bsiterations,
        gene_knn_neighbours=gene_knn_neighbours,
        n_test_pls_comp=n_test_pls_comp,
        fixed_n_components=fixed_n_components,
        projid_prune_list=projid_prune_list
        )

    return f"Completed PLS for {set_id}"


def main():
    
    # Argument parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=int, default=0, help="Which row of the CSV to use (0-based).")
    args = parser.parse_args()

    # Define main paths 
    base_main_dir = "/proj/nobackup/sens2023026/wharf/gabriv/gabriv-sens2023026/MCcharacterization/data/Results/"
    pheno_path = os.path.join(base_main_dir, "other_data", "pheno_filt_touse.xlsx")
    param_csv = "/proj/nobackup/sens2023026/wharf/gabriv/gabriv-sens2023026/MCcharacterization/data/Results/other_data/parameters_combos_expanded.csv"


    # Load pheno file and define filtering
    pheno_df = pd.read_excel(pheno_path)
    male_projids   = pheno_df.loc[pheno_df["msex"] == 1, "projid"].astype(int).tolist()
    female_projids = pheno_df.loc[pheno_df["msex"] == 0, "projid"].astype(int).tolist()


    # Settings
    pheno_cols_drop = ['ad_reagan', 'ci_num2_gct', 'niareagansc', 'braaksc', 'ceradsc', 'nft']

    pls_results_folder = 'PLS_results_adjR2'
    fixed_n_components = None
    print(female_projids)
    print(len(female_projids))
    #projid_prune_list = female_projids
    projid_prune_list = None

    
    # Load parameters file and validate index
    df_params = pd.read_csv(param_csv)
    if args.index < 0 or args.index >= len(df_params):
        print(f"Index {args.index} is out of range (0 to {len(df_params)-1})")
        sys.exit(1)

    row = df_params.iloc[args.index]
    sid = row["SetID"]
    print(f"Selected param set: {sid} (index={args.index})")

    result = run_pls_for_paramrow(row, 
                                  base_main_dir, 
                                  pheno_path, 
                                  pheno_cols_drop, 
                                  pls_results_folder,
                                  fixed_n_components=fixed_n_components,
                                  projid_prune_list=projid_prune_list
                                  )
    print(result)


if __name__ == "__main__":
    main()
