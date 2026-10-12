import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import subprocess
import multiprocessing

sys.path.insert(0, "/proj/nobackup/sens2023026/wharf/gabriv/gabriv-sens2023026/MCcharacterization/Scripts/Pipeline_code_nosexlmer/")

from Steps.i_summarize_files import summarize_lmer_files
from Mains.main_param_pls import run_pls_pipe

plt.switch_backend("Agg")


def run_major_r_and_summary(main_dir, major, k, r_script_lmer_path, r_script_lm_path):
    """For one major cell type: run R (LMER or LM) then summarize."""
    if isinstance(k, int) or (isinstance(k, float) and k.is_integer()):
        print(f"[{major}] Running LMER script")
        subprocess.run(["Rscript", r_script_lmer_path, main_dir, major], check=True)
    elif k == "pseudobulk":
        print(f"[{major}] Running LM script (pseudobulk)")
        subprocess.run(["Rscript", r_script_lm_path, main_dir, major], check=True)
    else:
        raise ValueError(f"Invalid value for k: {k}. Must be an integer or 'pseudobulk'.")

    print(f"[{major}] Summarizing")
    summarize_lmer_files(main_dir=main_dir, cell_type=major)

    return major


def master_pipeline_summary_pls(
    main_dir,
    pheno_path,
    pheno_cols_drop,
    major_celltypes,
    k,
    r_script_lmer_path,
    r_script_lm_path,
    gene_filter_threshold=0.7,
    n_bsiterations=1000,
    gene_knn_neighbours=5,
    n_test_pls_comp=10,
    pls_results_folder="PLS_res_nosexlmer",
    fixed_n_components=None,
    max_workers_inner=4,
):
    """
    Truncated pipeline for one parameter set:
      1) Run R per major cell type (LMER or LM depending on k)
      2) Summarize per major cell type
      3) Run PLS once
    """
    plt.switch_backend("Agg")

    print(f"\nRunning majors in parallel. main_dir={main_dir}")

    effective_workers = min(len(major_celltypes), max_workers_inner)

    with ProcessPoolExecutor(max_workers=effective_workers) as executor:
        future_to_major = {
            executor.submit(
                run_major_r_and_summary,
                main_dir,
                major,
                k,
                r_script_lmer_path,
                r_script_lm_path,
            ): major
            for major in major_celltypes
        }

        for fut in as_completed(future_to_major):
            major = future_to_major[fut]
            # crash if something fails (no try/except)
            fut.result()
            print(f"[{major}] DONE")

    print("\nAll majors done. Now running PLS")
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

    print("Summary and PLS finished for set.")


def run_pipeline_for_paramrow(
    row,
    base_main_dir,
    pheno_path,
    r_script_lmer_path,
    r_script_lm_path,
    pheno_cols_drop,
    pls_results_folder,
    major_celltypes,
    fixed_n_components=None,
    max_workers_inner=4,
):
    set_id = row["SetID"]
    set_main_dir = os.path.join(base_main_dir, set_id)
    os.makedirs(set_main_dir, exist_ok=True)

    # parse k from CSV
    k_raw = row["k"]
    k_val = "pseudobulk" if str(k_raw).lower() == "pseudobulk" else int(k_raw)

    master_pipeline_summary_pls(
        main_dir=set_main_dir,
        pheno_path=pheno_path,
        pheno_cols_drop=pheno_cols_drop,
        major_celltypes=major_celltypes,
        k=k_val,
        r_script_lmer_path=r_script_lmer_path,
        r_script_lm_path=r_script_lm_path,
        pls_results_folder=pls_results_folder,
        fixed_n_components=fixed_n_components,
        max_workers_inner=max_workers_inner,
    )

    return f"Completed rerun (R->summary->PLS) for param set {set_id}."


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", type=int, default=0, help="Which row of the CSV to use (0-based).")
    args = parser.parse_args()

    base_main_dir = "/proj/nobackup/sens2023026/wharf/gabriv/gabriv-sens2023026/MCcharacterization/data/Results/"
    pheno_path = os.path.join(base_main_dir, "other_data", "pheno_filt_touse.xlsx")

    r_script_lmer_path = "/proj/nobackup/sens2023026/wharf/gabriv/gabriv-sens2023026/MCcharacterization/Scripts/Pipeline_code_nosexlmer/Steps/h_train_val_lmer_log_summary_nosex.R"
    r_script_lm_path   = "/proj/nobackup/sens2023026/wharf/gabriv/gabriv-sens2023026/MCcharacterization/Scripts/Pipeline_code_nosexlmer/Steps/h_train_val_lm_log_summary_nosex.R"

    param_csv = os.path.join(base_main_dir, "other_data", "parameters_combos_expanded.csv")

    pheno_cols_drop = ["ad_reagan", "ci_num2_gct", "niareagansc", "braaksc", "ceradsc", "nft"]
    pls_results_folder = "PLS_res_nosexlmer"
    fixed_n_components = None

    major_celltypes = ["Inhibitory", "Excitatory", "Astrocytes", "OPCs",
                       "Oligodendrocytes", "Immune", "Vasculature"]

    df_params = pd.read_csv(param_csv)

    if args.index < 0 or args.index >= len(df_params):
        print(f"Index {args.index} is out of range (0 to {len(df_params)-1}).")
        sys.exit(1)

    row = df_params.iloc[args.index]
    sid = row["SetID"]
    print(f"Selected param set: {sid} (index={args.index})")

    # choose parallelism per set
    num_cores = multiprocessing.cpu_count()
    max_workers_inner = min(7, num_cores)  # 7 majors max

    result_msg = run_pipeline_for_paramrow(
        row=row,
        base_main_dir=base_main_dir,
        pheno_path=pheno_path,
        r_script_lmer_path=r_script_lmer_path,
        r_script_lm_path=r_script_lm_path,
        pheno_cols_drop=pheno_cols_drop,
        pls_results_folder=pls_results_folder,
        major_celltypes=major_celltypes,
        fixed_n_components=fixed_n_components,
        max_workers_inner=max_workers_inner,
    )
    print(result_msg)


if __name__ == "__main__":
    main()
