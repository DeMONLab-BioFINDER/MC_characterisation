import os
import sys
import matplotlib.pyplot as plt

from Steps.j_merge_sign_clusters import extract_significant_clusters, aggregate_significant_clusters, clean_impute_minmax_gene_data, impute_minmax_pheno
from Steps.k_PLS import pls_choose_comp, eval_pls_mod, pls_mod_comp_loads


def run_pls_pipe (main_dir, 
                  pheno_path, 
                  pheno_cols_drop, 
                  gene_filter_threshold,
                  pls_results_folder,
                  n_bsiterations=1000, 
                  gene_knn_neighbours = 5, 
                  n_test_pls_comp = 10,
                  fixed_n_components=None,
                  projid_prune_list=None):
    
    val_pval_clusters, val_adj_pval_clusters = extract_significant_clusters(main_dir)

    final_data_df = aggregate_significant_clusters(main_dir, val_pval_clusters)
    
    if final_data_df is None:
        msg = "No significant clusters found. Skipping PLS for this parameter set."
        print(msg)
        with open(os.path.join(main_dir, "no_clusters_message.txt"), "w") as f:
            f.write(msg + "\n")
        return  # do not proceed with PLS
    
    if projid_prune_list is not None:
        final_data_df = final_data_df[~final_data_df.index.isin(projid_prune_list)]
        print(f"Pruned {len(projid_prune_list)} projids. Remaining: {final_data_df.shape[0]}")

    imputed_minmax_genes = clean_impute_minmax_gene_data(final_df=final_data_df, 
                                                         filter_threshold=gene_filter_threshold, 
                                                         knn_neighbors=gene_knn_neighbours
                                                         )
    
    
    pls_results_dir = os.path.join(main_dir, pls_results_folder)
    os.makedirs(pls_results_dir, exist_ok=True)
    imputed_minmax_genes.to_csv(os.path.join(pls_results_dir, "Clusters_train_PLS_bypval_proc.csv"))

    scaled_pheno = impute_minmax_pheno(pheno_file=pheno_path, 
                                    imputed_genes=imputed_minmax_genes,
                                    pheno_cols_drop=pheno_cols_drop)

    num_clusters = imputed_minmax_genes.shape[1]  # number of cluster columns
    
    adjusted_ncs = min(n_test_pls_comp, num_clusters)

    if num_clusters < n_test_pls_comp:
        print(f"Warning: You requested up to {n_test_pls_comp} components, "
            f"but only {num_clusters} clusters are available. "
            f"Adjusting down to {num_clusters} components.")
    if num_clusters < 3:
        # If there's fewer gene clusters than 3 - skip
        msg = f"Only {num_clusters} gene clusters found, skipping PLS."
        print(msg)
        with open(os.path.join(main_dir, "too_few_clusters_message.txt"), "w") as f:
            f.write(msg + "\n")
        return # do not proceed with PLS

    # Join folder name to main_dir path
    pls_results_dir = os.path.join(main_dir, pls_results_folder)
    os.makedirs(pls_results_dir, exist_ok=True)


    if fixed_n_components is not None:
        optimal_components = fixed_n_components
        print(f"Using fixed number of PLS components: {fixed_n_components}")
    else:
        optimal_components = pls_choose_comp(gene_df=imputed_minmax_genes,
                                            pheno_df=scaled_pheno,
                                            main_dir=main_dir,
                                            ncs=adjusted_ncs,
                                            pls_results_dir=pls_results_dir)

    bsrx, bsry, allx_bs, ally_bs = eval_pls_mod(gene_df=imputed_minmax_genes,
                                                pheno_df=scaled_pheno,
                                                main_dir=main_dir,
                                                optimal_components=optimal_components,
                                                n_bsiterations=n_bsiterations,
                                                pls_results_dir=pls_results_dir)

    summary_df = pls_mod_comp_loads(main_dir=main_dir,
                                    gene_df=imputed_minmax_genes,
                                    pheno_df=scaled_pheno,
                                    bsrx=bsrx,
                                    bsry=bsry,
                                    allx_bs=allx_bs,
                                    ally_bs=ally_bs,
                                    n_bsiterations=n_bsiterations,
                                    pls_results_dir=pls_results_dir)
