import os
import sys

# Add the pipeline folder to Python path
sys.path.insert(0, '/proj/nobackup/sens2023026/wharf/gabriv/gabriv-sens2023026/MCcharacterization/Scripts/Pipeline_code/') 

import scanpy as sc
import matplotlib.pyplot as plt
from Steps.a_sc_preproc import preprocess_sc
from Steps.b_get_metacells import get_metacells
from Steps.b_get_pseudobulked import get_pseudobulk
from Steps.c_clean_metacells import preproc_metacell_file
from Steps.d_split_train_test import metacell_split_filter
from Steps.e_get_gene_silscore import calculate_sil_scores
from Steps.f_get_gene_clusters import cluster_genes_main
from Steps.g_validation_prep import val_set_prep
from Steps.i_summarize_files import summarize_lmer_files
import subprocess
from Scripts.cell_mappings import cell_subtype_mappings as cell_mapping

sc.settings.verbosity = 2             # verbosity: errors (0), warnings (1), info (2), hints (3)\n",
plt.switch_backend('Agg')


def run_celltype_pipeline(
    cell_type_file,
    main_dir,
    raw_h5ad_dir,
    pheno_path,
    projid_split_dir,
    r_script_lmer_path,
    r_script_lm_path,
    sparc_filt_thresh,
    hvg_min_mean,
    hvg_max_mean,
    hvg_min_disp,
    min_peak,
    max_peak,
    k,
    normz_target=1e4,
    k_max=26
):
    
    import tracemalloc

    # Start tracing
    tracemalloc.start()
    
    adata, base_file_name = preprocess_sc(data_file=cell_type_file, 
                                        input_dir=raw_h5ad_dir, 
                                        output_dir=main_dir, 
                                        hvg_min_mean = hvg_min_mean, 
                                        hvg_max_mean = hvg_max_mean, 
                                        hvg_min_disp = hvg_min_disp,
                                        normz_target=normz_target)


    if isinstance(k, int) or (isinstance(k, float) and k.is_integer()):
        # Call get_metacells when k is a number
        k = int(k)  # Convert to integer
        print(f"Running metacell generation with k={k}")
        raw_metacells = get_metacells(adata=adata, subtype_map=cell_mapping, main_dir=main_dir, base_file_name=base_file_name, k=k)
    elif k == 'pseudobulk':
        # Call get_pseudobulk when k is 'pseudobulk' string
        print("Running pseudobulk generation")
        raw_metacells = get_pseudobulk(adata=adata, main_dir=main_dir, base_file_name=base_file_name)
    else:
        raise ValueError(f"Invalid value for k: {k}. Must be an integer or 'pseudobulk'.")

    if "Immune" in base_file_name:
        raw_metacells = raw_metacells[~raw_metacells.index.str.endswith(("_Tcells", "_CAMs"))]

    preproc_metacells, major_cell_type = preproc_metacell_file(pheno_path = pheno_path, 
                                                            main_dir=main_dir, 
                                                            k=k,
                                                            cell_type=base_file_name,
                                                            metacell_file=raw_metacells)


    celltype_val_subset, celltype_train_subset = metacell_split_filter(metacells=preproc_metacells, 
                                                                            main_dir=main_dir, 
                                                                            cell_type=major_cell_type,
                                                                            sparc_filt_threshold=sparc_filt_thresh,  
                                                                            projid_split_list_dir=projid_split_dir)

    peaks = calculate_sil_scores(train_metacells=celltype_train_subset, 
                                main_dir=main_dir, 
                                cell_type=major_cell_type,
                                k_max= k_max,
                                min_peak= min_peak,
                                max_peak= max_peak
                                )

    metacell_clus_genes, geneclus_list= cluster_genes_main(main_dir=main_dir, 
                                                        train_metacells=celltype_train_subset, 
                                                        cell_type=major_cell_type, 
                                                        subtype_k=peaks)


    val_set_dict = val_set_prep(val_data=celltype_val_subset, gene_clus_list=geneclus_list, cell_type=major_cell_type, main_dir=main_dir)

    # Run the R script using subprocess for both train and validate sets, output pvals and coef
    if isinstance(k, int) or (isinstance(k, float) and k.is_integer()):
        # Call get_metacells when k is a number
        print(f"Running LMER script")
        subprocess.run(['Rscript', r_script_lmer_path, main_dir, major_cell_type], check=True)  
    elif k == 'pseudobulk':
        # Call get_pseudobulk when k is 'pseudobulk'
        print("Running running LM with pseudobulk")
        subprocess.run(['Rscript', r_script_lm_path, main_dir, major_cell_type], check=True)
    else:
        raise ValueError(f"Invalid value for k: {k}. Must be an integer or 'pseudobulk'.")

    summary = summarize_lmer_files(main_dir=main_dir, cell_type=major_cell_type)
    
    current, peak = tracemalloc.get_traced_memory()
    print(f"Current memory usage: {current / 10**6} MB")
    print(f"Peak memory usage: {peak / 10**6} MB")
    with open(os.path.join(main_dir, f"{base_file_name}_memory_usage.txt"), "w") as file:
        file.write(f"Current memory usage: {current / 10**6:.2f} MB\n")
        file.write(f"Peak memory usage: {peak / 10**6:.2f} MB\n")
    
    tracemalloc.stop()
    
    return summary
