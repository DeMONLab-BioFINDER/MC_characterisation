import os
import multiprocessing
from functools import partial
import sys
#sys.path.insert(0, '/Users/ga0028vi/Documents/MCcharactertization/ROSMAP/Pipeline_param_code/') 

from Steps.a_sc_preproc import preprocess_sc
from Steps.a_get_subtype_map import cluster_subtypes

# Define input and output directories
raw_h5ad_dir = '/Users/ga0028vi/Documents/MCcharactertization/ROSMAP/test/raw_h5ad/'
output_dir = '/Users/ga0028vi/Documents/MCcharactertization/ROSMAP/test/processed_output/'
# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Gather files to process
files_to_process = [f for f in os.listdir(raw_h5ad_dir) if f.endswith('.h5ad')]

# Parameters for the preprocessing and clustering
hvg_min_mean = 0.0125 
hvg_max_mean = 3
hvg_min_disp = 0.5
normz_target = 1e4
cell_count_threshold = 20000

def sc_proc_clus_subtypes(data_file, 
                          raw_h5ad_dir, 
                          output_dir, 
                          hvg_min_mean, hvg_max_mean, hvg_min_disp, normz_target, 
                          cell_count_threshold):

    """
    Preprocesses and clusters subtypes for a single h5ad file.
    """
    try:
        # Preprocess
        adata, base_file_name, preproc_output_dir = preprocess_sc(
            data_file=data_file,
            input_dir=raw_h5ad_dir,
            output_dir=output_dir,
            hvg_min_mean=hvg_min_mean,
            hvg_max_mean=hvg_max_mean,
            hvg_min_disp=hvg_min_disp,
            normz_target=normz_target,
            save_h5ad=False,
            save_figures=False,
            return_output_dir=True
        )

        # Cluster subtypes
        cluster_subtypes(
            adata=adata,
            base_file_name=base_file_name,
            preproc_output_dir=preproc_output_dir,
            cell_count_threshold=cell_count_threshold
        )
        print(f"Completed processing for: {data_file}")

    except Exception as e:
        print(f"Error processing file {data_file}: {e}")

if __name__ == '__main__':
    # Ensure multiprocessing compatibility
    num_processes = multiprocessing.cpu_count()

    # Use partial to prefill function arguments
    process_partial = partial(
        sc_proc_clus_subtypes,
        raw_h5ad_dir=raw_h5ad_dir,
        output_dir=output_dir,
        hvg_min_mean=hvg_min_mean,
        hvg_max_mean=hvg_max_mean,
        hvg_min_disp=hvg_min_disp,
        normz_target=normz_target,
        cell_count_threshold=cell_count_threshold
    )

    # Parallelize across all available cores
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.map(process_partial, files_to_process)