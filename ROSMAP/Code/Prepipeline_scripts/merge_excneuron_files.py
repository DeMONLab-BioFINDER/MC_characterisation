'''
Quick script to merge the Excitatory neurons raw h5ad files before the preprocessing.
'''
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np

def merge():
    #Read in each file and store the AnnData objects in a list
    files_to_merge = ['D:/Data_Projects/MCcharacterisation/ROSMAP/_core_data/Excitatory_neurons_set1.h5ad', 
                    'D:/Data_Projects/MCcharacterisation/ROSMAP/_core_data/Excitatory_neurons_set2.h5ad', 
                    'D:/Data_Projects/MCcharacterisation/ROSMAP/_core_data/Excitatory_neurons_set3.h5ad']
    adata_list = [sc.read_h5ad(file) for file in files_to_merge]
    # Concatenate the AnnData objects
    merged_adata = ad.concat(adata_list, join='outer', merge='unique')
    # Define the filename for the merged file
    merged_filename = 'D:/Data_Projects/MCcharacterisation/ROSMAP/_core_data/Excitatory_neurons.h5ad'
    #Save the merged AnnData object to a new .h5ad file
    merged_adata.write(merged_filename)
    print(f'Saved merged data to {merged_filename}')

if __name__ == '__main__':
    merge()