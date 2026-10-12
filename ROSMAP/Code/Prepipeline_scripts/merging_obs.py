'''
This used to be used in the single cell preprocessing script, but now we add pheno table at a later pipeline stage.
'''

import pandas as pd
import scanpy as sc

def merge_obs(meta_path, adata_path, merge_col):
    '''
    Merge the extra pheno data from csv file to AnnData objects .obs

    Parameters:
    meta_file_path (str): Path to the metadata Excel file.
    adata_file_path (str): Path to the input AnnData file.
    merge_col (str): column on which the meta and adata.obs will be merged on

    Returns:
    AnnData: with added pheno data

    Note: 
    Note: this doesn't save file anywhere
    '''
    adata = sc.read(adata_path)
    meta = pd.read_excel(meta_path)
    #Add cell IDs as a column to the original .obs DataFrame
    adata.obs['cell_id'] = adata.obs_names

    #create new df which would have the original metadata + meta from csv
    newobs = pd.merge(adata.obs, meta, on = merge_col)
    newobs.set_index('cell_id', inplace=True)

    #assign the new metadata to .obs
    adata.obs = newobs

    #creating extra column in .var with the gene names
    adata.var['Genes'] = adata.var_names

    return adata
