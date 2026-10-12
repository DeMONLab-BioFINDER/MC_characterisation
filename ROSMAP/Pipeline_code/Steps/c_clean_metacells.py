import os
import re
import pandas as pd
#from Scripts.metacell_proc_utils import clean_cell_names

def replace_characters(val):
    if isinstance(val, str):
        return val.replace('/', '').replace('.', '').replace(' ', '').replace('-', '')
    return val


def preproc_metacell_file(pheno_path, main_dir, k, cell_type = None, filepath=None, metacell_file = None):
    '''
    **STEP 3**
    
    Processes an individual (major cell type) metacell file by cleaning up metacell names, extracting project IDs, cell subtypes,
    and metacell classes from the metacell names, and merging phenotype data into the metacell data.

    Args:
        filename (str): The path to the CSV file containing metacell data.
        pheno_path (str): The path to the Excel file containing phenotype data that needs to be merged with metacell data.
        main_dir (str): The directory path where the output file will be saved.

    Returns:
        tuple: A tuple containing:
        - metacell_file (DataFrame): The processed DataFrame with metacell data including merged phenotype information.
        - cell_type (str): The major cell type extracted from the filename, used in naming the output files.
    '''
    # Ensure either adata or h5ad_path is provided
    if metacell_file is None and filepath is None:
        raise ValueError("Either `metacell_df` or `filepath` must be provided.")
    
    # If adata is None, read from file and derive base_file_name
    if metacell_file is None:
        metacell_file = pd.read_csv(filepath, index_col=0)
        if cell_type is None:
            cell_type = os.path.splitext(os.path.basename(filepath))[0]

    metacell_file.index.name = "metacell"

    #read-in pheno data and merge to the metacell file
    pheno = pd.read_excel(pheno_path)
    print(pheno['projid'].dtype)
    pheno['projid'] = pheno['projid'].astype(str)
    
    
    # Clean up the metacell names using predefined function
    metacell_file = metacell_file.applymap(replace_characters)
    
    
    # Extract metadata based on `k`
    if isinstance(k, int) or (isinstance(k, float) and k.is_integer()):
        k = int(k)  # Convert to integer to ensure consistency
        metacell_file.reset_index(inplace=True)  # Moves 'metacell' from index to a column
        metacell_file[['projid', 'cell_subtype', 'metacell_class']] = metacell_file['metacell'].str.split('_', n=2, expand=True)
        metacell_file = pd.merge(metacell_file, pheno, on='projid', how='inner')
        # Set 'metacell' as the index of the DataFrame
        metacell_file.set_index('metacell', inplace=True)
        
        
    elif k == 'pseudobulk':
        metacell_file.reset_index(inplace=True)  # Moves 'metacell' from index to a column
        metacell_file[['projid', 'cell_subtype']] = metacell_file['metacell'].str.split('_', n=1, expand=True)
        print('Pseudobulk data doesnt require spltting metacell column')

        metacell_file['projid'] = metacell_file['projid'].astype(str)
        metacell_file = pd.merge(metacell_file, pheno, on='projid', how='inner')
        # Set 'metacell' as the index of the DataFrame
        metacell_file.set_index('metacell', inplace=True)
    else:
        raise ValueError(f"Invalid value for k: {k}. Must be an integer or 'pseudobulk'.")
    
    # Define output directory for Preprocessed_metacell files
    preprocessed_dir = os.path.join(main_dir, 'Preprocessed_metacell')
    os.makedirs(preprocessed_dir, exist_ok=True)
    
    output_filename = f"{cell_type}_preproc_metacells.csv"
    output_path = os.path.join(preprocessed_dir, output_filename)
    metacell_file.to_csv(output_path)
    print(f'finished preproc_metacells {cell_type}')
    return metacell_file, cell_type
