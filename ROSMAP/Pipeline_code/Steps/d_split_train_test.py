import os
import multiprocessing
import pandas as pd

#to be continued
def metacell_split_filter(metacells, main_dir, cell_type, sparc_filt_threshold, projid_split_list_dir):
    """
    **STEP 4**
    
    Processes metacell data by filtering, splitting by project ID, and adjusting for sparsity.
    
    This function:
    1. Filters out entries tagged as 'small garbage cells'.
    2. Splits the metacells dataset into training, validation, and test subsets based on project IDs.
    3. Applies a sparsity filter to each subset to exclude entries with a sparsity ratio above 0.5.
    4. Saves the new filtered datasets into respective directories within the main directory.
    
    Parameters:
        metacells (DataFrame): The metacell data to process.
        main_dir (str): The main directory path where output subdirectories will be created.
        cell_type (str): The cell type label, used to name output files.
        projid_split_list_dir (str): The directory containing CSV files with project IDs for splitting the dataset.
    
    Returns:
        tuple: A tuple containing:
        - filtered_val_subset (DataFrame): The filtered validation dataset.
        - filtered_train_subset (DataFrame): The filtered training dataset.
        - filtered_test_subset (DataFrame): The filtered test dataset.
        - matrix_train_count (DataFrame): A pivot table with cell counts per donor per cell type for the training set.
    
    Outputs:
    - CSV files saved in specified subdirectories for training, validation, and test datasets, both before and after sparsity filtering.
    """
    # Create output directories if they don't exist
    train_output_dir = os.path.join(main_dir, 'Training_set')
    val_output_dir = os.path.join(main_dir, 'Validation_set')
    test_output_dir = os.path.join(main_dir, 'Test_set')
    os.makedirs(train_output_dir, exist_ok=True)
    os.makedirs(val_output_dir, exist_ok=True)
    os.makedirs(test_output_dir, exist_ok=True)
    
    print(f'start metacell_split_filter {cell_type}')

    # SPLIT BY TRAIN TEST VALIDATE
    #read in the train test validate lists
    trainset_path = os.path.join(projid_split_list_dir, 'ROSMAP_train_projid.csv')
    trainset = pd.read_csv(trainset_path, index_col=0)
    testset_path = os.path.join(projid_split_list_dir, 'ROSMAP_test_projid.csv')
    testset = pd.read_csv(testset_path, index_col=0)
    validate_path = os.path.join(projid_split_list_dir, 'ROSMAP_validate_projid.csv')
    validateset = pd.read_csv(validate_path, index_col=0)

    #extract projid lists
    train_projids = trainset['projid'].tolist()
    val_projids = validateset['projid'].tolist()
    test_projids = testset['projid'].tolist()
    
    #convert to string
    train_projids = [str(projid) for projid in train_projids]
    val_projids = [str(projid) for projid in val_projids]
    test_projids = [str(projid) for projid in test_projids]
    
    metacells.index = metacells.index.astype(str)


    #subset the RNAseq data based on projid lists
    train_subset = metacells[metacells['projid'].isin(train_projids)]
    val_subset = metacells[metacells['projid'].isin(val_projids)]
    test_subset = metacells[metacells['projid'].isin(test_projids)]

    # TEST SET FILTER SPARSITY and save
    #filter out high sparsity metacells
    filtered_test_subset = test_subset[test_subset['sparsity'] <= sparc_filt_threshold]
    filtered_test_subset.to_csv(
        os.path.join(test_output_dir,f'{cell_type}_metacells_test_subset_sparcfilt.csv'))

    # VALIDATION SET FILTER SPARSITY and save
    #filter out high sparsity metacells
    filtered_val_subset = val_subset[val_subset['sparsity'] <= sparc_filt_threshold]
    filtered_val_subset.to_csv(
        os.path.join(val_output_dir,f'{cell_type}_metacells_val_subset_sparcfilt.csv'))
    
    
    # TRAINING SET FILTER SPARSITY and save
    #filter out high sparsity metacells
    filtered_train_subset = train_subset[train_subset['sparsity'] <= sparc_filt_threshold]
    filtered_train_subset.to_csv(
        os.path.join(train_output_dir, f'{cell_type}_metacells_train_subset_sparcfilt.csv'))
    
    
    # Training set cell counts per donor per cell type
    matrix_train_count = filtered_train_subset.pivot_table(index='projid', 
                                               columns='cell_subtype', aggfunc='size', fill_value=0)
    matrix_train_count.to_csv(os.path.join(train_output_dir, f'{cell_type}_filtered_metacell_counts.csv'))

    print(f'finished metacell_split_filter {cell_type}')
    return filtered_val_subset, filtered_train_subset