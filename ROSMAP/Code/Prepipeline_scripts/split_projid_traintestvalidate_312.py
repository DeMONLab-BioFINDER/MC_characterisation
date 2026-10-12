'''
This script is used to:
1. Take Exc neurons metacells (file that has cells from all donors)
2. Stratify projids into train test validate
3. Save the LISTs of projids for each data set
'''
import pandas as pd
import os
from sklearn.model_selection import train_test_split

main_dir = '/Users/ga0028vi/Documents/MCcharactertization/ROSMAP/Metacells15Jul/Pipeline_build/'
raw_metacell_dir = os.path.join(main_dir, 'Preprocessed_metacell') 

# Load data
df = pd.read_csv(os.path.join(raw_metacell_dir, 'Excitatory_neurons_metacell_expr_preprocessed.csv'), index_col=0)
#get projid, msex cogdx into new dataframe, only unique projid
projid_df = df[['projid', 'msex', 'cogdx']].copy().drop_duplicates(subset='projid')
projid_df

projid_df.info

projid_df['stratify_col'] = projid_df['msex'].astype(str) + projid_df['cogdx'].astype(str)  # create stratification column

# Grouping by the stratification column
grouped = projid_df.groupby('stratify_col')

# Separating groups based on size
small_groups = []
other_groups = []
for name, group in grouped:
    if len(group) <= 3:
        small_groups.append(group)
    else:
        other_groups.append(group)

# Manual assignment for small groups
for group in small_groups:
    group_size = len(group)
    if group_size == 1:
        group['set'] = 'train'  # If only one element, assign it to train
    elif group_size == 2:
        group['set'] = ['train', 'test']
    elif group_size == 3:
        group['set'] = ['train', 'test', 'validate']

# Concatenating small groups back to a DataFrame
small_groups_df = pd.concat(small_groups)

# Handling other groups with train_test_split
train_df = pd.DataFrame()
test_validate_df = pd.DataFrame()
for group in other_groups:
    train_part, test_validate_part = train_test_split(group, test_size=0.66, stratify=group['stratify_col'])
    train_df = pd.concat([train_df, train_part])
    test_validate_df = pd.concat([test_validate_df, test_validate_part])

# Splitting the remaining data into test and validate
test_df, validate_df = train_test_split(test_validate_df, test_size=0.5, stratify=test_validate_df['stratify_col'])

# Combining all splits
train_df = pd.concat([train_df, small_groups_df[small_groups_df['set'] == 'train']])
test_df = pd.concat([test_df, small_groups_df[small_groups_df['set'] == 'test']])
validate_df = pd.concat([validate_df, small_groups_df[small_groups_df['set'] == 'validate']])

# Ensure that indices are reset if necessary
train_df.reset_index(drop=True, inplace=True)
test_df.reset_index(drop=True, inplace=True)
validate_df.reset_index(drop=True, inplace=True)

# Displaying the resulting sets
print("Train Set:", train_df)
print("Test Set:", test_df)
print("Validate Set:", validate_df)


output_dir = '/Users/ga0028vi/Documents/MCcharactertization/ROSMAP/_core_data/Split_33_33_33/'
train_df.to_csv(os.path.join(output_dir, 'ROSMAP_train_projid.csv'))
test_df.to_csv(os.path.join(output_dir, 'ROSMAP_test_projid.csv'))
validate_df.to_csv(os.path.join(output_dir, 'ROSMAP_validate_projid.csv'))