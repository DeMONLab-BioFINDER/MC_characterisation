import itertools
import pandas as pd
import os
# Define all parameters to test

#HVG selection
min_mean_max_mean_pairs = [
    (0.0125, 3), #current
    (0.025, 3.5),
    (0.05, 4)
    ]
min_disp = [0.5, 0.75, 0.4]

# Metacell aggregation
k = [25, 50, 100, 'pseudobulk']

# Metacell split filter
sparc_filt_thresh = [0.5, 0.4, 0.7]

# Calculate sil_scores
choose_peak_range_pairs = [
    (10, 21), (5, 16)]

#filter_metacells = [
#    [],  # No filtering
#    ['garbage_cell_leftovers_small', 'garbage_cell_small'],  # Current filters you used
#    ['garbage_cell_leftovers_small']  # Another variant
#]
# PLS filtering
#filter_PLS_donors = [0.7, 0.6] #currently have 0.3 as default


# Param dict
param_dict = {
    "min_mean_max_mean": min_mean_max_mean_pairs,
    "min_disp": min_disp,
    "k": k,
    "sparsity_filt":     sparc_filt_thresh,
    #"filter_metacells": filter_metacells,
    "choose_peak_range": choose_peak_range_pairs
    #"filter_PLS_donors": filter_PLS_donors
}

# Get all possible combos
all_combinations = list(itertools.product(*param_dict.values()))

df = pd.DataFrame(all_combinations, columns=param_dict.keys())

# Split min_mean and max_mean into separate cols
df[['min_mean', 'max_mean']] = df['min_mean_max_mean'].apply(pd.Series)
df.drop('min_mean_max_mean', axis=1, inplace=True)

# Split min_mean and max_mean into separate cols
df[['min_peak', 'max_peak']] = df['choose_peak_range'].apply(pd.Series)
df.drop('choose_peak_range', axis=1, inplace=True)

df.insert(0, 'SetID', [f'Set{i+1}' for i in range(len(df))])

output_dir = '/Users/ga0028vi/Documents/MCcharactertization/ROSMAP/Parameters_table'
#output_dir = 'ROSMAP/Parameters_table'
df.to_csv(os.path.join(output_dir, 'parameters_combos.csv'), index=False)