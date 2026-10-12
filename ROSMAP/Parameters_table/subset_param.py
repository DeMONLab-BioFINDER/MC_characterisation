import os
import pandas as pd

param = pd.read_csv('/Users/ga0028vi/Downloads/Parameters_table/parameters_combos.csv')

#pseudo = param[param['k']=='pseudobulk']
# Subset the DataFrame using proper parentheses for multiple conditions
comboA_params = param[(param['min_disp'] == 0.4) & (param['min_mean'] == 0.05) & (param['max_mean'] == 4)]
comboB_params = param[(param['min_disp'] == 0.5) & (param['min_mean'] == 0.05) & (param['max_mean'] == 4)]
comboC_params = param[(param['min_disp'] == 0.75) & (param['min_mean'] == 0.05) & (param['max_mean'] == 4)]

comboD_params = param[(param['min_disp'] == 0.4) & (param['min_mean'] == 0.025) & (param['max_mean'] == 3.5)]
comboE_params = param[(param['min_disp'] == 0.5) & (param['min_mean'] == 0.025) & (param['max_mean'] == 3.5)]
comboF_params = param[(param['min_disp'] == 0.75) & (param['min_mean'] == 0.025) & (param['max_mean'] == 3.5)]

comboG_params = param[(param['min_disp'] == 0.4) & (param['min_mean'] == 0.0125) & (param['max_mean'] == 3)]
comboH_params = param[(param['min_disp'] == 0.5) & (param['min_mean'] == 0.0125) & (param['max_mean'] == 3)]
comboI_params = param[(param['min_disp'] == 0.75) & (param['min_mean'] == 0.0125) & (param['max_mean'] == 3)]


comboA_params.to_csv('comboA_params.csv', index=False)
comboB_params.to_csv('comboB_params.csv', index=False)
comboC_params.to_csv('comboC_params.csv', index=False)
comboD_params.to_csv('comboD_params.csv', index=False)
comboE_params.to_csv('comboE_params.csv', index=False)
comboF_params.to_csv('comboF_params.csv', index=False)
comboG_params.to_csv('comboG_params.csv', index=False)
comboH_params.to_csv('comboH_params.csv', index=False)
comboI_params.to_csv('comboI_params.csv', index=False)
