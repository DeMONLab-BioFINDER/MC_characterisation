import os
import sys
import pandas as pd

params = pd.read_csv("/proj/nobackup/sens2023026/wharf/gabriv/gabriv-sens2023026/MCcharacterization/data/Results/other_data/parameters_combos_expanded.csv", index_col=0)
sets_directory = "/proj/nobackup/sens2023026/wharf/gabriv/gabriv-sens2023026/MCcharacterization/data/Results/" 

# Initialize a list for extracted values
results = []

# Iterate through the Set directories
for folder in os.listdir(sets_directory):
    if folder.startswith("Set"):  # Ensure it's a Set folder
        set_id = folder  # Store SetID
        metrics_file = os.path.join(sets_directory, folder, "PLS_results_logt_filtpheno", "final_pls_model_metrics.csv")
        loadings_file = os.path.join(sets_directory, folder, "PLS_results_logt_filtpheno", "bootstrap_loadings_summary.csv")

        mean_r2 = None
        total_x_feats = None

        # Extract mean_r2 from final_pls_model_metrics.csv
        if os.path.exists(metrics_file):
            try:
                df = pd.read_csv(metrics_file, sep=",")
                mean_r2 = df.loc[0, "mean_r2"]  # Extract first row's mean_r2
            except Exception as e:
                print(f"Error reading {metrics_file}: {e}")

        # Extract total_x_feats from bootstrap_loadings_summary.csv
        if os.path.exists(loadings_file):
            try:
                loadings_df = pd.read_csv(loadings_file, sep=",")
                total_x_feats = loadings_df["total_x_feats"].iloc[0]  # Extract first row's value
            except Exception as e:
                print(f"Error reading {loadings_file}: {e}")

        # Store extracted values
        results.append({"SetID": set_id, "mean_r2": mean_r2, "total_x_feats": total_x_feats})

# Convert extracted values to DataFrame
extracted_df = pd.DataFrame(results)

# Set SetID as index to match params_df
extracted_df.set_index("SetID", inplace=True)

# Merge extracted values into params_df
merged_df = params.merge(extracted_df, left_index=True, right_index=True, how="left")

merged_df.to_csv(os.path.join(sets_directory, 'pipeline_PLS_summary.csv'))
