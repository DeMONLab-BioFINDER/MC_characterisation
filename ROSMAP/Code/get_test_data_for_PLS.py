import os
import pandas as pd
import numpy as np

main_dir = '/Users/ga0028vi/Documents/MCcharactertization/ROSMAP/Top_sets/Set1'
cluster_folder = os.path.join(main_dir, 'Gene_clusters')

# Load the existing train+val cluster file to extract columns
train_val_df = pd.read_csv(os.path.join(main_dir, "Clusters_train_PLS_bypval_nosexlmer.csv"), index_col=0)
cluster_columns = train_val_df.columns.tolist()

final_test_data = []

for col in cluster_columns:
    # Parse cluster column name: e.g., Astrocytes_AstCHI3L1_Cluster_6
    parts = col.split('_')
    major_type = parts[0]
    subtype = parts[1]
    cluster = '_'.join(parts[2:])  # e.g., Cluster_6
    subtype_name = f"{major_type}_{subtype}"

    test_file_path = os.path.join(cluster_folder, f'{subtype_name}_test_gene_clust_metacell_zscorefilt.csv')
    
    if not os.path.exists(test_file_path):
        print(f"[Missing File] Skipping: {test_file_path}")
        continue

    test_df = pd.read_csv(test_file_path, index_col=0)

    if cluster not in test_df.columns:
        print(f"[Missing Cluster] {cluster} not in {test_file_path}")
        continue

    # Group and average per donor
    cluster_data = test_df[['projid', cluster]]
    avg_cluster_data = cluster_data.groupby('projid').mean().reset_index()
    avg_cluster_data = avg_cluster_data.rename(columns={cluster: col})
    final_test_data.append(avg_cluster_data.set_index('projid'))

# Merge all into one aligned test matrix
final_test_df = pd.concat(final_test_data, axis=1)

# Reindex to match column order exactly
final_test_df = final_test_df[cluster_columns]

# Save output
final_test_df.to_csv(os.path.join(main_dir, "Clusters_test_PLS_bypval_nosexlmer.csv"), index=True)

print("Test set saved with aligned cluster columns.")
