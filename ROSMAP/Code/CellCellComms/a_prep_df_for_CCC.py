import pandas as pd
import os
from glob import glob

# Path to main files
main_dir = '/Users/ga0028vi/Documents/MCcharactertization/ROSMAP/Top_sets/Set309'

pls_clusters_file = os.path.join(main_dir, "_set_data/Clusters_train_PLS_bypval.csv")

gene_cluster_folder = os.path.join(main_dir, 'Gene_clusters')
output_dir = os.path.join(main_dir, 'CCC_genes')
os.makedirs(output_dir, exist_ok=True)    

# Load PLS clusters and get headers
pls_df = pd.read_csv(pls_clusters_file)
cluster_names = pls_df.columns.tolist()

# Parse cluster names
parsed_info = []
for name in cluster_names:
    if '_Cluster_' in name:
        parts = name.split('_Cluster_')
        base = parts[0]
        clust_num = parts[1]
        parsed_info.append((name, base, clust_num))

# Process each unique celltype_subtype base
gene_dict = {} # dict to store results

for full_name, base, cluster_num in parsed_info:
    file_pattern = os.path.join(gene_cluster_folder, f"{base}_gene_clustering_list.csv")
    matching_files = glob(file_pattern)
    
    if not matching_files:
        print(f"File not found for base: {base}")
        continue
    gene_df = pd.read_csv(matching_files[0])
    if 'Cluster_k' not in gene_df.columns or 'Gene' not in gene_df.columns:
        print(f"Skipping file {matching_files[0]} due to missing columns")
        continue
    gene_list = gene_df.loc[gene_df['Cluster_k'] == int(cluster_num), 'Gene'].tolist()
    if not gene_list:
        print(f"No genes found for {full_name} (Cluster {cluster_num})")
    gene_dict[full_name] = gene_list

# Convert to df with nans and save
gene_df_padded = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in gene_dict.items()]))
gene_df_padded.to_csv(os.path.join(output_dir, "PLS_significant_gene_lists_padded.csv"), index=False)
print("Gene lists extracted and saved.")