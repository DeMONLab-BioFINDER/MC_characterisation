import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def cluster_subtypes(adata, base_file_name, preproc_output_dir, cell_count_threshold=30000):
    """
    Perform hierarchical clustering on subtype-averaged gene expression profiles.
    Highlights subtypes with fewer than 'cell_count_threshold' cells by 
    coloring their labels in red.

    Args:
        adata (AnnData): Preprocessed AnnData object with cell x gene matrix.
        base_file_name (str): The base name of the file/cell type for output naming.
        main_dir (str, optional): Directory where output figures will be stored (not mandatory if preproc_output_dir is provided).
        preproc_output_dir (str, optional): Directory where the preprocessed files and figures are stored.
        cell_count_threshold (int, optional): Threshold for highlighting low-cell subtypes. Default is 30000.

    Returns:
        None: Saves a figure of the clustermap in `hierarch_cell_map` directory.
    """

    cluster_dir = os.path.join(preproc_output_dir, 'hierarch_cell_map')
    os.makedirs(cluster_dir, exist_ok=True)

    if 'cell_subtype' not in adata.obs:
        raise ValueError(f"The AnnData object for '{base_file_name}' is missing the 'cell_subtype' column in adata.obs")
    
    # Compute cell counts per subtype
    subtype_counts = adata.obs['cell_subtype'].value_counts()

    # Aggregate mean expression per subtype
    expr_df = pd.DataFrame(
        adata.X.A if hasattr(adata.X, 'A') else adata.X,
        index=adata.obs_names,
        columns=adata.var_names
    )
    subtype_avg = expr_df.groupby(adata.obs['cell_subtype']).mean()

    # Create the clustermap
    cg = sns.clustermap(
        subtype_avg,
        method='ward', 
        metric='euclidean', 
        cmap='viridis',
        figsize=(10, 12),
        dendrogram_ratio=(0.2, 0.2),
        yticklabels=True
    )

    # Update y-axis labels
    yticklabels = cg.ax_heatmap.get_yticklabels()
    subtype_count_dict = subtype_counts.to_dict()

    for label in yticklabels:
        subtype_name = label.get_text()
        count = subtype_count_dict.get(subtype_name, 0)
        new_text = f"{subtype_name} ({count})"
        label.set_text(new_text)
        if count < cell_count_threshold:
            label.set_color('red')

    cg.ax_heatmap.set_yticklabels(yticklabels)
    cg.ax_heatmap.set_title(f'Hierarchical Clustering of {base_file_name} Subtypes')

    fig_path = os.path.join(cluster_dir, f'{base_file_name}_subtype_clustermap.png')
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()