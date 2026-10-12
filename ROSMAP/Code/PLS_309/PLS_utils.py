import pandas as pd
import os
import glob
import numpy as np
from sklearn.impute import KNNImputer
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import zscore
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats



def impute_minmax_genes(train_data, gene_filter_threshold=0.7, knn_neighbors=5):
    """
    Filter, impute, and MinMax-scale gene cluster data.

    Steps:
      1. Drop rows (donors) with too many NaNs.
      2. Drop columns (clusters) with too many NaNs.
      3. KNN-impute remaining NaNs.
      4. MinMax-scale to [0, 1].

    Args:
        train_data (pd.DataFrame): Gene/cluster data (rows = donors, cols = clusters).
        gene_filter_threshold (float): Fraction of non-NaN values required to keep
                                       a row/column (e.g. 0.7 = keep if ≥70% non-NaN).
        knn_neighbors (int): Number of neighbors for KNN imputation.

    Returns:
        train_minmax (pd.DataFrame): Scaled, imputed gene data.
        gene_scaler (MinMaxScaler): Fitted scaler.
        gene_imputer (KNNImputer): Fitted KNN imputer.
    """
    # Filtering by number of NaNs row-wise
    threshold = int(len(train_data.columns) * gene_filter_threshold)
    train_filtered = train_data.dropna(thresh=threshold)

    # Sparsity clusters columnswise 
    train_missing = train_filtered.isna()
    train_num_missing = train_missing.sum()
    train_sparc = train_num_missing / len(train_filtered)

    # Sparsity clusters columns wise
    valid_col_thresh = int(train_filtered.shape[0] * gene_filter_threshold)
    train_filtered = train_filtered.dropna(thresh=valid_col_thresh, axis=1)

    # IMPUTATION OF Genes
    gene_imputer = KNNImputer(n_neighbors=knn_neighbors)
    gene_imputer.fit(train_filtered)
    imputed_array = gene_imputer.transform(train_filtered)

    # convert to pd
    imputed_train = pd.DataFrame(
        imputed_array,
        index=train_filtered.index,
        columns=train_filtered.columns
    )

    # Scaling
    gene_scaler = MinMaxScaler()
    gene_scaler.fit(imputed_train)
    train_minmax = pd.DataFrame(
        gene_scaler.transform(imputed_train),
        columns=imputed_train.columns,
        index=imputed_train.index
    )

    return train_minmax, gene_scaler, gene_imputer

def impute_minmax_pheno(pheno, 
                        imputed_genes, 
                        pheno_cols_drop=None, 
                        rounding_cols=None, 
                        knn_neighbors=5, 
                        showheatmap = True):
    """
    Processes and imputes phenotype data, then scales it. HAVE TO MANUALLY CHANGE ROUNDED UP COLS

    Args:
        pheno_file (str): Path to the phenotype file.
        imputed_genes (pd.DataFrame): DataFrame of imputed gene data to align indices.
        pheno_cols_drop (list): List of columns to drop from the phenotype data (default None).
        knn_neighbors (int): Number of neighbors for KNN imputation (default 5).

    Returns:
        pd.DataFrame: A DataFrame with scaled and imputed phenotype data.
    """
    # Drop unnecessary columns
    if pheno_cols_drop:
        pheno.drop(columns=pheno_cols_drop, inplace=True, errors='ignore')

    # Align indices with gene data
    pheno = pheno[pheno.index.isin(imputed_genes.index)]

    # Impute missing values using KNN
    imputer = KNNImputer(n_neighbors=knn_neighbors)
    imputer = imputer.fit(pheno)
    imputed_pharray = imputer.transform(pheno)

    imputed_pheno = pd.DataFrame(imputed_pharray, index=pheno.index, columns=pheno.columns)
    imputed_pheno['caa_4gp'].unique()
    imputed_pheno['caa_4gp'] = imputed_pheno['caa_4gp'].apply(lambda x: round(x, 0))
    imputed_pheno['tdp_st4'].unique()
    imputed_pheno['tdp_st4'] = imputed_pheno['tdp_st4'].apply(lambda x: round(x, 0))
    

    #Round the cols
    if rounding_cols:
        for col in rounding_cols:
            if col in imputed_pheno.columns:
                imputed_pheno[col] = imputed_pheno[col].apply(lambda x: round(x, 0))

    # Log transform continuous pheno variables
    log_cols = ['tangles', 'amyloid', 'gpath']
    for col in log_cols:
        imputed_pheno[col] = np.log1p(imputed_pheno[col])

    # Scale the phenotype data
    scaler = MinMaxScaler()
    scaler.fit(imputed_pheno)
    pheno_minmax = pd.DataFrame(scaler.transform(imputed_pheno), columns=imputed_pheno.columns, index=imputed_pheno.index)

    if showheatmap:
        plt.figure(figsize=(6,5))
        sns.heatmap(pheno_minmax, cmap="viridis", vmin=0, vmax=1)
        plt.show()
        
    return pheno_minmax, scaler, imputer    


def plot_plsx_load(
    all_results_x,       # allx_bs from bootstrap
    bs_stats_x,          # bsrx
    x_feature_names,     # x_genes_train.columns
    component=0,         # 0 = C1
    invert=False,
    save=False,
    out_dir=None,
    file_name=None,
    pretty_labels=False,
    cluster_label_map=None, 
):
    """
    Gene-loadings bootstrap barplot for one PLS component.
    """

    # Get bootstrap df
    df = all_results_x[component] # (n_boot x n_features)
    n_iterations = df.shape[0]

    # Flatten + assign cluster names
    plotr = pd.concat([
        pd.Series(df.values.flatten()),
        pd.Series(list(x_feature_names) * n_iterations)
    ], axis=1)
    plotr.columns = ["Loading", "Gene_clus"]

    # invert sign if needed
    if invert:
        plotr["Loading"] *= -1

    # OPTIONAL pretty mapping
    if pretty_labels:
        if cluster_label_map is None:
            raise ValueError("pretty_labels=True but cluster_label_map is None.")
        plotr["Gene_clus"] = plotr["Gene_clus"].map(cluster_label_map).fillna(plotr["Gene_clus"])



    # Significance from bsrx
    sigs = [1 if x < 0.05 else 0 for x in bs_stats_x["p"].loc[component].astype(float)]

    # Colour by major type
    celltype_colors = {
        "Astrocytes": "tab:blue",
        "Excitatory": "tab:orange",
        "Immune": "tab:green",
        "Inhibitory": "tab:red",
        "Oligodendrocytes": "saddlebrown",
        "OPC": "tab:purple",
        "OPCs": "tab:purple",
    }
    def major_type(label: str) -> str:
        head = label.split()[0]

        pretty_map = {
            "Ast": "Astrocytes",
            "Exc": "Excitatory",
            "Inh": "Inhibitory",
            "Oli": "Oligodendrocytes",
            "OPC": "OPC",
            "Microglia": "Immune",   # keep palette called w/ Immune
        }
        if head in pretty_map:
            return pretty_map[head]

        # raw label format
        head_raw = label.split("_", 1)[0]
        return "Immune" if head_raw == "Immune" else head_raw

    # enforce feature order so sigs[i] matches the same cluster
    # What label does each p-value correspond to?
    raw_labels = list(x_feature_names)
    if pretty_labels and cluster_label_map is not None:
        raw_labels = [cluster_label_map.get(c, c) for c in raw_labels]

    sig_lookup = {}
    for lab, sig in zip(raw_labels, sigs):
        sig_lookup[lab] = max(sig_lookup.get(lab, 0), sig)

    # Alphabetical order of the labels that are actually plotted
    clus_order = sorted(plotr["Gene_clus"].unique())

    # Colors in that same alphabetical order
    colors = [
        celltype_colors[major_type(c)] if sig_lookup.get(c, 0) else "lightgrey"
        for c in clus_order
    ]

    # plot
    sns.set(font_scale=2)
    sns.set_style("ticks")

    plt.subplots(figsize=(9, 4))
    g = sns.barplot(x="Gene_clus", y="Loading", data=plotr, order=clus_order, errorbar=None)

    for bar, color in zip(g.patches, colors):
        bar.set_color(color)

    # CI lines - dist based on bootstrap count
    for clus in clus_order:
        vals = np.sort(plotr.loc[plotr["Gene_clus"] == clus, "Loading"].values)
        q = int(0.05 * len(vals))
        cil, ciu = vals[q], vals[-q]
        plt.plot([clus, clus], [cil, ciu], c="k", lw=0.5)

    # DECORATIONS
    vlines = [4.5, 9.5, 14.5, 19.5]
    for x in vlines:
        plt.axvline(x=x, color="black", linestyle="--", linewidth=1.2)

    plt.xticks(rotation=90, fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel("Gene Cluster", fontsize=14)
    plt.ylabel(f"C{component+1} Loading", fontsize=14)

    if save:
        if out_dir is None:
            out_dir = os.getcwd()
        os.makedirs(out_dir, exist_ok=True)
        if file_name is None:
            file_name = f"C{component+1}_gene_loadings.png"
        plt.savefig(os.path.join(out_dir, file_name), dpi=300, bbox_inches="tight")
        plt.show()
    else:
        plt.show()


def plot_plsy_load(
    all_results_y,      # ally_bs
    bs_stats_y,         # bsry
    y_feature_names,    # pheno_minmax.columns
    component=0,        # 0 = C1
    invert=False,
    save=False,
    out_dir=None,
    file_name=None,
):
    """
    Y-loadings bootstrap barplot for one PLS component.
    """
    # Extract bootstrap DF for component
    df = all_results_y[component] # (n_boot x n_features)
    n_iterations = df.shape[0]
    n_features = df.shape[1]

    # Flatten + phenotype labels
    plotr = pd.concat([
        pd.Series(df.values.flatten()), # all loadings
        pd.Series(list(y_feature_names) * n_iterations)
    ], axis=1)
    plotr.columns = ['Loading', 'Phenotypes']

    # Invert if needed
    if invert:
        plotr["Loading"] *= -1

    # Clean phenotype labels
    label_mapping = {
        'age_death': 'Age at Death',
        'msex': 'Sex: Men',
        'braaksc': 'Braak score',
        'ceradsc': 'CERAD score',
        'gpath': 'Global pathology',
        'niareagansc': 'NIA-Reagan score',
        'amyloid': 'Amyloid',
        'tangles': 'Tangles',
        'tdp_st4': 'TDP-43 load stages',
        'arteriol_scler': 'Arteriolosclerosis',
        'caa_4gp': 'Cerebral Amyloid Angiopathy (CAA)',
        'cvda_4gp2': 'Cerebral Atherosclerosis Rating',
        'ci_num2_mct': 'Chronic Microinfarcts',
        'E2': 'Apoe E2',
        'E4': 'Apoe E4',
        'new_cogdx': 'Clinical diagnosis'
    }

    plotr['Phenotypes'] = plotr['Phenotypes'].map(label_mapping).fillna(plotr['Phenotypes'])

    # Significance colours: bsry['p'].loc[component]
    sigs = [1 if x < 0.05 else 0 for x in bs_stats_y['p'].loc[component].astype(float)]
    colors = ['tab:red' if x else 'grey' for x in sigs]

    # plot
    plt.figure(figsize=(6, 4))
    g = sns.barplot(x='Phenotypes', y='Loading', data=plotr, errorbar=None)

    # set bar colours
    for bar, color in zip(g.patches, colors):
        bar.set_color(color)

    # CIs from dist[50] / dist[-50]
    for pheno in plotr.Phenotypes.unique():
        vals = np.sort(plotr[plotr.Phenotypes == pheno].Loading.values)
        ci_low = vals[int(0.05 * len(vals))]
        ci_high = vals[-int(0.05 * len(vals))]
        plt.plot([pheno, pheno], [ci_low, ci_high], c='k', lw=0.5)

    # DECORATIONS
    for x in [4.5, 9.5]:
        plt.axvline(x=x, color="black", linestyle="--", linewidth=1.2)
    plt.xticks(rotation=90, fontsize=14)
    plt.xlabel("Phenotypes", fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylabel(f"C{component+1} Loading", fontsize=14)

    if save:
        if out_dir is None:
            out_dir = os.getcwd()
        os.makedirs(out_dir, exist_ok=True)
        if file_name is None:
            file_name = f"C{component+1}_yloadings.png"
        plt.savefig(os.path.join(out_dir, file_name), dpi=300, bbox_inches="tight")
        plt.show()
    else:
        plt.show()

    
def plot_pls_multicCloads(
    allx_bs,
    ally_bs,
    bsrx,
    bsry,
    x_genes_train,
    y_pheno_train,
    components=(1, 2, 3),
    invert=True,
    save=False,
    fig_dir=None,
    file_name=None,
    pretty_labels=True,
    gene_cluster_mapping=None,
    phenotype_mapping=None,
):
    """
    Big 2-col subplot: gene clusters (left) + phenotypes (right) for chosen components.
    """

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    import os

    # OPTIONAL mappings (external)
    if gene_cluster_mapping is None:
        gene_cluster_mapping = {}   # default do nothing
    if phenotype_mapping is None:
        phenotype_mapping = {}      # default do nothing

    # Seaborn style
    sns.set(font_scale=1.8)
    sns.set_style("ticks")

    # infer n_iterations
    first_comp = components[0]
    n_iterations = allx_bs[first_comp].shape[0]

    # Create subplots
    n_rows = len(components)
    fig, axes = plt.subplots(
        nrows=n_rows,
        ncols=2,
        figsize=(25, 6.5 * n_rows),
        sharex=False,
        sharey=False,
        gridspec_kw={"width_ratios": [2.5, 1.5]},
    )
    if n_rows == 1:
        axes = np.array(axes).reshape(1, 2)


    # Helper to draw one barplot
    def make_barplot(ax, data, x_col, y_col, colors, cil_ciu_dict,
                     row_idx, component_label, is_gene, vlines=None):
        sns.barplot(ax=ax, x=x_col, y=y_col, data=data, hue=x_col,
                    palette=colors, errorbar=None)

        # CI lines
        for i, category in enumerate(data[x_col].unique()):
            cil, ciu = cil_ciu_dict[category]
            ax.plot([i, i], [cil, ciu], c="k", lw=1)

        ax.set_ylabel(component_label, fontsize=28)

        if vlines:
            for vline in vlines:
                ax.axvline(x=vline, color="#708090", linestyle="--", lw=1.5)

        if row_idx < n_rows - 1:
            ax.set_xticklabels([])
            ax.set_xlabel("")
        else:
            ax.set_xlabel("Gene Cluster" if is_gene else "Phenotypes", fontsize=28)
            xtick_pos = range(len(data[x_col].unique()))
            ax.set_xticks(xtick_pos)
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=26)

    # Loop over components - rows
    for row_idx, comp in enumerate(components):

        #  GENES (left)
        df_x = allx_bs[comp]
        plotr_genes = pd.concat([
            pd.Series(df_x.values.flatten()),
            pd.Series(list(x_genes_train.columns) * n_iterations)
        ], axis=1)
        plotr_genes.columns = ["Loading", "Gene_clus"]

        if invert:
            plotr_genes["Loading"] *= -1

        if pretty_labels and gene_cluster_mapping:
            plotr_genes["Gene_clus"] = plotr_genes["Gene_clus"].map(
                gene_cluster_mapping
            ).fillna(plotr_genes["Gene_clus"])

        # significance for genes
        sigs_genes = [1 if x < 0.05 else 0 for x in bsrx["p"].loc[comp].astype(float)]

        # align sigs to feature order, not .unique() order
        if pretty_labels and gene_cluster_mapping:
            gene_order = [gene_cluster_mapping.get(c, c) for c in list(x_genes_train.columns)]
        else:
            gene_order = list(x_genes_train.columns)

        colors_genes = dict(zip(gene_order, ["red" if sig else "grey" for sig in sigs_genes]))

        # CI for genes
        cil_ciu_genes = {}
        for cat in gene_order:
            dist = np.sort(plotr_genes.loc[plotr_genes["Gene_clus"] == cat, "Loading"].values)
            q = max(1, int(0.05 * len(dist)))
            cil_ciu_genes[cat] = (dist[q - 1], dist[-q])

        # alphabetical sort
        sorted_gene_clusters = sorted(plotr_genes["Gene_clus"].unique())
        plotr_genes["Gene_clus"] = pd.Categorical(
            plotr_genes["Gene_clus"],
            categories=sorted_gene_clusters,
            ordered=True
        )
        plotr_genes = plotr_genes.sort_values("Gene_clus")

        make_barplot(
            ax=axes[row_idx, 0],
            data=plotr_genes,
            x_col="Gene_clus",
            y_col="Loading",
            colors=colors_genes,
            cil_ciu_dict=cil_ciu_genes,
            row_idx=row_idx,
            component_label=f"C{comp+1} Loading",
            is_gene=True,
            vlines=[4.5, 9.5, 14.5, 19.5, 24.5, 29.5],
        )

        #  PHENOTYPES (right)
        df_y = ally_bs[comp]
        plotr_pheno = pd.concat([
            pd.Series(df_y.values.flatten()),
            pd.Series(list(y_pheno_train.columns) * n_iterations)
        ], axis=1)
        plotr_pheno.columns = ["Loading", "Phenotypes"]

        if invert:
            plotr_pheno["Loading"] *= -1

        if pretty_labels and phenotype_mapping:
            plotr_pheno["Phenotypes"] = plotr_pheno["Phenotypes"].map(
                phenotype_mapping
            ).fillna(plotr_pheno["Phenotypes"])

        sigs_pheno = [1 if x < 0.05 else 0 for x in bsry["p"].loc[comp].astype(float)]

        if pretty_labels and phenotype_mapping:
            pheno_order = [phenotype_mapping.get(c, c) for c in list(y_pheno_train.columns)]
        else:
            pheno_order = list(y_pheno_train.columns)

        colors_pheno = dict(zip(pheno_order, ["red" if sig else "grey" for sig in sigs_pheno]))

        cil_ciu_pheno = {}
        for cat in pheno_order:
            dist = np.sort(plotr_pheno.loc[plotr_pheno["Phenotypes"] == cat, "Loading"].values)
            q = max(1, int(0.05 * len(dist)))
            cil_ciu_pheno[cat] = (dist[q - 1], dist[-q])

        make_barplot(
            ax=axes[row_idx, 1],
            data=plotr_pheno,
            x_col="Phenotypes",
            y_col="Loading",
            colors=colors_pheno,
            cil_ciu_dict=cil_ciu_pheno,
            row_idx=row_idx,
            component_label=f"C{comp+1} Loading",
            is_gene=False,
            vlines=[4.5, 9.5],
        )

    plt.tight_layout()

    if save:
        if fig_dir is None:
            fig_dir = os.getcwd()
        os.makedirs(fig_dir, exist_ok=True)
        if file_name is None:
            file_name = "AllComps_loadings_subplots.png"
        plt.savefig(os.path.join(fig_dir, file_name), dpi=300, bbox_inches="tight")
        plt.show()
    else:
        plt.show()
