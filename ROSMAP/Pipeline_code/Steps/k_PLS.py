import os
import sys
sys.path.insert(0, '/Users/ga0028vi/Documents/MCcharactertization/ROSMAP/Pipeline_param_code/') 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests
from sklearn.cross_decomposition import PLSRegression
from Scripts.Jakesdecomp import evaluate_components, bootstrap_features
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.metrics import r2_score, mean_absolute_error
from sklearn.preprocessing import MinMaxScaler
from matplotlib.colors import Normalize
from scipy import stats
import numpy as np
import matplotlib
matplotlib.use('Agg')



def pls_choose_comp (gene_df, pheno_df, main_dir, pls_results_dir, ncs=10, save_nc_chooser=True):
    """
    Runs Partial Least Squares (PLS) analysis and generates results and visualizations.

    Args:
        gene_df (pd.DataFrame): Scaled gene data.
        pheno_df (pd.DataFrame): Scaled phenotype data.
        main_dir (str): Main directory where results will be saved.
        n_components (int): Maximum number of components to test (default 10).
        save_nc_chooser (bool): Save nc_chooser csv output in PLS output dir or not

    Returns:
        optimal_components (int): optimal number of components to be used in a model. Selected with Adjusted CV R2
    """
    
    #make sure the data indexes are alligned
    x_genes = gene_df.copy()
    y_pheno = pheno_df.copy()
    x_genes, y_pheno = x_genes.align(y_pheno, join='inner', axis=0)

    n = x_genes.shape[0]  # number of samples

    # Create results DataFrame
    nc_chooser = pd.DataFrame(index=range(ncs * 10), columns=['nc', 'r2', 'adj_r2', 'mae', 'iter'])
    i = 0
    
    # Loop over component counts for PLSR only
    for nc in range(1, ncs + 1):
        print(f"Testing {nc} components...")
        
        estimator = PLSRegression(n_components=nc)
        
        for it in range(10):  # 10 iterations of cross-validation
            cv = KFold(5, shuffle=True, random_state=123 + it)
            
            # cross-validation prediction
            pred = cross_val_predict(estimator, x_genes, y_pheno, cv=cv)
            
            # calculate normal CV R2
            r2 = r2_score(y_true=y_pheno, y_pred=pred)
            # calculate adjusted CV R2
            adj_r2 = 1 - ((1 - r2) * (n - 1) / (n - nc - 1))
            
            # Store performance metrics
            nc_chooser.loc[i, 'r2'] = r2
            nc_chooser.loc[i, 'adj_r2'] = adj_r2
            nc_chooser.loc[i, 'mae'] = mean_absolute_error(y_true=y_pheno, y_pred=pred)
            nc_chooser.loc[i, 'nc'] = nc
            nc_chooser.loc[i, 'iter'] = it
            i += 1
    
    # Plot MAE by number of components
    sns.catplot(x='nc', y='mae', data=nc_chooser, aspect=1, kind='point')
    plt.ylabel('Mean Absolute Error')
    plt.xlabel('# of Components')
    plt.savefig(os.path.join(pls_results_dir, 'mae_by_components.png'))


    # Plot R2 by number of components
    sns.catplot(x='nc', y='r2', data=nc_chooser, aspect=1, kind='point')
    plt.ylabel(r'$R^2$')
    plt.xlabel('# of Components')
    plt.savefig(os.path.join(pls_results_dir, 'r2_by_components.png'))
    
    # Plot R2 by number of components
    sns.catplot(x='nc', y='adj_r2', data=nc_chooser, aspect=1, kind='point')
    plt.ylabel(r'Adjusted $R^2$')
    plt.xlabel('# of Components')
    plt.savefig(os.path.join(pls_results_dir, 'adjr2_by_components.png'))
    

    optimal_components = nc_chooser.groupby('nc')['adj_r2'].mean().idxmax()
    
    sub = nc_chooser[nc_chooser['nc'] == optimal_components]
    mean_r2 = sub['r2'].mean()
    mean_adjr2 = sub['adj_r2'].mean()
    mean_mae = sub['mae'].mean()
    
    final_metrics = pd.DataFrame({
        'optimal_components': [optimal_components],
        'mean_r2': [mean_r2],
        'mean_adjr2': [mean_adjr2],
        'mean_mae': [mean_mae]
    })
    final_metrics.to_csv(os.path.join(pls_results_dir, 'final_pls_model_metrics.csv'), index=False)
    
    if save_nc_chooser == True:
        nc_chooser.to_csv(os.path.join(pls_results_dir, 'nc_chooser_model_metrics.csv'), index=False)
    return optimal_components

def eval_pls_mod(gene_df, pheno_df, main_dir, optimal_components, pls_results_dir, n_bsiterations = 1000):
    x_genes = gene_df.copy()
    y_pheno = pheno_df.copy()
    x_genes, y_pheno = x_genes.align(y_pheno, join='inner', axis=0)
    
    pls_mod = PLSRegression(n_components=optimal_components)
    pls_mod.fit(x_genes,y_pheno)

    
    # Model evaluation
    results, comp_sig, fig = evaluate_components(pls_mod, x_genes, y_pheno)
    if fig is not None:
        fig.savefig(os.path.join(pls_results_dir, 'pls_evaluation_plot.png'))
    comp_sig.to_csv(os.path.join(pls_results_dir, 'component_significance.csv'), index=False)
    
    
    # Bootstrapping
    x_genes = gene_df.copy()
    y_pheno = pheno_df.copy()
    x_genes, y_pheno = x_genes.align(y_pheno, join='inner', axis=0)
    
    # The optimal model needs to be rerun before bootstrapping
    pls_mod = PLSRegression(n_components=optimal_components)
    pls_mod.fit(x_genes,y_pheno)
    bsrx, bsry, allx_bs, ally_bs = bootstrap_features(PLSRegression(n_components=optimal_components), 
                                                    pls_mod, 
                                                    x_genes, 
                                                    y_pheno, 
                                                    n_iterations=n_bsiterations,
                                                    mode='rotations')
    return bsrx, bsry, allx_bs, ally_bs


def pls_mod_comp_loads(
    main_dir, 
    gene_df, 
    pheno_df, 
    bsrx, 
    bsry, 
    allx_bs, 
    ally_bs, 
    pls_results_dir,
    n_bsiterations=1000
):
    """
    Explore each PLS component's loadings in X and Y via bootstrap data.
    For each component, compute how many features are 'significant' in X and Y
    using uncorrected p < 0.05. Also plot barplots for both sides.

    Args:
        main_dir (str): Main directory to save outputs (if plotting).
        gene_df (pd.DataFrame): The scaled X data used in the PLS model.
        pheno_df (pd.DataFrame): The scaled Y data used in the PLS model.
        bsrx (dict): Bootstrap ratio results (like 'mean', 'se', 'bsr', 'p') for X side.
        bsry (dict): Same as bsrx but for Y side.
        allx_bs (dict): Dict of DataFrames (n_bsiterations x n_features) for each component in X.
        ally_bs (dict): Dict of DataFrames for Y side.
        n_bsiterations (int): Number of bootstrap iterations used.

    Returns:
        summary_df (pd.DataFrame): Table summarizing how many X and Y loadings
                                   are significant (p<0.05, uncorrected) per component.
    """
    
    # Align data
    x_genes = gene_df.copy()
    y_pheno = pheno_df.copy()
    x_genes, y_pheno = x_genes.align(y_pheno, join='inner', axis=0)

    summary_df = pd.DataFrame(
        columns=['component', 'num_sig_x', 'num_sig_y', 
                 'total_x_feats', 'total_y_feats']
    )

    # allx_bs and ally_bs have the same set of component keys
    n_components = len(allx_bs.keys())
    total_x_feats = x_genes.shape[1]
    total_y_feats = y_pheno.shape[1]

    # Loop over each component
    for c in range(n_components):
        #  X side
        x_comp_df = pd.concat([
            pd.Series(allx_bs[c].values.flatten(), name='Loading'),
            pd.Series(list(x_genes.columns) * n_bsiterations, name='Feature')
        ], axis=1)
        
        # Uncorrected p < 0.05
        pvals_x = bsrx['p'].loc[c].astype(float)
        sigs_x_bool = [1 if p < 0.05 else 0 for p in pvals_x]
        num_sig_x = sum(sigs_x_bool)
        
        # Plot X side barplot
        plt.figure(figsize=(12, 8))
        sns.set(font_scale=1.3)
        sns.set_style('ticks')

        colors_x = ['Red' if is_sig else 'Grey' for is_sig in sigs_x_bool]
        sns.barplot(x='Feature', y='Loading', data=x_comp_df, palette=colors_x, ci=None)
        
        # 5th/95th percentile lines
        for i_feat, feat_name in enumerate(x_genes.columns):
            dist = x_comp_df.loc[x_comp_df.Feature == feat_name, 
                                 'Loading'].sort_values().values
            cil, ciu = dist[int(0.05 * len(dist))], dist[int(0.95 * len(dist))]
            plt.plot([i_feat, i_feat], [cil, ciu], c='k', lw=1)
        
        plt.xticks(rotation=90)
        plt.ylabel(f'C{c+1} Loading')
        plt.title(f'X side: Component {c+1}')
        plt.tight_layout()
        plt.savefig(os.path.join(pls_results_dir, f'C{c+1}_X_loadings.png'), dpi=300)
        plt.close()

        # Y side
        y_comp_df = pd.concat([
            pd.Series(ally_bs[c].values.flatten(), name='Loading'),
            pd.Series(list(y_pheno.columns) * n_bsiterations, name='Feature')
        ], axis=1)
        
        pvals_y = bsry['p'].loc[c].astype(float)
        sigs_y_bool = [1 if p < 0.05 else 0 for p in pvals_y]
        num_sig_y = sum(sigs_y_bool)

        # Plot Y side barplot
        plt.figure(figsize=(12, 8))
        sns.set(font_scale=1.3)
        sns.set_style('ticks')

        colors_y = ['Red' if is_sig else 'Grey' for is_sig in sigs_y_bool]
        sns.barplot(x='Feature', y='Loading', data=y_comp_df, palette=colors_y, ci=None)
        
        # 5th/95th percentile lines
        for i_feat, feat_name in enumerate(y_pheno.columns):
            dist = y_comp_df.loc[y_comp_df.Feature == feat_name, 
                                 'Loading'].sort_values().values
            cil, ciu = dist[int(0.05 * len(dist))], dist[int(0.95 * len(dist))]
            plt.plot([i_feat, i_feat], [cil, ciu], c='k', lw=1)
        
        plt.xticks(rotation=90)
        plt.ylabel(f'C{c+1} Loading')
        plt.title(f'Y side: Component {c+1}')
        plt.tight_layout()
        plt.savefig(os.path.join(pls_results_dir, f'C{c+1}_Y_loadings.png'), dpi=300)
        plt.close()

        # Store info in summary_df
        summary_df.loc[c, 'component'] = c + 1
        summary_df.loc[c, 'num_sig_x'] = num_sig_x
        summary_df.loc[c, 'num_sig_y'] = num_sig_y
        summary_df.loc[c, 'total_x_feats'] = total_x_feats
        summary_df.loc[c, 'total_y_feats'] = total_y_feats

    # Save summary_df
    out_fn = os.path.join(pls_results_dir, 'bootstrap_loadings_summary.csv')
    summary_df.to_csv(out_fn, index=False)

    return summary_df

