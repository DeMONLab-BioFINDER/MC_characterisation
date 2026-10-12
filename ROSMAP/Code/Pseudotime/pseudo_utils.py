import matplotlib.pyplot as plt 
import numpy as np
import os

def plot_gam_sexcurves(
    cluster,
    df,
    fitted_models,
    cluster_col="Cluster",
    x_col="gpath",
    y_col="Expression",
    sex_col="msex",
    x_label="Global Pathology (log-scaled)",
    y_label="Gene expr. (z-scored)",
    save_dir=None,
    save_prefix="GAM_sex",
    save=True,
    figsize=(6, 5),
    n_grid=200,
):
    sub = df[df[cluster_col] == cluster]
    gam = fitted_models[cluster]
    g = np.linspace(sub[x_col].min(), sub[x_col].max(), n_grid)

    plt.figure(figsize=figsize)

    sex_info = {
        0: {"color": "#B66638", "label": "Female"},
        1: {"color": "#2B518D", "label": "Male"}
        }

    for sex, info in sex_info.items():
        Xp = np.c_[g, np.full_like(g, sex)]

        y = gam.predict(Xp)
        ci = gam.confidence_intervals(Xp)

        plt.plot(g, y, color=info["color"], lw=2, label=info["label"])
        plt.fill_between(g, ci[:, 0], ci[:, 1], color=info["color"], alpha=0.15)

        sex_sub = sub[sub[sex_col] == sex]

        plt.scatter(sex_sub[x_col], sex_sub[y_col], color=info["color"], alpha=0.9,s=15)

    cluster_label = "_".join(cluster.split("_")[1:])
    cluster_label = cluster_label.replace("_", " ")

    plt.title(cluster_label, fontsize=18)
    plt.xlabel(x_label, fontsize=16)
    plt.ylabel(y_label, fontsize=16)
    plt.xticks([0, 0.5, 1.0, 1.5], fontsize=10)
    plt.yticks(fontsize=14)
    plt.legend(loc="upper right", bbox_to_anchor=(1.02, 1), fontsize=22)
    plt.grid(alpha=0.3)
    plt.legend(
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        borderaxespad=0,
        fontsize=18
    )

    plt.tight_layout()

    if save and save_dir is not None:
        name = (cluster_label.replace(" ", "_"))
        save_path = os.path.join(save_dir, f"{save_prefix}_{name}.pdf")

        plt.savefig(save_path, dpi=300, bbox_inches="tight", transparent=True)
        
    plt.show()



def plot_cluster_overall(cluster, df, fitted_models, df_over, save_dir=None):
    """Plot the fitted GAM curve WITHOUT sex (one curve for everyone), with its
    95% confidence band, on top of the measured data points for one cluster.
    The FDR-corrected q-value for gpath is shown in the title.

    df            : long table with columns Cluster, gpath, Expression
    fitted_models : dict, cluster name -> GAM fitted on gpath only
    df_over       : results table with columns Cluster and q_gpath
    If save_dir is given, the figure is also saved there as a PNG.
    """
    # rows of the long table that belong to this cluster
    sub = df[df["Cluster"] == cluster]
    gam = fitted_models[cluster]

    # 200 evenly spaced gpath values covering this cluster's data;
    # the fitted curve is drawn through these points
    x_grid = np.linspace(sub["gpath"].min(), sub["gpath"].max(), 200)

    # this model has only one predictor (gpath), so the grid can go in directly
    y_pred = gam.predict(x_grid)
    ci = gam.confidence_intervals(x_grid)  # column 0 = lower, column 1 = upper

    plt.figure(figsize=(6, 4))
    plt.scatter(sub["gpath"], sub["Expression"], alpha=0.3, s=15, color="grey")
    plt.plot(x_grid, y_pred, color="black", lw=2)
    plt.fill_between(x_grid, ci[:, 0], ci[:, 1], color="grey", alpha=0.25)

    # look up this cluster's q-value in the results table
    result_row = df_over[df_over["Cluster"] == cluster]
    q = result_row["q_gpath"].values[0]

    plt.title(f"{cluster}\nq_gpath={q:.3g}")
    plt.xlabel("gpath")
    plt.ylabel("Z-scored Expression")
    plt.grid(alpha=0.3)
    plt.tight_layout()

    if save_dir is not None:
        plt.savefig(os.path.join(save_dir, f"{cluster}.png"), dpi=300)

    plt.show()

