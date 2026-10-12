import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
from matplotlib.patches import FancyArrowPatch
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
import hashlib
import re
# ====================================================================
#  Helper Functions
# ====================================================================

TYPE_COLS = {
    "Astrocytes": "#1f77b4",
    "Excitatory": "#ff7f0e",
    "Immune": "#2ca02c",
    "Inhibitory": "#d62728",
    "OPCs": "#9467bd",
    "Oligodendrocytes": "#8c564b",
}

SIGN_COLS = {
    "pos": "#e41a1c",
    "neg": "#377eb8",
    "ns": "#bdbdbd",
}

#Shorten the cell subtype names
def subtype_label(name: str) -> str:
    """Return the subtype part of a name like 'Astrocytes_AstDPP10' -> 'AstDPP10'."""
    return str(name).split("_", 1)[1] if "_" in str(name) else str(name)

# Count matches per AGM pair
def make_matrices(records, clusters, order, directed):
    """
    Count CellPhoneDB interactions between every pair of AGMs.

    For each AGM pair, checks whether the genes required for each interaction
    partner are present in the corresponding AGM gene sets. All genes belonging
    to a protein complex must be present for the interaction to match.

    Directed interactions are checked in the recorded partner order. Undirected
    interactions are checked in both orientations and copied across the matrix
    diagonal.

    Parameters
    ----------
    records : list
        Resolved CellPhoneDB interaction records for one directionality class.
    clusters : dict
        Mapping of AGM names to their gene sets.
    order : list
        AGM names in the order used for the matrix rows and columns.
    directed : bool
        Whether interaction direction should be preserved.

    Returns
    -------
    binary_mat : numpy.ndarray
        Matrix showing whether at least one interaction exists between each
        AGM pair.
    count_mat : numpy.ndarray
        Matrix containing the number of matched interactions between each
        AGM pair.
    """
    n = len(order)
    count_mat = np.zeros((n, n), dtype=int)

    for i, cluster_a in enumerate(order):
        genes_a = clusters[cluster_a]
        for j in range(n) if directed else range(i, n):
            genes_b = clusters[order[j]]

            for _, _, partner_1, partner_2, _, _ in records:
                forward = partner_1 <= genes_a and partner_2 <= genes_b
                reverse = partner_2 <= genes_a and partner_1 <= genes_b
                if forward or (not directed and reverse):
                    count_mat[i, j] += 1

            if not directed:
                count_mat[j, i] = count_mat[i, j]

    binary_mat = (count_mat > 0).astype(int)
    return binary_mat, count_mat

def draw_one_null(sizes_by_subtype, hvgs_by_major, rng):
    """

    Args:
        sizes_by_subtype (_type_): _description_
        hvgs_by_major (_type_): _description_
        rng (_type_): _description_

    Raises:
        ValueError: _description_
        ValueError: _description_

    Returns:
        _type_: _description_
        
    Returns {cluster_id: set(genes)}.
    For each (major, subtype), takes a shuffled copy of that major's HVG pool,
    slices disjoint chunks for that subtype's module sizes, assigns to cluster_ids.
    Genes CAN repeat across different subtypes (by design).
    """
    null_modules = {}
    for (major, subtype), items in sizes_by_subtype.items():
        pool = hvgs_by_major.get(major, []) #hvg pool for that specific cell type
        pool_copy = list(pool) # make a copy so other subtype lists wouldnt be reshuffled yet
        rng.shuffle(pool_copy)
        pos = 0
        for cluster_id, sz in items:
            genes = pool_copy[pos:pos+sz]
            pos += sz
            null_modules[cluster_id] = set(genes)
    return null_modules


def norm_cluster_name(x: str) -> str:
    """Clean cluster names"""
    x = str(x).strip()
    x = re.sub(r"\s+", "", x)
    x = re.sub(r"_Cluster_(\d+)\.0$", r"_Cluster_\1", x)
    return x


#pick all edges within N hops of a centre node, for one subtype interactions building
def edges_of_interest(
    G,
    center: str,
    mode: str = "out",      # "out" | "in" | "both"
    depth: int = 1,         # how many hops away to keep
    undirected: bool = True # True → return edges as sorted tuples
):
    """
    Return (keep_edges, keep_nodes) where
        keep_edges  is a *set* of edge keys (u,v)
        keep_nodes  is a *set* of node names

    mode="out"   → only edges leaving `center`
    mode="in"    → only edges entering `center`
    mode="both"  → all edges touching the frontier

    If undirected=True, every edge is stored as tuple(sorted([u,v]))
    so ('Ast','OPC') and ('OPC','Ast') are the same key.  That matches
    how we stored contact lines in edge_artist.
    """
    if center not in G:
        raise ValueError(f"{center!r} not in graph")

    frontier = {center}
    visited  = {center}
    keep_edges = set()

    for _ in range(depth):
        nxt = set()

        for n in frontier:
            if mode in ("out", "both"):
                for _, v in G.out_edges(n):
                    key = tuple(sorted([n, v])) if undirected else (n, v)
                    keep_edges.add(key)
                    if v not in visited:
                        nxt.add(v)

            if mode in ("in", "both"):
                for u, _ in G.in_edges(n):
                    key = tuple(sorted([u, n])) if undirected else (u, n)
                    keep_edges.add(key)
                    if u not in visited:
                        nxt.add(u)

        visited |= nxt
        frontier = nxt

    return keep_edges, visited


def plot_ccc_subtypenetwork(
    df,
    title,
    minlr=1,
    mincnt=1,
    spring_k=20.0,
    spring_iter=500,
    curve_rad=0.2,
    arrow_size=22,
    node_size=1000,
    nodefonts=12,
    figsize=(9, 5),
    offset_frac=0.01,
    black_frac=0.75,
    width_min=1,
    width_max=13,
    palette=None,
    savepath=None,
    legendfont=12,
):
    """
    Plot subtype-level CCC network.

    LR interactions are drawn as directed arrows. Non-LR contact-like interactions
    are drawn as separate black undirected lines, so they cannot overwrite LR arrows.
    """
    if palette is None:
        palette = sns.color_palette("tab10")

    df = df.copy()
    subtype = lambda n: "_".join(str(n).split("_")[:2])
    ctype = lambda s: str(s).split("_", 1)[0]

    df["s_sub"] = df["sender_cluster"].apply(subtype)
    df["r_sub"] = df["receiver_cluster"].apply(subtype)

    subtypes = sorted(set(df["s_sub"]).union(df["r_sub"]))
    ctypes = sorted({ctype(s) for s in subtypes})
    type2col = {t: TYPE_COLS.get(t, palette[i % len(palette)]) for i, t in enumerate(ctypes)}

    G = nx.DiGraph()
    for s in subtypes:
        G.add_node(s, color=type2col[ctype(s)], ctype=ctype(s))

    lr = (
        df[df["directionality"].eq("Ligand-Receptor")]
        .groupby(["s_sub", "r_sub"])
        .size()
        .reset_index(name="n")
    )

    for _, r in lr.iterrows():
        if r["n"] >= minlr and r["s_sub"] != r["r_sub"]:
            G.add_edge(r["s_sub"], r["r_sub"], cls="LR", weight=int(r["n"]))

    contacts = (
        df[~df["directionality"].eq("Ligand-Receptor")]
        .groupby(["s_sub", "r_sub"])
        .size()
        .reset_index(name="n")
    )

    contact_edges = {}
    for _, r in contacts.iterrows():
        if r["n"] < mincnt or r["s_sub"] == r["r_sub"]:
            continue
        u, v = sorted([r["s_sub"], r["r_sub"]])
        contact_edges[(u, v)] = contact_edges.get((u, v), 0) + int(r["n"])

    pos = nx.spring_layout(G, seed=42, k=spring_k, iterations=spring_iter)
    fig, ax = plt.subplots(figsize=figsize)

    nx.draw_networkx_nodes(
        G,
        pos,
        node_color=[G.nodes[n]["color"] for n in G],
        node_size=node_size,
        edgecolors="black",
        linewidths=0.9,
        ax=ax,
    )
    ax.margins(x=0.17, y=0.17)

    labels = {n: subtype_label(n) for n in G.nodes()}
    label_artists = nx.draw_networkx_labels(
        G, pos, labels=labels, font_size=nodefonts, font_weight="bold", ax=ax
    )
    for txt in label_artists.values():
        txt.set_zorder(6)
        txt.set_path_effects([pe.withStroke(linewidth=3, foreground="white")])

    lr_wts = [d["weight"] for _, _, d in G.edges(data=True) if d.get("cls") == "LR"]
    wmin, wmax = (min(lr_wts), max(lr_wts)) if lr_wts else (1, 1)

    def scale(w):
        if wmin == wmax:
            return width_min
        return width_min + (w - wmin) / (wmax - wmin) * (width_max - width_min)

    def offset(p, q, frac):
        dx, dy = q[0] - p[0], q[1] - p[1]
        length = math.hypot(dx, dy)
        if length == 0:
            return p, q
        ox, oy = -dy / length * frac, dx / length * frac
        return (p[0] + ox, p[1] + oy), (q[0] + ox, q[1] + oy)

    for (u, v), n in contact_edges.items():
        width = scale(n) * black_frac
        p1o, p2o = offset(pos[u], pos[v], offset_frac)
        ax.plot(
            [p1o[0], p2o[0]],
            [p1o[1], p2o[1]],
            color="black",
            linewidth=width,
            alpha=0.45,
            zorder=1,
        )

    pair_rad = {}
    for u, v, d in G.edges(data=True):
        if d.get("cls") != "LR":
            continue

        key = tuple(sorted([u, v]))
        if key not in pair_rad:
            h = int(hashlib.md5(("|".join(key)).encode()).hexdigest(), 16)
            pair_rad[key] = curve_rad if (h % 2 == 0) else -curve_rad

        arrow = FancyArrowPatch(
            pos[u],
            pos[v],
            connectionstyle=f"arc3,rad={pair_rad[key]}",
            arrowstyle="-|>",
            mutation_scale=arrow_size,
            linewidth=scale(d["weight"]),
            color=G.nodes[u]["color"],
            alpha=1,
            zorder=3,
            shrinkA=15,
            shrinkB=24,
        )
        ax.add_patch(arrow)

    ticks = sorted(set(np.linspace(wmin, wmax, 3, dtype=int)))
    edge_handles = [
        Line2D([0], [0], lw=scale(t), color="black", label=f"{int(t)} LR pairs")
        for t in ticks
    ]
    node_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=type2col[t],
            markeredgecolor="black",
            markersize=10,
            label=t,
        )
        for t in ctypes
    ]

    if edge_handles:
        leg1 = ax.legend(
            handles=edge_handles,
            fontsize=legendfont,
            loc="upper left",
            bbox_to_anchor=(1.02, 1),
            borderaxespad=0,
            frameon=False,
        )
        ax.add_artist(leg1)

    ax.legend(
        handles=node_handles,
        fontsize=legendfont,
        loc="lower left",
        bbox_to_anchor=(1.02, 0),
        borderaxespad=0,
        frameon=False,
    )

    plt.subplots_adjust(right=0.75)
    ax.set_title(title, fontsize=17, pad=20)
    ax.axis("off")

    if savepath is not None:
        fig.savefig(savepath, dpi=300, bbox_inches="tight", pad_inches=0.50, transparent=True)

    return fig, G, contact_edges, pos


def plot_ccc_spotlight(
    G,
    contact_edges,
    pos,
    focus_node,
    title,
    mode="both",
    depth=1,
    curve_rad=0.2,
    arrow_size=22,
    offset_frac=0.01,
    black_frac=0.75,
    width_min=1,
    width_max=13,
    savepath=None,
):
    """Plot only edges and nodes within N hops of a focus node."""
    if focus_node not in G:
        raise ValueError(f"focus_node {focus_node!r} not in G.nodes()")

    fig, ax = plt.subplots(figsize=(15, 11))
    labels = {n: subtype_label(n) for n in G.nodes()}

    lr_wts = [d["weight"] for _, _, d in G.edges(data=True) if d.get("cls") == "LR"]
    wmin, wmax = (min(lr_wts), max(lr_wts)) if lr_wts else (1, 1)

    def scale(w):
        if wmin == wmax:
            return width_min
        return width_min + (w - wmin) / (wmax - wmin) * (width_max - width_min)

    def offset(p, q, frac):
        dx, dy = q[0] - p[0], q[1] - p[1]
        length = math.hypot(dx, dy)
        if length == 0:
            return p, q
        ox, oy = -dy / length * frac, dx / length * frac
        return (p[0] + ox, p[1] + oy), (q[0] + ox, q[1] + oy)

    frontier = {focus_node}
    keep_nodes = {focus_node}
    keep_lr = set()
    keep_contact = set()

    for _ in range(depth):
        nxt = set()
        for n in frontier:
            if mode in ("out", "both"):
                for _, v in G.out_edges(n):
                    keep_lr.add((n, v))
                    nxt.add(v)
            if mode in ("in", "both"):
                for u, _ in G.in_edges(n):
                    keep_lr.add((u, n))
                    nxt.add(u)

            for a, b in contact_edges:
                if a == n or b == n:
                    keep_contact.add((a, b))
                    nxt.add(b if a == n else a)

        keep_nodes |= nxt
        frontier = nxt

    for n in keep_nodes:
        ax.scatter(
            pos[n][0],
            pos[n][1],
            s=2200,
            facecolor=G.nodes[n]["color"],
            edgecolor="black",
            linewidth=0.9,
            zorder=2,
        )

        txt = ax.text(
            pos[n][0],
            pos[n][1],
            labels[n],
            ha="center",
            va="center",
            fontsize=15,
            fontweight="bold",
            zorder=6,
        )
        txt.set_path_effects([pe.withStroke(linewidth=3, foreground="white")])

    for (u, v), n in contact_edges.items():
        if (u, v) not in keep_contact:
            continue
        width = scale(n) * black_frac
        p1o, p2o = offset(pos[u], pos[v], offset_frac)
        ax.plot(
            [p1o[0], p2o[0]],
            [p1o[1], p2o[1]],
            color="black",
            linewidth=width,
            alpha=0.45,
            zorder=1,
        )

    pair_rad = {}
    for u, v, d in G.edges(data=True):
        if d.get("cls") != "LR" or (u, v) not in keep_lr:
            continue

        key = tuple(sorted([u, v]))
        if key not in pair_rad:
            h = int(hashlib.md5(("|".join(key)).encode()).hexdigest(), 16)
            pair_rad[key] = curve_rad if (h % 2 == 0) else -curve_rad

        arrow = FancyArrowPatch(
            pos[u],
            pos[v],
            connectionstyle=f"arc3,rad={pair_rad[key]}",
            arrowstyle="-|>",
            mutation_scale=arrow_size,
            linewidth=scale(d["weight"]),
            color=G.nodes[u]["color"],
            alpha=0.8,
            zorder=3,
            shrinkA=15,
            shrinkB=24,
        )
        ax.add_patch(arrow)

    ctypes = sorted({G.nodes[n]["ctype"] for n in G.nodes()})
    type2col = {G.nodes[n]["ctype"]: G.nodes[n]["color"] for n in G.nodes()}

    ticks = sorted(set(np.linspace(wmin, wmax, 3, dtype=int)))
    edge_handles = [
        Line2D([0], [0], lw=scale(t), color="black", label=f"{int(t)} LR pairs")
        for t in ticks
    ]
    node_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=type2col[t],
            markeredgecolor="black",
            markersize=10,
            label=t,
        )
        for t in ctypes
    ]

    if edge_handles:
        leg1 = ax.legend(handles=edge_handles, title="Arrow width = LR count", fontsize=17, loc="upper left", frameon=False)
        ax.add_artist(leg1)

    ax.legend(handles=node_handles, title="Cell type", fontsize=17, loc="upper right", frameon=False)
    ax.set_title(title, fontsize=17, pad=20)
    ax.axis("off")

    if savepath is not None:
        fig.savefig(savepath, dpi=300, bbox_inches="tight", pad_inches=0.50)

    return fig, ax

def plot_ccc_cluster_network(
    df,
    title="Cluster-level CCC network",
    # thresholds & styling
    minlr=1,
    mincnt=1,
    arrow_size=18,
    width_min=0.8,
    width_max=6.0,
    black_frac=0.7,
    curve_rad=0.18,
    seed=42,
    spring_k=1.1,
    spring_iter=500,
    # cluster scattering around subtype anchors
    r_base=0.001,
    r_per_node=0.008,
    angle_jitter=0.06,
    r_jitter=0.08,
    # coloring
    color_by="celltype", # "celltype" or "pls"
    pls_df=None, # dataframe OR path to CSV from PLS export
    pls_component=1, # 1=C1, 2=C2, ...
    pls_sig_col_template="C{c}_sig",
    pls_val_col_template="C{c}_val",
    pls_alpha=None, # optional: if sig col missing, use p<alpha with pls_alpha
    pls_p_col_template="C{c}_p",
    sign_cols=None, # dict pos/neg/ns -> color
    type_cols=None,  # dict major -> color
    arrow_color="sender", # "sender" or "black"
    node_size=520,
    label_fontsize=13,
    figsize=(8, 6),
    savepath=None,
    return_objects=False,
):
    """
    Plot cluster-level CCC network.

    LR arrows are stored in G. Non-LR contact edges are stored separately so they
    cannot overwrite LR edges between the same node pair.
    """
    if type_cols is None:
        type_cols = TYPE_COLS
    if sign_cols is None:
        sign_cols = SIGN_COLS

    def split_label(cluster_name: str):
        cluster_name = str(cluster_name)
        if "_Cluster_" in cluster_name:
            left, cid = cluster_name.rsplit("_Cluster_", 1)
            parts = left.split("_", 1)
            major = parts[0]
            subtype = parts[1] if len(parts) == 2 else ""
            return major, subtype, cid, f"{major}_{subtype}"

        parts = cluster_name.split("_")
        major = parts[0]
        subtype = parts[1] if len(parts) > 1 else ""
        cid = parts[-1] if parts else ""
        return major, subtype, cid, f"{major}_{subtype}"

    def cluster_label(fullname):
        major, subtype, cid, _ = split_label(fullname)
        return f"{subtype}\nC{cid}" if subtype else f"{major}\nC{cid}"

    def scale_width(value, vmin, vmax):
        if vmax <= vmin:
            return (width_min + width_max) / 2.0
        return width_min + (value - vmin) * (width_max - width_min) / (vmax - vmin)

    df = df.copy()

    for col, newcols in [
        ("sender_cluster", ("s_major", "s_subtype", "s_cid", "s_key")),
        ("receiver_cluster", ("r_major", "r_subtype", "r_cid", "r_key")),
    ]:
        parsed = df[col].apply(split_label)
        for i, nc in enumerate(newcols):
            df[nc] = parsed.apply(lambda x: x[i])

    cluster_nodes = sorted(set(df["sender_cluster"]).union(df["receiver_cluster"]))
    node2key = {n: split_label(n)[3] for n in cluster_nodes}
    node2type = {n: split_label(n)[0] for n in cluster_nodes}

    # Layout subtype anchors first, then scatter clusters around their subtype anchor.
    SG = nx.Graph()
    for key in sorted(set(df["s_key"]).union(df["r_key"])):
        SG.add_node(key, major=key.split("_")[0])

    for _, row in df.iterrows():
        u, v = row["s_key"], row["r_key"]
        if u == v:
            continue
        if SG.has_edge(u, v):
            SG[u][v]["weight"] += 1
        else:
            SG.add_edge(u, v, weight=1)

    pos_super = nx.spring_layout(SG, seed=seed, k=spring_k, iterations=spring_iter)

    subtype2clusters = {}
    for n in cluster_nodes:
        subtype2clusters.setdefault(node2key[n], []).append(n)

    rng = np.random.default_rng(seed)
    pos = {}
    for subkey, nodes in subtype2clusters.items():
        cx, cy = pos_super[subkey]
        r = r_base + r_per_node * len(nodes)
        angles = np.linspace(0, 2 * np.pi, len(nodes), endpoint=False)
        for i, node in enumerate(nodes):
            a = angles[i] + rng.uniform(-angle_jitter, angle_jitter)
            rr = r * (1 - r_jitter / 2 + r_jitter * rng.random())
            pos[node] = (cx + rr * np.cos(a), cy + rr * np.sin(a))
    
    # ---------------- build directed graph ----------------
    G = nx.DiGraph()
    for n in cluster_nodes:
        G.add_node(n, major=node2type[n])

    is_lr = df["directionality"].str.lower().eq("ligand-receptor")
    lr_counts = (
        df[is_lr]
        .groupby(["sender_cluster", "receiver_cluster"])
        .size()
        .reset_index(name="n")
    )
    contact_counts = (
        df[~is_lr]
        .groupby(["sender_cluster", "receiver_cluster"])
        .size()
        .reset_index(name="n")
    )

    for _, r in lr_counts.iterrows():
        u, v, n = r["sender_cluster"], r["receiver_cluster"], int(r["n"])
        if n >= minlr and u != v:
            G.add_edge(u, v, cls="LR", weight=n)

    contact_agg = {}
    for _, r in contact_counts.iterrows():
        u, v = r["sender_cluster"], r["receiver_cluster"]
        if u == v:
            continue
        key = tuple(sorted([u, v]))
        contact_agg[key] = contact_agg.get(key, 0) + int(r["n"])

    contact_edges = {key: n for key, n in contact_agg.items() if n >= mincnt}

    lr_w = [d["weight"] for _, _, d in G.edges(data=True) if d.get("cls") == "LR"]
    wmin, wmax = (min(lr_w), max(lr_w)) if lr_w else (1, 1)

    if color_by.lower() == "pls":
        if pls_df is None:
            raise ValueError("color_by='pls' requires pls_df (DataFrame or CSV path).")
        if isinstance(pls_df, (str, bytes)):
            pls_df = pd.read_csv(pls_df)

        name_col = "gene_cluster" if "gene_cluster" in pls_df.columns else "Gene_clus" if "Gene_clus" in pls_df.columns else None
        if name_col is None:
            raise ValueError("pls_df must have a 'gene_cluster' or 'Gene_clus' column with cluster names.")

        c = pls_component
        sig_col = pls_sig_col_template.format(c=c)
        val_col = pls_val_col_template.format(c=c)
        p_col = pls_p_col_template.format(c=c)
        pls_df = pls_df.copy()

        if sig_col in pls_df.columns:
            sig = pls_df[sig_col].astype(bool).values
        elif pls_alpha is not None and p_col in pls_df.columns:
            sig = pls_df[p_col].astype(float).values < float(pls_alpha)
        else:
            raise ValueError(f"pls_df is missing '{sig_col}'. Provide pls_alpha and '{p_col}', or include '{sig_col}'.")

        if val_col not in pls_df.columns:
            raise ValueError(f"pls_df is missing '{val_col}' for component C{c}.")

        cluster2sign = {}
        for nm, s, v in zip(pls_df[name_col].astype(str), sig, pls_df[val_col].astype(float)):
            cluster2sign[nm] = "ns" if (not s or v == 0) else "pos" if v > 0 else "neg"

        def node_color(n):
            return sign_cols.get(cluster2sign.get(n, "ns"), sign_cols["ns"])

        legend_mode = "pls"
    else:
        def node_color(n):
            return type_cols.get(G.nodes[n]["major"], "#7f7f7f")

        legend_mode = "celltype"

    def arrow_col(u):
        return "black" if arrow_color == "black" else node_color(u)

    fig, ax = plt.subplots(figsize=figsize)

    nx.draw_networkx_nodes(
        G,
        pos,
        node_color=[node_color(n) for n in G.nodes()],
        node_size=node_size,
        edgecolors="black",
        linewidths=0.7,
        ax=ax,
    )

    text_art = nx.draw_networkx_labels(
        G,
        pos,
        labels={n: cluster_label(n) for n in G.nodes()},
        font_size=label_fontsize,
        font_weight="bold",
        ax=ax,
    )
    for txt in text_art.values():
        txt.set_zorder(6)
        txt.set_path_effects([pe.withStroke(linewidth=2.5, foreground="white")])

    # Draw contact edges from a separate dictionary so they cannot overwrite LR arrows.
    for (u, v), n in contact_edges.items():
        width = scale_width(n, wmin, wmax) * black_frac
        x1, y1 = pos[u]
        x2, y2 = pos[v]
        ax.plot([x1, x2], [y1, y2], color="black", alpha=0.4, linewidth=width, zorder=1)

    pair_rad = {}
    for u, v, d in G.edges(data=True):
        if d.get("cls") != "LR":
            continue

        key = tuple(sorted([u, v]))
        if key not in pair_rad:
            h = int(hashlib.md5(("|".join(key)).encode()).hexdigest(), 16)
            pair_rad[key] = curve_rad if (h % 2 == 0) else -curve_rad

        arrow = FancyArrowPatch(
            pos[u],
            pos[v],
            connectionstyle=f"arc3,rad={pair_rad[key]}",
            arrowstyle="-|>",
            mutation_scale=arrow_size,
            linewidth=scale_width(d["weight"], wmin, wmax),
            color=arrow_col(u),
            alpha=0.9,
            shrinkA=8,
            shrinkB=10,
            zorder=3,
        )
        ax.add_patch(arrow)

    ticks = sorted(set(np.linspace(wmin, wmax, 3, dtype=int)))
    edge_handles = [
        Line2D([0], [0], lw=scale_width(t, wmin, wmax), color="black", label=f"{int(t)} LR pairs")
        for t in ticks
    ]
    if edge_handles:
        leg1 = ax.legend(handles=edge_handles, title="Arrow width = LR count", loc="upper left", frameon=False)
        ax.add_artist(leg1)

    if legend_mode == "celltype":
        majors = []
        for n in G.nodes():
            m = G.nodes[n]["major"]
            if m not in majors:
                majors.append(m)
        node_handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                markerfacecolor=type_cols.get(m, "#7f7f7f"),
                markeredgecolor="black",
                markersize=8,
                label=m,
            )
            for m in majors
        ]
        ax.legend(handles=node_handles, title="Cell type", loc="upper right", frameon=False)
    else:
        dir_handles = [
            Line2D([0], [0], marker="o", linestyle="", markersize=8, markerfacecolor=sign_cols["pos"], markeredgecolor="black", label=f"C{pls_component} +"),
            Line2D([0], [0], marker="o", linestyle="", markersize=8, markerfacecolor=sign_cols["neg"], markeredgecolor="black", label=f"C{pls_component} −"),
            Line2D([0], [0], marker="o", linestyle="", markersize=8, markerfacecolor=sign_cols["ns"], markeredgecolor="black", label="Not significant"),
        ]
        ax.legend(handles=dir_handles, title=f"PLS C{pls_component} direction", loc="upper right", frameon=False)

    ax.set_title(title, pad=16)
    ax.axis("off")
    fig.tight_layout()

    if savepath is not None:
        fig.savefig(savepath, dpi=300, bbox_inches="tight")

    if return_objects:
        return fig, ax, G, pos, pos_super

    return fig, ax



def plot_ccc_subtypenetwork_sig_grey(
    df,
    title=None,
    sig_pairs=None,                 # set of unordered (subA, subB) pairs, e.g. {("Astrocytes_AstDPP10","OPCs_OPC"), ...}
    minlr=1,
    mincnt=1,
    spring_k=20.0,
    spring_iter=500,
    curve_rad=0.2,
    arrow_size=13,
    offset_frac=0.01,
    black_frac=0.75,
    width_min=1,
    width_max=4,
    palette=None,
    grey="#8585855c",
    grey_alpha=0.8,
    nodefonts=8,
    legendfont=8,
    figsize = (9,5),
    savepath=None
):
    """
    Network plot where ALL edges are shown, but NON-significant subtype-pairs are greyed out.
      - LR edges: directed arrows, colored by sender cell type if significant else grey
      - Contact edges (undirected black lines): black if significant else grey
      - Nodes: colored by cell type if they appear in ANY significant pair, else grey
    """
    if palette is None:
        palette = sns.color_palette("tab10")

    sig_pairs = set() if sig_pairs is None else set(sig_pairs)

    df = df.copy()

    subtype = lambda n: "_".join(str(n).split("_")[:2])     # e.g. "Astrocytes_AstDPP10"
    ctype   = lambda s: str(s).split("_", 1)[0]            # e.g. "Astrocytes"

    # subtype columns
    df["s_sub"] = df["sender_cluster"].apply(subtype)
    df["r_sub"] = df["receiver_cluster"].apply(subtype)

    # universe of nodes
    subtypes = sorted(set(df["s_sub"]).union(df["r_sub"]))
    ctypes   = sorted({ctype(s) for s in subtypes})
    type2col = {t: palette[i % len(palette)] for i, t in enumerate(ctypes)}

    # ---- build graph ----
    G = nx.DiGraph()
    for s in subtypes:
        G.add_node(s, color=type2col[ctype(s)], ctype=ctype(s))

    # ---- LR edges (directed) ----
    lr = (df[df.directionality == "Ligand-Receptor"]
            .groupby(["s_sub", "r_sub"]).size().reset_index(name="n"))

    for _, r in lr.iterrows():
        if r["n"] < minlr or r["s_sub"] == r["r_sub"]:
            continue
        u, v = r["s_sub"], r["r_sub"]
        key = tuple(sorted([u, v]))  # undirected key for significance
        G.add_edge(u, v, cls="LR", weight=int(r["n"]), sig=(key in sig_pairs))

    # ---- contact edges (undirected) ----
    contact_classes = ["Adhesion-Adhesion", "Gap-Gap", "Ligand-Ligand"]
    cnt = (df[df.directionality.isin(contact_classes)]
            .groupby(["s_sub", "r_sub"]).size().reset_index(name="n"))

    # store as undirected pair -> {"n": total_count, "sig": bool}
    contact_edges = {}
    for _, r in cnt.iterrows():
        if r["n"] < mincnt or r["s_sub"] == r["r_sub"]:
            continue
        u, v = sorted([r["s_sub"], r["r_sub"]])
        key = (u, v)
        if key not in contact_edges:
            contact_edges[key] = {"n": 0, "sig": (key in sig_pairs)}
        contact_edges[key]["n"] += int(r["n"])
        # sig stays as originally computed from sig_pairs (same for whole pair)

    # ---- which nodes are in any significant pair? ----
    sig_nodes = set()
    for a, b in sig_pairs:
        sig_nodes.add(a)
        sig_nodes.add(b)

    # layout
    pos = nx.spring_layout(G, seed=42, k=spring_k, iterations=spring_iter)

    #  plot 
    fig, ax = plt.subplots(figsize=figsize)

    # nodes (grey if not in any significant pair)
    node_colors = [(G.nodes[n]["color"] if n in sig_nodes else grey) for n in G.nodes()]
    nx.draw_networkx_nodes(
        G, pos,
        node_color=node_colors,
        node_size=1000,
        edgecolors="black",
        linewidths=0.9,
        ax=ax
    )
    ax.margins(x=0.17, y=0.17)

    # labels (grey if node not significant)
    labels = {n: subtype_label(n) for n in G.nodes()}
    label_artists = nx.draw_networkx_labels(
        G, pos, labels=labels,
        font_size=nodefonts, font_weight="bold", ax=ax
    )
    for n, txt in label_artists.items():
        txt.set_zorder(6)
        txt.set_path_effects([pe.withStroke(linewidth=3, foreground="white")])
        if n not in sig_nodes:
            txt.set_color("#666666")

    # ---- scale widths based on LR weights (and re-used for contacts) ----
    lr_wts = [d["weight"] for _, _, d in G.edges(data=True) if d.get("cls") == "LR"]
    wmin, wmax = (min(lr_wts), max(lr_wts)) if lr_wts else (1, 1)

    def scale(w):
        if wmin == wmax:
            return width_min
        return width_min + (w - wmin) / (wmax - wmin) * (width_max - width_min)

    def offset(p, q, frac):
        dx, dy = q[0] - p[0], q[1] - p[1]
        length = math.hypot(dx, dy)
        if length == 0:
            return p, q
        ox, oy = -dy / length * frac, dx / length * frac
        return (p[0] + ox, p[1] + oy), (q[0] + ox, q[1] + oy)

    # ---- draw contact edges (black if sig else grey) ----
    # These are undirected; draw as straight line with slight offset.
    for (u, v), info in contact_edges.items():
        n = info["n"]
        is_sig = info["sig"]

        width = scale(n) * black_frac
        p1o, p2o = offset(pos[u], pos[v], offset_frac)

        ax.plot(
            [p1o[0], p2o[0]], [p1o[1], p2o[1]],
            color=("black" if is_sig else grey),
            linewidth=width,
            alpha=(0.45 if is_sig else grey_alpha),
            zorder=1
        )

    # ---- draw LR edges (colored if sig else grey) ----
    pair_rad = {}
    for u, v, d in G.edges(data=True):
        if d.get("cls") != "LR":
            continue

        key = tuple(sorted([u, v]))
        if key not in pair_rad:
            h = int(hashlib.md5(("|".join(key)).encode()).hexdigest(), 16)
            pair_rad[key] = curve_rad if (h % 2 == 0) else -curve_rad
        rad = pair_rad[key]

        width = scale(d["weight"])

        is_sig = d.get("sig", False)
        col = (G.nodes[u]["color"] if is_sig else grey)
        a = (0.8 if is_sig else grey_alpha)

        arrow = FancyArrowPatch(
            pos[u], pos[v],
            connectionstyle=f"arc3,rad={rad}",
            arrowstyle="-|>",
            mutation_scale=arrow_size,
            linewidth=width,
            color=col,
            alpha=a,
            zorder=3,
            shrinkA=15,
            shrinkB=24
        )
        ax.add_patch(arrow)

    # ---- legends ----
    ticks = sorted(set(np.linspace(wmin, wmax, 3, dtype=int)))
    edge_handles = [
        Line2D([0], [0], lw=scale(t), color="black", label=f"{int(t)} LR pairs")
        for t in ticks
    ]

    node_handles = [
        Line2D([0], [0], marker="o", linestyle="",
               markerfacecolor=type2col[t], markeredgecolor="black",
               markersize=10, label=t)
        for t in ctypes
    ]

    # add "greyed out" legend entries
    grey_handles = [
        Line2D([0], [0], marker="o", linestyle="",
               markerfacecolor=grey, markeredgecolor="black",
               markersize=10, label="Not significant node"),
        Line2D([0], [0], lw=3, color=grey, alpha=grey_alpha, label="Not significant edge"),
    ]

    if edge_handles:
        leg1 = ax.legend(
            handles=edge_handles,
            title="",
            fontsize=legendfont,
            loc="upper left",
            bbox_to_anchor=(1.02, 1),
            borderaxespad=0,
            frameon=False
        )
        ax.add_artist(leg1)

        leg2 = ax.legend(
            handles=node_handles + grey_handles,
            title="",
            fontsize=legendfont,
            loc="lower left",
            bbox_to_anchor=(1.02, 0),
            borderaxespad=0,
            frameon=False
            )
    
    plt.subplots_adjust(right=0.75)
    ax.set_title(title, fontsize=17, pad=20)
    ax.axis("off")

    if savepath is not None:
        fig.savefig(savepath, dpi=300, bbox_inches="tight", pad_inches=0.50)

    return fig

def plot_ccc_cluster_network_sig_grey(
    df,
    pls_df,
    sig_pairs=None,
    title="",
    pls_component=1,
    minlr=1,
    mincnt=1,
    seed=42,
    spring_k=1.1,
    curve_rad=0.18,
    arrow_size=18,
    width_max=6,
    node_size=520,
    label_fontsize=13,
    legendfont=16,
    figsize=(8, 6),
    r_base=0.001,
    savepath=None,
):
    """
    PLS-coloured cluster-level CCC network.

    Significant nodes:
        coloured, alpha = 0.5

    Non-significant nodes:
        grey, alpha = 1

    Significant LR edges:
        coloured by sender, alpha = 1

    Significant adhesion/gap edges:
        black, alpha = 1

    Non-significant edges:
        grey, alpha = 1
    """

    grey = "#c7c7c7"
    sign_cols = {"pos": "#e41a1c", "neg": "#377eb8", "ns": "#bdbdbd"}

    sig_pairs = {
        tuple(sorted(pair))
        for pair in (sig_pairs or [])
    }

    nodes = sorted(
        set(df["sender_cluster"]) |
        set(df["receiver_cluster"])
    )

    # Parse every node name once
    meta = pd.Series(nodes, index=nodes).str.extract(
        r"^([^_]+)_(.+)_Cluster_(.+)$"
    )
    meta.columns = ["major", "subtype", "cid"]
    meta["key"] = meta["major"] + "_" + meta["subtype"]

    df = df.assign(
        s_key=df["sender_cluster"].map(meta["key"]),
        r_key=df["receiver_cluster"].map(meta["key"]),
    )

    # ---------------------------------------------------------
    # Layout: first position subtypes, then clusters within them
    # ---------------------------------------------------------
    subtype_edges = df.loc[
        df["s_key"].ne(df["r_key"]),
        ["s_key", "r_key"],
    ].to_numpy()

    subtype_edges = pd.DataFrame(
        np.sort(subtype_edges, axis=1),
        columns=["u", "v"],
    )

    subtype_edges = (
        subtype_edges
        .value_counts(sort=False)
        .rename("weight")
        .reset_index()
    )

    SG = nx.Graph()
    SG.add_nodes_from(meta["key"].unique())
    SG.add_weighted_edges_from(
        subtype_edges.itertuples(index=False, name=None)
    )

    pos_super = nx.spring_layout(
        SG,
        seed=seed,
        k=spring_k,
        iterations=500,
    )

    rng = np.random.default_rng(seed)
    pos = {}

    for key, group in meta.groupby("key").groups.items():
        group = list(group)
        n = len(group)
        radius = r_base + 0.008 * n

        angles = (
            np.linspace(0, 2 * np.pi, n, endpoint=False)
            + rng.uniform(-0.06, 0.06, n)
        )

        radii = radius * (
            0.96 + 0.08 * rng.random(n)
        )

        centre = pos_super[key]

        coordinates = np.column_stack((
            centre[0] + radii * np.cos(angles),
            centre[1] + radii * np.sin(angles),
        ))

        pos.update(zip(group, map(tuple, coordinates)))

    # ---------------------------------------------------------
    # Count both communication classes in one groupby
    # ---------------------------------------------------------
    counts = (
        df.assign(
            is_lr=df["directionality"]
            .str.casefold()
            .eq("ligand-receptor")
        )
        .groupby(
            ["is_lr", "sender_cluster", "receiver_cluster"],
            sort=False,
        )
        .size()
        .rename("weight")
        .reset_index()
    )

    counts = counts.loc[
        counts["sender_cluster"].ne(counts["receiver_cluster"])
    ]

    lr = counts.loc[
        counts["is_lr"] &
        counts["weight"].ge(minlr)
    ]

    contact = counts.loc[
        ~counts["is_lr"] &
        counts["weight"].ge(mincnt)
    ]

    # Collapse adhesion/gap connections to unordered pairs
    contact_pairs = np.sort(
        contact[["sender_cluster", "receiver_cluster"]].to_numpy(),
        axis=1,
    )

    contact = (
        pd.DataFrame(contact_pairs, columns=["u", "v"])
        .assign(weight=contact["weight"].to_numpy())
        .groupby(["u", "v"], sort=False, as_index=False)["weight"]
        .sum()
    )

    # ---------------------------------------------------------
    # Directed LR graph
    # ---------------------------------------------------------
    G = nx.DiGraph()
    G.add_nodes_from(nodes)
    G.add_weighted_edges_from(
        lr[
            ["sender_cluster", "receiver_cluster", "weight"]
        ].itertuples(index=False, name=None)
    )

    sig_nodes = {
        node
        for pair in sig_pairs
        for node in pair
    }

    weights = pd.concat(
        [lr["weight"], contact["weight"]],
        ignore_index=True,
    )

    wmin, wmax = (
        (weights.min(), weights.max())
        if not weights.empty
        else (1, 1)
    )

    def edge_width(value):
        if wmin == wmax:
            return (0.8 + width_max) / 2

        return 0.8 + (
            (value - wmin) /
            (wmax - wmin)
        ) * (width_max - 0.8)

    # ---------------------------------------------------------
    # PLS colours
    # ---------------------------------------------------------
    if not isinstance(pls_df, pd.DataFrame):
        pls_df = pd.read_csv(pls_df)

    name_col = (
        "gene_cluster"
        if "gene_cluster" in pls_df.columns
        else "Gene_clus"
    )

    values = pd.to_numeric(
        pls_df[f"C{pls_component}_val"],
        errors="coerce",
    )

    cluster_sign = dict(zip(
        pls_df[name_col]
        .astype(str)
        .map(norm_cluster_name),

        np.select(
            [values.gt(0), values.lt(0)],
            ["pos", "neg"],
            default="ns",
        ),
    ))

    def node_colour(node):
        return sign_cols[
            cluster_sign.get(
                norm_cluster_name(node),
                "ns",
            )
        ]

    # ---------------------------------------------------------
    # Plot nodes and labels
    # ---------------------------------------------------------
    fig, ax = plt.subplots(figsize=figsize)

    nx.draw_networkx_nodes(
        G,
        pos,
        node_color=[
            node_colour(node)
            if node in sig_nodes
            else grey
            for node in nodes
        ],
        alpha=[
            0.5 if node in sig_nodes else 1
            for node in nodes
        ],
        node_size=node_size,
        edgecolors="black",
        linewidths=0.7,
        ax=ax,
    )

    labels = nx.draw_networkx_labels(
        G,
        pos,
        labels={
            node: f"{meta.at[node, 'subtype']}\nC{meta.at[node, 'cid']}"
            for node in nodes
        },
        font_size=label_fontsize,
        font_weight="bold",
        ax=ax,
    )

    for node, label in labels.items():
        label.set_zorder(6)
        label.set_path_effects([
            pe.withStroke(
                linewidth=2.5,
                foreground="white",
            )
        ])

        if node not in sig_nodes:
            label.set_color("#666666")

    # ---------------------------------------------------------
    # Undirected adhesion/gap edges
    # ---------------------------------------------------------
    for edge in contact.itertuples(index=False):
        ax.plot(
            [pos[edge.u][0], pos[edge.v][0]],
            [pos[edge.u][1], pos[edge.v][1]],
            color=(
                "black"
                if (edge.u, edge.v) in sig_pairs
                else grey
            ),
            linewidth=edge_width(edge.weight) * 0.7,
            alpha=1,
            zorder=1,
        )

    # ---------------------------------------------------------
    # Directed LR arrows
    # ---------------------------------------------------------
    for sender, receiver, edge in G.edges(data=True):
        pair = tuple(sorted((sender, receiver)))
        is_sig = pair in sig_pairs

        curve_direction = (
            1
            if int(
                hashlib.md5(
                    "|".join(pair).encode()
                ).hexdigest(),
                16,
            ) % 2 == 0
            else -1
        )

        ax.add_patch(
            FancyArrowPatch(
                pos[sender],
                pos[receiver],
                connectionstyle=(
                    f"arc3,rad={curve_rad * curve_direction}"
                ),
                arrowstyle="-|>",
                mutation_scale=arrow_size,
                linewidth=edge_width(edge["weight"]),
                color=node_colour(sender) if is_sig else grey,
                alpha=1,
                shrinkA=8,
                shrinkB=18,
                zorder=3,
            )
        )

    # ---------------------------------------------------------
    # Legends
    # ---------------------------------------------------------
    ticks = np.unique(
        np.linspace(wmin, wmax, 3).astype(int)
    )

    edge_legend = ax.legend(
        handles=[
            Line2D(
                [0],
                [0],
                color="black",
                linewidth=edge_width(tick),
                label=f"{tick} LR pairs",
            )
            for tick in ticks
        ],
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        frameon=False,
        fontsize=legendfont,
    )

    ax.add_artist(edge_legend)

    ax.legend(
        handles=[
            Line2D(
                [0], [0],
                marker="o",
                linestyle="",
                markerfacecolor=sign_cols["pos"],
                markeredgecolor="black",
                markersize=8,
                alpha=0.5,
                label=f"C{pls_component} +",
            ),
            Line2D(
                [0], [0],
                marker="o",
                linestyle="",
                markerfacecolor=sign_cols["neg"],
                markeredgecolor="black",
                markersize=8,
                alpha=0.5,
                label=f"C{pls_component} −",
            ),
            Line2D(
                [0], [0],
                marker="o",
                linestyle="",
                markerfacecolor=grey,
                markeredgecolor="black",
                markersize=8,
                label="Not significant",
            ),
        ],
        loc="lower left",
        bbox_to_anchor=(1.02, 0),
        frameon=False,
        fontsize=legendfont,
    )

    plt.subplots_adjust(right=0.75)

    ax.set_title(title, pad=16)
    ax.axis("off")

    if savepath is not None:
        fig.savefig(
            savepath,
            dpi=300,
            bbox_inches="tight",
            pad_inches=0.5,
            transparent=True,
        )

    return fig, ax