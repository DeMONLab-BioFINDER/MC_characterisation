import os
import time
import requests
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path
from io import StringIO
from matplotlib.lines import Line2D
import gseapy as gp


def query_string_network_with_channels(
    genes,
    species=9606,
    required_score=400,
    network_type="functional",
    channels=None,
    # caller_identity="gabriele_AGM_pipeline"
):
    """
    Query STRING network API with optional evidence-channel on/off controls.

    Parameters
    ----------
    genes : list
        Gene/protein names to query.

    species : int
        NCBI taxonomy ID. Human = 9606.

    required_score : int
        STRING confidence cutoff, from 0 to 999.

    network_type : str
        Usually "functional" or "physical".

    channels : dict or None
        Dictionary controlling which STRING evidence channels are on/off.

        Expected keys:
        - "neighborhood"
        - "fusion"
        - "cooccurrence"
        - "coexpression"
        - "experiments"
        - "database"
        - "textmining"

        Values should be True/False.

        Example:
        channels = {
            "neighborhood": False,
            "fusion": False,
            "cooccurrence": False,
            "coexpression": True,
            "experiments": True,
            "database": False,
            "textmining": False
        }

    Returns
    -------
    pd.DataFrame
        STRING network edge table.
    """

    string_api_url = "https://string-db.org/api"
    output_format = "tsv-no-header"
    method = "network"

    request_url = "/".join([string_api_url, output_format, method])

    genes_str = "\r".join(genes)

    params = {
        "identifiers": genes_str,
        "species": species,
        "required_score": required_score,
        "network_type": network_type,
        # "caller_identity": caller_identity
    }

    # STRING channel mapping
    channel_map = {
        "neighborhood": "channel1",   # nscore
        "fusion": "channel2",         # fscore
        "cooccurrence": "channel3",   # pscore
        "coexpression": "channel4",   # ascore
        "experiments": "channel5",    # escore
        "database": "channel6",       # dscore
        "textmining": "channel7"      # tscore
    }

    if channels is not None:
        for channel_name, channel_param in channel_map.items():
            use_channel = channels.get(channel_name, True)
            params[channel_param] = "on" if use_channel else "off"

    response = requests.post(request_url, data=params, timeout=30)
    response.raise_for_status()

    if response.text.strip() == "":
        return pd.DataFrame()

    network_df = pd.read_csv(StringIO(response.text), sep="\t", header=None)

    network_df.columns = [
        "stringId_A",
        "stringId_B",
        "preferredName_A",
        "preferredName_B",
        "ncbiTaxonId",
        "score",
        "nscore",  # neighborhood
        "fscore",  # fusion
        "pscore",  # co-occurrence
        "ascore",  # co-expression
        "escore",  # experimental
        "dscore",  # database
        "tscore"   # text-mining
    ]

    return network_df

def get_string_connected_groups(network_df, genes=None, include_singletons=True):
    """
    Find groups of genes that are connected to each other in a STRING network.

    Parameters
    ----------
    network_df : pd.DataFrame
        STRING network dataframe returned by your API query.

    genes : list or None
        Original AGM gene list. Needed if you want to include genes with no edges.

    include_singletons : bool
        If True, genes from the original AGM list that do not appear in any edge
        are included as isolated genes.

    Returns
    -------
    groups_df : pd.DataFrame
        Table showing which STRING connected group each gene belongs to.
    """

    G = nx.Graph()

    # Add edges from STRING table
    if not network_df.empty:
        for _, row in network_df.iterrows():
            G.add_edge(row["preferredName_A"], row["preferredName_B"])

    # Add all original genes so isolated genes are not lost
    if include_singletons and genes is not None:
        G.add_nodes_from(genes)

    rows = []

    components = list(nx.connected_components(G))
    components = sorted(components, key=len, reverse=True)

    for i, component in enumerate(components, start=1):
        group_id = f"Group_{i}"
        group_size = len(component)

        for gene in sorted(component):
            rows.append({
                "AGM_gene": gene,
                "STRING_group": group_id,
                "group_size": group_size,
                "is_singleton": group_size == 1
            })

    groups_df = pd.DataFrame(rows)

    return groups_df

def save_string_network_images(
    genes,
    out_prefix,
    species=9606,
    required_score=400,
    network_type="functional",
    network_flavor="evidence",
    channels=None,
    add_color_nodes=0,
    add_white_nodes=0,
    hide_disconnected_nodes=True,
    save_svg=True,
    save_png=True,
    save_highres_png=False,
    caller_identity=None,
    timeout=30
):
    """
    Query STRING and save the rendered network as SVG and/or PNG.

    Parameters
    ----------
    genes : list
        Gene/protein names.

    out_prefix : str or Path
        Output path WITHOUT file extension.
        Example:
        out_dir / "Astro_AGM1_STRING_network"

        This will create:
        - Astro_AGM1_STRING_network.svg
        - Astro_AGM1_STRING_network.png

    species : int
        NCBI taxonomy ID. Human = 9606.

    required_score : int
        STRING confidence cutoff, from 0 to 999.

    network_type : str
        "functional" or "physical".

    network_flavor : str
        "evidence", "confidence", or "actions".
        evidence = colored evidence-channel edges.
        confidence = edge thickness reflects confidence.

    channels : dict or None
        Optional evidence-channel on/off controls.

        Keys:
        - "neighborhood"
        - "fusion"
        - "cooccurrence"
        - "coexpression"
        - "experiments"
        - "database"
        - "textmining"

        Values:
        True  = use this evidence channel
        False = turn this evidence channel off

    add_color_nodes : int
        Adds extra colored STRING neighbour nodes.
        Use 0 if you only want your input genes.

    add_white_nodes : int
        Adds extra white STRING neighbour nodes.
        Use 0 if you only want your input genes.

    hide_disconnected_nodes : bool
        If True, STRING hides input genes/proteins that have no network edges.

    save_svg : bool
        If True, save SVG.

    save_png : bool
        If True, save normal PNG.

    save_highres_png : bool
        If True, save high-resolution PNG.

    caller_identity : str or None
        Optional identifier for STRING API usage.

    timeout : int
        Request timeout in seconds.
    """

    string_api_url = "https://string-db.org/api"
    method = "network"

    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    genes_str = "\r".join(genes)

    base_params = {
        "identifiers": genes_str,
        "species": species,
        "required_score": required_score,
        "network_type": network_type,
        "network_flavor": network_flavor,
        "add_color_nodes": add_color_nodes,
        "add_white_nodes": add_white_nodes,
        "hide_disconnected_nodes": 1 if hide_disconnected_nodes else 0,
    }

    if caller_identity is not None:
        base_params["caller_identity"] = caller_identity

    channel_map = {
        "neighborhood": "channel1",
        "fusion": "channel2",
        "cooccurrence": "channel3",
        "coexpression": "channel4",
        "experiments": "channel5",
        "database": "channel6",
        "textmining": "channel7",
    }

    if channels is not None:
        for channel_name, channel_param in channel_map.items():
            use_channel = channels.get(channel_name, True)
            base_params[channel_param] = "on" if use_channel else "off"

    output_jobs = []

    if save_svg:
        output_jobs.append(("svg", out_prefix.with_suffix(".svg")))

    if save_png:
        output_jobs.append(("image", out_prefix.with_suffix(".png")))

    if save_highres_png:
        output_jobs.append(("highres_image", Path(str(out_prefix) + "_highres.png")))

    if len(output_jobs) == 0:
        print("Nothing saved because save_svg, save_png, and save_highres_png are all False.")
        return

    for output_format, out_file in output_jobs:
        request_url = "/".join([string_api_url, output_format, method])

        response = requests.post(request_url, data=base_params, timeout=timeout)
        response.raise_for_status()

        with open(out_file, "wb") as f:
            f.write(response.content)

        print(f"Saved STRING {output_format} to: {out_file}")

        time.sleep(1)

def query_string_geneset_description(
    genes,
    species=9606,
    caller_identity="gabriele_AGM_pipeline",
    timeout=30
):
    """
    Query STRING geneset_description API for a set of genes/proteins.

    Returns
    -------
    pd.DataFrame
        Usually one row with:
        primary_description, secondary_description, tertiary_description, proteins
    """

    if len(genes) < 2:
        return pd.DataFrame([{
            "primary_description": "-",
            "secondary_description": "-",
            "tertiary_description": "-",
            "proteins": ",".join(genes)
        }])

    string_api_url = "https://version-12-0.string-db.org/api"
    output_format = "json"
    method = "geneset_description"

    request_url = "/".join([string_api_url, output_format, method])

    params = {
        "identifiers": "\r".join(genes),
        "species": species,
        "caller_identity": caller_identity
    }

    response = requests.post(request_url, data=params, timeout=timeout)
    response.raise_for_status()

    data = response.json()

    if len(data) == 0:
        return pd.DataFrame([{
            "primary_description": "-",
            "secondary_description": "-",
            "tertiary_description": "-",
            "proteins": ",".join(genes)
        }])

    return pd.DataFrame(data)


def annotate_string_connected_groups(
    groups_df,
    module_name=None,
    species=9606,
    caller_identity="gabriele_AGM_pipeline",
    sleep_time=1
):
    """
    Annotate each connected STRING group using STRING geneset_description API.

    Parameters
    ----------
    groups_df : pd.DataFrame
        Output from get_string_connected_groups().
        Must contain:
        - AGM_gene
        - STRING_group
        - group_size

    module_name : str or None
        Optional AGM/module name to add to output.

    Returns
    -------
    group_descriptions_df : pd.DataFrame
        One row per connected group.
    """

    rows = []

    for group_id, sub_df in groups_df.groupby("STRING_group"):
        group_genes = (
            sub_df["AGM_gene"]
            .dropna()
            .astype(str)
            .unique()
            .tolist()
        )

        desc_df = query_string_geneset_description(
            genes=group_genes,
            species=species,
            caller_identity=caller_identity
        )

        desc_row = desc_df.iloc[0].to_dict()

        rows.append({
            "AGM_ID": module_name,
            "STRING_group": group_id,
            "group_size": len(group_genes),
            "genes": ";".join(sorted(group_genes)),
            "primary_description": desc_row.get("primary_description", "-"),
            "secondary_description": desc_row.get("secondary_description", "-"),
            "tertiary_description": desc_row.get("tertiary_description", "-"),
            "STRING_returned_proteins": desc_row.get("proteins", "-")
        })

        time.sleep(sleep_time)

    group_descriptions_df = pd.DataFrame(rows)

    return group_descriptions_df

def run_enrichr_for_unannotated_groups(
    unannotated_groups,
    module_name,
    ref_genes,
    libraries,
    top_n_to_save=20,
    FDR_filter=0.05,
    sleep_time=2
):
    all_enrichr_results = []

    for _, row in unannotated_groups.iterrows():

        target_group = row["STRING_group"]

        # Split genes like: "MYBPC1;MYLPF;PDLIM3"
        group_genes = str(row["genes"]).split(";")

        print(
            f"Running {module_name} | {target_group} | "
            f"{len(group_genes)} group genes | "
            f"{len(ref_genes)} background genes"
        )

        # Skip tiny groups
        if len(group_genes) < 3:
            print(f"Skipping {target_group}: fewer than 3 genes")
            continue

        try:
            enr = gp.enrichr(
                gene_list=group_genes,
                gene_sets=libraries,
                organism="human",
                background=ref_genes,
                outdir=None,
                cutoff=FDR_filter
            )

            res = enr.results.copy()

            # Keep only significant terms
            res = res.loc[res["Adjusted P-value"] < FDR_filter].copy()

            if res.empty:
                print(f"No significant terms for {module_name} | {target_group}")
                continue

            # Keep top N terms per database for this group
            res = (
                res
                .sort_values(["Gene_set", "Adjusted P-value", "P-value"])
                .groupby("Gene_set", as_index=False)
                .head(top_n_to_save)
            )

            # Add metadata
            res.insert(0, "AGM_ID", module_name)
            res.insert(1, "STRING_group", target_group)
            res.insert(2, "group_size", row["group_size"])
            res.insert(5, "n_background_genes", len(ref_genes))

            all_enrichr_results.append(res)

        except Exception as e:
            print(f"Failed {module_name} | {target_group}: {e}")

        time.sleep(sleep_time)

    # Combine results
    if all_enrichr_results:
        unannotated_group_enrichr_df = pd.concat(
            all_enrichr_results,
            ignore_index=True
        )
    else:
        unannotated_group_enrichr_df = pd.DataFrame()

    return unannotated_group_enrichr_df