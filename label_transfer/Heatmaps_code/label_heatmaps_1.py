import os 
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy as sp
plt.switch_backend('Agg')  # Non-interactive backend safe for parallel processing


rosmap_path = '/Users/ga0028vi/Documents/MCcharactertization/ROSMAP/_core_data/scanpy/Astrocytes.h5ad'
seaad_path = '/Users/ga0028vi/Documents/MCcharactertization/SEAAD_PFC/_core_data/SEAAD_PFC_finalnuclei_Astrocytes.h5ad'
out_dir = '/Users/ga0028vi/Documents/MCcharactertization/label_transfer/plots/'
os.makedirs(out_dir, exist_ok=True)

# Deal with rosmap
rosmap = sc.read_h5ad(rosmap_path)
rosmap.obs['projid'] = rosmap.obs['projid'].astype(str)
rosmap.obs['newdons'] = rosmap.obs['projid']
rosmap.obs['cellsubtypes'] = rosmap.obs['cell_type_high_resolution']
rosmap.obs['cellsubtypes'] = rosmap.obs['cellsubtypes'].str.replace('_', '', regex=True)
rosmap.X = rosmap.layers['counts']


#deal with seaad
seaad = sc.read_h5ad(seaad_path)

### SEAAD cleanup
cols_to_drop = ['Race (choice=White)',
 'Race (choice=Black/ African American)',
 'Race (choice=Asian)',
 'Race (choice=American Indian/ Alaska Native)',
 'Race (choice=Native Hawaiian or Pacific Islander)',
 'Race (choice=Unknown or unreported)',
 'Race (choice=Other)',
 'specify other race',
 'Hispanic/Latino',
 'Highest level of education',
 'Years of education',
 'PMI',
 'Fresh Brain Weight',
 'Brain pH',
 'Overall AD neuropathological Change',
 'Thal',
 'Braak',
 'CERAD score',
 'Overall CAA Score',
 'Highest Lewy Body Disease',
 'Total Microinfarcts (not observed grossly)',
 'Total microinfarcts in screening sections',
 'Atherosclerosis',
 'Arteriolosclerosis',
 'LATE',
 'Cognitive Status',
 'Last CASI Score',
 'Interval from last CASI in months',
 'Last MMSE Score',
 'Interval from last MMSE in months',
 'Last MOCA Score',
 'Interval from last MOCA in months',
 'APOE Genotype',
 'Primary Study Name',
 'Secondary Study Name',
 'NeuN positive fraction on FANS',
 'RIN',
 'cell_prep_type',
 'facs_population_plan',
 'rna_amplification',
 'sample_name',
 'sample_quantity_count',
 'expc_cell_capture',
 'method',
 'pcr_cycles',
 'percent_cdna_longer_than_400bp',
 'rna_amplification_pass_fail',
 'amplified_quantity_ng',
 'load_name',
 'library_prep',
 'library_input_ng',
 'r1_index',
 'avg_size_bp',
 'quantification_fmol',
 'library_prep_pass_fail',
 'exp_component_vendor_name',
 'batch_vendor_name',
 'experiment_component_failed',
 'alignment',
 'Genome',
 'ar_id',
 'bc',
 'GEX_Estimated_number_of_cells',
 'GEX_number_of_reads',
 'GEX_sequencing_saturation',
 'GEX_Mean_raw_reads_per_cell',
 'GEX_Q30_bases_in_barcode',
 'GEX_Q30_bases_in_read_2',
 'GEX_Q30_bases_in_UMI',
 'GEX_Percent_duplicates',
 'GEX_Q30_bases_in_sample_index_i1',
 'GEX_Q30_bases_in_sample_index_i2',
 'GEX_Reads_with_TSO',
 'GEX_Sequenced_read_pairs',
 'GEX_Valid_UMIs',
 'GEX_Valid_barcodes',
 'GEX_Reads_mapped_to_genome',
 'GEX_Reads_mapped_confidently_to_genome',
 'GEX_Reads_mapped_confidently_to_intergenic_regions',
 'GEX_Reads_mapped_confidently_to_intronic_regions',
 'GEX_Reads_mapped_confidently_to_exonic_regions',
 'GEX_Reads_mapped_confidently_to_transcriptome',
 'GEX_Reads_mapped_antisense_to_gene',
 'GEX_Fraction_of_transcriptomic_reads_in_cells',
 'GEX_Total_genes_detected',
 'GEX_Median_UMI_counts_per_cell',
 'GEX_Median_genes_per_cell',
 'Multiome_Feature_linkages_detected',
 'Multiome_Linked_genes',
 'Multiome_Linked_peaks',
 'ATAC_Confidently_mapped_read_pairs',
 'ATAC_Fraction_of_genome_in_peaks',
 'ATAC_Fraction_of_high_quality_fragments_in_cells',
 'ATAC_Fraction_of_high_quality_fragments_overlapping_TSS',
 'ATAC_Fraction_of_high_quality_fragments_overlapping_peaks',
 'ATAC_Fraction_of_transposition_events_in_peaks_in_cells',
 'ATAC_Mean_raw_read_pairs_per_cell',
 'ATAC_Median_high_quality_fragments_per_cell',
 'ATAC_Non-nuclear_read_pairs',
 'ATAC_Number_of_peaks',
 'ATAC_Percent_duplicates',
 'ATAC_Q30_bases_in_barcode',
 'ATAC_Q30_bases_in_read_1',
 'ATAC_Q30_bases_in_read_2',
 'ATAC_Q30_bases_in_sample_index_i1',
 'ATAC_Sequenced_read_pairs',
 'ATAC_TSS_enrichment_score',
 'ATAC_Unmapped_read_pairs',
 'ATAC_Valid_barcodes',
 'Number of mapped reads',
 'Number of unmapped reads',
 'Number of multimapped reads',
 'Number of reads',
 'Number of UMIs',
 'Genes detected',
 'Doublet score',
 'Fraction mitochondrial UMIs',
 'Used in analysis',
 'Class confidence',
 'Class',
 'Subclass confidence','Supertype confidence',
 'Supertype (non-expanded)','Continuous Pseudo-progression Score','Organism','Neurotypical reference',
 'Brain Region', 'Sex', 'Gender', 'Age at Death',
 'Severely Affected Donor']
seaad.obs = seaad.obs.drop(columns=cols_to_drop, errors="ignore")
del seaad.uns
del seaad.obsm
del seaad.obsp

seaad.X = seaad.layers['UMIs']
#seaad.obs["major_cell_type"] = "OPC"
seaad.obs["cellsubtypes"] = seaad.obs['Supertype']
seaad.obs['newdons'] = seaad.obs['Donor ID']


# common genes
genes = rosmap.var_names.intersection(seaad.var_names)
seaad = seaad[:, genes]
rosmap = rosmap[:, genes]



##### LOOK FOR MARKER GENES ######
# Normalizing to median total counts
sc.pp.normalize_total(rosmap)
# Logarithmize the data
sc.pp.log1p(rosmap)
sc.tl.rank_genes_groups(rosmap, groupby="cellsubtypes", method="wilcoxon")

sc.tl.filter_rank_genes_groups(
    rosmap,
    key="rank_genes_groups",
    key_added="rank_genes_groups_filtered"
)

sc.pl.rank_genes_groups_dotplot(
    rosmap,
    groupby="cellsubtypes",
    key="rank_genes_groups_filtered",
    n_genes=15,
    standard_scale="var"
)



#### SEAAD
# Normalizing to median total counts
sc.pp.normalize_total(seaad)
# Logarithmize the data
sc.pp.log1p(seaad)

sc.tl.rank_genes_groups(seaad, groupby="cellsubtypes", method="wilcoxon")

sc.tl.filter_rank_genes_groups(
    seaad,
    key="rank_genes_groups",
    key_added="rank_genes_groups_filtered")

sc.pl.rank_genes_groups_dotplot(
    seaad,
    groupby="cellsubtypes",
    key="rank_genes_groups_filtered",
    n_genes=15,
    standard_scale="var"
)


#### PLOT TOP 10 genes overlap
n_top=10
groups = (
    list(rosmap.obs['cellsubtypes'].cat.categories)
)

rosmap_markers = {}
for g in groups:
    df = sc.get.rank_genes_groups_df(rosmap, group=g, key='rank_genes_groups_filtered')

    # filter_rank_genes_groups keeps rows but sets names to NA -> drop them
    df = df.dropna(subset=["names"]).head(n_top).copy()

    rosmap_markers[g] = df

groups = (
    list(seaad.obs['cellsubtypes'].cat.categories)
)

seaad_markers = {}
for g in groups:
    df = sc.get.rank_genes_groups_df(seaad, group=g, key='rank_genes_groups_filtered')

    # filter_rank_genes_groups keeps rows but sets names to NA -> drop them
    df = df.dropna(subset=["names"]).head(n_top).copy()

    seaad_markers[g] = df
    
# make sure your keys are in a stable order
ros_groups = list(rosmap_markers.keys())
sea_groups = list(seaad_markers.keys())

overlap = pd.DataFrame(index=ros_groups, columns=sea_groups, dtype=int)

for r_sub in ros_groups:
    r_genes = set(rosmap_markers[r_sub]["names"])
    for s_sub in sea_groups:
        s_genes = set(seaad_markers[s_sub]["names"])
        overlap.loc[r_sub, s_sub] = len(r_genes & s_genes)

overlap

overlap_norm_ntop = overlap / n_top


plt.figure(figsize=(14,10))
sns.heatmap(overlap_norm_ntop, cmap="viridis", vmin=0, vmax=1, annot=True)
plt.title("Exc  ROS seaAD detected marker genes overlap n=10")
plt.xlabel("SEAAD")
plt.ylabel("ROSMAP")
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "markergenes10_corr.png"))




#### PLOT TOP 50 genes overlap
n_top=50
groups = (
    list(rosmap.obs['cellsubtypes'].cat.categories)
)

rosmap_markers = {}
for g in groups:
    df = sc.get.rank_genes_groups_df(rosmap, group=g, key='rank_genes_groups_filtered')

    # filter_rank_genes_groups keeps rows but sets names to NA -> drop them
    df = df.dropna(subset=["names"]).head(n_top).copy()

    rosmap_markers[g] = df

groups = (
    list(seaad.obs['cellsubtypes'].cat.categories)
)

seaad_markers = {}
for g in groups:
    df = sc.get.rank_genes_groups_df(seaad, group=g, key='rank_genes_groups_filtered')

    # filter_rank_genes_groups keeps rows but sets names to NA -> drop them
    df = df.dropna(subset=["names"]).head(n_top).copy()

    seaad_markers[g] = df
    
# make sure your keys are in a stable order
ros_groups = list(rosmap_markers.keys())
sea_groups = list(seaad_markers.keys())

overlap = pd.DataFrame(index=ros_groups, columns=sea_groups, dtype=int)

for r_sub in ros_groups:
    r_genes = set(rosmap_markers[r_sub]["names"])
    for s_sub in sea_groups:
        s_genes = set(seaad_markers[s_sub]["names"])
        overlap.loc[r_sub, s_sub] = len(r_genes & s_genes)

overlap

overlap_norm_ntop = overlap / n_top


plt.figure(figsize=(14,10))
sns.heatmap(overlap_norm_ntop, cmap="viridis", vmin=0, vmax=1, annot=True)
plt.title("Exc  ROS seaAD detected marker genes overlap n=10")
plt.xlabel("SEAAD")
plt.ylabel("ROSMAP")
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "markergenes50_corr.png"))
    


#### PLOT TOP 100 genes overlap
n_top=100
groups = (
    list(rosmap.obs['cellsubtypes'].cat.categories)
)

rosmap_markers = {}
for g in groups:
    df = sc.get.rank_genes_groups_df(rosmap, group=g, key='rank_genes_groups_filtered')

    # filter_rank_genes_groups keeps rows but sets names to NA -> drop them
    df = df.dropna(subset=["names"]).head(n_top).copy()

    rosmap_markers[g] = df

groups = (
    list(seaad.obs['cellsubtypes'].cat.categories)
)

seaad_markers = {}
for g in groups:
    df = sc.get.rank_genes_groups_df(seaad, group=g, key='rank_genes_groups_filtered')

    # filter_rank_genes_groups keeps rows but sets names to NA -> drop them
    df = df.dropna(subset=["names"]).head(n_top).copy()

    seaad_markers[g] = df
    
# make sure your keys are in a stable order
ros_groups = list(rosmap_markers.keys())
sea_groups = list(seaad_markers.keys())

overlap = pd.DataFrame(index=ros_groups, columns=sea_groups, dtype=int)

for r_sub in ros_groups:
    r_genes = set(rosmap_markers[r_sub]["names"])
    for s_sub in sea_groups:
        s_genes = set(seaad_markers[s_sub]["names"])
        overlap.loc[r_sub, s_sub] = len(r_genes & s_genes)

overlap

overlap_norm_ntop = overlap / n_top


plt.figure(figsize=(14,10))
sns.heatmap(overlap_norm_ntop, cmap="viridis", vmin=0, vmax=1, annot=True)
plt.title("Exc  ROS seaAD detected marker genes overlap n=10")
plt.xlabel("SEAAD")
plt.ylabel("ROSMAP")
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "markergenes100_corr.png"))


