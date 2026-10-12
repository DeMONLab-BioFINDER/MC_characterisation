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


#rosmap = sc.read_h5ad(r"C:\Users\gabri\Documents\Data_projects\MCcharactertization\ROSMAP\_core_data\scanpy\Astrocytes.h5ad")
rosmap = sc.read_h5ad('/Users/ga0028vi/Documents/MCcharactertization/ROSMAP/_core_data/scanpy/Astrocytes.h5ad')

rosmap.obs['projid'] = rosmap.obs['projid'].astype(str)
rosmap.obs['newdons'] = rosmap.obs['projid']
rosmap.obs['cellsubtypes'] = rosmap.obs['cell_type_high_resolution']
rosmap.obs['cellsubtypes'] = rosmap.obs['cellsubtypes'].str.replace('_', '', regex=True)
rosmap.X = rosmap.layers['counts']

#seaad = sc.read_h5ad(r"C:\Users\gabri\Documents\Data_projects\MCcharactertization\SEAAD_PFC\_core_data\SEAAD_PFC_finalnuclei_Astrocytes.h5ad")
seaad = sc.read_h5ad('/Users/ga0028vi/Documents/MCcharactertization/SEAAD_PFC/_core_data/SEAAD_PFC_finalnuclei_Astrocytes.h5ad')

seaad.obs['Donor ID'] = seaad.obs['Donor ID'].astype(str)
seaad.obs['Donor ID'] = seaad.obs['Donor ID'].str.replace('.', '')

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

sc.pp.filter_cells(rosmap, min_genes=100)
sc.pp.filter_genes(rosmap, min_cells=3) 
sc.pp.filter_cells(seaad, min_genes=100)
sc.pp.filter_genes(seaad, min_cells=3)

# common genes
genes = rosmap.var_names.intersection(seaad.var_names)
seaad = seaad[:, genes]
rosmap = rosmap[:, genes]

# Normalizing and logtransform 
sc.pp.normalize_total(rosmap)
sc.pp.log1p(rosmap)
sc.pp.normalize_total(seaad)
sc.pp.log1p(seaad)

rosmap.obs['cellsubtypes'] = rosmap.obs['cellsubtypes'].str.replace(' ', '', regex=True)

pseudobulk_means_list = []

groups = rosmap.obs.groupby(["newdons", "cellsubtypes"]).groups  # (donor, subtype) -> row indices

##### PB
### ROSMAP pseudobulk
X = rosmap.X
if not sp.sparse.issparse(rosmap.X):
    rosmap.X = sp.csr_matrix(rosmap.X)

# group codes for subtypes
rosmap.obs['cellsubtypes'] = rosmap.obs['cellsubtypes'].str.replace(' ', '', regex=True)

sub = rosmap.obs["cellsubtypes"].astype("category")
codes = sub.cat.codes.values
levels = list(sub.cat.categories)

n_cells = rosmap.n_obs
n_groups = len(levels)

# indicator matrix: cell -> group
G = sp.sparse.csr_matrix(
    (np.ones(n_cells, dtype=np.float32), (codes, np.arange(n_cells))),
    shape=(n_groups, n_cells)
)

# sums per group (groups x genes), still sparse-ish
sums = G @ rosmap.X

# divide by group sizes to get means
sizes = np.asarray(G.sum(axis=1)).ravel()
means = sums.multiply(1.0 / sizes[:, None])

# convert to a DataFrame (this WILL be dense, but now it's only n_groups x n_genes)
ros_pb = pd.DataFrame(means.toarray(), index=levels, columns=genes)


### SEAAD pseudobulk

X = seaad.X
if not sp.sparse.issparse(X):
    X = sp.csr_matrix(X)

# group codes for subtypes
sub = seaad.obs["cellsubtypes"].astype("category")
codes = sub.cat.codes.values
levels = list(sub.cat.categories)

n_cells = seaad.n_obs
n_groups = len(levels)

# indicator matrix: cell -> group
G = sp.sparse.csr_matrix(
    (np.ones(n_cells, dtype=np.float32), (codes, np.arange(n_cells))),
    shape=(n_groups, n_cells)
)

# sums per group (groups x genes), still sparse-ish
sums = G @ X

# divide by group sizes to get means
sizes = np.asarray(G.sum(axis=1)).ravel()
means = sums.multiply(1.0 / sizes[:, None])

# convert to a DataFrame (this WILL be dense, but now it's only n_groups x n_genes)
sea_pb = pd.DataFrame(means.toarray(), index=levels, columns=genes)



### ALL GENES corr
from sklearn.preprocessing import StandardScaler

ros_z = pd.DataFrame(
    StandardScaler().fit_transform(ros_pb),
    index=ros_pb.index,
    columns=ros_pb.columns,
)

sea_z = pd.DataFrame(
    StandardScaler().fit_transform(sea_pb),
    index=sea_pb.index,
    columns=sea_pb.columns,
)

expr_sim = pd.DataFrame(
    index=ros_z.index,
    columns=sea_z.index,
    dtype=float
)

for r in ros_z.index:
    for s in sea_z.index:
        expr_sim.loc[r, s] = np.corrcoef(
            ros_z.loc[r].values,
            sea_z.loc[s].values
        )[0, 1]
        
plt.figure(figsize=(14,10))
sns.heatmap(expr_sim, vmin=-1, vmax=1, cmap="coolwarm", annot=True)
plt.xlabel("SEAAD subtypes")
plt.ylabel("ROSMAP subtypes")
plt.title("Subtype expression similarity (pseudobulk r val)")
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "all_corr_ps.png"))


### Select and sub  HVGs
sc.pp.highly_variable_genes(rosmap, n_top_genes = 10)
### Select and sub  HVGs
sc.pp.highly_variable_genes(seaad, n_top_genes = 10)

# Subset BOTH rosmap and seaad HVGs
hvg = rosmap.var_names[rosmap.var["highly_variable"]].union(
      seaad.var_names[seaad.var["highly_variable"]]
)
rosmap = rosmap[:, hvg].copy()
seaad  = seaad[:, hvg].copy()

# - OR - 
# Subset just to rosmap or just SEAAD HVGs
#hvg_genes = rosmap.var_names[rosmap.var["highly_variable"]].tolist()
#rosmap  = rosmap[:, hvg_genes].copy()
#seaad  = seaad[:, hvg_genes].copy()



### ROSMAP pseudobulk

genes = rosmap.var_names

X = rosmap.X
if not sp.sparse.issparse(X):
    X = sp.csr_matrix(X)

# group codes for subtypes
sub = rosmap.obs["cellsubtypes"].astype("category")
codes = sub.cat.codes.values
levels = list(sub.cat.categories)

n_cells = rosmap.n_obs
n_groups = len(levels)

# indicator matrix: cell -> group
G = sp.sparse.csr_matrix(
    (np.ones(n_cells, dtype=np.float32), (codes, np.arange(n_cells))),
    shape=(n_groups, n_cells)
)

# sums per group (groups x genes), still sparse-ish
sums = G @ X

# divide by group sizes to get means
sizes = np.asarray(G.sum(axis=1)).ravel()
means = sums.multiply(1.0 / sizes[:, None])

# convert to a DataFrame (this WILL be dense, but now it's only n_groups x n_genes)
ros_pb = pd.DataFrame(means.toarray(), index=levels, columns=genes)


### SEAAD pseudobulk
genes = seaad.var_names

X = seaad.X
if not sp.sparse.issparse(X):
    X = sp.csr_matrix(X)

# group codes for subtypes
sub = seaad.obs["cellsubtypes"].astype("category")
codes = sub.cat.codes.values
levels = list(sub.cat.categories)

n_cells = seaad.n_obs
n_groups = len(levels)

# indicator matrix: cell -> group
G = sp.sparse.csr_matrix(
    (np.ones(n_cells, dtype=np.float32), (codes, np.arange(n_cells))),
    shape=(n_groups, n_cells)
)

# sums per group (groups x genes), still sparse-ish
sums = G @ X

# divide by group sizes to get means
sizes = np.asarray(G.sum(axis=1)).ravel()
means = sums.multiply(1.0 / sizes[:, None])

# convert to a DataFrame (this WILL be dense, but now it's only n_groups x n_genes)
sea_pb = pd.DataFrame(means.toarray(), index=levels, columns=genes)


### HVG GENES corr
from sklearn.preprocessing import StandardScaler

ros_z = pd.DataFrame(
    StandardScaler().fit_transform(ros_pb),
    index=ros_pb.index,
    columns=ros_pb.columns,
)

sea_z = pd.DataFrame(
    StandardScaler().fit_transform(sea_pb),
    index=sea_pb.index,
    columns=sea_pb.columns,
)

expr_sim = pd.DataFrame(
    index=ros_z.index,
    columns=sea_z.index,
    dtype=float
)

for r in ros_z.index:
    for s in sea_z.index:
        expr_sim.loc[r, s] = np.corrcoef(
            ros_z.loc[r].values,
            sea_z.loc[s].values
        )[0, 1]
        
plt.figure(figsize=(14,10))
sns.heatmap(expr_sim, vmin=-1, vmax=1, cmap="coolwarm", annot=True)
plt.xlabel("SEAAD subtypes")
plt.ylabel("ROSMAP subtypes")
plt.title("Subtype expression similarity (pseudobulk r val)")
plt.tight_layout()
plt.savefig(os.path.join(out_dir, "HVGs_corr_ps.png"))