import os
import scanpy as sc

file_dir  = "/proj/nobackup/sens2023026/wharf/gabriv/gabriv-sens2023026/MCcharacterization/SEAAD_MTG/"
file_name = "SEAAD_MTG.h5ad"
in_path = os.path.join(file_dir, file_name)

print("Reading in the file")
adata = sc.read_h5ad(in_path)
print("File loaded")


# define masks
print("Creating masks")
mask_exc = (adata.obs["Class"] == "Neuronal: Glutamatergic")
mask_inh = (adata.obs["Class"] == "Neuronal: GABAergic")

mask_astro = (adata.obs["Subclass"] == "Astrocyte")
mask_opc   = (adata.obs["Subclass"] == "OPC")
mask_oligo = (adata.obs["Subclass"] == "Oligodendrocyte")

mask_immune = adata.obs["Subclass"].isin(["Microglia-PVM"])

mask_vasc = adata.obs["Subclass"].isin([
    "Endothelial",
    "VLMC",
])

print("Maskss done")


base = os.path.splitext(file_name)[0]
print("Saving")
adata[mask_exc].copy().write_h5ad(os.path.join(file_dir, f"{base}_Excitatories.h5ad"))
print("Exc saved")

adata[mask_inh].copy().write_h5ad(os.path.join(file_dir, f"{base}_Inhibitories.h5ad"))
print("Inh saved")

adata[mask_astro].copy().write_h5ad(os.path.join(file_dir, f"{base}_Astrocytes.h5ad"))
print("Ast saved")

adata[mask_opc].copy().write_h5ad(os.path.join(file_dir, f"{base}_OPCs.h5ad"))
print("OPC saved")

adata[mask_oligo].copy().write_h5ad(os.path.join(file_dir, f"{base}_Oligodendrocytes.h5ad"))
print("Oli saved")

adata[mask_immune].copy().write_h5ad(os.path.join(file_dir, f"{base}_Immune.h5ad"))
print("Immun saved")

adata[mask_vasc].copy().write_h5ad(os.path.join(file_dir, f"{base}_Vasculature.h5ad"))
print("Vasc saved")

