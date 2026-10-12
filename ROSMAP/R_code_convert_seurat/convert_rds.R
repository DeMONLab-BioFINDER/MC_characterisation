library("tidyverse")
library(SeuratObject)
library(Seurat)
#library(SeuratDisk)
#library(SeuratData)
library(MuDataSeurat)
#library(SingleCellExperiment)

data <- read_rds("C:/Users/gabri/Documents/Data_projects/MCcharactertization/ROSMAP_MTG/RDS/Midtemporal_cortex.rds")
MuDataSeurat::WriteH5AD(data, "C:/Users/gabri/Documents/Data_projects/MCcharactertization/ROSMAP_MTG/_core_data/Midtemporal_cortex.h5ad", assay="RNA")
rm(data)

