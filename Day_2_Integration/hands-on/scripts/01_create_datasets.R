#--------------------------Create data sets script-----------------------------#
# Author: António Sousa (e-mail: aggode@utu.fi)
# Description: create data sets for the course: 
#(1) w/o shared cell types
#(2) w/ shared cell types
#(3) w/ some shared cell types
#(4) reference-mapping task
# Date: 22/06/2024
# Last update: 30/06/20224
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## Import packages & set env (R v.4.1.0)

# Packages
library("Seurat") # v.5.1.0
if (!"SeuratData" %in% installed.packages()) devtools::install_github("satijalab/seurat-data") # if isn't installed, install 'SeuratData' 
library("SeuratData") # v.0.2.2.9001
library("Azimuth")

# Set WD
setwd("scripts")

# Set seed
set.seed(1024)

# Define output directory
data.dir <- "../data"

# Download options
options(timeout = 60*10)
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## (1) w/o shared cell types: pbmc3k x panc8

# Create directory to save files
output <- file.path(data.dir, "pbmc3k_panc8.rds")

# Install data
InstallData("pbmc3k")
InstallData("panc8")

# Load data & subset data sets 
pbmc3k <- LoadData("pbmc3k")
panc8 <- LoadData("panc8")
panc8 <- subset(panc8, cells = colnames(panc8)[panc8$dataset=="indrop1"]) # select only batch sample from 'indrop1' tech
gc(); 

# Merge them
pbmc3k$batch <- "pbmc"
colnames(pbmc3k@meta.data)[4] <- "cell_type"
panc8$batch <- "pancreas"
colnames(panc8@meta.data)[7] <- "cell_type"
panc8@meta.data <- panc8@meta.data[,c(1:3, 7, 9)]
table(row.names(panc8) %in% row.names(pbmc3k)) # check no. of genes shared between data sets
seu <- merge(x = pbmc3k, y = panc8, add.cell.ids = c("pbmc", "pancreas"), 
             project = "pbmc3k_panc8")
rm(list = c("pbmc3k", "panc8")); gc(); 
seu[['RNA']]  <- JoinLayers(object = seu[['RNA']])

# Save object
saveRDS(object=seu, file=output)

# Clean env vars
rm(list = ls()[ls() != "data.dir"]); gc(); 
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## (2) w/ shared cell types: ifnb

# Create directory to save files
output <- file.path(data.dir, "ifnb.rds")

# Install data
InstallData("ifnb") # run it only if not installed above

# Load data
ifnb <- LoadData("ifnb")

# Save object
saveRDS(object=ifnb, file=output)

# Clean env vars
rm(list = ls()[ls() != "data.dir"]); gc(); 
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## (3) w/ some shared cell types: jurkat

# Create file name to save Seurat file 
output <- file.path(data.dir, "jurkat.rds")

# List of data sets to download
data2down <- list(
  "t293" = "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/293t/293t_filtered_gene_bc_matrices.tar.gz", 
  "jurkat" = "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/jurkat/jurkat_filtered_gene_bc_matrices.tar.gz", 
  "jurkat_t293_50:50" = "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/jurkat:293t_50:50/jurkat:293t_50:50_filtered_gene_bc_matrices.tar.gz"
)

# Download & import data sets 1-by-1
seu.list <- list()
for (f in names(data2down)) {
  f.path <- gsub(".tar.gz", "", file.path(data.dir, basename(data2down[[f]])))
  untar.path <- file.path(data.dir, "filtered_matrices_mex/hg19")
  cat("Downloading data set", f, "to", paste0(f.path, ".tar.gz"), "...\n")
  download.file(url=data2down[[f]], destfile=paste0(f.path, ".tar.gz"))
  cat("Untar compressed data set", f, "to", untar.path, "...\n")
  untar(tarfile=paste0(f.path, ".tar.gz"), exdir=data.dir)
  cat("Importing gene expression matrix for dat set", f,"from", 
      untar.path,"...\n")
  seu.list[[f]] <- CreateSeuratObject(counts = Read10X(data.dir = untar.path))
  cat("Cleaning folder...\n")
  file.remove(paste0(f.path, ".tar.gz")); unlink(dirname(untar.path), recursive=TRUE);
  cat("Finished!\n\n")
}

# Check if gene names are exactly the same & add 'batch' labels to 'meta.data'
g.list <- lapply(seu.list, row.names)
for (g in length(g.list)) {
  for (i in g.list[-g]) {
    stopifnot(all(g.list[[g]]==i))
  } 
}
for (d in names(seu.list)) {
  seu.list[[d]]$batch <- d
}

# Merge seurat list
seu <- merge(x = seu.list$t293, y = list(seu.list$jurkat, seu.list$`jurkat_t293_50:50`), 
             add.cell.ids = c("t293", "jurkat", "jurkat_t293_50:50"), project = "jurkat")
rm(seu.list); gc(); 
seu[['RNA']]  <- JoinLayers(object = seu[['RNA']])

# Cell type
counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
seu$jurkat.pos.CD3D <- (counts["CD3D",]>0)
seu$t293.pos.XIST <- (counts["XIST",]>0)
seu$cell_type <- ifelse(seu$jurkat.pos.CD3D & !seu$t293.pos.XIST, "Jurkat", 
                        ifelse(!seu$jurkat.pos.CD3D & seu$t293.pos.XIST, "293T", 
                               NA))	
rm(counts); gc(); 
seu <- subset(seu, cells = colnames(seu)[!is.na(seu$cell_type)])

# Save object
saveRDS(object=seu, file=output)

# Clean env vars
rm(list = ls()[ls() != "data.dir"]); gc(); 
#
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
#
## (4) reference-mapping task

# Create file name to save Seurat file 
output.query <- file.path(data.dir, "covid.rds")
output.ref <- file.path(data.dir, "pbmcref.rds")

## Query
# Download file from https://cellxgene.cziscience.com - study from Guo et al. 
data2down <- "https://datasets.cellxgene.cziscience.com/c6371ff2-d6b8-42e4-8f5d-099e25b4de8e.rds"
download.file(url = data2down, destfile = output.query)
query <- readRDS(file = output.query)
stopifnot(all(row.names(query) == row.names(query@assays$RNA@meta.features)))
counts <- query@assays$RNA@counts
data <- query@assays$RNA@data
row.names(counts) <- row.names(data) <- query@assays$RNA@meta.features$feature_name
query.new <- CreateSeuratObject(counts = counts, data = data, meta.data = query@meta.data)
query.new@reductions <- query@reductions
saveRDS(object = query.new, file = output.query)

## Reference
InstallData("pbmcref")
ref <- LoadData("pbmcref", type = "azimuth")
saveRDS(object = ref$map, file = file.path(data.dir, "pbmcref.rds"))

# Clean env vars
rm(list = ls()[ls() != "data.dir"]); gc(); 
#
#------------------------------------------------------------------------------#