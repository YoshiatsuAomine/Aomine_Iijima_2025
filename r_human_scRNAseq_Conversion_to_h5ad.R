
# ---------- Install the packages ----------
install.packages("promises")
install.packages("data.table") 
install.packages("hdf5r")      
install.packages("Rtsne")      
install.packages("ggplot2")     


install.packages("Seurat")
packageVersion("Seurat")

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("ncborcherding/scRepertoire")

devtools::install_github("10XGenomics/loupeR")

install.packages("R.utils")

if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
  install.packages("remotes")
  remotes::install_github("mojaveazure/seurat-disk")
}

# ---------- Run ----------
# Load the packages
library(Seurat)
library(scRepertoire)
library(loupeR)
loupeR::setup()
library(Matrix)
library(dplyr)
library(readr)
library(hdf5r)
library(data.table)
library(Rtsne)
library(ggplot2)
library(R.utils)
library(SeuratDisk)


setwd("D:/aomine/Single_cell_Iijima/scRNA-seq_Neuron2022")

meta_file <- "GSE197289_snRNA-seq_human_barcode_meta.csv.gz"
counts_file <- "GSE197289_snRNA-seq_human_raw_counts.RDS.gz"

meta_data <- read.csv(gzfile(meta_file), stringsAsFactors = FALSE)

tmp_file <- tempfile(fileext = ".rds")
gunzip(counts_file, destname = tmp_file, overwrite = TRUE, remove = FALSE)
counts_matrix <- readRDS(tmp_file)

rownames(meta_data) <- meta_data$V1

if (!all(meta_data$V1 %in% colnames(counts_matrix))) {
  stop("meta_data$V1 does not match the cell names in counts_matrix.")
}

rownames(meta_data) <- meta_data$V1
meta_data <- meta_data[colnames(counts_matrix), , drop = FALSE]

seurat_obj <- CreateSeuratObject(counts = counts_matrix, meta.data = meta_data)

convert_to_legacy_assay <- function(seurat_obj) {
  assay_data <- GetAssayData(seurat_obj, layer = "counts")
  
  if ("data" %in% Layers(seurat_obj[["RNA"]])) {
    norm_data <- GetAssayData(seurat_obj, layer = "data")
    if (all(norm_data == 0)) {
      message("The ‘data’ layer is empty. Only counts will be used.")
      legacy_assay <- CreateAssayObject(counts = assay_data)
    } else {
      legacy_assay <- CreateAssayObject(counts = assay_data, data = norm_data)
    }
  } else {
    message("The ‘data’ layer does not exist. Only ‘counts’ will be used.")
    legacy_assay <- CreateAssayObject(counts = assay_data)
  }
  
  seurat_obj[["RNA"]] <- legacy_assay
  return(seurat_obj)
}

seurat_obj <- convert_to_legacy_assay(seurat_obj)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
DefaultAssay(seurat_obj) <- "RNA"

cat("Seurat Object Metadata Column:\n")
print(colnames(seurat_obj@meta.data))
cat("\The first five lines:\n")
print(head(seurat_obj@meta.data))

SaveH5Seurat(seurat_obj, filename = "Neuron2022_Human_TG.h5Seurat", overwrite = TRUE)

Convert("Neuron2022_Human_TG.h5Seurat", dest = "h5ad", overwrite = TRUE)

cat("Conversion to h5ad complete: Neuron2022_Human_TG.h5ad'\n")


