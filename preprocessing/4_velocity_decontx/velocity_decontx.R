#!/usr/bin/env Rscript

# decontX on the stacked spliced/unspliced/ambiguous counts of one sample, as in
# nf-core/scdownstream's CELDA_DECONTX module (modules/local/celda/decontx), with
# the sample's empty droplets as background. Run by velocity_decontx.py in the
# decontx container (containers/decontx), which writes the inputs.
#
# Usage: Rscript velocity_decontx.R <cells h5ad> <background h5ad> <output h5ad>
#
# Inputs have the stacked counts in X (droplets x features). The output has the
# corrected counts in X, with the same cells and features as <cells h5ad>, and
# the contamination per cell in obs$decontX_contamination. decontX clusters
# the cells itself.

suppressPackageStartupMessages({
    library(anndataR)
    library(celda)
    library(SingleCellExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)

cells <- read_h5ad(args[1])$as_SingleCellExperiment(x_mapping = "counts", assays_mapping = FALSE)
background <- read_h5ad(args[2])$as_SingleCellExperiment(x_mapping = "counts", assays_mapping = FALSE)
stopifnot(identical(rownames(cells), rownames(background)))

corrected <- decontX(cells, assayName = "counts", background = background)

output <- AnnData(
    X = Matrix::t(assay(corrected, "decontXcounts")),
    obs = data.frame(decontX_contamination = corrected$decontX_contamination,
                     row.names = colnames(corrected)),
    var = data.frame(row.names = rownames(corrected))
)
write_h5ad(output, args[3], compression = "gzip")

cat("cells:", ncol(corrected), "| background droplets:", ncol(background),
    "| median contamination:", median(corrected$decontX_contamination),
    "| celda", as.character(packageVersion("celda")), "\n")
