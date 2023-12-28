## Installation
# install.packages("remotes")
# library(remotes)
# 
# remotes::install_github("carmonalab/UCell")
# remotes::install_github("carmonalab/scGate")
# 
# remotes::install_github("carmonalab/ProjecTILs")

library(ggplot2)
library(ProjecTILs)
library(gridExtra)

table(query.object$Time)

ref <- readRDS("/Users/bianjh/Documents/R files/Kanakry/Manuscript_Revision/Project_TILs/ref_TILAtlas_mouse_v1.rds")

refCols <- c("#edbe2a", "#A58AFF", "#53B400", "#F8766D", "#00B6EB", "#d1cfcc", "#FF0000", "#87f6a5", "#e812dd")
DimPlot(ref,label = T, cols = refCols)

# Load query sample from package
querydata <- ProjecTILs::query_example_seurat

# Load query data from GEO
library(GEOquery)
geo_acc <- "GSE86028"
getGEOSuppFiles(geo_acc)

fname3 <- sprintf("%s/GSE86028_TILs_sc_wt_mtko.tpm.log2.txt.gz", geo_acc)
querydata3 <- read.sc.query(fname3, type = "raw.log2")

# Kanakry's CART mice data
kanakry_CART <- SO

# Project query onto reference data
query.projected <- make.projection(kanakry_CART, ref=ref)

plot.projection(ref, query.projected)

# Predict counts based on nearest neighbor
query.projected <- cellstate.predict(ref=ref, query=query.projected)
table(query.projected$functional.cluster)

kanakry_CART@meta.data$ProjectTILS_calls <- query.projected$functional.cluster[match(rownames(kanakry_CART@meta.data),names(query.projected$functional.cluster))]

write.csv(kanakry_CART@meta.data,"ProjectTIL_metadata.csv")
