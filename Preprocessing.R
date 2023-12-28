
library(Seurat)
library(reshape2)
library(tidyverse)
library(gridExtra)
library(RColorBrewer)
library(stringr)
library(ggplot2)

localFilePaths <- paste("/rstudio-files/ccbr-data/users/Jing/CCBR_1171/data/",
                        list.files("/rstudio-files/ccbr-data/users/Jing/CCBR_1171/data"), sep = "")

#obj.list <- lapply(localFilePaths, function(x) { return(c(Read10X_h5(x, use.names=FALSE), Read10X_h5(x, use.names=TRUE))) })
obj.list <- lapply(localFilePaths, function(x) { return(Read10X_h5(x, use.names=TRUE)) })

names(obj.list) <- c("Day0_CART_Ctrl","Day0_CART_on_D3","Day5_CART_on_D8","Day0_CART_on_D8","Day5_CART_on_D13","Day0_CART_on_D16","Day5_CART_on_D21","Day0_CART_Ctrl_rep2","Day0_CART_on_D3_rep2","Day5_CART_on_D8_rep2","Day0_CART_on_D8_rep2","Day5_CART_on_D13_rep2","Day0_CART_on_D16_rep2","Day5_CART_on_D21_rep2")
obj.list <- obj.list[sort(names(obj.list))]

#obj.list <- CS028876_Kanakry_Manuscript_Revision$value

mincells = 3
mingenes = 200

mitoch = "^mt-"
cc.genes$g2m.genes= str_to_title(cc.genes$g2m.genes)
cc.genes$s.genes = str_to_title(cc.genes$s.genes)

seurat_object <- function(i) {
  if (class(obj.list[[i]]) == "dgCMatrix"){
    so.nf <- CreateSeuratObject(counts = obj.list[[i]], assay = "RNA", project=names(obj.list)[[i]], min.cells = 1, min.features = 0)
  } else {
    so.nf <- CreateSeuratObject(counts = obj.list[[i]][1]$`Gene Expression`, assay = "RNA", project=names(obj.list)[[i]], min.cells = 1, min.features = 0)
  }
  so.nf <- NormalizeData(so.nf, normalization.method = "LogNormalize", scale.factor = 10000)
  so.nf[["percent.mt"]] <- PercentageFeatureSet(object = so.nf, pattern = mitoch)
  so.nf$log10GenesPerUMI <- log10(so.nf$nFeature_RNA) / log10(so.nf$nCount_RNA)
  
  #Filtered Seurat Object:
  if (class(obj.list[[i]]) == "dgCMatrix"){
    so <- CreateSeuratObject(counts = obj.list[[i]], assay = "RNA", project=names(obj.list)[[i]], min.cells = mincells, min.features = mingenes)
  } else {
    so <- CreateSeuratObject(counts = obj.list[[i]][1]$`Gene Expression`, assay = "RNA", project=names(obj.list)[[i]], min.cells = mincells, min.features = mingenes)
  }
  so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
  so[["percent.mt"]] <- PercentageFeatureSet(object = so, pattern = mitoch)
  # so <- CellCycleScoring(object = so, g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)
  so$log10GenesPerUMI <- log10(so$nFeature_RNA) / log10(so$nCount_RNA)
  
  so.origcount = dim(so.nf)[2]

  #Start with filtering here:
  maxgenes = 2500
  complexity = 0.8
  minUMI = 500
  MAD_gene <- TRUE
  ngenestdev <- mad(so@meta.data$nFeature_RNA)
  ngenemed <- median(so@meta.data$nFeature_RNA)
  ngenemaxlim <- ngenemed+(3*ngenestdev)
  gl = format(round(ngenemaxlim,0),nsmall=0)
  
  maxmitoch = 15
  
  MAD_mitoch <- FALSE
  mitostdev <- mad(so@meta.data$percent.mt)
  mitomed <- median(so@meta.data$percent.mt)
  mitomaxlim <- mitomed+(3*mitostdev)
  ml = format(round(mitomaxlim,2),nsmall=2)
  
  if (MAD_gene == TRUE & MAD_mitoch == TRUE)       {
    cat(paste0("Gene Count Filter = low:",mingenes," high:",gl),"\n")
    cat(paste0("Mitochondrial Percentage Filter =",ml,"\n"))
    cat(paste0("Complexity Filter =",complexity,"\n"))
    so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < ngenemaxlim & so@meta.data$percent.mt <= mitomaxlim & so@meta.data$log10GenesPerUMI > complexity), ]))
    perc.remain = (dim(so)[2]/so.origcount)*100
    perc.remain=formatC(perc.remain,format = "g",digits=3)
    cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
    cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
  }
  else if (MAD_gene == FALSE & MAD_mitoch == TRUE) {
    cat(paste0("Gene Count Filter = low:", mingenes," high:", maxgenes),"\n")
    cat(paste0("Mitochondrial Percentage Filter =",ml,"\n"))
    so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < maxgenes & so@meta.data$percent.mt <= mitomaxlim & so@meta.data$log10GenesPerUMI > complexity & so@meta.data$nCount_RNA > minUMI), ]))
    perc.remain = (dim(so)[2]/so.origcount)*100
    perc.remain=formatC(perc.remain,format = "g",digits=3)
    cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
    cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
  }
  else if (MAD_gene == TRUE & MAD_mitoch == FALSE){
    cat(paste0("Gene Count Filter = low:",mingenes," high:",gl),"\n")
    cat(paste0("Mitochondrial Percentage Filter =", maxmitoch,"\n"))
    so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < ngenemaxlim & so@meta.data$nFeature_RNA > mingenes & so@meta.data$percent.mt < maxmitoch & so@meta.data$log10GenesPerUMI > complexity & so@meta.data$nCount_RNA > minUMI), ]))
    perc.remain = (dim(so)[2]/so.origcount)*100
    perc.remain=formatC(perc.remain,format = "g",digits=3)
    cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
    cat(paste0("Percent Remaining=" ,perc.remain,"%\n\n"))
  }
  else {
    cat(paste0("Gene Count Filter = low:", mingenes," high:", maxgenes),"\n")
    cat(paste0("Mitochondrial Percentage Filter =", maxmitoch,"\n"))
    so <- subset(so, cells = rownames(so@meta.data[which(so@meta.data$nFeature_RNA < maxgenes & so@meta.data$nFeature_RNA > mingenes & so@meta.data$percent.mt < maxmitoch & so@meta.data$log10GenesPerUMI > complexity & so@meta.data$nCount_RNA > minUMI), ]))
    perc.remain = (dim(so)[2]/so.origcount)*100
    perc.remain=formatC(perc.remain,format = "g",digits=3)
    cat(paste0("Filtered Cell Count=" ,dim(so)[2]),"\n")
    cat(paste0("Percent Remaining=" ,perc.remain),"\n\n")
  }
  
  df.m <- melt(so@meta.data)
  df.m$filt <- "filt"
  df.m$filt <- as.factor(df.m$filt)
  df2.m <- melt(so.nf@meta.data)
  df2.m$filt <- "raw"
  df2.m$filt <- as.factor(df2.m$filt)
  
  v <- unique(df.m$variable)
  so2.list <- list(so,so.nf)
  
  return(so2.list)
}

so.list <- lapply(seq_along(obj.list), seurat_object)

SO <- lapply(so.list,function(x) x[[1]])
names(SO) <- sapply(names(obj.list), function(x) gsub("_filtered.h5", "", x))



#initialize Citeseq functionality as false, 
#later the template will check for a Protein assay and run if it finds it
dat = vector()

for(i in 2:length(SO)){dat=c(dat,SO[[i]]) }
SO_merge <- merge(SO[[1]], y = dat, add.cell.ids = names(SO), project = "scRNAProject", merge.data = TRUE)
allgenes <- rownames(SO_merge)
SO_merge <- ScaleData(SO_merge, assay = "RNA", features=allgenes)

npcs = 20
Do_SCTransform = TRUE
vars_to_regress = c("nCount_RNA")

SO_merge <- SCTransform(SO_merge,do.correct.umi = TRUE,vars.to.regress=vars_to_regress,return.only.var.genes = FALSE)
 
all_features <- lapply(SO, row.names) %>% Reduce(intersect, .)

SO_merge <- FindVariableFeatures(object = SO_merge, nfeatures = 3000, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, 100000), selection.method = "vst", verbose = FALSE)
SO_merge <- RunPCA(object = SO_merge, npcs = npcs, verbose = FALSE,seed.use = 42)
SO_merge <- RunUMAP(object = SO_merge, reduction = "pca", dims = 1:npcs, seed.use=42)
SO_merge <- RunTSNE(object = SO_merge, reduction = "pca", dim.embed = 2, dims = 1:npcs, seed.use = 1)
SO_merge <- FindNeighbors(SO_merge, dims = 1:npcs)

for (i in seq(0.2,1.2,0.2)){
  SO_merge <- FindClusters(SO_merge, resolution = i, algorithm = 1)
}

return(SO_merge)
