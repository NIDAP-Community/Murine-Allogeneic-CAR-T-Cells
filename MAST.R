MAST <- function(object,
                 contrasts_to_run = c("Day0_CART_on_D3-Day5_CART_on_D8","Day0_CART_on_D8-Day5_CART_on_D13","Day0_CART_on_D16-Day5_CART_on_D21")){
  
  library(Seurat)
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
  library(tidyverse)
  library(ggrepel)
  library(gdata)
  library(reshape2)
  library(tools)
  library(grid)
  library(gridBase)
  library(gridExtra)
  library(parallel)
  library(MAST)
  
  SO = object
  
  #collect parameters here:
  contrasts <- contrasts_to_run
  samples = unique(SO@meta.data$orig.ident)
  
  if (length(samples) == 0) {
    samples = unique(SO@meta.data$sample_name)
  }
  
  colnames(SO@meta.data) <- gsub("orig_ident","orig.ident",colnames(SO@meta.data))
  if("active.ident" %in% slotNames(SO)){
    sample_name = as.factor(SO@meta.data$orig.ident)
    names(sample_name)=names(SO@active.ident)
    SO@active.ident <- as.factor(vector())
    SO@active.ident <- sample_name
    SO.sub = subset(SO, ident = samples)
  } else {
    sample_name = as.factor(SO@meta.data$orig.ident)
    names(sample_name)=names(SO@active.ident)
    SO@active.ident <- as.factor(vector())
    SO@active.ident <- sample_name
    SO.sub = subset(SO, ident = samples)
  }
  
  print("selected samples:")
  print(SO.sub)
  
  meta.df <- SO.sub@meta.data
  colnames(SO.sub@meta.data) = gsub("\\.","_",colnames(SO.sub@meta.data))
  
  #define contrasts
  newcont <- list()
  for (i in 1:length(contrasts)){
    newcont[[i]] <- c(paste(unlist(strsplit(contrasts[i],"-"))))
  }
  contrasts <- newcont
  
  #ERROR CATCHING
  #collect valid names of valid columns
  validColumns <- character()
  for (i in colnames(meta.df)) {
    if (!any(is.na(meta.df[[i]]))) {
      validColumns <-c(validColumns,i)
    }
  }
  
  param2test <- "merged_samples"
  
  if (param2test =="") {
    mcols = colnames(SO.sub@meta.data)
    param2test <-mcols[grepl("RNA_snn",mcols)][[1]]
    print(paste("No parameter selected, defaulting to",param2test))
  }
  
  contrastTarget <- SO.sub@meta.data[[param2test]]
  contrastType <- param2test
  contrastCounts = as.data.frame(table(contrastTarget))
  validContrasts = subset(contrastCounts, Freq>2)[[1]]
  
  #catch malformed contrasts
  for (i in contrasts) {
    if (!(i[[1]] %in% contrastTarget)) {
      print(paste(i[[1]],"is not a valid contrast for contrast type:", contrastType))
      print("Please see below for an example of valid contrasts for your selected contrast type.")
      print(validContrasts)
      stop("You have entered an invalid group to contrast against.")
    } else if (!(i[[2]] %in% contrastTarget) & (i[[2]] != "all")) {
      print(paste(i[[2]],"is not a valid contrast for contrast type:", contrastType))
      print("Please see below for an example of valid contrasts for your selected contrast type.")
      print(validContrasts)
      stop("You have entered an invalid group to contrast against.")
    } else if (length(i)>2) {
      print("Contrasts are as follows..")
      print(i)
      stop("The console says there are too many inputs in your contrasts. A contrast should only contain Group1-Group2, but the console thinks you have inputed Group1-Group2-Group3")
    } else if (!(i[[2]] %in% validContrasts) & (i[[2]] != "all")) {
      print(paste(i[[2]],"has two few values (less than 3 cells) to contrast against. Please see below for contrasts with enough cells:", validContrasts))
      stop("You have entered an invalid group to contrast against.")
    } else if (!(i[[1]] %in% validContrasts)) {
      print(paste(i[[1]],"has two few values (less than 3 cells) to contrast against. Please see below for contrasts with enough cells:", validContrasts))
      stop("You have entered an invalid group to contrast against.")
    }
  }
  
  #print out contrast cell contrastCounts
  for (i in seq_along(contrasts)) {
    firstGroup <- contrasts[[i]][[1]]
    firstGroupCount <- subset(contrastCounts, contrastTarget == firstGroup)$Freq
    if  (contrasts[[i]][[2]]!= "all") {
      secondGroup <-contrasts[[i]][[2]]
      secondGroupCount <-subset(contrastCounts, contrastTarget == secondGroup)$Freq      
      print(paste("Contrast No.",i,"contrasts cluster",firstGroup,"with",firstGroupCount,"cells vs. cluster",secondGroup,"with",secondGroupCount,"cells."))
    } else {
      secondGroupCount <-ncol(SO.sub)-firstGroupCount
      print(paste("Contrast No.",i,"contrasts cluster",firstGroup,"with",firstGroupCount,"cells vs. all other clusters, totalling",secondGroupCount,"cells."))
    } 
  }
  
  #define and call function for running DEG
  get_deg_table <- function(n) {
    library(Seurat)
    
    firstCluster <-n[[1]]
    secondCluster <- n[[2]]
    
    if (n[[2]]=='all') {
      secondCluster <- NULL
    }
    
    Idents(SO.sub) <- param2test
    markers = FindMarkers(SO.sub, ident.1 = firstCluster, ident.2 = secondCluster, test.use = "MAST", logfc.threshold = 0, verbose=FALSE, assay = "SCT")
    colnames(markers) <- chartr(old=" ",new="_",paste(colnames(markers), n[[1]],"vs",n[[2]],sep = "_"))
    return(markers)
  }
  
  deg_tables <- lapply(contrasts, get_deg_table) 

  for(i in seq_along(deg_tables)){
    degtab <- deg_tables[[i]]
    degtab %>% dplyr::filter(.[[1]] < 0.05) %>% dplyr::filter(.[[2]] > 0) %>% dim() -> pos 
    degtab %>% dplyr::filter(.[[1]] < 0.05) %>% dplyr::filter(.[[2]] < 0) %>% dim() -> neg
    print(paste0("The number of upregulated genes at p<0.05 in contrast number ", i, " is:"))
    print(pos[1])
    print(paste0("The number of downregulated genes at p<0.05 in contrast number ", i, " is:"))
    print(neg[1]) 
  }
  
  #Merge the deg tables together
  out_df <- NULL
  for (i in deg_tables) {
    if (is.null(out_df)) {
      out_df <- deg_tables[1]
      out_df <- as.data.frame(out_df)
    } else {
      out_df <- merge(out_df, i, by="row.names", all=TRUE)
      rownames(out_df) <- out_df$Row.names #set the rownames
      out_df$Row.names <- NULL #drop the row.names columns which we no longer need
    }
  }
  
  out_df$Gene <- rownames(out_df)
  out_df$Row.names <- NULL
  out_df <- out_df %>% dplyr::select(Gene, everything())
  return(out_df)
  
}
