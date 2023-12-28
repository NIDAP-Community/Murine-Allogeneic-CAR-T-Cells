# Preprocessing h5 Data ~ required for downstream steps
so <- source("Preprocessing.R")

# garbage collection
rm(list = setdiff(ls(), "so"))
gc()

so <- so$value

# Run ProjecTILs for cell annotation and merge replicates
so$Barcode <- colnames(so)

so_meta <- so@meta.data

so_meta <- so_meta %>% mutate(merged_samples = case_when(
  grepl("Day0_CART_Ctrl",orig.ident) ~ "Day0_CART_Ctrl",
  grepl("Day0_CART_on_D16",orig.ident) ~ "Day0_CART_on_D16",
  grepl("Day0_CART_on_D3",orig.ident) ~ "Day0_CART_on_D3",
  grepl("Day0_CART_on_D8",orig.ident) ~ "Day0_CART_on_D8",
  grepl("Day5_CART_on_D13",orig.ident) ~ "Day5_CART_on_D13",
  grepl("Day5_CART_on_D21",orig.ident) ~ "Day5_CART_on_D21",
  grepl("Day5_CART_on_D8",orig.ident) ~ "Day5_CART_on_D8",
))

so@meta.data$merged_samples <- so_meta$merged_samples[match(rownames(so@meta.data),so_meta$Barcode)]

meta <- read.csv("files_for_analysis/ProjectTIL_metadata.csv")
so@meta.data$ProjectTILS_calls <- meta$ProjectTILS_calls[match(meta$X, so@meta.data$Barcode)]

# Filter out NA and Day0 Ctrl
Idents(so) <- so$ProjectTILS_calls
keep_celltype <- unique(Idents(so))[!unique(Idents(so)) %in% "NA"]
so <- subset(so, idents = keep_celltype)
Idents(so) <- so$merged_samples
keep_timepoint <- unique(Idents(so))[!unique(Idents(so)) %in% "Day0_CART_Ctrl"]
so <- subset(so, idents = keep_timepoint)

# order by timepoint
so$merged_samples <- factor(so$merged_samples, levels = c(
  "Day0_CART_on_D3",
  "Day5_CART_on_D8",
  "Day0_CART_on_D8",
  "Day5_CART_on_D13",
  "Day0_CART_on_D16",
  "Day5_CART_on_D21"
))

# Figures 6A and B
source("plotMeta.R")

# tSNE shape and Number of clusters (FindClusters) may vary across package versions 
plotMetadata(
    object = so,
    metadata.to.plot = c("merged_samples","SCT_snn_res_0_4"),
    columns.to.summarize = "c()",
    summarization.cut.off = 5,
    reduction.type = "tsne",
    use.cite.seq = FALSE,
    show.labels = FALSE,
    legend.text.size = 1,
    legend.position = "right",
    dot.size = 0.01
)

# Figures 6C and D
source("Composition_Barplot_Sample.R")
source("Composition_Barplot_Celltype.R")

# Figure 7 DEG and Volcanoes / Figure S14 Pathway Plots 
comp_celltype <- read.csv("files_for_analysis/Compiled_Celltype_meta.csv")
so@meta.data$Likely_Celltype <- comp_celltype$Likely_CellType[match(so$Barcode, comp_celltype$Barcode)]

Idents(so) <- so$Likely_Celltype
cd4_treg <- subset(so, idents = c("CD4","Tregs"))
cd8 <- subset(so, idents = c("CD8"))

source("MAST.R")

cd8_deg <- MAST(object = cd8)
cd4_treg_deg <- MAST(object = cd4_treg)

source("Volcano.R")
source("l2p_Pathway.R")

# Cd8+ cells
volcanoPlot(deg_tab = cd8_deg,
            significance_col = "p_val_Day0_CART_on_D3_vs_Day5_CART_on_D8",
            log_fc_col = "avg_logFC_Day0_CART_on_D3_vs_Day5_CART_on_D8")

l2p_Pathway(cd8_deg,
  significance_col = "p_val_Day0_CART_on_D3_vs_Day5_CART_on_D8",
  log_fc_col = "avg_logFC_Day0_CART_on_D3_vs_Day5_CART_on_D8",
  upregulated = T)

l2p_Pathway(cd8_deg,
            significance_col = "p_val_Day0_CART_on_D3_vs_Day5_CART_on_D8",
            log_fc_col = "avg_logFC_Day0_CART_on_D3_vs_Day5_CART_on_D8",
            upregulated = F)

volcanoPlot(deg_tab = cd8_deg,
            significance_col = "p_val_Day0_CART_on_D8_vs_Day5_CART_on_D13",
            log_fc_col = "avg_logFC_Day0_CART_on_D8_vs_Day5_CART_on_D13")

l2p_Pathway(cd8_deg,
            significance_col = "p_val_Day0_CART_on_D8_vs_Day5_CART_on_D13",
            log_fc_col = "avg_logFC_Day0_CART_on_D8_vs_Day5_CART_on_D13",
            upregulated = T)

l2p_Pathway(cd8_deg,
            significance_col = "p_val_Day0_CART_on_D8_vs_Day5_CART_on_D13",
            log_fc_col = "avg_logFC_Day0_CART_on_D8_vs_Day5_CART_on_D13",
            upregulated = F)

volcanoPlot(deg_tab = cd8_deg,
            significance_col = "p_val_Day0_CART_on_D16_vs_Day5_CART_on_D21",
            log_fc_col = "avg_logFC_Day0_CART_on_D16_vs_Day5_CART_on_D21")

l2p_Pathway(cd8_deg,
            significance_col = "p_val_Day0_CART_on_D16_vs_Day5_CART_on_D21",
            log_fc_col = "avg_logFC_Day0_CART_on_D16_vs_Day5_CART_on_D21",
            upregulated = T)

l2p_Pathway(cd8_deg,
            significance_col = "p_val_Day0_CART_on_D16_vs_Day5_CART_on_D21",
            log_fc_col = "avg_logFC_Day0_CART_on_D16_vs_Day5_CART_on_D21",
            upregulated = F)
# Cd4+ cells
volcanoPlot(deg_tab = cd4_treg_deg,
            significance_col = "p_val_Day0_CART_on_D3_vs_Day5_CART_on_D8",
            log_fc_col = "avg_logFC_Day0_CART_on_D3_vs_Day5_CART_on_D8")

l2p_Pathway(cd4_treg_deg,
            significance_col = "p_val_Day0_CART_on_D3_vs_Day5_CART_on_D8",
            log_fc_col = "avg_logFC_Day0_CART_on_D3_vs_Day5_CART_on_D8",
            upregulated = T)

l2p_Pathway(cd4_treg_deg,
            significance_col = "p_val_Day0_CART_on_D3_vs_Day5_CART_on_D8",
            log_fc_col = "avg_logFC_Day0_CART_on_D3_vs_Day5_CART_on_D8",
            upregulated = F)

volcanoPlot(deg_tab = cd4_treg_deg,
            significance_col = "p_val_Day0_CART_on_D8_vs_Day5_CART_on_D13",
            log_fc_col = "avg_logFC_Day0_CART_on_D8_vs_Day5_CART_on_D13")

l2p_Pathway(cd4_treg_deg,
            significance_col = "p_val_Day0_CART_on_D8_vs_Day5_CART_on_D13",
            log_fc_col = "avg_logFC_Day0_CART_on_D8_vs_Day5_CART_on_D13",
            upregulated = T)

l2p_Pathway(cd4_treg_deg,
            significance_col = "p_val_Day0_CART_on_D8_vs_Day5_CART_on_D13",
            log_fc_col = "avg_logFC_Day0_CART_on_D8_vs_Day5_CART_on_D13",
            upregulated = F)

volcanoPlot(deg_tab = cd4_treg_deg,
            significance_col = "p_val_Day0_CART_on_D16_vs_Day5_CART_on_D21",
            log_fc_col = "avg_logFC_Day0_CART_on_D16_vs_Day5_CART_on_D21")

l2p_Pathway(cd4_treg_deg,
            significance_col = "p_val_Day0_CART_on_D16_vs_Day5_CART_on_D21",
            log_fc_col = "avg_logFC_Day0_CART_on_D16_vs_Day5_CART_on_D21",
            upregulated = T)

l2p_Pathway(cd4_treg_deg,
            significance_col = "p_val_Day0_CART_on_D16_vs_Day5_CART_on_D21",
            log_fc_col = "avg_logFC_Day0_CART_on_D16_vs_Day5_CART_on_D21",
            upregulated = F)

# Figure S12A (Dimension Plot by Gene Expression) and S12B (Violin)
source("Color_by_Gene.R")

colorByGene(object = so,
            gene = c("Pdcd1","Lag3","Ctla4","Havcr2","Tigit","Il2ra"),
            reduction.type = "tsne")

source("Violin_Plots_by_Metadata.R")

violinPlot(object = so, 
           assay = "SCT", 
           slot = "scale.data", 
           group.by = "merged_samples", 
           group.subset = c(), 
           genes.of.interest = c("Pdcd1","Lag3","Ctla4","Havcr2","Tigit","Il2ra"),
           filter.outliers = T)

# Figure S15
violinPlot(object = cd8, 
           assay = "SCT", 
           slot = "scale.data", 
           group.by = "merged_samples", 
           group.subset = c(), 
           genes.of.interest = c("Gzma"),
           filter.outliers = T)

cd8_deg_clust4v3 <- MAST(object = cd8,
                contrasts_to_run = c("4-3"))

volcanoPlot(deg_tab = cd8_deg_clust4v3,
            significance_col = "p_val_4_vs_3",
            log_fc_col = "avg_logFC_4_vs_3")


