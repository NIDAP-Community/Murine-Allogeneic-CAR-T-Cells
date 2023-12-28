 
bar_table <- table(so@meta.data$merged_samples, so@meta.data$SCT_snn_res.0.4)

# For aggregated samples
counts_percentage_aggregate_df <- table(so@meta.data$merged_samples, so@meta.data$SCT_snn_res.0.4)
colsums_df <- colSums(bar_table)

col_norm_df <- list()
for (col_name in colnames(counts_percentage_aggregate_df)){
  col_norm_df[[col_name]] <- counts_percentage_aggregate_df[,col_name]/colsums_df[col_name]
}

row_col_norm_df <- as.data.frame(do.call(cbind, col_norm_df))

# row_col_norm_df <- list()
# rowsums_df <- rowSums(col_norm_df)
# for (row_name in rownames(col_norm_df)){
#   row_col_norm_df[[row_name]] <- col_norm_df[row_name,]/rowsums_df[row_name]
# }

#row_col_norm_df <- do.call(rbind,row_col_norm_df)
#row_col_norm_df <- as.data.frame(row_col_norm_df)

row_col_norm_df$CellType <- rownames(row_col_norm_df)

counts_percentage_aggregate_df <- pivot_longer(row_col_norm_df, colnames(row_col_norm_df)[colnames(row_col_norm_df) != "CellType"], names_to = "Sample", values_to = "Proportion")

counts_percentage_aggregate_df <- na.omit(counts_percentage_aggregate_df)
#factor(counts_percentage_aggregate_df$Cluster, levels = sort(unique(as.numeric(counts_percentage_aggregate_df$Cluster))))

# From early timepoint to endpoint
#counts_percentage_aggregate_df$Cluster <- factor(counts_percentage_aggregate_df$Cluster, levels = c("0","8","10","11","12","1","2","4","5","6","7","9","13","14","3"))

scale_colour_discrete <- function(...) {
  scale_colour_brewer(..., palette="cols")
}

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#e6beff","#9A6324","#800000","#aaffc3","#808000","yellow")

g <- ggplot(counts_percentage_aggregate_df, aes(fill=CellType, y=`Proportion`, x=Sample)) + ggtitle("Sample Proportion in each Sample") + theme_classic() +
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values=cbPalette) + theme(axis.text.x = element_text(angle = 90)) + theme(plot.title = element_text(lineheight=.8, face="bold", hjust = 0.5), legend.title = element_blank())

#ggsave(plot = g, filename = paste0(celltype, ".png"))
print(g)

