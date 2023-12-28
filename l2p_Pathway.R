l2p_Pathway <- function(deg_tab,
                        significance_col,
                        log_fc_col,
                        upregulated = T){

  library(l2p)
  library(magrittr)
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
  library(RCurl)
  library(RColorBrewer)
  
  fold_change_column <- significance_col  ##Add Fold Change Column
  p_value_column <- log_fc_col ###Add pvalue column
  
  deg_tab %>% SparkR::select("Gene", fold_change_column , p_value_column) %>% SparkR::collect() -> genesmat
  genes_universe = as.vector(unique(unlist(genesmat["Gene"])))
  
  if(upregulated){
    genesmat %>% filter(.data[[fold_change_column]] < -0.2 & .data[[p_value_column]] < 0.05) %>% pull(Gene) -> genes_to_include
  } else {
    genesmat %>% filter(.data[[fold_change_column]] > 0.2 & .data[[p_value_column]] < 0.05) %>% pull(Gene) -> genes_to_include
    }
  print(length(genes_to_include))
  
  genes_to_include <- as.vector(unique(unlist(genes_to_include)))
  gene_set_sources_to_include = c("GO","REACTOME","KEGG")
  categories_string <- paste(gene_set_sources_to_include, collapse=",")
  
  organism = "Mouse"
  if (organism != "Human") {
    organism_of_interest = "Mouse"
    orthology_table = SparkR::filter(Ortholog_Map_for_RNA_Seq,     
                                     Ortholog_Map_for_RNA_Seq$Organism==organism_of_interest)
    
    if ("from human"=="from human") {
      orthology_reference_column = "Human_gene_name"
      orthology_conversion_column = "Nonhuman_gene_name"
    } else {
      orthology_reference_column = "Nonhuman_gene_name"
      orthology_conversion_column = "Human_gene_name"
    }
    orthology_table %>% SparkR::withColumnRenamed(orthology_reference_column,"orthology_reference") %>%
      SparkR::withColumnRenamed(orthology_conversion_column, "orthology_conversion") %>% SparkR::select("orthology_reference", "orthology_conversion") -> orthology_table
    orthology_table <- SparkR::collect(orthology_table)
    orthology_table %>% dplyr::filter(orthology_conversion %in% genes_to_include) %>% dplyr::select(orthology_reference) -> genes_to_include
    genes_to_include <- as.vector(genes_to_include$orthology_reference)
    genes_to_include <- unique(genes_to_include)
    orthology_table %>% dplyr::filter(orthology_conversion %in% genes_universe) %>% dplyr::select(orthology_reference) -> genes_universe
    genes_universe <- as.vector(genes_universe$orthology_reference)
    genes_universe <- unique(genes_universe)
  } else {
    genes_to_include = genes_to_include
    genes_universe = genes_universe
  }
  
  use_built_in_gene_universe = FALSE
  if (use_built_in_gene_universe) {
    x <- l2pwcats(genes_to_include, categories_string)
    print("Using built-in gene universe.")
  } else {
    x <- l2puwcats(genes_to_include, genes_universe, categories_string)
    print("Using all genes in differential expression analysis as gene universe.")
  }
  x %>%
    dplyr::arrange(pval) %>% dplyr::mutate(hitsPerc=(pwhitcount/(pwnohitcount+pwhitcount))*100) %>% dplyr::mutate(pathtotal=pwnohitcount+pwhitcount) %>% dplyr::filter(ratio >= 0) %>%
    dplyr::select("pathwayname", "category", "pathwayaccessionidentifier", "pval", "fdr", "pwhitcount", "genesinpathway", "pwnohitcount","pathtotal","hitsPerc", "inputcount", "pwuniverseminuslist","ratio")  %>% dplyr::rename(diff_ratio = ratio) %>% dplyr::filter(pval < 0.05)%>% dplyr::filter(pwhitcount >= 5) -> x
  
  print(paste0("Total number of pathways: ", nrow(x)))
  goResults <- x
  goResults %>% top_n(10, wt=-log(pval)) %>%
    dplyr::arrange(-log(pval)) -> goResults
  minp = min(goResults$pval) - 0.1*min(goResults$pval)
  maxp = max(goResults$pval) + 0.1*max(goResults$pval)
  sizemax = ceiling(max(goResults$pwhitcount)/10)*10  
  goResults %>% dplyr::mutate(pathwayname2 = stringr::str_replace_all(pathwayname, "_", " ")) -> goResults
  goResults %>% dplyr::mutate(pathwayname2 = stringr::str_wrap(pathwayname2,30)) -> goResults
  
  if (FALSE){
    goResults %>% dplyr::mutate(percorder = order(goResults$pval)) -> goResults
    goResults$pathwayname2 <- factor(goResults$pathwayname2, levels = goResults$pathwayname2[goResults$percorder])
    xmin = floor(min(goResults$pval))
    xmax = max(goResults$pval) 
    gplot <- goResults %>% 
      ggplot(aes(x=pval,
                 y=pathwayname2, 
                 colour=hitsPerc, 
                 size=pwhitcount)) +
      geom_point() +
      theme(text = element_text(size=20), legend.position = "right", legend.key.height = unit(1, "cm"),
            axis.title.x = element_text(size = rel(1.2)),axis.title.y = element_text(size = rel(1.2))) +
      xlim(xmin,xmax) +
      expand_limits(colour = seq(minp, maxp, by = 10),
                    size = seq(0, sizemax,by=10)) +
      labs(x="p value", y="GO term", colour="Hits (%)", size="Count") 
    print(gplot)
  } else {
    goResults %>% dplyr::mutate(percorder = order(goResults$hitsPerc)) -> goResults
    goResults$pathwayname2 <- factor(goResults$pathwayname2, levels = goResults$pathwayname2[goResults$percorder])
    xmin = floor(min(goResults$hitsPerc)-5)
    xmax = ceiling(max(goResults$hitsPerc)+5) 
    gplot <- goResults %>% 
      ggplot(aes(x=hitsPerc,
                 y=pathwayname2, 
                 colour=pval, 
                 size=pwhitcount)) +
      geom_point() +
      theme_classic() +
      theme(text = element_text(size=30), legend.position = "right", legend.key.height = unit(1, "cm"),
            axis.title.x = element_text(size = rel(1.2)),axis.title.y = element_text(size = rel(1.2))) +
      xlim(xmin,xmax) +
      expand_limits(colour = seq(minp, maxp, by = 10),
                    size = seq(0, sizemax,by=10)) +
      labs(x="Hits (%)", y="Pathway", colour="p value", size="Count") +
      guides(colour = guide_colourbar(order = 1), size = guide_legend(order=2))
    print(gplot)
  }
  return(x)
}
