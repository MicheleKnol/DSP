###############################################################
# 18_ORA.R
# Over-Representation Analysis on only the significantly differentially expressed genes
# Author: Michèle Knol
# 11.02.2026
###############################################################

library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(readr)

# ----------------------------
# User settings
# ----------------------------
dge_dir    <- "results/DGE"
output_dir <- "results/ORA"

# ----------------------------
# Define contrasts
# ----------------------------
comparisons <- c("CONV - Control", "FLASH - Control", "FLASH - CONV")

# ----------------------------
# ORA function
# ----------------------------
ora <- function(dge, region, segment, output_dir) {
  
  ontologies <- c("BP", "MF", "CC")
  
  if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
  }
  
  for(comp in comparisons){
    
    dge_sub <- dge %>%
      filter(Contrast == comp)
    
    sig_genes <- dge_sub %>%
      filter(Pr_t < 0.05) %>%
      filter(abs(Estimate) > 1) %>%
      pull(Gene) %>%
      unique()
    
    gene_df <- bitr(sig_genes,
                    fromType = "SYMBOL",
                    toType   = "ENTREZID",
                    OrgDb    = org.Mm.eg.db)
    
    background_genes <- dge_sub$Gene %>% unique()
    
    bg_df <- bitr(background_genes,
                  fromType = "SYMBOL",
                  toType   = "ENTREZID",
                  OrgDb    = org.Mm.eg.db)
    
    combined_results <- list()
    
    # GO enrichment
    for(ont in ontologies){
      
      ego <- enrichGO(
        gene          = gene_df$ENTREZID,
        universe      = bg_df$ENTREZID,
        OrgDb         = org.Mm.eg.db,
        keyType       = "ENTREZID",
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.5,
        readable      = TRUE
      )
      
      if(!is.null(ego) && nrow(as.data.frame(ego)) > 0){
        df <- as.data.frame(ego)
        df$Ontology <- paste0("GO_", ont)
        combined_results[[length(combined_results) + 1]] <- df
      }
    }
    
    # KEGG enrichment
    ekegg <- enrichKEGG(
      gene          = gene_df$ENTREZID,
      universe      = bg_df$ENTREZID,
      organism      = "mmu",
      pvalueCutoff  = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 0.5
    )
    
    if(!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0){
      ekegg <- setReadable(ekegg,
                           OrgDb = org.Mm.eg.db,
                           keyType = "ENTREZID")
      df <- as.data.frame(ekegg)
      df$Ontology <- "KEGG"
      combined_results[[length(combined_results) + 1]] <- df
    }
    
    if(length(combined_results) == 0){
      message("No enrichment for ",
              region, " ", segment, " | ", comp)
      next
    }
    
    # Combine all enrichment results
    combined_df <- bind_rows(combined_results)
    
    # Convert GeneRatio to numeric
    combined_df$GeneRatioNumeric <- sapply(
      combined_df$GeneRatio,
      function(x) eval(parse(text = x))
    )
    
    # Select top 20 pathways across ALL ontologies
    combined_df <- combined_df %>%
      arrange(p.adjust) %>%
      slice_head(n = 20)
    
    comp_clean <- gsub(" ", "_", comp)
    comp_clean <- gsub("-", "vs", comp_clean)
    
    # Single unified dotplot
    p <- ggplot(combined_df,
                aes(x = GeneRatioNumeric,
                    y = reorder(Description, GeneRatioNumeric),
                    size = Count,
                    color = p.adjust)) +
      geom_point() +
      scale_size_continuous(
        name = "Gene Count",
        limits = c(min(combined_df$Count),
                   max(combined_df$Count)),
        range = c(2, 10)
      ) +
      labs(title = paste("ORA for ", region, " | ", segment, " | ", comp),
           x = "Gene Ratio",
           y = NULL) +
      theme_minimal()
    
    pdf(file.path(output_dir,
                  paste0("ORA_combined_dotplot_",
                         region, "_",
                         segment, "_",
                         comp_clean, ".pdf")),
        width = 6, height = 8, useDingbats = FALSE)
    
    print(p)
    dev.off()
    
    message("Finished combined ORA for ",
            region, " ", segment, " | ", comp)
  }
}



# ----------------------------
# Run over all DGE files
# ----------------------------
regions  <- c("CA1", "CA3", "DG")
segments <- c("Neurons", "Astrocytes")

for(r in regions){
  for(s in segments){
    
    dge_file <- file.path(dge_dir, paste0("DGE_", r, "_", s, ".csv"))
    
    if(!file.exists(dge_file)){
      message("File not found: ", dge_file)
      next
    }
    
    dge <- read_csv(dge_file, show_col_types = FALSE)
    
    ora(dge, region = r, segment = s, output_dir = output_dir)
  }
}
