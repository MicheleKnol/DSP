###############################################################
# 12_Reactome.R - 
# Michèle Knol
# 25-09-2025
###############################################################

library(ReactomePA)
library(clusterProfiler)
library(org.Mm.eg.db)
library(readr)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(tidyr)
library(pheatmap)

# Parameters
base_path <- "C:/Users/knol/OneDrive - Université de Genève/Documents/DSP Paediatrics/Nanostring"

# Define regions and cell types
pairs <- list(
  # c("Glial cells", "VTA"),
  c("Glial cells", "DG")
  # c("Glial cells", "Cortex"),
  # c("Glial cells", "Cerebellum"),
  # c("Glial cells", "CA1"),
  # c("Glial cells", "CA3"),
  # c("Glial cells", "all_r"),
  # c("Glial cells", "Hippo"),
  # c("Neurons", "VTA"),
  # c("Neurons", "DG"),
  # c("Neurons", "Cortex"),
  # c("Neurons", "Cerebellum"),
  # c("Neurons", "CA1"),
  # c("Neurons", "CA3"),
  # c("Neurons", "all_r"),
  # c("Neurons", "Hippo")
)

all_results <- list()

for (pair in pairs) {
  celltype <- pair[1]
  region <- pair[2]
  label <- paste(celltype, region)
  dge_file <- file.path(base_path, label, "DGE.csv")
  output_dir <- file.path(base_path, label)
  
  if (!file.exists(dge_file)) {
    warning(paste("Missing file:", dge_file))
    next
  }
  
  message("Processing: ", label)
  dge_data <- read_csv(dge_file, show_col_types = FALSE)
  dge_data$Estimate <- -dge_data$Estimate  # Flip sign
  
  if (!"Gene" %in% names(dge_data)) {
    if ("SYMBOL" %in% names(dge_data)) {
      names(dge_data)[names(dge_data) == "SYMBOL"] <- "Gene"
    } else {
      warning("Missing 'Gene' or 'SYMBOL' column in: ", label)
      next
    }
  }
  
  gene_map <- bitr(dge_data$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db) %>% na.omit()
  dge_mapped <- inner_join(dge_data, gene_map, by = c("Gene" = "SYMBOL"))
  
  comparisons <- unique(dge_mapped$Contrast)
  for (comp in comparisons) {
    comp_dge <- dge_mapped %>% filter(Contrast == comp)
    sig_genes <- comp_dge %>% filter(abs(Estimate) > 0.5, Pr_t < 0.05)
    
    sig_entrez <- unique(na.omit(sig_genes$ENTREZID))
    universe_entrez <- unique(na.omit(comp_dge$ENTREZID))
    
    if (length(sig_entrez) == 0) {
      message("No significant genes for ", label, " - ", comp)
      next
    }
    
    reactome_result <- tryCatch({
      enrichPathway(gene = sig_entrez,
                    universe = universe_entrez,
                    organism = "mouse",
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    readable = TRUE)
    }, error = function(e) NULL)
    
    if (!is.null(reactome_result) && nrow(reactome_result@result) > 0) {
      comp_id <- gsub(" ", "_", comp)
      out_file <- file.path(output_dir, paste0("Reactome_ORA_", comp_id, ".xlsx"))
      write.xlsx(reactome_result@result, out_file)
      
      # Save individual barplot
      top_res <- head(reactome_result@result[order(reactome_result@result$p.adjust), ], 10)
      top_res$GeneRatio_num <- sapply(top_res$GeneRatio, function(x) {
        parts <- unlist(strsplit(x, "/"))
        as.numeric(parts[1]) / as.numeric(parts[2])
      })
      
      p <- ggplot(top_res, aes(x = reorder(Description, GeneRatio_num),
                               y = GeneRatio_num,
                               fill = p.adjust)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        scale_fill_gradient(low = "red", high = "blue", name = "adj p-value") +
        labs(title = paste("Reactome ORA:", label, comp),
             x = "Pathway", y = "Gene Ratio") +
        theme_minimal(base_size = 12)
      
      ggsave(file.path(output_dir, paste0("Reactome_ORA_", comp_id, "_barplot.png")), p, width = 10, height = 5)
      
      # Store for summary heatmap
      temp_df <- reactome_result@result %>%
        mutate(comparison = paste(label, comp, sep = "_")) %>%
        select(comparison, Description, p.adjust)
      all_results[[paste(label, comp, sep = "_")]] <- temp_df
    } else {
      message("No Reactome terms enriched for ", label, " - ", comp)
    }
    
    # Optional: plot top pathway per contrast as before (comment out if not needed)
    if (!is.null(reactome_result) && nrow(reactome_result@result) > 0) {
      sig_terms <- reactome_result@result %>% filter(p.adjust < 0.05)
      
      if (nrow(sig_terms) > 0) {
        top_desc <- sig_terms$Description[1]
        top_pathway_id <- sig_terms$ID[1]
        
        message("Plotting pathway: ", top_desc)
        
        # Clean pathway name
        clean_path_name <- gsub("[^[:alnum:] /\\-]", "", top_desc)
        
        # Prepare foldChange vector
        fold_change_vec <- setNames(as.numeric(sig_genes$Estimate), as.character(sig_genes$ENTREZID))
        fold_change_vec <- fold_change_vec[!is.na(fold_change_vec)]
        fold_change_vec <- fold_change_vec[!duplicated(names(fold_change_vec))]
        
        tryCatch({
          png(filename = file.path(output_dir, paste0("Reactome_Pathway_", comp_id, "_", 
                                                      sub("[^A-Za-z0-9]", "", top_desc), ".png")),
              width = 1200, height = 800)
          
          viewPathway(
            pathName = clean_path_name,
            organism = "mouse",
            readable = TRUE,
            foldChange = fold_change_vec,
            keyType = "ENTREZID"
          )
          
          dev.off()
        }, error = function(e) {
          message("Error plotting pathway: ", e$message)
          dev.off()
        })
      }
    }
  }
  
  # Combine all results for summary heatmap
  combined <- bind_rows(all_results)
  
  if (nrow(combined) > 0) {
    combined <- combined %>%
      mutate(negLogAdjP = -log10(p.adjust),
             Pathway = Description)
    
    # Select top 20 pathways by minimum p.adjust across all comparisons
    top_pathways <- combined %>%
      group_by(Pathway) %>%
      summarize(min_p = min(p.adjust)) %>%
      arrange(min_p) %>%
      slice(1:20) %>%
      pull(Pathway)
    
    heatmap_data <- combined %>%
      filter(Pathway %in% top_pathways) %>%
      select(comparison, Pathway, negLogAdjP) %>%
      pivot_wider(names_from = comparison, values_from = negLogAdjP, values_fill = 0)
    
    heatmap_mat <- as.matrix(heatmap_data[,-1])
    rownames(heatmap_mat) <- heatmap_data$Pathway
    
    pheatmap(heatmap_mat,
             main = "Reactome ORA summary: top pathways",
             fontsize_row = 8,
             fontsize_col = 10,
             filename = file.path(output_dir, "Reactome_ORA_summary_heatmap.png"),
             width = 10,
             height = 8)
  }
  
  ##Plot user-specified Reactome pathways with fold changes
  
  pathways_of_interest <- c("Antigen processing-Cross presentation")
  
  fold_change_vec_all <- setNames(as.numeric(comp_dge$Estimate), as.character(comp_dge$ENTREZID))
  fold_change_vec_all <- fold_change_vec_all[!is.na(fold_change_vec_all)]
  fold_change_vec_all <- fold_change_vec_all[!duplicated(names(fold_change_vec_all))]
  
  for (pw in pathways_of_interest) {
    message("Plotting user-selected pathway: ", pw)
    
    # Clean pathway name to avoid special characters (for filename)
    pw_clean <- gsub("[^A-Za-z0-9]", "", pw)
    
    # Plot pathway with fold changes
    tryCatch({
      png(filename = file.path(output_dir, paste0("Reactome_CustomPathway_", pw_clean, "_", comp_id, ".png")),
          width = 1200, height = 800)
      
      viewPathway(
        pathName = pw,
        organism = "mouse",
        readable = TRUE,
        foldChange = fold_change_vec_all,
        keyType = "ENTREZID"
      )
      
      dev.off()
    }, error = function(e) {
      message("Error plotting pathway ", pw, ": ", e$message)
      dev.off()
    })
  }
  
}