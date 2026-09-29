###############################################################
# 08_GSEA.R - Gene Set Enrichment Analysis
# Michèle Knol
# 25-09-2025
###############################################################
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(readr)
library(GOSemSim)

# ----------------------------
# User settings
# ----------------------------
dge_dir    <- "results/DGE"
output_dir <- "results/GSEA"

# ----------------------------
# Define contrasts
# ----------------------------
comparisons <- c("CONV - Control", "FLASH - Control", "FLASH - CONV")

# ----------------------------
# Simplify function
# ----------------------------
# ----------------------------
# Precompute semantic similarity data (done once, reused across all calls)
# ----------------------------
semdata_list <- list(
  BP = godata(OrgDb = org.Mm.eg.db, keytype = "SYMBOL", ont = "BP"),
  MF = godata(OrgDb = org.Mm.eg.db, keytype = "SYMBOL", ont = "MF"),
  CC = godata(OrgDb = org.Mm.eg.db, keytype = "SYMBOL", ont = "CC")
)

# ----------------------------
# Simplify function
# ----------------------------
simplify_gsea <- function(gsea_go_obj, semdata, sig_cutoff = 0.1,
                          cutoff = 0.7, by = "p.adjust", select_fun = min) {
  if (is.null(gsea_go_obj) || nrow(gsea_go_obj@result) == 0) {
    message("  No GO results to simplify.")
    return(NULL)
  }
  
  # Filter to significant terms only before simplifying
  n_total <- nrow(gsea_go_obj@result)
  gsea_go_obj@result <- gsea_go_obj@result %>% filter(p.adjust < sig_cutoff)
  n_sig <- nrow(gsea_go_obj@result)
  
  if (n_sig == 0) {
    message("  No significant GO terms (p.adjust < ", sig_cutoff, ") to simplify.")
    return(NULL)
  }
  
  message("  Simplifying ", n_sig, " / ", n_total, " significant GO terms...")
  
  simplified <- tryCatch({
    simplify(gsea_go_obj,
             cutoff     = cutoff,
             by         = by,
             select_fun = select_fun,
             semData    = semdata)
  }, error = function(e) {
    message("  simplify() failed: ", e$message)
    NULL
  })
  
  if (is.null(simplified) || nrow(simplified@result) == 0) return(NULL)
  return(simplified@result)
}

# ----------------------------
# GSEA function
# ----------------------------
gsea <- function(dge_file, output_dir) {
  label <- tools::file_path_sans_ext(basename(dge_file))
  message("Processing file: ", label)
  
  # Load DGE
  dge <- read_csv(dge_file, show_col_types = FALSE)
  
  # Map to Entrez IDs
  dge$ENTREZID <- mapIds(org.Mm.eg.db,
                         keys      = dge$Gene,
                         column    = "ENTREZID",
                         keytype   = "SYMBOL",
                         multiVals = "first")
  
  for (comp in comparisons) {
    message("Contrast: ", comp)
    comp_id  <- gsub(" - ", "_", comp)
    comp_dge <- dge %>% filter(Contrast == comp)
    
    # Ranking metric
    comp_dge <- comp_dge %>%
      mutate(
        FCsign  = sign(Estimate),
        logPval = -log10(Pr_t),
        metric  = FCsign * logPval
      )
    
    ranked_genes_symbol <- comp_dge$metric
    names(ranked_genes_symbol) <- comp_dge$Gene
    ranked_genes_symbol <- sort(ranked_genes_symbol, decreasing = TRUE)
    
    ranked_genes_entrez <- comp_dge$metric
    names(ranked_genes_entrez) <- comp_dge$ENTREZID
    ranked_genes_entrez <- ranked_genes_entrez[!is.na(names(ranked_genes_entrez))]
    ranked_genes_entrez <- sort(ranked_genes_entrez, decreasing = TRUE)
    
    # Run GSEA
    gsea_results          <- list()
    gsea_results_simplified <- list()
    
    # --- GO ---
    # --- GO (run separately per ontology so simplify() works correctly) ---
    go_results_full       <- list()
    go_results_simplified <- list()
    
    for (ont in c("BP", "MF", "CC")) {
      message("  Running gseGO for ", ont, "...")
      
      gsea_go_ont <- tryCatch({
        gseGO(geneList     = ranked_genes_symbol,
              OrgDb        = org.Mm.eg.db,
              keyType      = "SYMBOL",
              ont          = ont,
              pvalueCutoff = 1,
              verbose      = FALSE)
      }, error = function(e) NULL)
      
      if (is.null(gsea_go_ont) || nrow(gsea_go_ont@result) == 0) next
      
      # Collect full results
      ont_res          <- gsea_go_ont@result
      ont_res$ONTOLOGY <- ont
      go_results_full[[ont]] <- ont_res
      
      # Simplify
      ont_simplified <- simplify_gsea(gsea_go_ont, semdata = semdata_list[[ont]])
      if (!is.null(ont_simplified)) {
        ont_simplified$ONTOLOGY        <- ont
        go_results_simplified[[ont]]   <- ont_simplified
      }
    }
    
    if (length(go_results_full) > 0) {
      go_res          <- bind_rows(go_results_full)
      go_res$ONTOLOGY <- "GO"
      gsea_results[["GO"]] <- go_res
    }
    
    if (length(go_results_simplified) > 0) {
      go_simplified          <- bind_rows(go_results_simplified)
      go_simplified$ONTOLOGY <- "GO"
      gsea_results_simplified[["GO"]] <- go_simplified
    }
    
    # # --- KEGG ---
    # gsea_kegg <- tryCatch({
    #   gseKEGG(geneList     = ranked_genes_entrez,
    #           organism     = "mmu",
    #           keyType      = "ncbi-geneid",
    #           pvalueCutoff = 1,
    #           verbose      = FALSE)
    # }, error = function(e) NULL)
    # 
    # if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
    #   kegg_res          <- gsea_kegg@result
    #   kegg_res$ONTOLOGY <- "KEGG"
    #   gsea_results[["KEGG"]]            <- kegg_res
    # }
    # 
    # # --- WikiPathways ---
    # gsea_wp <- tryCatch({
    #   gseWP(geneList     = ranked_genes_entrez,
    #         organism     = "Mus musculus",
    #         pvalueCutoff = 1,
    #         verbose      = FALSE)
    # }, error = function(e) NULL)
    # 
    # if (!is.null(gsea_wp) && nrow(gsea_wp@result) > 0) {
    #   wp_res          <- gsea_wp@result
    #   wp_res$ONTOLOGY <- "WP"
    #   gsea_results[["WP"]]            <- wp_res
    # }
    
    # --- Save full results ---
    if (length(gsea_results) > 0) {
      all_res  <- bind_rows(gsea_results)
      out_file <- file.path(output_dir, paste0("Pathways_", label, "_", comp_id, ".csv"))
      write_csv(all_res, out_file)
      message("Saved: ", out_file)
    } else {
      message("No GSEA results for ", comp)
    }
    
    # --- Save simplified results ---
    if (length(gsea_results_simplified) > 0) {
      all_res_simplified  <- bind_rows(gsea_results_simplified)
      out_file_simplified <- file.path(output_dir, paste0("Pathways_simplified_", label, "_", comp_id, ".csv"))
      write_csv(all_res_simplified, out_file_simplified)
      message("Saved (simplified): ", out_file_simplified)
    } else {
      message("No simplified GSEA results for ", comp)
    }
  }
}

# ----------------------------
# Run over all DGE files
# ----------------------------
dge_files <- list.files(dge_dir, pattern = "\\.csv$", full.names = TRUE)
for (f in dge_files) {
  try(gsea(f, output_dir))
}