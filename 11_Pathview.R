###############################################################
# 11_Pathview.R - Show pathways in more detail
# Michèle Knol
# 22-10-2025
###############################################################
library(pathview)
library(readr)
library(dplyr)
library(tidyr)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(here)
library(KEGGREST)
library(purrr)

# -------------------
# Directories
# -------------------
setwd(here())
dge_path    <- "Results/DGE"
gsea_path   <- "Results/GSEA"
output_dir  <- "Results/Pathview"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -------------------
# Inputs
# -------------------
regions   <- c("VTA", "Cortex", "Cerebellum","CA1","CA3","DG")
celltypes <- c("Neurons","Astrocytes")
kegg_pathway <- "mmu04064"

gsea_all <- tibble()
dge_all <- tibble()

# -------------------
# Loop through all regions and cell types
# -------------------
for(region in regions){
  for(celltype in celltypes){
    
    dge_file <- file.path(dge_path, paste0("DGE_", region, "_", celltype, ".csv"))
    
    if(!file.exists(dge_file)){
      message("Skipping missing file: ", dge_file)
      next
    }
    
    dge <- read_csv(dge_file, col_names = TRUE)
    
    #Add Region and CellType info to DGE
    dge <- dge %>%
      mutate(Region = region, CellType = celltype)
    
    # Append to dge_all
    dge_all <- bind_rows(if (exists("dge_all")) dge_all else NULL, dge)
    
    
    #Get NES values for later
    gsea_file_conv  <- file.path(gsea_path, paste0("Pathways_", region, "_", celltype, "_CONV-Control.csv"))
    gsea_file_flash <- file.path(gsea_path, paste0("Pathways_", region, "_", celltype, "_FLASH-Control.csv"))
    
    for (f in c(gsea_file_conv, gsea_file_flash)) {
      if (!file.exists(f)) {
        message("Missing GSEA file: ", f)
        next
      }
      tmp <- read_csv(f, show_col_types = FALSE)
      tmp$Region <- region
      tmp$CellType <- celltype
      
      # Extract contrast name automatically from filename
      tmp$Contrast <- ifelse(grepl("FLASH-Control", f), "FLASH - Control", "CONV - Control")
      
      # Keep relevant columns (Pathway + NES)
      tmp <- tmp %>% dplyr::select(Region, CellType, Contrast, Pathway = ID, Description, NES, p.adjust)
      gsea_all <- bind_rows(gsea_all, tmp)
    }
    
    # Separate contrasts
    dge_conv  <- dge %>% filter(Contrast == "CONV - Control") %>% dplyr::select(Gene, Estimate)
    dge_flash <- dge %>% filter(Contrast == "FLASH - Control") %>% dplyr::select(Gene, Estimate)
    
    # Map gene symbols to Entrez IDs
    gene_conv_entrez  <- mapIds(org.Mm.eg.db, keys = dge_conv$Gene,  column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    gene_flash_entrez <- mapIds(org.Mm.eg.db, keys = dge_flash$Gene, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
    
    # Remove NAs
    dge_conv  <- dge_conv[!is.na(gene_conv_entrez), ]
    dge_flash <- dge_flash[!is.na(gene_flash_entrez), ]
    gene_conv_entrez  <- gene_conv_entrez[!is.na(gene_conv_entrez)]
    gene_flash_entrez <- gene_flash_entrez[!is.na(gene_flash_entrez)]
    
    # Create named vectors for Pathview
    conv_vec  <- setNames(dge_conv$Estimate, gene_conv_entrez)
    flash_vec <- setNames(dge_flash$Estimate, gene_flash_entrez)
    
    # Merge into matrix with two columns
    all_genes <- union(names(conv_vec), names(flash_vec))
    gene_matrix <- cbind(CONV  = conv_vec[all_genes],
                         FLASH = flash_vec[all_genes])
    gene_matrix[is.na(gene_matrix)] <- 0
    
    # -------------------
    # Run Pathview
    # -------------------
    # Save current working directory
    original_wd <- getwd()

    # Change to output directory
    setwd(output_dir)

    pathview(gene.data = gene_matrix,
             pathway.id = kegg_pathway,
             species = "mmu",
             kegg.native = TRUE,
             out.suffix = paste0(region, "_", celltype),
             same.layer = TRUE,
             gene.idtype = "entrez",
             low = "blue", high = "red", mid = "white",
             limit = c(-2,2),
             out.dir = normalizePath(output_dir))

    # Restore original working directory
    setwd(original_wd)
  }
}

# -------------------
# Plot Barplot for All Regions/CellTypes
# -------------------
pathway_info <- keggGet(kegg_pathway)[[1]]$NAME
pathway_info <- sub(" - .*", "", pathway_info)

ggplot(
  gsea_all %>% filter(Pathway == kegg_pathway),
  aes(x = Region, y = NES, fill = Contrast)
) +
  geom_col(position = "dodge") +
  facet_wrap(~CellType, scales = "free_y") +
  scale_fill_manual(values = c("CONV - Control" = "#9acd32", "FLASH - Control" = "#ff69b4")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 10) +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(
    title = paste0("KEGG ", kegg_pathway, " (", pathway_info, ")"),
    y = "Normalized Enrichment Score (NES)",
    x = "Region"
  )
ggsave(file.path(output_dir, paste0("Pathview_Barplot_", kegg_pathway, ".png")),
       width = 10, height = 6, dpi = 300)
message("All Pathview plots saved to ", output_dir)



# -------------------
# Plot key genes for the selected pathway
# -------------------
region = "Cortex"
celltype = "Astrocytes"

# Retrieve pathway gene list
kegg_data <- keggGet(kegg_pathway)[[1]]
pathway_genes <- kegg_data$GENE

# Extract the gene symbols (even positions)
pathway_symbols <- pathway_genes[seq(2, length(pathway_genes), 2)] %>% unique()
pathway_genes <- kegg_data$GENE
pathway_symbols <- pathway_genes[seq(2, length(pathway_genes), 2)] %>%
  gsub(";.*", "", .) %>%
  trimws() %>%
  unique()


# Summarize to get mean Estimate per gene/contrast/region/celltype
dge_plot <- dge_all %>%
  filter(Gene %in% pathway_symbols) %>%
  filter(Contrast %in% c("CONV - Control", "FLASH - Control")) %>%
  group_by(Gene, Contrast, Region, CellType) %>%
  summarise(Estimate = mean(Estimate, na.rm = TRUE), .groups = "drop")

# Function to plot each Region/CellType
plot_region_celltype <- function(region, celltype) {
  subset_df <- dge_plot %>%
    filter(Region == region, CellType == celltype)
  
  # Pivot wider to get estimates side by side
  wide_df <- subset_df %>%
    pivot_wider(names_from = Contrast, values_from = Estimate)
  
  # Calculate absolute difference and select top 10 genes
  top_genes <- wide_df %>%
    mutate(Diff = abs(`CONV - Control` - `FLASH - Control`)) %>%
    arrange(desc(Diff)) %>%
    slice_head(n = 10) %>%
    pull(Gene)
  
  # Filter original subset for top genes
  plot_df <- subset_df %>%
    filter(Gene %in% top_genes) %>%
    mutate(
      Direction = ifelse(Contrast == "CONV - Control", -1, 1),
      EstimateAbs = abs(Estimate),
      EstimatePlot = Direction * EstimateAbs,
      Regulation = Estimate > 0
    )
  
  # Plot
  p <- ggplot(plot_df, aes(x = EstimatePlot, y = Gene, fill = Regulation)) +
    geom_col(width = 0.7, alpha = 1, color = "black") +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
    annotate("text", x = -2, y = Inf, label = "CONV",
             vjust = 2, hjust = 0, color = "#9acd32", fontface = "bold", size = 4) +
    annotate("text", x =  2, y = Inf, label = "FLASH",
             vjust = 2, hjust = 1, color = "#ff69b4", fontface = "bold", size = 4) +
    scale_fill_manual(
      values = c("TRUE" = "#d81b60", "FALSE" = "#1565c0"),
      labels = c("FALSE" = "Downregulated", "TRUE" = "Upregulated"),
      name = "Regulation"
    ) +
    scale_x_continuous(
      limits = c(-2, 2),
      labels = abs
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom") +
    labs(
      title = paste0("Top pathway genes – ", pathway_info),
      subtitle = paste(region, "-", celltype),
      x = "Absolute mean logFC",
      y = "Gene"
    )
  
  # Save each plot
  out_name <- paste0("Pathway_", kegg_pathway, "_", region, "_", celltype, "_Top10.png")
  ggsave(file.path(output_dir, out_name), plot = p, width = 7, height = 5, dpi = 300)
  
  return(p)
}


# Generate plots for all Region/CellType pairs
region_celltype_pairs <- dge_plot %>%
  distinct(Region, CellType)

plots <- purrr::pmap(region_celltype_pairs, function(Region, CellType) {
  plot_region_celltype(Region, CellType)
})