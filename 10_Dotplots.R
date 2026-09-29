###############################################################
# 09_GSEA_Plots.R - Dotplots + UpSet plots for GSEA results
# Michèle Knol
# 26-09-2025
###############################################################

library(tidyverse)
library(plotly)
library(ComplexUpset)
library(htmlwidgets)
library(dplyr)
library(tidyverse)
library(plotly)
library(htmlwidgets)

# ----------------------------
# Define directories
# ----------------------------
gsea_dir   <- "Results/GSEA"
output_dir <- "Results/Dotplots"
output_dir_upset <- "Results/UpsetPathways"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir_upset, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Parameters
# ----------------------------
regions <- c("CA1", "CA3", "DG")
# regions <- c("VTA", "Cortex", "Cerebellum", "CA1", "CA3", "DG")
# segments <- c("Neurons", "Glial cells")
segments <- c("Neurons", "Astrocytes")
comparisons <- c("FLASH-Control", "CONV-Control", "FLASH-CONV")
# ranking_methods <- c("After FLASH", "After CONV")
ranking_methods <- c("After FLASH", "After CONV", "Biggest Difference", "After RT", "After FLASH, not CONV")
directions <- c("up", "down")

# -------------------
# Define colors
# -------------------
region_colors <- c("VTA" = "#B0F2B4", "Cortex" = "#F2E2BA", "Cerebellum" = "#F2BAC9", 
                   "CA1" = "#BAF2E9", "CA3" = "#BAD7F2", "DG" = "#D6C9DE")


# ----------------------------
# Function clean the data
# ----------------------------
clean <- function(df,region, segment) { 
  df$n_genes <- sapply(strsplit(as.character(df$core_enrichment), "/"), length)
  df$gene_ratio <- df$n_genes / df$setSize
  df <- df %>%
    # filter(ONTOLOGY == "KEGG") %>%
    mutate(
      Description = gsub(" - Mus musculus \\(house mouse\\)", "", Description),
      Pathway = paste(ONTOLOGY, Description),
      Pathway = str_wrap(Pathway, 50),
      Region = region,
      Segment = segment)
}

# ----------------------------
# Function to make the dotplots
# ----------------------------
dotplots <- function(region, segment, gsea_dir, output_dir, rank_by) {
  
  qval_cutoff <- 0.05
  top_n <- 10
  
  # Open GSEAs file
  flashfile <- file.path(gsea_dir, paste0("Pathways_simplified_DGE_", region, "_", segment, "_FLASH_Control.csv"))
  convfile <- file.path(gsea_dir, paste0("Pathways_simplified_DGE_", region, "_", segment, "_CONV_Control.csv"))
  
  if (file.exists(flashfile)) {
    flashfile <- read_csv(flashfile, show_col_types = FALSE)
    convfile <- read_csv(convfile, show_col_types = FALSE)
    message("Start creating plots for: ", region, "|", segment, ", ranked by ", rank_by)
  } else {
    message("Skipped (not found): ", flashfile)
    return(NULL)
  }
  
  # Clean dataframes and combine
  flashfile <- clean(flashfile, region, segment)
  convfile <- clean(convfile, region, segment)
  flashfile$Condition <- "FLASH"
  convfile$Condition <- "CONV"
  df <- bind_rows(flashfile, convfile) %>%
    dplyr::select(Condition, Pathway, NES, pvalue, qvalue, gene_ratio)
  
  # Prepare rankings (use NES = 0 only for ranking purposes, not for plotting)
  df_wide <- df %>%
    dplyr::select(Pathway, Condition, NES) %>%
    pivot_wider(names_from = Condition, values_from = NES)
  
  df_wide <- df_wide %>%
    mutate(
      FLASH_rank = ifelse(is.na(FLASH), 0, FLASH),
      CONV_rank  = ifelse(is.na(CONV), 0, CONV),
      difference = abs(FLASH_rank - CONV_rank),
      overall    = abs(FLASH_rank) + abs(CONV_rank)
    )
  
  if (rank_by == "After FLASH") {
    top_pathways <- df_wide %>%
      arrange(desc(FLASH_rank)) %>%
      slice_head(n = top_n) %>%
      pull(Pathway)
    
  } else if (rank_by == "After CONV") {
    top_pathways <- df_wide %>%
      arrange(desc(CONV_rank)) %>%
      slice_head(n = top_n) %>%
      pull(Pathway)
    
  } else if (rank_by == "Biggest Difference") {
    top_pathways <- df_wide %>%
      arrange(desc(difference)) %>%
      slice_head(n = top_n) %>%
      pull(Pathway)
    
  } else if (rank_by == "After RT") {
    top_pathways <- df_wide %>%
      arrange(desc(overall)) %>%
      slice_head(n = top_n) %>%
      pull(Pathway)
    
  } else if (rank_by == "After FLASH, not CONV") {
    top_pathways <- df_wide %>%
      arrange(desc(FLASH_rank)) %>%
      filter(CONV_rank <= 0 | is.na(CONV)) %>%
      slice_head(n = top_n) %>%
      pull(Pathway)
    
  } else {
    stop("Invalid ranking method specified.")
  }
  
  # Force a complete grid: every top pathway x every condition gets a row.
  # Missing combinations stay NA -> plotted as grey dots, distinct from a real NES = 0
  full_grid <- expand_grid(
    Pathway   = top_pathways,
    Condition = c("FLASH", "CONV")
  )
  
  df <- full_grid %>%
    left_join(df, by = c("Pathway", "Condition")) %>%
    mutate(
      Missing   = is.na(NES),
      plot_size = ifelse(is.na(qvalue), 0, -log10(qvalue))  # smallest size if missing
    ) %>%
    mutate(Pathway = factor(Pathway, levels = rev(top_pathways)))
  
  dotplot <- ggplot(df, aes(x = Condition, y = Pathway, size = plot_size, color = NES)) +
    geom_point() +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                          limits = c(-2.5, 2.5), name = "NES", na.value = "grey80") +
    scale_size_continuous(name = "-log10(q)", range = c(2, 12), limits = c(0, 8), breaks = c(2, 4, 6, 8)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) +
    ggtitle(paste("Top Pathways (", rank_by, ")", "\n",  region, "|", segment))
  
  # Save PNG
  ggsave(file.path(output_dir, paste0("Dotplot_", rank_by, "_", region, "_", segment, ".png")),
         dotplot, width = 6, height = 6, dpi = 300)
  
  pdf(file.path(output_dir,paste0("Dotplot_", rank_by, "_", region, "_", segment, ".pdf")),
      width = 6,height = 6,useDingbats = FALSE)
  print(dotplot)
  dev.off()
  
  message("Plots saved for ", region, "/", segment, "/", rank_by)
}

# ----------------------------
# Function to make upset plots
# ----------------------------
make_upset <- function(gsea_dir, segment, contrast, direction, outdir, qval_cutoff) {
  all_list <- list()
  
  for (r in regions) {
    fname <- file.path(gsea_dir, paste0("Pathways_", r, "_", segment, "_", contrast, ".csv"))
    
    if (!file.exists(fname)) {
      message("Skipped (not found): ", fname)
      next
      }
    else {
       df <- read.csv(fname, check.names = FALSE, stringsAsFactors = FALSE)
      }
    
    df_clean <- clean(df, r, segment)
    
    keep_cols <- c("Pathway", "NES", "qvalue", "pvalue", "core_enrichment", "gene_ratio")
    keep_cols <- keep_cols[keep_cols %in% colnames(df_clean)]
    
    df_clean <- df_clean %>%
      dplyr::select(all_of(keep_cols)) %>%
      mutate(Region = r, Condition = contrast)
    
    all_list[[length(all_list) + 1]] <- df_clean
  }
  
  combined <- bind_rows(all_list)
  
  # significance filter
  if (direction == "up") {
    sub <- combined %>% filter(NES > 0, qvalue < qval_cutoff)
  } else {
    sub <- combined %>% filter(NES < 0, qvalue < qval_cutoff)
  }

  # ensure unique pathway-region pairs
  sub <- sub %>% distinct(Pathway, Region, .keep_all = TRUE)
  
  # build wide binary matrix: Pathway x Region
  mat <- sub %>%
    mutate(flag = 1L) %>%
    dplyr::select(Pathway, Region, flag) %>%
    tidyr::pivot_wider(names_from = Region, values_from = flag, values_fill = 0)
  
  # make sure all regions are included
  for (r in regions) {
    if (!(r %in% colnames(mat))) {
      mat[[r]] <- 0L
    }
  }
  
  upset_df <- as.data.frame(mat)
  set_cols <- regions
  
  # construct upset plot
  upsetplot <- ComplexUpset::upset(
    upset_df,
    intersect = set_cols,
    n_intersections = 20,
    keep_empty_groups = TRUE,
    queries = if (!is.null(region_colors)) {
      lapply(regions, function(r) {
        ComplexUpset::upset_query(
          set = r,
          fill = region_colors[r]
        )
      })
    } else NULL,
    base_annotations = list(
      'Intersection size' = (
        intersection_size(counts = TRUE) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = 'black')
          )
      )
    ),
    matrix = intersection_matrix(
      geom = geom_point(
        shape = 'circle filled',
        size = 3.5,
        stroke = 0.45
      )
    ),
    set_sizes = (
      upset_set_size(geom = geom_bar(width = 0.4)) +
        {if (!is.null(region_colors)) scale_fill_manual(values = region_colors) else NULL} +
        theme(
          axis.line.x = element_line(colour = 'black'),
          axis.ticks.x = element_line()
        )
    ),
    stripes = upset_stripes(
      geom = geom_segment(size = 12),
      colors = c('grey95', 'white')
    ),
    sort_sets = "ascending",
    sort_intersections = "descending"
  ) +
    ggtitle(paste0(ifelse(direction == "up", "Higher Expressed", "Lower Expressed"),
                   " Pathways — ", contrast, " | ", segment)) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )
  
  # Save
  outfile <- file.path(outdir, paste0("UpSet_", direction, "_", contrast, "_", segment, ".png"))
  ggsave(outfile, upsetplot, width = 12, height = 8, dpi = 300)
  
  # Save Excel with binary presence/absence matrix
  out_xlsx <- file.path(outdir, paste0("UpSet_", direction, "_", contrast, "_", segment, ".xlsx"))
  openxlsx::write.xlsx(upset_df, out_xlsx, rowNames = FALSE)
  
  message("Saved Upset plot for ", segment, " | ", contrast)
}

# ----------------------------
# Loop over all regions, celltypes, and comparisons
# ----------------------------
for (r in regions) {
  for (s in segments) {
    for (rank in ranking_methods) {
      dotplots(r, s, gsea_dir, output_dir, rank)
    }
  }
}

# for (s in segments) {
#   for (c in comparisons) {
#     for (d in directions) {
#       make_upset(gsea_dir, s, c, d, output_dir_upset, 0.1)
#     }
#   }
# }

# To check a single plot
# make_upset(gsea_dir, "Neurons", "FLASH-Control", "up", output_dir_upset, 0.1)
# dotplots("VTA", "Neurons", gsea_dir, output_dir, "Biggest Difference")



# ----------------------------
# Search for neuroinflammation / microglia / astrocyte terms
# ----------------------------
search_terms_gsea <- function(region, segment, gsea_dir,
                              keywords = c("inflamm", "microglia", "microglial",
                                           "astrocyte", "astroglia", "cytokine",
                                           "complement", "interleukin", "TNF",
                                           "interferon", "glial", "gliosis",
                                           "phagocyt", "innate immune", "immune response", "lactate"),
                              qval_cutoff = 0.05) {
  
  flashfile <- file.path(gsea_dir, paste0("Pathways_simplified_DGE_", region, "_", segment, "_FLASH_Control.csv"))
  convfile  <- file.path(gsea_dir, paste0("Pathways_simplified_DGE_", region, "_", segment, "_CONV_Control.csv"))
  
  if (!file.exists(flashfile) || !file.exists(convfile)) {
    message("Skipped (not found): ", region, " | ", segment)
    return(NULL)
  }
  
  flash <- read_csv(flashfile, show_col_types = FALSE) %>% mutate(Condition = "FLASH")
  conv  <- read_csv(convfile,  show_col_types = FALSE) %>% mutate(Condition = "CONV")
  
  df <- bind_rows(flash, conv)
  
  # Build a single regex pattern from all keywords (case-insensitive)
  pattern <- paste(keywords, collapse = "|")
  
  hits <- df %>%
    filter(grepl(pattern, Description, ignore.case = TRUE)) %>%
    mutate(
      Region   = region,
      Segment  = segment,
      Significant = qvalue < qval_cutoff
    ) %>%
    dplyr::select(Region, Segment, Condition, ONTOLOGY, Description, NES, pvalue, qvalue, Significant) %>%
    arrange(Description, Condition)
  
  if (nrow(hits) == 0) {
    message("No matching terms found for ", region, " | ", segment)
    return(NULL)
  }
  
  return(hits)
}

# ----------------------------
# Run across all region/segment combinations of interest
# ----------------------------
regions   <- c("DG", "CA3", "CA1")   # adjust to match your actual region names
segments  <- c("Neurons")           # adjust to match your actual segment names

all_hits <- list()
for (r in regions) {
  for (s in segments) {
    res <- search_terms_gsea(r, s, gsea_dir = "results/GSEA")
    if (!is.null(res)) {
      all_hits[[paste(r, s)]] <- res
    }
  }
}

all_hits_df <- bind_rows(all_hits)

# View full table
print(all_hits_df, n = Inf)

# Optionally save
write_csv(all_hits_df, "results/GSEA/Neuroinflammation_terms_summary.csv")

# Quick summary: which terms are significant, and in which condition, sign of NES
all_hits_df %>%
  filter(Significant) %>%
  arrange(Region, Description, Condition) %>%
  dplyr::select(Region, Segment, Condition, Description, NES, qvalue)