###############################################################
# 09_GSEAplots.R - Visualization of the GSEA
# Michèle Knol
# 25-09-2025
###############################################################

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# ----------------------------
# Define directories
# ----------------------------
gsea_dir   <- "results/GSEA"
output_dir   <- "results/GSEA_plots"

# ----------------------------
# Define thresholds
# ----------------------------
qval_cutoff <- 0.1
top_n       <- 15

# ----------------------------
# Load all GSEA results
# ----------------------------
gsea_files <- list.files(gsea_dir, pattern = "^Pathways_simplified.*\\.csv$", full.names = TRUE)

all_gsea <- lapply(gsea_files, function(f) {
  df <- read_csv(f, show_col_types = FALSE)
  df$file <- basename(f)
  df
}) %>% bind_rows()

# ----------------------------
# Extract metadata from filenames
# ----------------------------
all_gsea <- all_gsea %>%
  mutate(file_noext = str_remove(file, "\\.csv$")) %>%
  mutate(parts = str_split(file_noext, "_")) %>%
  rowwise() %>%
  mutate(
    Region   = parts[[2]],  # 3rd element after "Pathways"
    Celltype = parts[[3]],  # 4th element
    Contrast = parts[[4]]   # 5th element
  ) %>%
  ungroup() %>%
  select(-parts, -file_noext)

# ----------------------------
# Filter significant results
# ----------------------------
sig_gsea <- all_gsea %>%
  filter(qvalue < qval_cutoff)

# ----------------------------
# Function to plot one subset
# ----------------------------
plot_top_pathways <- function(df, region, celltype, contrast, output, top_n = 10,
                              qval_limits = c(0, 0.1)) {
  df_sub <- df %>%
    filter(Region == region,
           Celltype == celltype,
           Contrast == contrast) %>%
    arrange(qvalue) %>%
    slice_head(n = top_n) %>%
    mutate(Label = paste0("[", ONTOLOGY, "] ", Description))
  
  if (nrow(df_sub) == 0) return(NULL)
  
  p <- ggplot(df_sub,
              aes(x = NES,
                  y = reorder(Label, NES),
                  fill = qvalue)) +
    geom_col() +
    scale_fill_gradient(
      low = "red",
      high = "blue",
      limits = qval_limits,
      oob = scales::squish,
      name = "q-value"
    ) +
    labs(
      title = paste("Top GSEA pathways\n", region, "|", celltype, "-", contrast),
      x = "Normalized Enrichment Score (NES)",
      y = "Pathway"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  
  ggsave(file.path(output_dir,
                   paste0("TopPathways_", region, "_", celltype, "_", contrast, ".png")),
         p, width = 8, height = 6)
  message("Saved plot: ", region, "|", celltype, "|", contrast)
}


# ----------------------------
# Loop over all regions, celltypes, and comparisons
# ----------------------------
for (r in unique(sig_gsea$Region)) {
  for (ct in unique(sig_gsea$Celltype)) {
    for (c in unique(sig_gsea$Contrast)) {
      plot_top_pathways(sig_gsea, r, ct, c, plot_dir)
    }
  }
}