###############################################################
# 15_UpsetPathways.R - Upset plots of the pathways
# Michèle Knol
# 22-10-2025
###############################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(ComplexUpset)
library(purrr)
library(stringr)

#-----------------------------
# Directories
#-----------------------------
gsea_dir <- "Results/GSEA"
output_dir <- "Results/UpsetPathways"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -----------------------------
# Regions, Cell Types, Contrasts
# -----------------------------
regions <- c("Cortex", "Cerebellum", "CA1", "CA3", "DG")
segments <- c("Neurons", "Astrocytes")
contrasts <- c("FLASH-CONV", "CONV-Control", "FLASH-Control")

# -----------------------------
# Load all pathways
# -----------------------------
all_gsea <- list()

for (r in regions) {
  for (s in segments) {
    for (c in contrasts) {
      
      gsea_file <- file.path(gsea_dir, paste0("Pathways_", r, "_", s, "_", c, ".csv"))
      
      if (!file.exists(gsea_file)) {
        message("Missing file: ", gsea_file)
        next
      }
      
      gsea <- read_csv(gsea_file, show_col_types = FALSE) %>%
        mutate(
          Region = r,
          CellType = s,
          Contrast = c,
          Group = paste(r, s, sep = "_")
        )
      
      all_gsea[[paste(r, s, c, sep = "_")]] <- gsea
    }
  }
}

gsea_df <- bind_rows(all_gsea)

# -----------------------------
# Loop over contrasts and generate UpSet plots
# -----------------------------
contrasts = c("CONV-Control")
for (c in contrasts) {
  
  # Filter for the comparison
  df_c <- gsea_df %>% 
    filter(Contrast == c)
  
  #Filter for significant pathways
  sig_gsea <- df_c %>%
    filter(qvalue < 0.1) %>%
    mutate(Direction = ifelse(NES > 0, "Up", "Down"))
  
  # Build UpSet input
  upset_df <- sig_gsea %>%
    select(Description, Group, Direction) %>%
    distinct() %>%
    mutate(value = 1) %>%
    pivot_wider(
      names_from = Group,
      values_from = value,
      values_fill = 0
    )
  
  # --------------------------
  # UP pathways for this contrast
  # --------------------------
  up_df <- upset_df %>% 
    inner_join(sig_gsea %>% filter(Direction == "Up") %>% select(Description) %>% distinct(),
               by = "Description")
  
  png(file.path(output_dir, paste0("Upregulated_", c, "_UpSet.png")),
      width = 2400, height = 2000, res = 300)
  ComplexUpset::upset(
    up_df,
    intersect = unique(sig_gsea$Group),
    name = paste("Upregulated -", c),
    width_ratio = 0.15,
    n_intersections = 20
  )
  dev.off()
  
  # --------------------------
  # DOWN pathways for this contrast
  # --------------------------
  down_df <- upset_df %>% 
    inner_join(sig_gsea %>% filter(Direction == "Down") %>% select(Description) %>% distinct(),
               by = "Description")
  
  png(file.path(output_dir, paste0("Downregulated_", c, "_UpSet.png")),
      width = 2400, height = 2000, res = 300)
  ComplexUpset::upset(
    down_df,
    intersect = unique(sig_gsea$Group),
    name = paste("Downregulated -", c),
    width_ratio = 0.15,
    n_intersections = 20
  )
  dev.off()
}

# -----------------------------
# Barplot for shared pathways (ALL contrasts)
# -----------------------------
sig_gsea <- gsea_df %>%
  filter(qvalue < 0.25) %>%
  mutate(Direction = ifelse(NES > 0, "Up", "Down"))

upset_df <- sig_gsea %>%
  select(Description, Group, Direction) %>%
  distinct() %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = Group,
    values_from = value,
    values_fill = 0
  )

group_cols <- setdiff(colnames(upset_df), c("Description", "Direction"))

count_df <- upset_df %>%
  mutate(n_groups = rowSums(across(all_of(group_cols)))) %>%
  inner_join(sig_gsea %>% select(Description, Direction) %>% distinct(), by = "Description")

png(file.path(output_dir, "Shared_Pathways_Barplot.png"), width = 2000, height = 1600, res = 300)
ggplot(count_df, aes(x = reorder(Description, n_groups), y = n_groups, fill = Direction)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Shared Enriched Pathways Across Regions / Cell Types / Contrasts",
    x = "Pathway",
    y = "Number of Groups"
  ) +
  theme_bw()
dev.off()


