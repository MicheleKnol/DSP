###############################################################
# 19_DETable.R
# Creates a master table of all significant DE genes
# Author: Michèle Knol
# 11.02.2026
###############################################################

library(dplyr)
library(readr)

# ----------------------------
# User settings
# ----------------------------
dge_dir <- "Results/DGE"
output_dir <- "Results"

regions <- c("CA1", "CA3", "DG")
celltypes <- c("Neurons", "Astrocytes")
comparisons <- c("CONV - Control", "FLASH - Control")

# ----------------------------
# Combine all significant DE genes
# ----------------------------
all_de_genes <- list()

for(r in regions){
  for(c in celltypes){
    dge_file <- file.path(dge_dir, paste0("DGE_", r, "_", c, ".csv"))
    
    if(!file.exists(dge_file)){
      message("File not found: ", dge_file)
      next
    }
    
    dge <- read_csv(dge_file, show_col_types = FALSE)
    
    for(comp in comparisons){
      
      sig <- dge %>%
        filter(Contrast == comp) %>%
        filter(Pr_t < 0.05 & abs(Estimate) > 1) %>%
        mutate(
          Region = r,
          CellType = c,
          Direction = ifelse(Estimate > 0, "Up", "Down")
        ) %>%
        select(Region, CellType, Contrast, Gene, Estimate, Pr_t, Direction)
      
      if(nrow(sig) > 0){
        all_de_genes[[length(all_de_genes) + 1]] <- sig
      }
    }
  }
}

# Combine all into one master table
combined_de <- bind_rows(all_de_genes)

# Order for clarity: Region → CellType → Contrast → p-value
combined_de <- combined_de %>%
  arrange(Region, CellType, Contrast, Pr_t)

# Save master table
write_csv(combined_de, file.path(output_dir, "DE_genes_all_regions_celltypes.csv"))
message("Master DE gene table saved.")

# ----------------------------
# Summary table: counts per region/celltype/contrast
# ----------------------------
summary_table <- combined_de %>%
  group_by(Region, CellType, Contrast) %>%
  summarise(
    N_Significant = n(),
    N_Up = sum(Direction == "Up"),
    N_Down = sum(Direction == "Down"),
    .groups = "drop"
  )

# Save summary table
write_csv(summary_table, file.path(output_dir, "DE_genes_summary.csv"))
message("Summary table saved.")