###############################################################
# 06_VenndiagramsGenes.R - Venndiagram plots of higher/lower expressed genes
# Michèle Knol
# 18-09-2025
###############################################################

library(ggVennDiagram)
library(ggplotify)
library(openxlsx)
library(gt)
library(dplyr)
library(tidyr)
library(readr)
library(grid)
library(ggplot2)
library(UpSetR)
library(ComplexUpset)
library(patchwork)
library(purrr)
library(msigdbr)

# -------------------
# Define directories
# -------------------
dge_dir <- "Results/DGE"
output_dir <- "Results/VennGenes"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Created directory: ", output_dir)
}
# -------------------
# Define colors
# -------------------
region_colors <- c("VTA" = "#B0F2B4", "Cortex" = "#F2E2BA", "Cerebellum" = "#F2BAC9", 
                   "CA1" = "#BAF2E9", "CA3" = "#BAD7F2", "DG" = "#D6C9DE")

# -------------------
# Define cutoff values
# -------------------
p_cutoff <- 0.05
fc_cutoff <- 1

# -------------------
# Load preprocessed object
# -------------------
Data <- readRDS("preprocessed_data.rds")

# -------------------
# Load MSigDB gene to pathway mapping
# -------------------
genesets <- msigdbr(species = "Mus musculus") %>%
  select(gs_name, gene_symbol)

# -------------------
# Loop over Region and Segment to create a Venn diagram
# -------------------
# regions <- unique(Data$Region)
regions <- c("CA1", "CA3", "DG", "Cerebellum")
segments <- unique(Data$Segment)

for(r in regions){
  for(s in segments){
    message(paste("Creating the Venn diagram + Excel tables for:", r, "|", s))
    
    dge_file <- file.path(dge_dir, paste0("DGE_", r, "_", s, ".csv"))
    if (!file.exists(dge_file)) {
      message("  Missing file, skipping: ", basename(dge_file))
      next
    }
    dge <- read_csv(dge_file, show_col_types = FALSE)
    
    # Initialize Excel workbook
    results <- createWorkbook()
    
    # Filter the data
    flash_data <- dge %>% filter(Contrast == "FLASH - Control", Pr_t < p_cutoff)
    conv_data  <- dge %>% filter(Contrast == "CONV - Control",  Pr_t < p_cutoff)
    
    # Significant genes
    flash_high <- flash_data %>% filter(Estimate >  fc_cutoff) %>% pull(Gene) %>% unique()
    flash_low  <- flash_data %>% filter(Estimate < -fc_cutoff) %>% pull(Gene) %>% unique()
    conv_high  <- conv_data  %>% filter(Estimate >  fc_cutoff) %>% pull(Gene) %>% unique()
    conv_low   <- conv_data  %>% filter(Estimate < -fc_cutoff) %>% pull(Gene) %>% unique()
    
    # Create Venn diagrams
    venn_up <- ggVennDiagram(
      list(FLASH = flash_high, CONV = conv_high),
      set_size = 4,
      label = "count",
      label_alpha = 0,
      label_size = 10
    ) +
      # scale_fill_gradient(low = "#fce4ec", high = "#d81b60", limits = c(0, 950)) +
      scale_fill_gradient(low = "#fce4ec", high = "#d81b60") +
      ggtitle(paste0("Higher Expressed Genes in ", r, " | ", s)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      coord_fixed(ratio = 1)
    ggsave(
      filename = file.path(output_dir, paste0(r, "_", s, "_VennHigh.png")),
      plot = venn_up, width = 8, height = 6, dpi = 300
    )
    
    venn_down <- ggVennDiagram(
      list(FLASH = flash_low, CONV = conv_low),
      set_size = 4,
      label = "count",
      label_alpha = 0,
      label_size = 10
    ) +
      # scale_fill_gradient(low = "#e3f2fd", high = "#1565c0", limits = c(0, 950)) +
      scale_fill_gradient(low = "#e3f2fd", high = "#1565c0") +
      ggtitle(paste0("Lower Expressed Genes in ", r, " | ", s)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      coord_fixed(ratio = 1)
    ggsave(
      filename = file.path(output_dir, paste0(r, "_", s, "_VennLow.png")),
      plot = venn_down, width = 8, height = 6, dpi = 300
    )
    
    # Add pathway annotation
    annotate_pathways <- function(df) {
      df %>%
        left_join(genesets, by = c("Gene" = "gene_symbol")) %>%
        group_by(Gene, Estimate, Pr_t, adjp) %>%
        summarise(
          Pathways = paste(unique(gs_name), collapse = "; "),
          .groups = "drop"
        )
    }
    
    # Write Excel sheets
    addWorksheet(results, paste0(r, "_", s, "_FLASH_H"))
    writeData(results, paste0(r, "_", s, "_FLASH_H"),
              annotate_pathways(flash_data %>% filter(Estimate > fc_cutoff)))
    addWorksheet(results, paste0(r, "_", s, "_CONV_H"))
    writeData(results, paste0(r, "_", s, "_CONV_H"),
              annotate_pathways(conv_data %>% filter(Estimate > fc_cutoff)))
    addWorksheet(results, paste0(r, "_", s, "_FLASH_L"))
    writeData(results, paste0(r, "_", s, "_FLASH_L"),
              annotate_pathways(flash_data %>% filter(Estimate < -fc_cutoff)))
    addWorksheet(results, paste0(r, "_", s, "_CONV_L"))
    writeData(results, paste0(r, "_", s, "_CONV_L"),
              annotate_pathways(conv_data %>% filter(Estimate < -fc_cutoff)))
    
    # Save Excel workbook
    saveWorkbook(results,
                 file.path(output_dir, paste0(r, "_", s, "_Venn_GeneLists.xlsx")),
                 overwrite = TRUE)
  }
}

#########Faceted venndiagrams
all_venns <- list()

for(r in regions){
  for(s in segments){

    dge_file <- file.path(dge_dir, paste0("DGE_", r, "_", s, ".csv"))
    if (!file.exists(dge_file)) next
    dge <- read_csv(dge_file, show_col_types = FALSE)

    flash_data <- dge %>% filter(Contrast == "FLASH - Control" & Pr_t < p_cutoff)
    conv_data  <- dge %>% filter(Contrast == "CONV - Control" & Pr_t < p_cutoff)

    venn_list <- list(
      FLASH_Up   = flash_data %>% filter(Estimate > fc_cutoff) %>% pull(Gene) %>% unique(),
      FLASH_Down = flash_data %>% filter(Estimate < -fc_cutoff) %>% pull(Gene) %>% unique(),
      CONV_Up    = conv_data  %>% filter(Estimate > fc_cutoff) %>% pull(Gene) %>% unique(),
      CONV_Down  = conv_data  %>% filter(Estimate < -fc_cutoff) %>% pull(Gene) %>% unique()
    )

    if(all(sapply(venn_list, length) == 0)) next

    # Directly create the ggplot object
    p <- ggVennDiagram(venn_list, label_alpha = 0, edge_size = 0.5, label= "count") +
      scale_fill_gradient(low = "white", high = "coral") +
      ggtitle(paste(r, s, sep = " | ")) +
      theme(plot.title = element_text(hjust = 0.5, size = 10),
            legend.position = "none")

    all_venns[[paste(r,s,sep="_")]] <- p
  }
}

# Combine all plots in a grid: rows = Segment, cols = Region
venn_grid <- wrap_plots(all_venns, ncol = length(regions)) +
  plot_annotation(title = "Venn diagrams across regions and segments")

ggsave(file.path(output_dir, "AllRegions_AllSegments_4setVenn.png"),
       venn_grid, width = 16, height = 10, dpi = 300)


############# Less upsetplots ###################
contrasts <- c("FLASH - CONV", "FLASH - Control", "CONV - Control")

for (contrast in contrasts) {
  message("Processing contrast: ", contrast)
  
  # collect rows (one row per significant gene occurrence)
  rows <- list()
  
  for (s in segments) {
    for (r in regions) {
      dge_file <- file.path(dge_dir, paste0("DGE_", r, "_", s, ".csv"))
      if (!file.exists(dge_file)) {
        message("  Missing file, skipping: ", basename(dge_file))
        next
      }
      
      dge <- read_csv(dge_file, show_col_types = FALSE)
      dge <- dge %>% filter(Contrast == !!contrast)
      if (nrow(dge) == 0) next
      
      # extract only significant genes
      dge_sig <- dge %>%
        mutate(
          Direction = case_when(
            Estimate >  0.5 & Pr_t < 0.05  ~ "Up",
            Estimate < -0.5 & Pr_t < 0.05  ~ "Down",
            TRUE                            ~ NA_character_
          )
        ) %>%
        filter(!is.na(Direction)) %>%
        select(Gene, Estimate, Pr_t, Direction)
      
      if (nrow(dge_sig) == 0) next
      
      # append region & segment meta for each row
      dge_sig <- dge_sig %>%
        mutate(Region = r, Segment = s)
      
      rows[[length(rows) + 1]] <- dge_sig
    }
  }
  
  # bind all rows
  if (length(rows) == 0) {
    message("No significant genes found for contrast: ", contrast)
    next
  }
  df_combined <- bind_rows(rows)
  
  # Make sure Gene rows are unique per Region x Segment
  df_combined <- df_combined %>% distinct(Gene, Region, Segment, Direction, .keep_all = TRUE)
  
  # Build wide `upset_df` with boolean columns for each region
  upset_df <- df_combined %>%
    mutate(Value = TRUE) %>%
    pivot_wider(names_from = Region, values_from = Value, values_fill = FALSE)
  
  # Keep Direction & Segment columns for later
  upset_df$Direction <- factor(upset_df$Direction, levels = c("Up", "Down"))
  upset_df$Segment   <- factor(upset_df$Segment, levels = segments)
  head(upset_df)
  
  # Custom barplot: for each intersection, show counts by Segment and Direction
  # Create an "intersection signature" string for every gene (which regions it belongs to)
  region_cols <- regions
  missing_cols <- setdiff(region_cols, colnames(upset_df))
  if (length(missing_cols) > 0) {
    for (mc in missing_cols) upset_df[[mc]] <- FALSE
  }
  
  # compute intersection label: comma-separated included regions
  upset_df <- upset_df %>%
    rowwise() %>%
    mutate(
      intersection_label = paste(sort(region_cols[which(c_across(all_of(region_cols)))]), collapse = ";")
    ) %>%
    ungroup() %>%
    mutate(
      intersection_label = ifelse(intersection_label == "", "None", intersection_label)
    )
  
  # aggregate counts per intersection x Segment x Direction
  agg <- upset_df %>%
    group_by(intersection_label, Segment, Direction) %>%
    summarise(n = n(), .groups = "drop")
  
  agg$intersection_label <- factor(agg$intersection_label, levels = c("Cerebellum", "DG", "CA1", "CA3"))
  
  # Plot: facet per Segment, bars stacked by Direction (or dodge if you prefer)
  p_bar <- ggplot(agg, aes(x = intersection_label, y = n, fill = Direction)) +
    geom_col(position = position_dodge()) +                      # stacked bars: position_stack() side-by-side:position_dodge() 
    facet_wrap(~Segment, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c("Up" = "#d81b60", "Down" = "#1565c0")) +
    labs(x = "Regions", y = "Number of significant genes",
         title = paste("Significantly Differentially Expressed Genes |", contrast)) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          panel.grid.minor = element_blank())
  
  ggsave(file.path(output_dir, paste0("IntersectionBySegment_", gsub(" ", "_", contrast), ".png")),
         plot = p_bar, width = 16, height = 5, dpi = 300)
  message(" Saved intersection composition plot: ", paste0("IntersectionBySegment_", gsub(" ", "_", contrast), ".png"))
}



#######Upsetplot############
regions <- c("Cerebellum", "DG", "CA1", "CA3")

make_upset_plot <- function(upset_df, regions, contrast, celltype) {
  upset(
    upset_df,
    regions,
    n_intersections = 40,
    sort_intersections = "descending",
    name = paste("Lower Expressed Genes", contrast, "|", celltype),
    base_annotations = list(
      "Intersection size" = intersection_size(
        text = list(size = 3),
        mapping = aes(fill = Direction)
      )
    ),
    set_sizes = upset_set_size(
      mapping = aes(fill = Region)
    )
  ) +
    scale_fill_manual(
      values = c(region_colors, Up = "#d81b60", Down = "#1565c0")
    ) +
    theme(
      text = element_text(size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}


contrasts <- c("FLASH - CONV", "FLASH - Control", "CONV - Control")
celltypes <- c("Neurons", "Astrocytes", "Microglia")
for (ct in celltypes) {
  for (contrast in contrasts) {
    
    upset_df_ct <- upset_df_all %>%
      filter(Celltype == ct, Contrast == contrast, Direction == "Down")
    
    if (nrow(upset_df_ct) == 0) next
    
    p <- make_upset_plot(
      upset_df_ct,
      regions,
      contrast,
      ct
    )
    
    ggsave(
      file.path(
        output_dir,
        paste0("Upset_", ct, "_", gsub(" ", "_", contrast), ".png")
      ),
      p,
      width = 14,
      height = 6,
      dpi = 300
    )
  }
}

ggsave(file.path(output_dir, "AllRegionsUpSet.png"),
       upsetplotall, width = 10, height = 6, dpi = 600)


#######All regions and segments together (barplot)############

# Collect results
bar_data <- list()

for(r in regions){
  for(s in segments){
    dge_file <- file.path(dge_dir, paste0("DGE_", r, "_", s, ".csv"))
    if (!file.exists(dge_file)) next
    dge <- read_csv(dge_file, show_col_types = FALSE)
    
    # FLASH vs Control
    flash_up <- dge %>% filter(Contrast == "FLASH - Control", Pr_t < p_cutoff, Estimate >  fc_cutoff)
    flash_dn <- dge %>% filter(Contrast == "FLASH - Control", Pr_t < p_cutoff, Estimate < -fc_cutoff)
    
    # CONV vs Control
    conv_up  <- dge %>% filter(Contrast == "CONV - Control", Pr_t < p_cutoff, Estimate >  fc_cutoff)
    conv_dn  <- dge %>% filter(Contrast == "CONV - Control", Pr_t < p_cutoff, Estimate < -fc_cutoff)
    
    # Combine counts
    bar_data[[paste(r, s, sep="_")]] <- tibble(
      Region   = r,
      Segment  = s,
      Set      = c("FLASH_Up", "FLASH_Down", "CONV_Up", "CONV_Down"),
      Count    = c(nrow(flash_up), nrow(flash_dn), nrow(conv_up), nrow(conv_dn))
    )
  }
}

# Merge all into one dataframe
bar_df <- bind_rows(bar_data)
bar_df$Segment <- factor(bar_df$Segment, levels = c("Neurons", "Astrocytes", "Microglia"))
bar_df$Region <- factor(bar_df$Region, levels = c("Cerebellum", "DG", "CA1", "CA3"))

# Plot bar chart
bar_plot <- ggplot(bar_df, aes(x = Set, y = Count, fill = Set)) +
  geom_bar(stat = "identity") +
  facet_grid(Segment ~ Region, scales = "fixed") +
  scale_fill_manual(values = c(
    "FLASH_Up"   = "#ffb3d9",
    "FLASH_Down" = "#b3d9ff",
    "CONV_Up"    = "#9c0038",
    "CONV_Down"  = "#0f3a7a"
  )) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text = element_text(size = 20),
    legend.position = "bottom",
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    axis.title = element_text(size = 18)
  ) +
  labs(
    title = "",
    x = "",
    y = "Number of significant genes",
    fill = "Set"
  )

# Save plot and table
ggsave(file.path(output_dir, "AllRegions_AllSegments_BarCounts.png"),
       bar_plot, width = 14, height = 10, dpi = 600)

write_csv(
  bar_df,
  file.path(output_dir, "AllRegions_AllSegments_BarCounts_Table.csv")
)

wide_table <- bar_df %>%
  pivot_wider(
    id_cols = c(Region, Segment),
    names_from = Set,
    values_from = Count
  ) %>%
  arrange(Region, Segment)

table <- wide_table %>%
  gt() %>%
  tab_header(
    title = "Significant Gene Counts Across Regions and Cell Types",
    subtitle = "FLASH and CONV compared to Control"
  ) %>%
  cols_label(
    Segment    = "Cell Type",
    FLASH_Up   = "FLASH Up",
    FLASH_Down = "FLASH Down",
    CONV_Up    = "CONV Up",
    CONV_Down  = "CONV Down"
  ) %>%
  fmt_number(
    columns = vars(FLASH_Up, FLASH_Down, CONV_Up, CONV_Down),
    decimals = 0
  ) %>%
  tab_spanner(
    label = "FLASH",
    columns = c(FLASH_Up, FLASH_Down)
  ) %>%
  tab_spanner(
    label = "CONV",
    columns = c(CONV_Up, CONV_Down)
  ) %>%
  tab_style(
    style = cell_fill(color = "#f8f8f8"),
    locations = cells_body(
      rows = seq(2, nrow(wide_table), 2)
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  opt_align_table_header(align = "center") %>%
  tab_options(
    table.font.size = 12
  )

gtsave(table,filename = file.path(output_dir, "PathwayCounts_Table.png"))


message("Saved bar plot for significant gene counts.")

