###############################################################
# 06_Venn_UpSet_GenePlots.R
# Differential expression overlap plots
# Author: Michèle Knol
###############################################################

############################
# Libraries
############################
library(tidyverse)
library(ggVennDiagram)
library(ComplexUpset)
library(openxlsx)
library(msigdbr)
library(patchwork)
library(gt)
library(writexl)

############################
# Directories
############################
dge_dir    <- "Results/DGE"
output_dir <- "Results/VennGenes"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

############################
# Global parameters
############################
regions   <- c("DG", "CA1", "CA3")
segments  <- c("Neurons", "Astrocytes", "Microglia")
contrasts <- c("FLASH - Control", "CONV - Control", "FLASH - CONV")

p_cutoff  <- 0.05
fc_cutoff <- 1

region_colors <- c(
  # "Cerebellum" = "#F2BAC9",
  "DG"         = "#D6C9DE",
  "CA1"        = "#BAF2E9",
  "CA3"        = "#BAD7F2"
)

direction_colors <- c(
  "Up"   = "#d81b60",
  "Down" = "#1565c0")

############################
# MSigDB annotation
############################
genesets <- msigdbr(species = "Mus musculus") |>
  select(gs_name, gene_symbol)

annotate_pathways <- function(df) {
  df |>
    left_join(genesets, by = c("Gene" = "gene_symbol")) |>
    group_by(Gene, Estimate, Pr_t) |>
    summarise(
      Pathways = paste(unique(gs_name), collapse = "; "),
      .groups = "drop"
    )
}

############################
# Load DGE files
############################
load_dge <- function(region, segment) {
  f <- file.path(dge_dir, paste0("DGE_", region, "_", segment, ".csv"))
  if (!file.exists(f)) return(NULL)
  read_csv(f, show_col_types = FALSE)
}

############################
# Venn diagrams
############################
for (r in regions) {
  for (s in segments) {
    
    dge <- load_dge(r, s)
    if (is.null(dge)) next
    
    sig <- dge |> filter(Pr_t < p_cutoff)
    
    split_sets <- function(contrast, direction) {
      sig |>
        filter(
          Contrast == contrast,
          if (direction == "up")   Estimate >  fc_cutoff else Estimate < -fc_cutoff
        ) |>
        pull(Gene) |>
        unique()
    }
    
    venn_sets_up <- list(
      FLASH = split_sets("FLASH - Control", "up"),
      CONV  = split_sets("CONV - Control",  "up")
    )
    
    venn_sets_dn <- list(
      FLASH = split_sets("FLASH - Control", "down"),
      CONV  = split_sets("CONV - Control",  "down")
    )
    
    plot_venn <- function(sets, title, fill) {
      ggVennDiagram(sets, label = "count", label_alpha = 0) +
        scale_fill_gradient(low = "white", high = fill) +
        ggtitle(title) +
        coord_fixed() +
        theme(plot.title = element_text(hjust = 0.5))
    }
    
    ggsave(
      file.path(output_dir, paste0(r, "_", s, "_Venn_Up.png")),
      plot_venn(venn_sets_up, paste("Upregulated |", r, s), "#d81b60"),
      width = 7, height = 5, dpi = 300
    )
    pdf(file.path(output_dir,paste0(r, "_", s, "_Venn_Up.pdf")),width = 14,height = 10,useDingbats = FALSE)
    print(plot_venn(venn_sets_up, paste("Upregulated |", r, s), "#d81b60"))
    dev.off()
    
    ggsave(
      file.path(output_dir, paste0(r, "_", s, "_Venn_Down.png")),
      plot_venn(venn_sets_dn, paste("Downregulated |", r, s), "#1565c0"),
      width = 7, height = 5, dpi = 300)
    pdf(file.path(output_dir,paste0(r, "_", s, "_Venn_Down.pdf")),width = 14,height = 10,useDingbats = FALSE)
    print(plot_venn(venn_sets_dn, paste("Downregulated |", r, s), "#1565c0"))
    dev.off()
    
    # wb <- createWorkbook()
    # 
    # write_sheet <- function(name, df) {
    #   addWorksheet(wb, name)
    #   writeData(wb, name, annotate_pathways(df))
    # }
    # 
    # write_sheet("FLASH_Up",   sig |> filter(Contrast=="FLASH - Control", Estimate >  fc_cutoff))
    # write_sheet("FLASH_Down", sig |> filter(Contrast=="FLASH - Control", Estimate < -fc_cutoff))
    # write_sheet("CONV_Up",    sig |> filter(Contrast=="CONV - Control",  Estimate >  fc_cutoff))
    # write_sheet("CONV_Down",  sig |> filter(Contrast=="CONV - Control",  Estimate < -fc_cutoff))
    # 
    # saveWorkbook(
    #   wb,
    #   file.path(output_dir, paste0(r, "_", s, "_GeneLists.xlsx")),
    #   overwrite = TRUE
    # )
    
    #Four way venn
    venn_sets_4 <- list(
      FLASH_Up   = split_sets("FLASH - Control", "up"),
      FLASH_Down = split_sets("FLASH - Control", "down"),
      CONV_Up    = split_sets("CONV - Control",  "up"),
      CONV_Down  = split_sets("CONV - Control",  "down")
    )
    plot_venn_4 <- function(sets, title) {
      ggVennDiagram(sets, label = "count", label_alpha = 0) +
        scale_fill_gradient(low = "white", high = "coral") +
        ggtitle(title) +
        coord_fixed() +
        theme(plot.title = element_text(hjust = 0.5))
    }
    ggsave(
      file.path(output_dir, paste0(r, "_", s, "_Venn_UpDown_FLASH_CONV.png")),
      plot_venn_4(venn_sets_4, paste("FLASH vs CONV | Up & Down |", r, s)),
      width = 7, height = 6, dpi = 300
    )
    
    pdf(
      file.path(output_dir, paste0(r, "_", s, "_Venn_UpDown_FLASH_CONV.pdf")),
      width = 14, height = 12, useDingbats = FALSE
    )
    print(plot_venn_4(venn_sets_4, paste("FLASH vs CONV | Up & Down |", r, s)))
    dev.off()
    
  }
}

############################
# Combined Venn (faceted)
############################
venn_plots <- list()

for (r in regions) {
  for (s in segments) {
    
    dge <- load_dge(r, s)
    if (is.null(dge)) next
    
    sig <- dge |> filter(Pr_t < p_cutoff)
    
    venn_sets <- list(
      FLASH_Up   = sig |> filter(Contrast=="FLASH - Control", Estimate >  fc_cutoff) |> pull(Gene),
      FLASH_Down = sig |> filter(Contrast=="FLASH - Control", Estimate < -fc_cutoff) |> pull(Gene),
      CONV_Up    = sig |> filter(Contrast=="CONV - Control",  Estimate >  fc_cutoff) |> pull(Gene),
      CONV_Down  = sig |> filter(Contrast=="CONV - Control",  Estimate < -fc_cutoff) |> pull(Gene)
    )
    
    if (all(lengths(venn_sets) == 0)) next
    
    venn_plots[[paste(r, s, sep = "_")]] <-
      ggVennDiagram(venn_sets, label = "count", label_alpha = 0) +
      scale_fill_gradient(low = "white", high = "coral") +
      ggtitle(paste(r, s, sep = " | ")) +
      theme(plot.title = element_text(hjust = 0.5, size = 9))
  }
}

ggsave(
  file.path(output_dir, "AllRegions_AllSegments_4SetVenn.png"),
  wrap_plots(venn_plots, ncol = length(regions)),
  width = 16, height = 10, dpi = 300
)

############################
# UpSet plots
############################
upset_df_all <- map_dfr(segments, function(s) {
  map_dfr(regions, function(r) {

    dge <- load_dge(r, s)
    if (is.null(dge)) return(NULL)

    dge %>%
      filter(Pr_t < p_cutoff) %>%
      mutate(
        Direction = case_when(
          Estimate >  fc_cutoff ~ "Up",
          Estimate < -fc_cutoff ~ "Down",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(Direction)) %>%
      mutate(
        Region  = r,
        Segment = s,
        Value   = TRUE
      ) %>%
      select(Gene, Contrast, Segment, Direction, Region, Value)
  })
})

upset_df_all <- upset_df_all %>%
  pivot_wider(
    names_from  = Region,
    values_from = Value,
    values_fill = FALSE
  ) %>%
  distinct(Gene, Contrast, Segment, Direction, .keep_all = TRUE) %>%
  mutate(Direction = factor(Direction, levels = c("Up", "Down"))
  )



for (ct in segments) {
  for (con in contrasts) {
    df_ct <- upset_df_all |>
      filter(
        Segment  == ct,
        Contrast == con)
    overlap_df <- df_ct %>%
      select(Gene, Direction, all_of(regions)) %>%
      mutate(
        OverlapRegions = apply(
          select(., all_of(regions)),
          1,
          function(x) paste(regions[x], collapse = ", ")
        ),
        N_regions = rowSums(select(., all_of(regions)))
      ) %>%
      arrange(desc(N_regions), Gene)
    write_xlsx(
      overlap_df,
      path = file.path(
        output_dir,
        paste0("UpSet_overlap_", ct, "_", gsub(" ", "_", con), ".xlsx")
      )
    )
    
    p <- ComplexUpset::upset(
      df_ct,
      intersect = regions,
      set_sizes = upset_set_size(
        geom = geom_bar(width = 0.5, fill = region_colors),
        position = "left"),
      base_annotations = list(
        "Intersection size" = intersection_size(
          text = list(),
          counts = TRUE,
          mapping = aes(fill = Direction)))) +
      scale_fill_manual(
        values = direction_colors) +
      labs(title = paste("Differentially expressed genes |", con, "|", ct)) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(hjust = 0.5))
    
    ggsave(file.path(output_dir,paste0("UpSet_", ct, "_", gsub(" ", "_", con), ".png")),p,width = 8,height = 6,dpi = 600)
    pdf(file.path(output_dir,paste0("UpSet_", ct, "_", gsub(" ", "_", con), ".pdf")),width = 8,height = 6,useDingbats = FALSE)
    print(p)
    dev.off()
    }}


############################
# Barplot summary
############################
bar_df <- upset_df_all |>
  count(Region, Segment, Contrast, Direction) |>
  pivot_wider(names_from = Direction, values_from = n, values_fill = 0)

ggsave(
  file.path(output_dir, "AllRegions_AllSegments_BarCounts.png"),
  ggplot(bar_df, aes(x = Contrast, y = Up + Down, fill = Contrast)) +
    geom_col() +
    facet_grid(Segment ~ Region) +
    theme_bw(),
  width = 14, height = 10, dpi = 600
)

write_csv(bar_df, file.path(output_dir, "AllRegions_AllSegments_BarCounts_Table.csv"))




###### Old function, to implement in this script
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
  facet_grid(Region ~ Segment, scales = "fixed") +
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
       bar_plot, width = 10, height = 14, dpi = 600)

pdf(file.path(output_dir,paste0("AllRegions_AllSegments_BarCounts.pdf")),width = 10,height = 14,useDingbats = FALSE)
print(bar_plot)
dev.off()




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







message("All plots successfully generated.")
