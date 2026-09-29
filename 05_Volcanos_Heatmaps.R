###############################################################
# 05_Volcanos_Heatmaps.R - Plot volcanoplots and heatmaps of the DGE results
# Michèle Knol
# 17-09-2025
###############################################################

library(EnhancedVolcano)
library(pheatmap)
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(textshape)
library(tibble)
library(ComplexHeatmap)
library(circlize)

# -------------------
# Define directories
# -------------------
output_dir_vol <- "Results/VolcanoPlots"
output_dir_heat <- "Results/Heatmaps_DGE"
dge_dir <- "Results/DGE"

# -------------------
# Load preprocessed object
# -------------------
Data <- readRDS("preprocessed_data.rds")

# -------------------
# Define colors
# -------------------
group_colors <- c("Control" = "#87cefa", "CONV" = "#9acd32", "FLASH" = "#ff69b4")
region_colors <- c("VTA" = "#B0F2B4", "Cortex" = "#F2E2BA", "Cerebellum" = "#F2BAC9", 
                   "CA1" = "#BAF2E9", "CA3" = "#BAD7F2", "DG" = "#D6C9DE")
segment_colors <- c("Neurons" = "#E58EAA", "Astrocytes" = "#7A96C3", Microglia = "#F2C94C")
# segment_colors <- c("Neurons" = "#E58EAA", "Glial cells" = "#8F8A93")

# -------------------
# Loop over Region and Segment
# -------------------
regions <- c("Cerebellum", "CA1", "CA3", "DG")
segments <- c("Microglia", "Astrocytes", "Neurons")

for(r in regions){
  for(s in segments){
    
    dge_file <- file.path(dge_dir, paste0("DGE_", r, "_", s, ".csv"))
    
    if (file.exists(dge_file)) {
      dge <- read_csv(dge_file, show_col_types = FALSE)
      message("Loaded: ", dge_file)
    } else {
      message("Skipped (not found): ", dge_file)
      next
    }
    
    for(comp in unique(dge$Contrast)){
      comp_data <- dge %>% filter(Contrast == comp)
      comp_data$Estimate <- as.numeric(comp_data$Estimate)
      comp_data$Pr_t <- as.numeric(comp_data$Pr_t)
      
      # # Volcano Plot
      # volcano <- EnhancedVolcano(
      #   comp_data,
      #   lab = comp_data$Gene,
      #   x = 'Estimate',
      #   y = 'Pr_t',
      #   xlab = bquote(~Log[2]~ 'fold change'),
      #   ylab = bquote(~-Log[10]~ 'p-value'),
      #   title = paste0(r, " | ", s),
      #   subtitle = comp,
      #   pCutoff = 0.01,
      #   FCcutoff = 1,
      #   pointSize = 1,
      #   labSize = 3,
      #   xlim = c(-3, 3),
      #   ylim = c(0, 6),
      #   col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
      #   legendPosition = 'right',
      #   legendLabSize = 10,
      #   legendIconSize = 3.0,
      #   drawConnectors = TRUE,
      #   widthConnectors = 0.5
      # )
      # 
      # ggsave(file.path(output_dir_vol, paste0(r, "_", s, "_", gsub(" ", "_", comp), "_Volcano.png")), volcano, width = 12, height = 6)
      
      ## -------------------
      ## Heatmap for top genes
      ## -------------------
      #per region and segment and contrast
      topgenes <- comp_data %>%
        filter(Pr_t < 0.05, abs(Estimate) > 1) %>%
        mutate(signal = abs(Estimate) * -log10(Pr_t)) %>%
        arrange(desc(signal)) %>%
        slice_head(n = 25) %>%
        pull(Gene)
      
      if(length(topgenes) > 1){
        mat <- assayDataElement(Data, "log_q")[topgenes, ]
        ann <- pData(Data)[, c("Group","Region","Segment"), drop = FALSE]
        
        # Subset only matching region + segment
        keep_samples <- ann$Region == r & ann$Segment == s & ann$Group %in% c("CONV", "FLASH")
        mat <- mat[, keep_samples, drop = FALSE]
        ann <- ann[keep_samples, , drop = FALSE]
        mat <- mat[, order(ann$Group)]
        ann <- ann[order(ann$Group), ]
        
        ann$Group <- factor(ann$Group, levels = names(group_colors))
        ann$Region <- factor(ann$Region, levels = names(region_colors))
        ann$Segment <- factor(ann$Segment, levels = names(segment_colors))
        
        if (ncol(mat) > 1) {
          
          png(file.path(output_dir_heat, paste0(r, "_", s, "_", comp, "_Heatmap.png")),
              width = 3000, height = 3000, res = 300)
          
          ht <- pheatmap(
            mat,
            annotation_col = ann,
            annotation_colors = list(
              Group = group_colors,
              Region = region_colors,
              Segment = segment_colors),
            scale = 'row',
            cluster_rows = TRUE,
            cluster_cols = FALSE,
            show_rownames = TRUE,
            show_colnames = TRUE,
            cellwidth = 20,
            cellheight = 15,
            clustering_distance_rows = "correlation",
            fontsize = 8,
            main = paste0("Top 25 genes in ", r, " | ", s)
          )
          draw(ht)
          dev.off()
          
          pdf(file.path(output_dir_heat,paste0(r, "_", s, "_", comp, "_Heatmap.pdf")),width = 6,height = 6,useDingbats = FALSE)
          print(ht)
          dev.off()
        }
      }
    }
  }
}

# -------------------
# Faceted overview volcano grid across all Regions & Segments
# -------------------

# Gather all DGE results into a single dataframe
all_dge_files <- list.files(dge_dir, pattern = "^DGE_.*\\.csv$", full.names = TRUE)

dge_all <- lapply(all_dge_files, function(f){
  df <- read_csv(f, show_col_types = FALSE)
  
  # Extract region and segment from filename (DGE_Region_Segment.csv)
  parts <- strsplit(basename(f), "_|\\.")[[1]]
  region <- parts[2]
  segment <- parts[3]
  
  df <- df %>%
    mutate(Region = region,
           Segment = segment)
  return(df)
}) %>% bind_rows()

# Filter and clean
dge_all <- dge_all %>%
  mutate(Estimate = as.numeric(Estimate),
         Pr_t = as.numeric(Pr_t),
         Contrast = factor(Contrast)) %>%
  filter(!is.na(Estimate), !is.na(Pr_t)) %>%
  filter(Region %in% c("Cerebellum", "CA1", "CA3", "DG"))

dge_all$Segment <- factor(dge_all$Segment,levels = c("Neurons", "Astrocytes", "Microglia"))
dge_all$Region <- factor(dge_all$Region,levels = c("Cerebellum", "DG", "CA3", "CA1"))

# Add significance category
dge_all <- dge_all %>%
  mutate(
    log10p = -log10(Pr_t),
    sig = case_when(
      Pr_t < 0.05 & Estimate > 1  ~ "Up",
      Pr_t < 0.05 & Estimate < -1 ~ "Down",
      TRUE                        ~ "NS"
    )
  )

# ---- Faceted grid volcano ----
volcano_grid <- ggplot(dge_all, aes(x = Estimate, y = log10p, color = sig)) +
  geom_point(alpha = 0.5, size = 0.9) +
  facet_grid(Region ~ Segment + Contrast, scales = "fixed") +
  scale_color_manual(values = c("Down" = "blue", "Up" = "red", "NS" = "grey70")) +
  theme_bw(base_size = 8) +
  theme(
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    axis.title = element_text(size = 16)
  ) +
  labs(
    x = expression(Log[2]~fold~change),
    y = expression(-Log[10]~p~value),
    color = "Direction")

# Save
ggsave(file.path(output_dir_vol, "AllRegions_AllSegments_FacetedVolcanoGrid.png"),
       volcano_grid, width = 14, height = 10, dpi = 600)

pdf(file.path(output_dir_vol,paste0("AllRegions_AllSegments_FacetedVolcanoGrid.pdf")),width = 14,height = 10,useDingbats = FALSE)
print(volcano_grid)
dev.off()

message("Saved faceted volcano grid")



# -------------------
# Heatmap summary of DGE counts
# -------------------
metrics <- c("Up", "Down", "Total")
contrasts <- unique(summary_df$Contrast)

# Static pheatmap version
for (m in metrics) {
  for (c in contrasts) {
    sub <- summary_df %>% filter(Contrast == c)
    mat_sub <- sub %>%
      select(Region, Segment, !!sym(m)) %>%
      tidyr::pivot_wider(names_from = Segment, values_from = !!sym(m), values_fill = 0) %>%
      column_to_rownames(loc=1)

    # Metric-specific colors
    col_fun <- switch(m,
                      "Up" = colorRampPalette(c("#fce4ec", "#d81b60"))(100),
                      "Down" = colorRampPalette(c("#e3f2fd", "#1565c0"))(100),
                      "Total" = colorRampPalette(c("white", "orange"))(100)
    )

    png(file.path(output_dir_heat,
                  paste0("SummaryHeatmap_", c, "_", m, ".png")),
        width = 1000, height = 1000, res = 300)

    pheatmap(as.matrix(mat_sub),
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             color = col_fun,
             main = paste0(m, " genes (", c, ")"),
             fontsize = 10)

    dev.off()
  }
}
message("Saved individual heatmaps for each contrast and metric.")



# Faceted ggplot version
summary_long <- summary_df %>%
  pivot_longer(cols = all_of(metrics),
               names_to = "Metric",
               values_to = "Value")

# Set color scales for ggplot
get_palette <- function(metric) {
  if (metric == "Up") scale_fill_gradient(low = "#fce4ec", high = "#d81b60")
  else if (metric == "Down") scale_fill_gradient(low = "#e3f2fd", high = "#1565c0")
  else scale_fill_gradient2(low = "white", high = "orange")
}

# Faceted plot per metric
for (m in metrics) {
  data_m <- summary_long %>% filter(Metric == m)

  p <- ggplot(data_m, aes(x = Segment, y = Region, fill = Value)) +
    geom_tile(color = "grey80") +
    facet_wrap(~Contrast, nrow = 1) +
    get_palette(m) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    labs(title = paste0("Summary of ", m, " genes across regions and segments"),
         x = "Segment", y = "Region")

  ggsave(file.path(output_dir_heat,
                   paste0("Faceted_", m, "_Heatmap.png")),
         p, width = 14, height = 5, dpi = 300)
}

message("Saved 3 faceted plots (one per metric, faceted by contrast).")


# Combined plot
p_all <- ggplot(summary_long, aes(x = Segment, y = Region, fill = Value)) +
  geom_tile(color = "grey80") +
  facet_grid(Metric ~ Contrast, scales = "free") +
  scale_fill_gradient2(low = "deepskyblue", mid = "white", high = "darksalmon") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(title = "Summary of Up, Down, and Total DEG counts across all contrasts",
       x = "Segment", y = "Region")

ggsave(file.path(output_dir_heat, "Faceted_AllMetrics_AllContrasts.png"),
       p_all, width = 14, height = 10, dpi = 300)

message("Saved one combined faceted plot (all metrics × all contrasts).")
