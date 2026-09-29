###############################################################
# 02_QCPlots.R - Quality Control Plots
# Michèle Knol
# 14-09-2025
###############################################################

library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(dplyr)
library(ggplot2)
library(ggforce)
library(networkD3)
library(knitr)
library(scales)
library(cowplot)
library(patchwork)
library(rlang)
library(uwot)
library(Rtsne)
library(tidyr)
library(dplyr)
library(ggforce)

# -------------------
# Define output directories
# -------------------
output_dir  <- "Results/QC"

# -------------------
# Load preprocessed object
# -------------------
Data <- readRDS("preprocessed_data.rds")

# -------------------
# Define colors
# -------------------
group_colors <- c("Control" = "#87cefa", "CONV" = "#9acd32", "FLASH" = "#ff69b4")
# region_colors <- c("VTA" = "#B0F2B4", "Cortex" = "#F2E2BA", "Cerebellum" = "#F2BAC9", "CA1" = "#BAF2E9", "CA3" = "#BAD7F2", "DG" = "#D6C9DE")
region_colors <- c( "Cerebellum" = "#F2BAC9", "CA1" = "#F2E2BA", "CA3" = "#BAD7F2", "DG" = "#B0F2B4")
segment_colors <- c("Neurons" = "#E58EAA", "Astrocytes" = "#7A96C3", Microglia = "#F2C94C")
# segment_colors <- c("Neurons" = "#E58EAA", "Glial cells" = "#8F8A93")

# ------------------------------------------------------------
# Experimental design visualization (Sankey)
# ------------------------------------------------------------
pData_df <- as.data.frame(pData(Data))

make_sankey <- function(df) {
  sankeyCols <- c("source", "target", "value")
  
  # Create links
  link1 <- dplyr::count(df, Group, Mm_ID)
  link2 <- dplyr::count(df, Mm_ID, Region)
  link3 <- dplyr::count(df, Region, Segment)
  colnames(link1) <- sankeyCols
  colnames(link2) <- sankeyCols
  colnames(link3) <- sankeyCols
  links <- rbind(link1, link2, link3)
  
  # Standardize node names
  links$source <- trimws(links$source)
  links$target <- trimws(links$target)
  
  # Create nodes
  nodes <- data.frame(name = unique(c(links$source, links$target)), stringsAsFactors = FALSE)

  # Assign node colors
  nodes$color <- sapply(nodes$name, function(x) {
    if (x %in% names(group_colors)) return(group_colors[[x]])
    if (x %in% names(region_colors)) return(region_colors[[x]])
    if (x %in% names(segment_colors)) return(segment_colors[[x]])
    return("#CCCCCC")
  })
  
  # Map indices
  links$source <- match(links$source, nodes$name) - 1
  links$target <- match(links$target, nodes$name) - 1
  
  # Build JS color scale
  color_scale <- sprintf(
    'd3.scaleOrdinal().domain(["%s"]).range(["%s"])',
    paste(nodes$name, collapse = '","'),
    paste(nodes$color, collapse = '","')
  )
  
  # Create Sankey network
  sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    colourScale = JS(color_scale),
    fontSize = 14,
    nodeWidth = 30
  )
}

sankey <- make_sankey(pData_df)
saveNetwork(sankey, file.path(output_dir, "sankey.html"))

# ------------------------------------------------------------
# QC Histograms (sequencing, alignment, saturation, morphology)
# ------------------------------------------------------------
QC_histogram <- function(df, var, fill_by = "Region", thr = NULL, scale_trans = NULL, title = NULL) {
  
  # Ensure numeric column
  x <- df[[var]]
  if (is.data.frame(x) || is.list(x)) x <- unlist(x)
  
  df_plot <- df[!is.na(x), , drop = FALSE]
  x <- x[!is.na(x)]
  
  plt <- ggplot(df_plot, aes(x = x, fill = .data[[fill_by]])) +
    geom_histogram(bins = 50, color = "white") +
    theme_light(base_size = 14) +
    labs(x = var,
         y = "# Segments",
         title = ifelse(is.null(title), var, title),
         fill = fill_by) +
    theme(legend.position = "bottom",
          panel.grid.major = element_line(color = "grey90"),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "grey98"))
  
  # Apply custom colors
  if (fill_by == "Region") plt <- plt + scale_fill_manual(values = region_colors)
  if (fill_by == "Group") plt <- plt + scale_fill_manual(values = group_colors)
  if (fill_by == "Segment") plt <- plt + scale_fill_manual(values = segment_colors)
  
  # Threshold line
  if (!is.null(thr)) plt <- plt + geom_vline(xintercept = thr, linetype = "dashed", color = "red")
  
  # Scale transformation
  if (!is.null(scale_trans)) {
    if (scale_trans == "log10") plt <- plt + scale_x_continuous(trans = "log10")
    else if (scale_trans == "sqrt") plt <- plt + scale_x_continuous(trans = "sqrt")
    else plt <- plt + scale_x_continuous(trans = scale_trans)
  }
  
  return(plt)
}

# Generate QC plots
qc_plots <- list(
  raw       = QC_histogram(sData(Data), "Raw", thr = 1e5, title = "Raw Sequencing Reads"),
  trimmed   = QC_histogram(sData(Data), "Trimmed (%)", thr = 80),
  stitched  = QC_histogram(sData(Data), "Stitched (%)", thr = 80),
  aligned   = QC_histogram(sData(Data), "Aligned (%)", thr = 75),
  saturated = QC_histogram(sData(Data), "Saturated (%)", thr = 50, title = "Sequencing Saturation (%)"),
  area      = QC_histogram(sData(Data), "Area", thr = 1000, scale_trans = "log10"),
  nuclei    = QC_histogram(sData(Data), "Nuclei", thr = 25)
)

# Save plots
for (n in names(qc_plots)) ggsave(file.path(output_dir, paste0(n, "_qc.png")), qc_plots[[n]], dpi = 300)

# ------------------------------------------------------------
# Plot with number of nuclei per condition per celltype per region
# ------------------------------------------------------------
#Cerebellum
nuclei_plot_cere <- ggplot(
  pData(Data) %>% filter(Region %in% c("Cerebellum")) %>% filter(Segment %in% c("Microglia")),
  aes(x = Group, y = Nuclei, fill = Group)
) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = Segment), size = 1,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  scale_fill_manual(values = segment_colors) +
  scale_color_manual(values = segment_colors) +
  facet_grid(Segment ~ Region, scales = "free") +
  theme_bw(base_size = 14) +
  labs(title = "Number of Nuclei per Region / Cell Type / Group",
       x = "Group", y = "Number of Nuclei", fill = "Cell Type", color = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(output_dir, "Nuclei_Cere.png"), nuclei_plot_cere, width = 10, height = 6)
pdf(file.path(output_dir, "Nuclei_Cere.pdf"), width = 10, height = 6, useDingbats = FALSE)
print(nuclei_plot_cere)
dev.off()

#Hippocampus
nuclei_plot_hippo <- ggplot(
  pData(Data) %>% filter(Region %in% c("CA1", "CA3", "DG")) %>% filter(Segment %in% c("Neurons", "Astrocytes")),
  aes(x = Group, y = Nuclei, fill = Segment)
) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = Segment), size = 1,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  scale_fill_manual(values = segment_colors) +
  scale_color_manual(values = segment_colors) +
  facet_grid(Segment ~ Region, scales = "free") +
  theme_bw(base_size = 14) +
  labs(title = "Number of Nuclei per Region / Cell Type / Group",
       x = "Group", y = "Number of Nuclei", fill = "Cell Type", color = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(output_dir, "Nuclei_Hippo.png"), nuclei_plot_hippo, width = 10, height = 6)
pdf(file.path(output_dir, "Nuclei_Hippo.pdf"), width = 10, height = 6, useDingbats = FALSE)
print(nuclei_plot_hippo)
dev.off()


# ------------------------------------------------------------
# QC Flags summary (PASS / WARNING)
# ------------------------------------------------------------
QCresults <- protocolData(Data)[["QCFlags"]]
QCresults$QCStatus <- ifelse(rowSums(QCresults) == 0, "PASS", "WARNING")

QC_Summary <- data.frame(
  Pass = sum(QCresults$QCStatus == "PASS"),
  Warning = sum(QCresults$QCStatus == "WARNING")
)

kable(QC_Summary, caption = "QC Summary Table")
write.csv(QC_Summary, file.path(output_dir, "QC_Summary.csv"))

# ------------------------------------------------------------
# LOQ & gene detection rate
# ------------------------------------------------------------
# Add a factor column for nicer plotting
pData(Data)$DetectionThreshold <- cut(
  pData(Data)$GeneDetectionRate,
  breaks = c(0, 0.01, 0.05, 0.1, 0.15, 0.3, 1),
  labels = c("<1%","1-5%", "5-10%", "10-15%", "15-30%", ">30%")
)

# Gene detection rate plots per Region
plot_detection_rate <- function(df, x_var = "DetectionThreshold", fill_var = "Region", facet_vars = NULL, filename) {
  
  plt <- ggplot(df, aes(x = .data[[x_var]], fill = .data[[fill_var]])) +
    geom_bar(position = "stack") +
    theme_bw(base_size = 14) +
    labs(x = "Gene Detection Rate", y = "# Segments", fill = fill_var) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Apply custom colors
  if (fill_var == "Region") plt <- plt + scale_fill_manual(values = region_colors)
  if (fill_var == "Group") plt <- plt + scale_fill_manual(values = group_colors)
  if (fill_var == "Segment") plt <- plt + scale_fill_manual(values = segment_colors)
  
  if (!is.null(facet_vars)) plt <- plt + facet_wrap(facet_vars, scales = "fixed")
  
  ggsave(file.path(output_dir, filename), plt, width = 8, height = 6)
  return(plt)
}

# Plots
detection_region <- plot_detection_rate(pData(Data), facet_vars = NULL,
                                        filename = "GeneDetection_byRegion.png")

detection_segment <- plot_detection_rate(pData(Data), facet_vars = "Segment",
                                         filename = "GeneDetection_bySegment.png")

detection_group <- plot_detection_rate(pData(Data), facet_vars = "Group",
                                       filename = "GeneDetection_byGroup.png")


pdf(file.path(output_dir,paste0("GeneDetection_bySegment.pdf")),width = 8,height = 6,useDingbats = FALSE)
print(detection_segment)
dev.off()

# Plot: stacked bar with facets
gene_detection_plot <- ggplot(pData(Data), 
                              aes(x = DetectionThreshold, fill = Region)) +
  geom_bar(position = "stack") +
  theme_bw(base_size = 14) +
  scale_y_continuous(expand = expansion(mult = c(0,0.1))) +
  scale_fill_manual(values = region_colors) +  # <-- Apply your region_colors here
  labs(x = "Gene Detection Rate", y = "# Segments", fill = "Region") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(Region ~ Segment, scales = "fixed")

ggsave(file.path(output_dir, "GeneDetection_Rate_Faceted.png"), gene_detection_plot, width = 10, height = 10)

# Summarize gene detection per Group/Region/Segment
gene_detection_summary <- pData(Data) %>%
  group_by(Group, Region, Segment) %>%
  summarise(
    Mean = mean(GeneDetectionRate, na.rm = TRUE),
    SD = sd(GeneDetectionRate, na.rm = TRUE),
    .groups = "drop"
  )

# Violin + boxplot overlay
gene_detection_violin <- ggplot(pData(Data), 
                                aes(x = Region, y = GeneDetectionRate, fill = Group)) +
  geom_boxplot(width = 0.5, position = position_dodge(width = 0.9), outlier.size = 0.5) +
  geom_point(data = gene_detection_summary, 
             aes(x = Region, y = Mean, color = Group),
             position = position_dodge(width = 0.9), size = 2) +
  facet_wrap(~Segment, scales = "free_x") +
  theme_bw(base_size = 14) +
  labs(x = "Region", y = "Gene Detection Rate", fill = "Group", color = "Group") +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")


ggsave(file.path(output_dir, "GeneDetectionRate_Violin.png"), gene_detection_violin, width = 10, height = 6)

# ------------------------------------------------------------
# Dimensionality reduction (UMAP)
# ------------------------------------------------------------
expr_mat <- assayDataElement(Data, "q_norm")
log_expr <- log2(expr_mat + 1) # avoid log(0)

# UMAP
set.seed(24)
umap_out <- uwot::umap(t(log_expr), n_neighbors = 15, min_dist = 0.1, metric = "correlation")
pData(Data)$UMAP1 <- umap_out[,1]
pData(Data)$UMAP2 <- umap_out[,2]

# UMAP / tSNE plotting function
plot_embedding <- function(df, x, y, color_by = "Region", title = NULL, filename) {
  
  plt <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], color = .data[[color_by]])) +
    geom_point(size = 3) +
    theme_bw(base_size = 14) +
    labs(title = title, color = color_by) +
    theme(legend.position = "bottom")
  
  # Apply custom colors
  if (color_by == "Region") plt <- plt + scale_color_manual(values = region_colors)
  if (color_by == "Group") plt <- plt + scale_color_manual(values = group_colors)
  if (color_by == "Segment") plt <- plt + scale_color_manual(values = segment_colors)
  
  ggsave(file.path(output_dir, filename), plt, width = 8, height = 6)
  return(plt)
}

umap_celltype <- plot_embedding(pData(Data), "UMAP1", "UMAP2", color_by = "Segment",
                            title = "UMAP - Cell type", filename = "UMAP_byCelltype.png")

umap_plot <- plot_embedding(pData(Data), "UMAP1", "UMAP2", color_by = "Region",
                            title = "UMAP - Region", filename = "UMAP_byRegion.png")

umap_facet <- plot_embedding(pData(Data), "UMAP1", "UMAP2", color_by = "Group",
                             title = "UMAP - Group", filename = "UMAP_byGroup.png")

pdf(file.path(output_dir,paste0("UMAP_byCelltype.pdf")),width = 8,height = 6,useDingbats = FALSE)
print(umap_celltype)
dev.off()

pdf(file.path(output_dir,paste0("UMAP_byRegion.pdf")),width = 8,height = 6,useDingbats = FALSE)
print(umap_plot)
dev.off()

# Faceted by region, shape=segment, color=group
plot_embedding_facet <- function(df, x, y, color_by = "Group", shape_by = "Segment", facet_by = "Region", filename) {
  
  # Map colors and shapes
  color_values <- switch(color_by,
                         "Group" = group_colors,
                         "Region" = region_colors,
                         "Segment" = segment_colors,
                         NULL)
  
  plt <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], color = .data[[color_by]], shape = .data[[shape_by]])) +
    geom_point(size = 3) +
    theme_bw(base_size = 14) +
    labs(title = paste0(x, " vs ", y), color = color_by, shape = shape_by) +
    facet_wrap(as.formula(paste("~", facet_by)), scales = "fixed") +
    scale_color_manual(values = color_values) +
    theme(legend.position = "bottom",
          axis.text.x = element_text(hjust = 1))
  
  ggsave(file.path(output_dir, filename), plt, width = 10, height = 6)
  return(plt)
}

# UMAP plots
umap_facet <- plot_embedding_facet(pData(Data), 
                                   "UMAP1", "UMAP2", 
                                   color_by = "Group", shape_by = "Segment", facet_by = "Region", 
                                   filename = "UMAP_Group_Segment_Facet.png")


# Add ellipsoids for cell types (Segment)  
plt <- ggplot(pData(Data), aes(x = .data[["UMAP1"]], y = .data[["UMAP2"]],
                      color = .data[["Group"]],
                      shape = .data[["Segment"]])) +
  geom_point(size = 2, alpha = 1) +
  
  stat_ellipse(aes(group = .data[["Segment"]], color = .data[["Segment"]]),
               type = "norm", level = 0.95, size = 1, alpha = 0.6) +
  
  facet_wrap(as.formula(paste("~", "Region")), scales = "fixed") +
  scale_color_manual(values = c(group_colors, segment_colors)) +
  theme_bw(base_size = 14) +
  labs(title = paste0("UMAP1", " vs ", "UMAP2"),
       color = "Group", shape = "Segment") +
  theme(legend.position = "bottom",
        axis.text.x = element_text(hjust = 1))

ggsave(file.path(output_dir, "UMAP_Facet_Ellipsoids.png"), plt, width = 10, height = 6)

# ----------------------------
# Normalization Check Plots
# ----------------------------
# Before normalization (raw counts)
raw_expr <- assayDataElement(Data, "exprs")
raw_df <- as.data.frame(raw_expr)
raw_df$Gene <- rownames(raw_df)

raw_long <- raw_df %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Segment",
    values_to = "Expression"
  ) %>%
  mutate(Stage = "Raw")

# After normalization (q_norm)
norm_df <- as.data.frame(expr_mat)
norm_df$Gene <- rownames(norm_df)

norm_long <- norm_df %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Segment",
    values_to = "Expression"
  ) %>%
  mutate(Stage = "Normalized")

# Combine
combined_long <- bind_rows(raw_long, norm_long)

# Define colors for stages
stage_colors <- c("Raw" = "#ff9999", "Normalized" = "#87cefa")

# Plot
norm_plot <- ggplot(combined_long, aes(x = Stage, y = Expression, fill = Stage)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = stage_colors) +
  theme_bw(base_size = 14) +
  labs(title = "Expression Distribution Before/After Normalization",
       y = "Expression (log2 counts + 1)", x = "") +
  scale_y_continuous(trans = "log10") +
  theme(legend.position = "bottom")

# Save
ggsave(file.path(output_dir, "Normalization_Check.png"), norm_plot, width = 7, height = 5)

message("QC finished. All plots and filtered object saved.")
