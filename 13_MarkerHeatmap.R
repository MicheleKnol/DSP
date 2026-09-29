###############################################################
# 13_MarkerHeatmap.R - Heatmap of marker panel for the different celltypes
# Michèle Knol
# 30-09-2025
###############################################################

library(GeomxTools)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

#-----------------------------
# Directories
#-----------------------------
output_dir <- "Results/CellMarkers"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

#-----------------------------
# Load preprocessed data
#-----------------------------
data <- readRDS("preprocessed_data.rds")

# expr_mat <- assayDataElement(data, "exprs")  # raw aggregated counts
expr_mat <- assayDataElement(data, "q_norm") # normalized counts
meta     <- pData(data)

#-----------------------------
# Define marker genes
#-----------------------------
marker_genes <- c(
  "Rbfox3", "Syn1", "Map2", "Snap25", "Tubb3", "Syt1", "Nefl",  # Neurons
  "Gfap", "Aqp4", "S100b",  # Astrocytes
  "Aif1", "Cx3cr1", "Tmem119"   # Microglia
)

# In case one of the markers is not present in the dataset
marker_genes <- marker_genes[marker_genes %in% rownames(expr_mat)]

#-----------------------------
# Build expression matrix
#-----------------------------
heatmap_mat <- expr_mat[marker_genes, ]

segment_order <- c("Neurons", "Astrocytes", "Microglia")
region_order <- c("VTA", "Cortex", "Cerebellum", "CA1", "CA3", "DG")
meta$Segment <- factor(meta$Segment, levels = segment_order)
meta$Region <- factor(meta$Region, levels = region_order)

ord <- order(meta$Segment, meta$Region)
heatmap_mat <- heatmap_mat[, ord]
meta <- meta[ord, ]

#-----------------------------
# Annotation for heatmap
#-----------------------------
ann_col <- data.frame(
  Segment = factor(meta$Segment),
  Group   = factor(meta$Group),
  Region  = factor(meta$Region)
)
rownames(ann_col) <- colnames(heatmap_mat)

colors <- list(
  Segment = c("Neurons" = "#E58EAA", "Astrocytes" = "#7A96C3", Microglia = "#F2C94C"),
  Group   = c("Control" = "#87cefa", "CONV" = "#9acd32", "FLASH" = "#ff69b4"),
  Region  = c("VTA" = "#B0F2B4", "Cortex" = "#F2E2BA", "Cerebellum" = "#F2BAC9",
              "CA1" = "#BAF2E9", "CA3" = "#BAD7F2", "DG" = "#D6C9DE")
)

# Apply log10 scaling
heatmap_mat <- log10(heatmap_mat)

# -----------------------------
# Heatmap: individual samples
# -----------------------------
png(file.path(output_dir, "marker_panel_heatmap.png"), width = 2400, height = 2000, res = 300)

pheatmap::pheatmap(
  heatmap_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  annotation_col = ann_col,
  annotation_colors = colors,
  scale = "row",
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
  fontsize_row = 10,
  fontsize_col = 6,
  main = "Marker Panel Heatmap (Individual Samples)"
)

dev.off()

# -----------------------------
# Heatmap: summary (mean per cell type)
# -----------------------------
summary_mat <- sapply(segment_order, function(seg){
  cols <- which(meta$Segment == seg)
  rowMeans(heatmap_mat[, cols, drop = FALSE], na.rm = TRUE)
})
rownames(summary_mat) <- marker_genes

# Log10 transform AFTER averaging
summary_mat <- log10(summary_mat + 1)

png(file.path(output_dir, "marker_panel_heatmap_summary.png"), width = 2400, height = 2000, res = 300)

pheatmap::pheatmap(
  summary_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  annotation_col = NULL,
  scale = "row",
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
  fontsize_row = 10,
  fontsize_col = 10,
  main = "Marker Panel Heatmap (Mean per Celltype)"
)

dev.off()


#Save as pdf
pdf(file.path(output_dir, "marker_panel_heatmap_summary.pdf"),
    width = 12, height = 10, useDingbats = FALSE)

pheatmap::pheatmap(
  summary_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_colnames = TRUE,
  annotation_col = NULL,
  scale = "row",
  color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
  fontsize_row = 10,
  fontsize_col = 10,
  main = "Marker Panel Heatmap (Mean per Celltype)"
)

dev.off()

# #-----------------------------
# #Check if there are outliers
# #-----------------------------
# # Compute total expression per sample
# sample_sums <- colSums(heatmap_mat, na.rm = TRUE)
# 
# # Threshold: 1.5 * IQR above 75th percentile
# outlier_threshold <- quantile(sample_sums, 0.75) + 1 * IQR(sample_sums) #I put it at 1*IQR, maybe a bit low
# outliers <- names(sample_sums[sample_sums > outlier_threshold])
# meta[outliers, c("Region", "Segment", "Group")]
# 
# # Remove outliers from heatmap matrix and metadata
# heatmap_mat_filt <- heatmap_mat[, !(colnames(heatmap_mat) %in% outliers)]
# meta_filt <- meta[!(rownames(meta) %in% outliers), ]
# 
# # Re-plot heatmap without outliers
# png(file.path(output_dir, "marker_panel_heatmap.png"),  width = 2400, height = 2000, res = 300)
# 
# pheatmap(
#   heatmap_mat_filt,
#   cluster_rows = FALSE,
#   cluster_cols = FALSE,
#   show_colnames = FALSE,
#   annotation_col = ann_col,
#   annotation_colors = colors,
#   scale = "none",
#   color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
#   fontsize_row = 10,
#   fontsize_col = 6,
#   main = "Marker Panel Heatmap"
# )
# 
# dev.off()