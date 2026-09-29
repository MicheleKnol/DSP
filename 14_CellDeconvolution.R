###############################################################
# 14_CellDeconvolution.R - Cell Deconvolution Analysis
# Michèle Knol
# 22-10-2025
###############################################################

library(SpatialDecon)
library(GeomxTools)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(tidyr)

#-----------------------------
# Directories
#-----------------------------
output_dir <- "Results/CellDeconvolution"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -----------------------------
# Load data
# -----------------------------
target_data <- readRDS("unprocessed_data.rds")
data <- readRDS("preprocessed_data.rds")
featureType(data) <- "Target"

# -----------------------------
# Expression matrices
# -----------------------------
norm_mat <- assayDataElement(data, "q_norm")
raw_mat  <- exprs(data)

# Use raw_mat if norm_mat is too flat
use_raw <- FALSE
expr_mat <- if (use_raw) raw_mat else norm_mat
expr_mat <- as.matrix(expr_mat)
rownames(expr_mat) <- make.names(rownames(expr_mat), unique = TRUE)

# -----------------------------
# Background matrix
# -----------------------------
mean_neg <- assayDataApply(
  target_data[fData(target_data)$CodeClass == "Negative", ],
  MARGIN = 2,
  FUN = ngeoMean,
  elt = "exprs"
)
mean_neg <- mean_neg[colnames(expr_mat)]
bg_mat <- sweep(expr_mat * 0, 2, mean_neg, "+")

# -----------------------------
# Load cell profile matrix
# -----------------------------
profile_file <- "Brain_AllenBrainAtlas.RData"
env <- new.env()
load(profile_file, envir = env)

# Use the correct object
profile_mat <- env$profile_matrix
profile_mat <- as.matrix(profile_mat)

# Clean up rownames
rownames(profile_mat) <- make.names(rownames(profile_mat), unique = TRUE)
rownames(expr_mat) <- make.names(rownames(expr_mat), unique = TRUE)

# Check overlap
shared_genes <- intersect(rownames(expr_mat), rownames(profile_mat))
message("Shared genes: ", length(shared_genes), " / ", nrow(profile_mat))


expr_sub <- expr_mat[shared_genes, , drop = FALSE]
bg_sub   <- bg_mat[shared_genes, , drop = FALSE]
profile_sub <- profile_mat[shared_genes, , drop = FALSE]

# -----------------------------
# Run SpatialDecon
# -----------------------------
message("Running SpatialDecon...")
res <- spatialdecon(norm = expr_sub, bg = bg_sub, X = profile_sub)
saveRDS(res, file = file.path(output_dir, "spatialdecon_result.rds"))

# -----------------------------
# Save results
# -----------------------------
write.csv(res$beta, file = file.path(output_dir, "deconv_beta.csv"))
write.csv(res$prop_of_all, file = file.path(output_dir, "deconv_prop_of_all.csv"))
write.csv(res$p, file = file.path(output_dir, "deconv_pvals.csv"))

# -----------------------------
# Cell type mapping
# -----------------------------
celltype_map <- c(
  "Astro" = "Astrocytes", "Micro.PVM" = "Microglia", "Oligo" = "Oligodendrocytes",
  "Endo" = "Endothelial", "SMC.Peri" = "Vascular", "VLMC" = "Vascular", "V3d" = "Vascular",
  "CA1" = "Neurons", "CA2" = "Neurons", "CA3" = "Neurons", "DG" = "Neurons",
  "L5.PT.CTX" = "Neurons", "L5.IT.CTX" = "Neurons", "L5.PPP" = "Neurons",
  "L6.CT.CTX" = "Neurons", "L6b.CTX" = "Neurons", "L3.RSP.ACA" = "Neurons",
  "Sst" = "Neurons", "Sst.Chodl" = "Neurons", "Pvalb" = "Neurons", "NP.PPP" = "Neurons",
  "Car3" = "Neurons", "Sncg" = "Neurons", "CR" = "Neurons", "CT.SUB" = "Neurons",
  "Meis2" = "Neurons"
)

# -----------------------------
# Tidy proportions
# -----------------------------
prop_all <- res$prop_of_all
prop_all[!is.finite(prop_all)] <- 0
annot <- pData(data)
common_samples <- intersect(colnames(prop_all), rownames(annot))
prop_all <- prop_all[, common_samples]
annot <- annot[common_samples, , drop = FALSE]

prop_tidy <- as.data.frame(t(prop_all)) %>%
  tibble::rownames_to_column("Sample") %>%
  pivot_longer(-Sample, names_to = "CellType", values_to = "Prop") %>%
  left_join(annot %>% tibble::rownames_to_column("Sample"), by = "Sample") %>%
  mutate(BroadCellType = recode(CellType, !!!celltype_map, .default = "Other"))

# -----------------------------
# Summary and visualization
# -----------------------------
summary_df <- prop_tidy %>%
  group_by(Region, Segment, BroadCellType) %>%
  summarise(MeanProp = mean(Prop), .groups = "drop") %>%
  filter(MeanProp > 0) %>%
  group_by(Region, Segment) %>%
  mutate(NormMeanProp = MeanProp / sum(MeanProp)) %>%
  ungroup()

# Heatmap
broad_mat <- prop_tidy %>%
  group_by(BroadCellType, Sample) %>%
  summarise(Prop = mean(Prop), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Prop, values_fill = 0)

broad_mat_mtx <- as.matrix(broad_mat[,-1])
rownames(broad_mat_mtx) <- broad_mat$BroadCellType

png(file.path(output_dir, "broad_celltypes_heatmap.png"), width = 2400, height = 1400, res = 200)
pheatmap(broad_mat_mtx, cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Broad Cell-Type Proportions")
dev.off()

# Barplot
ggplot(summary_df, aes(x = BroadCellType, y = MeanProp, fill = BroadCellType)) +
  geom_col() +
  facet_grid(Segment ~ Region, scales = "free_y") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Mean Broad Cell-Type Proportions")
ggsave(file.path(output_dir, "broad_celltypes_barplot.png"), width = 14, height = 10, dpi = 300)

message("Deconvolution completed and visualizations saved.")