###############################################################
# 07_Geneplots.R - Gene expression boxplots
# Michèle Knol
# 24-09-2025
###############################################################

library(ggplot2)
library(GeomxTools)
library(ggpubr)
library(rstatix)
library(dplyr)
library(purrr)

#-----------------------------
# Directories
#-----------------------------
output_dir <- "Results/Geneplots"

#-----------------------------
# Gene(s) of interest and region/celltype
#-----------------------------
gene_of_interest <- "miR-124-2"
region   <- "CA1"
celltype <- "Astrocytes"

#-----------------------------
# Load preprocessed data
#-----------------------------
data <- readRDS("preprocessed_data.rds")

# -------------------
# Extract expression matrix
# -------------------
# expr_mat <- assayDataElement(data, "exprs")  # raw aggregated counts
expr_mat <- assayDataElement(data, "q_norm") # normalized counts
# expr_mat <- assayDataElement(data, "log_q")   # log2 normalized counts

 
if (!(gene_of_interest %in% rownames(expr_mat))) {
  stop("Gene ", gene_of_interest, " not found in dataset.")
}

expr_values <- expr_mat[gene_of_interest, ]
meta <- pData(data)

# combine metadata + expression
meta$Expression <- as.numeric(expr_values)
meta$CountsperNuclei <- meta$Expression / meta$Nuclei
meta$Group   <- factor(meta$Group, levels = c("Control","CONV","FLASH"))
meta$Segment <- factor(meta$Segment)
meta$Region  <- factor(meta$Region)

# -------------------
# Define colors
# -------------------
group_colors <- c("Control" = "#87cefa", "CONV" = "#9acd32", "FLASH" = "#ff69b4")
region_colors <- c("VTA" = "#B0F2B4", "Cortex" = "#F2E2BA", "Cerebellum" = "#F2BAC9", 
                   "CA1" = "#BAF2E9", "CA3" = "#BAD7F2", "DG" = "#D6C9DE")
segment_colors <- c("Neurons" = "#E58EAA", "Astrocytes" = "#7A96C3", Microglia = "#4AAE9B")
# segment_colors <- c("Neurons" = "#E58EAA", "Glial cells" = "#8F8A93")

# -------------------
# Define comparisons
# -------------------
comparisons <- list(c("Control","CONV"), c("Control","FLASH"), c("CONV","FLASH"))
segment_comp <- list(c("Astrocytes", "Microglia"), c("Microglia", "Neurons"), c("Astrocytes", "Neurons"))
# segment_comp <- list(c("Astrocytes", "Neurons"))

# -------------------
# Group × Segment faceted by Region
# -------------------
plot_g_r_s <- ggplot(meta, aes(x = Group, y = Expression, fill = Segment)) +
# plot_g_r_s <- ggplot(meta, aes(x = Group, y = CountsperNuclei, fill = Segment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6,
               position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = Segment), size = 1, alpha = 1,
              position = position_dodge(width = 0.8)) +
  facet_wrap(~Region, scales = "free_y") +
  scale_fill_manual(values = segment_colors) +
  scale_color_manual(values = segment_colors) +
  labs(
    title = gene_of_interest,
    y = paste0(gene_of_interest, " expression (counts)")
    # y = paste0(gene_of_interest, " expression (counts per nuclei)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title   = element_text(hjust = 0.5),
    plot.subtitle= element_text(hjust = 0.5, size = 10),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(output_dir, paste0(gene_of_interest, "_expression_group_segment_region.png")),
  plot = plot_g_r_s, width = 10, height = 6
)
print(plot_g_r_s)

# ------------------- 
# Region × Segment
# -------------------
plot_r_s <- ggplot(meta, aes(x = Region, y = Expression, fill = Segment)) +
# plot_r_s <- ggplot(meta, aes(x = Region, y = CountsperNuclei, fill = Segment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5,
               position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = Segment), size = 1, alpha = 1,
              position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = segment_colors) +
  scale_color_manual(values = segment_colors) +
  labs(
    title    = gene_of_interest,
    x        = "Region",
    y = paste0(gene_of_interest, " expression (counts)")
    # y = paste0(gene_of_interest, " expression (counts per nuclei)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text.x   = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(output_dir, paste0(gene_of_interest, "_expression_region_segment.png")),
  plot = plot_r_s, width = 10, height = 5
)
print(plot_r_s)

# ------------------- 
# Segment
# -------------------
# plot_s <- ggplot(meta, aes(x = Segment, y = Expression, fill = Segment)) +
# # plot_s <- ggplot(meta, aes(x = Segment, y = CountsperNuclei, fill = Segment)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.5,
#                position = position_dodge(width = 0.8)) +
#   geom_jitter(aes(color = Segment), size = 1, alpha = 1,
#               position = position_dodge(width = 0.8)) +
#   stat_compare_means(
#     aes(label = ..p.signif..),
#     method      = "wilcox.test",
#     hide.ns     = FALSE,
#     comparisons = segment_comp,
#     position    = position_dodge(width = 0.8)
#   ) +
#   scale_fill_manual(values = segment_colors) +
#   scale_color_manual(values = segment_colors) +
#   labs(
#     title    = gene_of_interest,
#     x        = "Celltype",
#     y = paste0(gene_of_interest, " expression (counts)")
#     # y = paste0(gene_of_interest, " expression (counts per nuclei)")
#   ) +
#   theme_bw(base_size = 12) +
#   theme(
#     plot.title    = element_text(hjust = 0.5),
#     plot.subtitle = element_text(hjust = 0.5, size = 10),
#     axis.text.x   = element_text(angle = 45, hjust = 1)
#   )
# 
# ggsave(
#   filename = file.path(output_dir, paste0(gene_of_interest, "_expression_segment.png")),
#   plot = plot_s, width = 10, height = 5
# )
# print(plot_s)

# -------------------
# Gene in one region + one celltype
# -------------------
meta_sub <- meta %>%
  filter(Region == region, Segment == celltype)

# plot_gene <- ggplot(meta_sub, aes(x = Group, y = Expression, fill = Group)) +
# # plot_gene <- ggplot(meta_sub, aes(x = Group, y = CountsperNuclei, fill = Group)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.5,
#                position = position_dodge(width = 0.8)) +
#   scale_fill_manual(values = group_colors) +
#   scale_color_manual(values = group_colors) +
#   geom_jitter(aes(color = Group), size = 1, alpha = 1,
#               position = position_dodge(width = 0.8)) +
#   # stat_compare_means(
#   #   aes(label = ..p.signif..),
#   #   method      = "wilcox.test",
#   #   hide.ns     = FALSE,
#   #   comparisons = comparisons,
#   #   position    = position_dodge(width = 0.8)
#   # ) +
#   labs(
#     title = paste(gene_of_interest, "in", celltype, "|", region),
#     x = "Condition",
#     y = paste0(gene_of_interest, " expression (counts)")
#     # y = paste0(gene_of_interest, " expression (counts per nuclei)")
#   ) +
#   theme_bw(base_size = 12) +
#   # scale_y_continuous(lim=c(0,max(meta_sub$Expression)*1.1)) +
#   theme(
#     plot.title    = element_text(hjust = 0.5),
#     plot.subtitle = element_text(hjust = 0.5, size = 10),
#     axis.text.x   = element_text(angle = 45, hjust = 1)
#   )


#Bars instead of boxplot
plot_gene <- ggplot(meta_sub, aes(x = Group, y = Expression, fill = Group)) +
  stat_summary(
    fun = mean,
    geom = "bar",
    position = position_dodge(width = 0.8),
    alpha = 0.75
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  geom_point(
    aes(color = Group),
    size = 1,
    alpha = 1
  ) +
  labs(
      title = paste(gene_of_interest, "in", celltype, "|", region),
      x = "Condition",
      y = paste0(gene_of_interest, " expression (counts)")
    ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  theme_minimal(base_size = 12)


ggsave(
  filename = file.path(output_dir, paste0(gene_of_interest, "_expression_", region, "_", celltype, ".png")),
  plot = plot_gene, width = 6, height = 5
)

pdf(file.path(output_dir,paste0(gene_of_interest, "_expression_", region, "_", celltype, ".pdf")),
    width = 6,height = 5,useDingbats = FALSE)
print(plot_gene)
dev.off()

print(plot_gene)

# Summary statistics
print(meta_sub %>%
        group_by(Group) %>%
        get_summary_stats(Expression, type = "mean_sd") )
