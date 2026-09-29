###############################################################
# 17_SmallHeatmaps.R - 
# Michèle Knol
# 10-02-2026
###############################################################

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
output_dir <- "Results/Heatmaps"

# -------------------
# Define genes and region and celltype
# -------------------
#cerebellum micro		c("Pwwp2a", "Ifna1", "Lnx1", "Fam187a", "Klf5")

#ca1 neurons		c("Dnajc21", "Dsg1a", "Gm52537")
# ca1 astro		c("Car15", "Clec2j", "Col1a1", "Dnal4", "Gm13420", "Gm35486", "Ints13", "Tram2")

#ca3 neurons		c("Clec2j", "Dock3", "Efna1", "Emx1", "Gm11554", "Smcr8", "V1ra8")
#CA3 astro		c("2300009A05Rik", "Acvr2a", "Agbl4", "Anapc7", "Arg2", "Atf3", "BC024139", "Bsn", "Calm1", 
#             "Camk2a","Ccl2", "Clec2j", "Cpe", "Eif4a2", "F8a", "Gabarap", "Ggt7","Git1", "Grik4", "Inpp5f", 
#             "Kiz", "Krtdap", "Morf4l2", "Ndufs6", "Ubb", "Ube2ql1", "Wfdc5", "Lipt2", "Nsdhl", "Olfr872")

#dg neurons		c("Clec2j", "H2-M1", "Olfr995")
#dg astro		c("1700093K21Rik", "Adck5", "Chrac1", "Chst7", "E130218I03Rik", "Entpd3", "Gm38499", "Hexim1", 
#            "Lamtor4", "Mlana", "Plp1", "Pus7l", "Spata9", "Vmn1r211", "Zfp414", "Arl13a", "Nacc1", 
#            "Ntsr2", "Nucb2", "Tfec")

#FLASH regions astro     c("Ampd2", "Atraid","Icosl","Olfr1113","Olfr713","Prss57","Ptbp2","Sdf2l1","Sobp",
#                            "AA792892","Atmin","Ces1d","Csf2rb","Gja5","Gm35486","Scrn3",
#                            "Clec2j","Eif3k","Erfl","Exd1","Lemd2","Oprd1","Pgrmc1","Ppp2r3a","R3hdml","Rims4","Sprr1b","Strc","Trim27")

#FLASH regions neurons  c("Clec2j", "Csf2rb2","Fev","Gm52537","Gm52799","Gm7827","Mib2","Olfr255","Olfr507","Prl7a1")

#CONV regions astro     c("Clec2j", "Trh")

#CONV regions neurons   c("Clec2j")


genes_of_interest <- c("Clec2j", "Trh")

celltype <- "Astrocytes"

region <- "CA1"

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

# -------------------
# Generate heatmap
# -------------------
mat <- assayDataElement(Data, "q_norm")[genes_of_interest, ]
ann <- pData(Data)[, c("Group","Region","Segment"), drop = FALSE]

# Subset only matching region + segment
keep_samples <- ann$Region == region & ann$Segment == celltype & ann$Group %in% c("Control", "CONV", "FLASH")
mat <- mat[, keep_samples, drop = FALSE]
ann <- ann[keep_samples, , drop = FALSE]
mat <- mat[, order(ann$Group)]
ann <- ann[order(ann$Group), ]

ann$Group <- factor(ann$Group, levels = names(group_colors))


# png(file.path(output_dir, paste0(region, "_", celltype, "_Heatmap.png")),
#     width = 3000, height = 1500, res = 300)
# 
# heatmap <- pheatmap(mat,
#                     annotation_col = ann,
#                     annotation_colors = list(Group = group_colors),
#                     cluster_rows = TRUE,
#                     cluster_cols = FALSE,
#                     show_colnames = TRUE,
#                     scale = "row",
#                     main = paste("Common Radiation-Induced Genes in", region,"|",celltype))
# 
# draw(heatmap)
# dev.off()
# 
# pdf(file.path(output_dir,paste0(region, "_", celltype, "_Heatmap.pdf")),width = 6,height = 3,useDingbats = FALSE)
# print(heatmap)
# dev.off()


# -------------------
# Heatmap grouped per group
# -------------------
mat <- as.matrix(mat)

expr_grouped <- mat %>%
  t() %>%
  as.data.frame() %>%
  mutate(Group = ann$Group) %>%
  group_by(Group) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  column_to_rownames("Group") %>%
  t()


# png(file.path(output_dir, paste0(region, "_", celltype, "_HeatmapGrouped.png")),
#     width = 1500, height = 1000, res = 300)
# 
# heatmap_group <- pheatmap(expr_grouped,
#                     cluster_rows = TRUE,
#                     cluster_cols = FALSE,
#                     scale = "row",
#                     main = paste0("Common Radiation-Induced Genes | ", region, " ", celltype ))
# 
# draw(heatmap_group)
# dev.off()
# 
# pdf(file.path(output_dir,paste0(region, "_", celltype, "_HeatmapGrouped.pdf")),width = 4,height = 5,useDingbats = FALSE)
# print(heatmap_group)
# dev.off()


# -------------------
# Heatmap grouped per region
# -------------------
regions  <- c("CA1", "CA3", "DG")
dge_dir  <- "Results/DGE"

sig_genes_df <- data.frame()

for (r in regions) {
  dge_file <- file.path(dge_dir, paste0("DGE_", r, "_", celltype, ".csv"))
  dge_data <- read_csv(dge_file, show_col_types = FALSE) %>%
    filter(
      Contrast == "CONV - Control",
      Pr_t < 0.05
    ) %>%
    mutate(Region = r) %>%
    select(Gene, Region, Pr_t, log2FC = Estimate)
  
  sig_genes_df <- bind_rows(sig_genes_df, dge_data)
}

log2fc_mat <- sig_genes_df %>%
  filter(Gene %in% genes_of_interest) %>%
  group_by(Gene, Region) %>%
  summarise(log2FC = mean(log2FC), .groups = "drop") %>%
  pivot_wider(
    names_from = Region,
    values_from = log2FC
  ) %>%
  column_to_rownames("Gene") %>%
  as.matrix()


png(file.path(output_dir, paste0(celltype, "_CONV_vs_Control_MeanPerRegion.png")),
    width = 1500, height = 1000, res = 300)

heatmap_region <- pheatmap(
  log2fc_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "none",
  na_col = "grey90",
  main = paste0("FLASH vs Control log2fc | Common Genes in Hippocampal Regions |",celltype))
draw(heatmap_region)
dev.off()

pdf(file.path(output_dir, paste0(celltype, "_CONV_vs_Control_MeanPerRegion.pdf")),width = 5, height = 4, useDingbats = FALSE)
print(heatmap_region)
dev.off()