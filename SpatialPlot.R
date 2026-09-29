library(GeoMxTools)
library(ggspavis)
library(ggplot2)
library(cowplot)
library(dplyr)

# --- Load your processed data ---
# Example: load a saved GeoMxSet or SpatialExperiment object
# Replace this with your own file:
load("your_processed_geomxset_or_spe_object.RData")  # e.g., spe <- your object



# For this example, assume you're working with SpatialExperiment
# and gene expression is in log2 scale

# Set your object
spe <- your_spe_object  # Replace this with actual name if needed

# --- Choose a gene to visualize ---
gene_to_plot <- "Cdkn1a"  # Change to any gene you like

# --- Get expression and spatial coordinates ---
expr_values <- assay(spe, "logcounts")[gene_to_plot, ]
coords <- spatialCoords(spe)
meta <- colData(spe)

plot_data <- data.frame(
  x = coords[, "x"],
  y = coords[, "y"],
  Expression = expr_values,
  Region = meta$region,           # Adjust depending on your metadata
  CellType = meta$cellType,       # Adjust depending on your metadata
  Sample = meta$slide_id          # Optional: your slide/sample name
)

# --- Plot spatial expression ---
p <- ggplot(plot_data, aes(x = x, y = y)) +
  geom_point(aes(color = Expression), size = 2) +
  scale_color_gradient(low = "gray90", high = "firebrick") +
  coord_fixed() +
  facet_wrap(~ Sample) +  # Remove if only one sample
  theme_void() +
  ggtitle(paste("Spatial Expression of", gene_to_plot))

# Save
ggsave(paste0("Spatial_", gene_to_plot, ".png"), p, width = 6, height = 6)
