###############################################################
# 04_CheckDGE.R - Some checks on the DGE results
# Michèle Knol
# 17-09-2025
###############################################################

library(dplyr)
library(readr)
library(ggplot2)

# -------------------
# Define directories
# -------------------
output_dir <- "Results/PvalHistograms"
dge_dir <- "Results/DGE"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Created directory: ", output_dir)
}

# -------------------
# Load preprocessed object
# -------------------
Data <- readRDS("preprocessed_data.rds")

# -------------------
# Get regions and segments
# -------------------
regions <- unique(Data$Region)
segments <- unique(Data$Segment)

# -------------------
# Check P-adj values in DGE results
# -------------------
# Loop over regions and segments
for(r in regions){
  for(s in segments){
    dge_file <- file.path(dge_dir, paste0("DGE_", r, "_", s, ".csv"))
    
    if (file.exists(dge_file)) {
      dge <- read.csv(dge_file, stringsAsFactors = FALSE)
      message("Loaded: ", dge_file)
    } else {
      message("Skipped (not found): ", dge_file)
    }
    
    # Basic summaries
    if(length(dge$Pr_t)==1){
      cat("min raw p:", min(dge$Pr_t, na.rm=TRUE), "\n")
      cat("count raw p < 0.05:", sum(dge$Pr_t < 0.05, na.rm=TRUE), "\n")
    }
    if(length(dge$adjp)==1){
      cat("min adj p:", min(dge$adjp, na.rm=TRUE), "\n")
      cat("count adj p < 0.05:", sum(dge$adjp < 0.05, na.rm=TRUE), "\n")
    }

    # P-value histogram (raw)
    png(file=file.path(output_dir, paste0(r, "_", s, "_p_raw.png")),
        width=600, height=350)
    hist(dge$Pr_t, breaks = 50, main = paste("P-value histogram", r, s), xlab="raw p")
    dev.off()
    # If raw p are uniform (flat), little signal. If spike near 0 but adj are large -> multiple testing/power issue.
    
    # P-value histogram (adjusted)
    png(file=file.path(output_dir, paste0(r, "_", s, "_p_adj.png")),
        width=600, height=350)
    hist(dge$adjp, breaks = 50, main = paste("P-value histogram", r, s), xlab="adjusted p")
    dev.off()
    
  }
}

# -------------------
# Check direction of estimate
# -------------------

# The first level of testGroup is the reference group, the Estimate (log2 fold change) is relative to that reference.
# 
# So, if you are looking at CONV vs CTRL:
#   Estimate > 0 → gene is upregulated in CONV compared to CTRL.
#   Estimate < 0 → gene is downregulated in CONV compared to CTRL.

gene <- "Smim22"
expr <- assayDataElement(Data, "log_q")[gene, ]

df <- data.frame(
  Expression = expr,
  Group = pData(Data)$Group,
  Slide = pData(Data)$Slide_Name,
  Region = pData(Data)$Region,
  Segment = pData(Data)$Segment
)

df <- subset(df, Region == "Cerebellum" & Segment == "Microglia")

boxplotgene <- ggplot(df, aes(x = Group, y = Expression, color = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2) +
  labs(title = paste("Expression of", gene))

print(boxplotgene)

#If sign fc is swapped, do this in the scripts:
# dge_data$Estimate <- -dge_data$Estimate

#Check distribution of Estimates in DGE files
all_dge <- dge_files %>%
  lapply(read_csv, show_col_types = FALSE) %>%
  bind_rows()

ggplot(all_dge, aes(x = Estimate)) +
  geom_histogram(bins = 100) +
  theme_bw() +
  labs(title = "Distribution of all log2 fold changes across all regions/celltypes")