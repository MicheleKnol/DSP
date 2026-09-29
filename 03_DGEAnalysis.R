###############################################################
# 03_dge_analysis.R - Differential gene expression (DGE)
# Michèle Knol
# 14-09-2025
###############################################################

library(limma)
library(edgeR)
library(dplyr)
library(DESeq2)
library(Biobase)

# -------------------
# Define output directories
# -------------------
output_dir  <- "Results/DGE"

# -------------------
# Load preprocessed object
# -------------------
Data <- readRDS("preprocessed_data.rds")

# -------------------
# Function to run DGE per region/segment
# -------------------
regions <- c("CA1", "CA3", "DG")
segments <- c("Neurons", "Astrocytes")

dge <- function(data, output_dir) {
  results <- list()
  
  for (r in regions) {
    for (s in segments) {
      
      message("Running mixed-model DGE for: ", r, " | ", s)
      
      # Subset for region and segment
      sub <- which(pData(data)$Region == r &
                         pData(data)$Segment == s)
      if (length(sub) < 3) next
      
      temp_Data <- data[, sub]
      pData(temp_Data)$slide     <- factor(pData(temp_Data)$Slide_Name)
      pData(temp_Data)$testGroup <- factor(pData(temp_Data)$Group,
                                           levels = c("FLASH","CONV","Control"))
      
      # Run LMM
      mixedOutmc <- mixedModelDE(
        temp_Data,
        elt = "log_q",
        modelFormula = ~ testGroup + (1|slide),
        groupVar = "testGroup",
        nCores = 1,
        multiCore = FALSE
      )
      
      # Format results to a df
      r_test <- do.call(rbind, mixedOutmc['lsmeans', ])
      r_test <- as.data.frame(r_test)
      r_test$Contrast <- rownames(r_test)
      
      # Multiple testing correction
      colnames(r_test)[which(names(r_test) == "Pr(>|t|)")] <- "Pr_t"
      r_test$adjp <- p.adjust(r_test$Pr_t, method = "fdr")
      
      # Add metadata
      r_test$Gene   <- unlist(lapply(colnames(mixedOutmc),
                                     rep, nrow(mixedOutmc["lsmeans", ][[1]])))
      r_test$Subset <- paste(r, s, sep = "_")
      r_test$MeanExp <- rowMeans(assayDataElement(temp_Data, "q_norm"))
      
      # Save all contrasts for this region+segment
      fname <- paste0(output_dir, "/DGE_", r, "_", s, ".csv")
      write.csv(r_test, file = fname, row.names = FALSE)
      
      results[[paste(r, s, sep = "_")]] <- r_test
      message("Finished mixed-model DGE for: ", r, " | ", s)
    }
  }
  
  return(results)
}


# -------------------
# Run DGE and format contrast column to correct names
# -------------------
dge_results <- dge(Data, output_dir)

for(r in regions){
  for(s in segments){
    dge_file <- file.path(output_dir, paste0("DGE_", r, "_", s, ".csv"))
    if (file.exists(dge_file)) {
      dge <- read.csv(dge_file, stringsAsFactors = FALSE)
      message("Loaded: ", dge_file)
    } else {
      message("Skipped (not found): ", dge_file)
      next
    }
    dge$Contrast <- gsub("\\.\\.\\.", " - ", dge$Contrast)
    dge$Contrast <- gsub("\\.\\d+$", "", dge$Contrast)
    write.csv(dge, dge_file, row.names = FALSE)
  }
}

message("DGE complete.")