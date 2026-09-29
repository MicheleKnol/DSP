###############################################################
# 01_Preprocessing.R - Loads in DSP data, QC checks, filtering, and normalization
# Michèle Knol
# 14-09-2025
###############################################################

library(Biobase)
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(dplyr)
library(ggplot2)
library(edgeR)
library(scales)
library(openxlsx)

# -------------------
# Define the paths
# -------------------
project_dir <- "C:/Users/knol/OneDrive - Université de Genève/Documents/DSP Paediatrics"
output_dir  <- file.path(project_dir, "DSP R Analysis/Results/QC")
datadir     <- file.path(project_dir, "Data")

dcc_files <- dir(file.path(datadir, "dccs"), pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
pkc_files <- dir(file.path(datadir, "pkcs"), pattern = ".pkc$", full.names = TRUE, recursive = TRUE)
sample_anno <- file.path(datadir, "data", "Annotation.xlsx")

# -------------------
# Load data
# -------------------
geomx_set <- readNanoStringGeoMxSet(
  dccFiles = dcc_files,
  pkcFiles = pkc_files,
  phenoDataFile = sample_anno,
  phenoDataSheet = "ModifiedSheet",
  phenoDataDccColName = "dcc"
)

saveRDS(geomx_set, file = file.path(output_dir, "unprocessed_data.rds"))

#Keep a log of when and why segments are removed
audit_log <- data.frame(
  SegmentID = rownames(pData(geomx_set)),
  Group     = pData(geomx_set)$Group,
  Region    = pData(geomx_set)$Region,
  Segment   = pData(geomx_set)$Segment,
  RemovedAt = NA_character_,
  Reason    = NA_character_,
  stringsAsFactors = FALSE
)


log_removals <- function(removed_ids, step, reasons = NULL) {
  if (length(removed_ids) == 0) return(invisible(NULL))
  
  audit_log <<- within(audit_log, {
    idx <- SegmentID %in% removed_ids & is.na(RemovedAt)
    RemovedAt[idx] <- step
    
    if (is.null(reasons)) {
      Reason[idx] <- NA_character_
    } else if (length(reasons) == 1) {
      Reason[idx] <- reasons
    } else {
      Reason[idx] <- reasons
    }
  })
}

# Exclude flagged ROIs
before <- rownames(pData(geomx_set))
geomx_set <- geomx_set[, pData(geomx_set)$Exclude == 0]
after <- rownames(pData(geomx_set))

removed <- setdiff(before, after)
log_removals(removed, "Exclude_flag", "Exclude == 1")


# Keep only selected brain regions
regions_keep <- c("DG", "CA1", "CA3", "Cerebellum")
geomx_set <- geomx_set[, pData(geomx_set)$Region %in% regions_keep]

# Shift counts (to avoid problems with count=0 later on)
# geomx_set <- shiftCountsOne(geomx_set, useDALogic = TRUE)

# -------------------
# Segment QC
# -------------------
qc_params <- list(
  minSegmentReads = 1000, # Minimum number of reads (1000)
  percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
  percentStitched = 80,   # Minimum % of reads stitched (80%)
  percentAligned = 80,    # Minimum % of reads aligned (80%) (75)
  percentSaturation = 50, # Minimum sequencing saturation (50%)
  minNegativeCount = 1,   # Minimum negative control counts (1-10)
  minNuclei = 25,         # Minimum # of nuclei estimated (100)
  minArea = 1000          # Minimum segment area (5000)
)

geomx_set <- setSegmentQCFlags(geomx_set, qcCutoffs = qc_params)

# Extract QC flags
qc_flags <- protocolData(geomx_set)[["QCFlags"]]

qc_flags$LowNuclei <- pData(geomx_set)$Nuclei < qc_params$minNuclei
qc_flags$LowArea   <- pData(geomx_set)$Area   < qc_params$minArea

protocolData(geomx_set)[["QCFlags"]] <- qc_flags

#Filter segments failing any QC
before <- rownames(pData(geomx_set))
keep_qc <- apply(qc_flags, 1, function(x) all(x == FALSE))
geomx_set_filt <- geomx_set[, keep_qc]
after <- rownames(pData(geomx_set_filt))
removed <- setdiff(before, after)


failed_flags <- apply(qc_flags[removed, , drop = FALSE], 1, function(x) {
  paste(names(x)[x], collapse = ",")
})
log_removals(removed, "Segment_QC", failed_flags)


# Save QC summary
qc_summary <- data.frame(
  Region  = pData(geomx_set)$Region,
  CellType = pData(geomx_set)$Segment,
  Group   = pData(geomx_set)$Group,
  Nuclei  = pData(geomx_set)$Nuclei,
  Area    = pData(geomx_set)$Area,
  qc_flags,
  check.names = FALSE
)

write.xlsx(qc_summary, "QC_summary.xlsx", rowNames = FALSE)

#Show how many segments were removed
message(ncol(geomx_set) - ncol(geomx_set_filt), " segments were removed during QC.")

# -------------------
# Probe QC
# -------------------
geomx_set_filt <- setBioProbeQCFlags(
  geomx_set_filt,
  qcCutoffs = list(minProbeRatio = 0.1, percentFailGrubbs = 20),
  removeLocalOutliers = TRUE
)

geomx_set_filt <- subset(
  geomx_set_filt,
  fData(geomx_set_filt)[["QCFlags"]][, "LowProbeRatio"] == FALSE &
    fData(geomx_set_filt)[["QCFlags"]][, "GlobalGrubbsOutlier"] == FALSE
)

# -------------------
# Aggregate counts per target (to correct for the multiple probes per gene)
# -------------------
target_data <- aggregateCounts(geomx_set_filt)

# -------------------
# LOQ
# -------------------
# Extract modules
modules <- gsub(".pkc", "", annotation(target_data))

# Calculate negative geometric means per module
negativeGeoMeans <- esBy(
  negativeControlSubset(geomx_set),
  GROUP = "Module",
  FUN = function(x) assayDataApply(x, 2, ngeoMean, elt = "exprs")
)

# Match filtered segments to full IDs
full_ids <- rownames(pData(geomx_set))
filtered_ids <- rownames(pData(target_data))
negativeGeoMeans_filtered <- negativeGeoMeans[match(filtered_ids, full_ids)]
names(negativeGeoMeans_filtered) <- filtered_ids

protocolData(target_data)[["NegGeoMean"]] <- negativeGeoMeans_filtered

# Copy to phenoData
negCols <- paste0("NegGeoMean_", modules)
pData(target_data)[, negCols] <- protocolData(target_data)[["NegGeoMean"]]

# LOQ calculation
minLOQ <- 2
cutoff <- 2
LOQ <- data.frame(matrix(NA, nrow = ncol(target_data), ncol = length(modules),
                         dimnames = list(colnames(target_data), modules)))

for (mod in modules) {
  neg_mean_col <- paste0("NegGeoMean_", mod)
  neg_sd_col   <- paste0("NegGeoSD_", mod)
  
  if (all(c(neg_mean_col, neg_sd_col) %in% colnames(pData(target_data)))) {
    LOQ[, mod] <- pmax(
      minLOQ,
      pData(target_data)[, neg_mean_col] * (pData(target_data)[, neg_sd_col] ^ cutoff)
    )
  }
}

pData(target_data)$LOQ <- LOQ

# Build LOQ matrix
LOQ_Mat <- matrix(FALSE, nrow = nrow(target_data), ncol = ncol(target_data),
                  dimnames = list(fData(target_data)$TargetName, colnames(target_data)))

for (mod in modules) {
  ind <- fData(target_data)$Module == mod
  if (sum(ind) == 0) next
  LOQ_Mat[ind, ] <- t(esApply(target_data[ind, ], 1, function(x) x > LOQ[, mod]))
}

# Per-segment detection
pData(target_data)$GenesDetected <- colSums(LOQ_Mat, na.rm = TRUE)
pData(target_data)$GeneDetectionRate <- pData(target_data)$GenesDetected / nrow(target_data)

# Filter segments with ≥5% genes detected
LOQ_FILTERING_THRESHOLD <- 0.04

# pData(target_data) %>%
#   dplyr::filter(Group == "FLASH",
#                 Region == "CA1",
#                 Segment == "Astrocytes") %>%
#   select(GenesDetected, GeneDetectionRate)

keep_segments <- pData(target_data)$GeneDetectionRate >= LOQ_FILTERING_THRESHOLD

before <- colnames(target_data)
after  <- colnames(target_data)[keep_segments]
removed <- setdiff(before, after)
rates_before <- pData(target_data)$GeneDetectionRate[match(removed, colnames(target_data))]

target_data <- target_data[, keep_segments]

reasons <- paste0("GeneDetectionRate = ", round(rates_before, 3)," < ", LOQ_FILTERING_THRESHOLD)
log_removals(removed_ids = removed,step = "LOQ_filter",reason = reasons)

# Per-gene detection
LOQ_Mat <- LOQ_Mat[, colnames(target_data)]
fData(target_data)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_data)$DetectionRate <- fData(target_data)$DetectedSegments / ncol(target_data)

message("LOQ / Detection Rate Filtering complete.")
message(ncol(geomx_set_filt) - ncol(target_data), " segments were removed during LOQ / detection rate filtering.")

# -------------------
# Normalization
# -------------------
counts <- assayDataElement(target_data, "exprs")
counts <- ceiling(counts)

dge <- DGEList(counts = counts, group = pData(target_data)$Group)
dge <- calcNormFactors(dge)

assayDataElement(target_data, "q_norm") <- cpm(dge, log = FALSE)

# Log2 transform
assayDataElement(target_data, "log_q") <- log2(assayDataElement(target_data, "q_norm") + 1)

message("Normalization complete.")


# -------------------
# Save processed object
# -------------------
saveRDS(target_data, file = file.path(output_dir, "preprocessed_data.rds"))
message("Preprocessing complete and saved.")


# pData(geomx_set) %>%
# dplyr::filter(Group == "FLASH") %>%
# dplyr::count(Region, Segment)

audit_log_out <- audit_log[!is.na(audit_log$RemovedAt), ]

write.table(
  audit_log_out,
  file = file.path(output_dir, "segment_removal_audit.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
