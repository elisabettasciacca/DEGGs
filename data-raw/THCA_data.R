## code to prepare `DATASET` dataset goes here
library(RTCGAToolbox)
thcaData <- getFirehoseData(dataset="THCA", clinical=TRUE,
                            RNASeq2GeneNorm = TRUE,
                            destdir = "C:/Users/Elisabetta/"
)
THCA_norm_counts <- as.data.frame(getData(thcaData,"RNASeq2GeneNorm")[[1]]@DataMatrix)

control_samples <- colnames(THCA_norm_counts)[
  startsWith(sapply(strsplit(colnames(THCA_norm_counts), "-", fixed = TRUE),
                    "[[", 4), "1") ]

tumor_samples <- colnames(THCA_norm_counts)[
  startsWith(sapply(strsplit(colnames(THCA_norm_counts), "-", fixed = TRUE),
                    "[[", 4), "0") ]

# let's prepare two subgroups of with same size
tumor_samples_subset <- sample(tumor_samples, 59)

THCA_norm_counts <- THCA_norm_counts[, c(control_samples, tumor_samples_subset)]
THCA_metadata <- data.frame('ID' = colnames(THCA_norm_counts), stringsAsFactors = F)
THCA_metadata$case.control <- ifelse(THCA_metadata$ID %in% control_samples, "control",
                                     "tumor")
rownames(THCA_metadata) <- THCA_metadata$ID


usethis::use_data(THCA_metadata, THCA_norm_counts)
