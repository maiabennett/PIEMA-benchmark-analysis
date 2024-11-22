# ============================================================
# File: process-inputs.R
# Description: This script processes input data for the PIEMA benchmark analysis.
# Author: Maia Bennett-Boehm
# University of Nebraska at Omaha
# ============================================================

source("../TCRPaired/util/filtering.R")
source("../TCRPaired/util/processing.R")
source("../TCRPaired/util/formatting.R")

# PIEMA benchmark data selection
# Import high confidence data
high.confidence <- read.csv("./data/high-confidence-paired-sequences.csv")

# Import all reference data (for negative data generation)
all.reference.distinct <- read.csv("./data/all-distinct-cdr-paired-sequences.csv")

# Compare column for combined CDR3a and CDR3b sequences
compare.cdrs <- high.confidence %>%
    mutate(Compare = paste0(CDR3a, CDR3b))

# Filter out sequences with > 90% CDR3 similarity
high.confidence.similarity.90 <- filterSequenceSimilarity(compare.cdrs, "Compare", 0.9)
data.frame(Distinct_values = sapply(high.confidence.similarity.90, function(x) length(unique(x))))

write.csv(high.confidence.similarity.90 %>% select(-Compare), "./data/high-confidence-similarity-90-paired-sequences.csv", row.names = FALSE)

# This list contains all epitopes with > 25 receptors in 90% CDR3 sequence similarity high confidence dataset (i.e., "run 3: larger pool benchmark")
# For reasons of simplifying the matter of cross-reactivity in logistic classifier training, the YVLDHLIVV epitope is excluded from this list
target.epitopes.info <- fetchEpitopes(high.confidence.similarity.90, column = "Epitope", threshold = 25)
target.epitopes <- target.epitopes.info %>% 
    filter(Epitope != "YVLDHLIVV") %>%
    pull(Epitope)
# Here, although there is a variable number of receptors per epitope, the number of negative receptors is 30, as the characteristics of the negative data are further constrained and not all will be successfully modeled by Rosetta
n.negative <- 30
exclude.source = "chan"
# Additionally, the CMV epitope KLG is excluded from the list of negative epitopes as it is disproportionally represented in the dataset, as are highly prevalent epitopes from EBV (in order to avoid overfitting and possible cross-reactivity problems)
all.reference.distinct <- all.reference.distinct %>% 
    filter(!(Epitope %in% c("KLGGALQAK",
        "AVFDRKSDAK",
        "RAKFKQLL",
        "IVTDFSVIK",
        "RLRAEAQVK"
    ))) 
negative.data <- generateNegatives(all.reference.distinct, n.negative, target.epitope = NULL, target.epitopes, exclude.source)
negative.data <- negative.data %>% 
    group_by(Epitope) %>%   
    mutate(clone.id = paste0("decoy_", Epitope, "_", row_number())) %>%
    ungroup()

write.csv(negative.data, "./data/piema-benchmark-negative-sequences.csv", row.names = FALSE)

# For the expanded benchmark dataset, the process differs a bit: 30 receptors are sampled from the high confidence dataset for each epitope with counts > 25 are included, as there may be a number of receptors that fail to be modeled by Rosetta
positive.data <- high.confidence.similarity.90 %>% 
    filter(Epitope %in% target.epitopes) %>%
    group_by(Epitope) %>% 
    sample_n(30) %>% 
    ungroup() 

write.csv(positive.data, "./data/piema-benchmark-positive-sequences.csv", row.names = FALSE)

# Expanded benchmark dataset
all.input <- rbind(
    positive.data %>% select(clone.id, AV, CDR3a, AJ, BV, CDR3b, BJ),
    negative.data %>% select(clone.id, AV, CDR3a, AJ, BV, CDR3b, BJ))

write.csv(all.input, "./data/piema-benchmark-input-sequences.csv", row.names = FALSE)
