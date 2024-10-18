# ============================================================
# File: process-inputs.R
# Description: This script processes input data for the PIEMA benchmark analysis.
# Author: Maia Bennett-Boehm
# University of Nebraska at Omaha
# ============================================================

source("../paired-tcr-data/util/filtering.R")
source("../paired-tcr-data/util/processing.R")
source("../paired-tcr-data/util/formatting.R")

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

# Filter out sequences with > 90% full sequence similarity
high.confidence.full.similarity.90 <- filterSequenceSimilarity(high.confidence, "full.seq", 0.9)
data.frame(Distinct_values = sapply(high.confidence.full.similarity.90, function(x) length(unique(x))))

write.csv(high.confidence.full.similarity.90, "./data/high-confidence-full-similarity-90-paired-sequences.csv", row.names = FALSE)

# Generate negative data with specified epitope list (generates evenly distributed data for epitopes)
# This list contains all epitopes with > 50 receptors in 90% CDR3 sequence similarity high confidence dataset (i.e., "run 1: CDR3 similarity benchmark")
target.epitopes.info <- fetchEpitopes(high.confidence.similarity.90, column = "Epitope", threshold = 50)
target.epitopes <- target.epitopes.info %>% pull(Epitope)
n.negative <- 50
exclude.source = "chan"
negative.data <- generateNegatives(all.reference.distinct, n.negative, target.epitope = NULL, target.epitopes, exclude.source)
negative.data <- negative.data %>% 
    group_by(Epitope) %>%   
    mutate(clone.id = paste0("decoy_", Epitope, "_", row_number())) %>%
    ungroup() 

write.csv(negative.data, "./data/piema-benchmark-negative-sequences.csv", row.names = FALSE)

# This list contains all epitopes with > 49 receptors in 90% full sequence similarity high confidence dataset; 49 is used to include YLQ, which is of interest and present in the first dataset (i.e., "run 2: full sequence similarity benchmark")
target.epitopes.info.full <- fetchEpitopes(high.confidence.full.similarity.90, column = "Epitope", threshold = 49)
target.epitopes.full <- target.epitopes.info.full %>% pull(Epitope)
n.negative <- 50
exclude.source = "chan"
negative.data.full <- generateNegatives(all.reference.distinct, n.negative, target.epitope = NULL, target.epitopes.full, exclude.source)
negative.data.full <- negative.data.full %>% 
    group_by(Epitope) %>%   
    mutate(clone.id = paste0("decoy_", Epitope, "_", row_number())) %>%
    ungroup()

write.csv(negative.data.full, "./data/piema-benchmark-negative-sequences-full.csv", row.names = FALSE)

# This list contains all epitopes with > 25 receptors in 90% CDR3 sequence similarity high confidence dataset (i.e., "run 3: larger pool benchmark")
# For reasons of simplifying the matter of cross-reactivity in logistic classifier training, the YVLDHLIVV epitope is excluded from this list
target.epitopes.info.large <- fetchEpitopes(high.confidence.similarity.90, column = "Epitope", threshold = 25)
target.epitopes.large <- target.epitopes.info.large %>% 
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
negative.data.large <- generateNegatives(all.reference.distinct, n.negative, target.epitope = NULL, target.epitopes.large, exclude.source)
negative.data.large <- negative.data.large %>% 
    group_by(Epitope) %>%   
    mutate(clone.id = paste0("decoy_", Epitope, "_", row_number())) %>%
    ungroup()

write.csv(negative.data.large, "./data/piema-benchmark-negative-sequences-expanded-benchmark.csv", row.names = FALSE)

# Randomly sample 50 receptors per epitope from high confidence dataset
# In this dataset, 36 are known to be cross-reactive between GIL and YVL epitopes
# Their inclusion is purposeful and may be used to evaluate the ability of a method to distinguish between cross-reactive epitopes or discarded later
# To minimize confusion, all cross-reactive receptors are assigned to the YVL epitope
positive.data <- high.confidence.similarity.90 %>% 
    mutate(Epitope = case_when(Epitope == "GILGFVFTL,YVLDHLIVV" ~ "YVLDHLIVV", TRUE ~ Epitope)) %>%
    filter(Epitope %in% target.epitopes) %>%
    group_by(Epitope) %>% 
    sample_n(50) %>% 
    ungroup()

write.csv(positive.data, "./data/piema-benchmark-positive-sequences.csv", row.names = FALSE)

# The same process is repeated for the full sequence similarity dataset
positive.data.full <- high.confidence.full.similarity.90 %>% 
    mutate(Epitope = case_when(Epitope == "GILGFVFTL,YVLDHLIVV" ~ "YVLDHLIVV", TRUE ~ Epitope)) %>%
    filter(Epitope %in% target.epitopes.full) %>%
    filter(Epitope != "YLQPRTFLL") %>%
    group_by(Epitope) %>% 
    sample_n(50) %>% 
    ungroup() %>%
    rbind(high.confidence.full.similarity.90 %>% filter(Epitope == "YLQPRTFLL"))

write.csv(positive.data.full, "./data/piema-benchmark-positive-sequences-full.csv", row.names = FALSE)

# For the expanded benchmark dataset, the process differs a bit: 30 receptors are sampled from the high confidence dataset for each epitope with counts > 25 are included, as there may be a number of receptors that fail to be modeled by Rosetta
positive.data.large <- high.confidence.similarity.90 %>% 
    filter(Epitope %in% target.epitopes.large) %>%
    group_by(Epitope) %>% 
    sample_n(30) %>% 
    ungroup() 

write.csv(positive.data.large, "./data/piema-benchmark-positive-sequences-expanded.csv", row.names = FALSE)

# Final data for benchmarking matrix: 
# True binders: 50 receptors per epitope
# "Unlikely" binders: the 250 true binding receptors not specific to epitope in question
# Negative binders: 50 randomly generated binding receptors not specific to epitope in question (drawn from full pool of receptors, not high confidence receptors)

all.input <- rbind(
    positive.data %>% select(clone.id, AV, CDR3a, AJ, BV, CDR3b, BJ),
    negative.data %>% select(clone.id, AV, CDR3a, AJ, BV, CDR3b, BJ))

write.csv(all.input, "./data/piema-benchmark-input-sequences.csv", row.names = FALSE)

# Full sequence similarity dataset
all.input.full <- rbind(
    positive.data.full %>% select(clone.id, AV, CDR3a, AJ, BV, CDR3b, BJ),
    negative.data.full %>% select(clone.id, AV, CDR3a, AJ, BV, CDR3b, BJ))

write.csv(all.input.full, "./data/piema-benchmark-input-sequences-full.csv", row.names = FALSE)

# Expanded benchmark dataset
all.input.large <- rbind(
    positive.data.large %>% select(clone.id, AV, CDR3a, AJ, BV, CDR3b, BJ),
    negative.data.large %>% select(clone.id, AV, CDR3a, AJ, BV, CDR3b, BJ))

write.csv(all.input.large, "./data/piema-benchmark-input-sequences-expanded.csv", row.names = FALSE)
