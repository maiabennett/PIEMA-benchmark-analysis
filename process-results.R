# ============================================================
# File: process-results.R
# Description: Script for processing and analyzing PIEMA benchmark results
# Author: Maia Bennett-Boehm
# University of Nebraska at Omaha
# ============================================================


# Load necessary libraries
library(tidyverse)

# Specify paths
negative.data.path <- "./results/negative/"
positive.data.path <- "./results/positive/"
out.path <- "./analysis/apbs/without-cross-reactives/CDR3-similarity/all/"
# out.path <- "./analysis/apbs/without-cross-reactives/full-similarity/all/"


# Make analysis directory
dir.create(out.path, showWarnings = FALSE, recursive = TRUE)


# Import data
negative.data <- data.frame()
positive.data <- read.csv(paste0(positive.data.path, "positive_final_results.csv"))

# Import all .csv files in the folder "./results/negative/apbs" and append them to negative.data
files <- list.files(path = negative.data.path, pattern = "*.csv", full.names = TRUE)
for (file in files) {
    temp.data <- read.csv(file)
    negative.data <- rbind(negative.data, temp.data)
}

# Clean data
# If '_rot' still present in ID column, remove all instances
negative.data <- negative.data %>%
    mutate(id = str_replace_all(id, "_rot", ""))

positive.data <- positive.data %>%
    mutate(id = str_replace_all(id, "_rot", ""))

# Break ID column into separate columns for each receptor
negative.data <- negative.data %>%
    separate(id, into = c("ref.id", "samp.id"), sep = "(?<=\\d)_", remove =FALSE)

positive.data <- positive.data %>%
    separate(id, into = c("ref.id", "samp.id"), sep = "(?<=\\d)_", remove =FALSE)

# Add columns for reference epitope and sample epitope
negative.data <- negative.data %>%
    mutate(ref.epitope = str_extract(ref.id, "(?<=_).*(?=_[^_]*$)"),
           samp.epitope = str_extract(samp.id, "(?<=_).*(?=_[^_]*$)"))

positive.data <- positive.data %>%
    mutate(ref.epitope = str_extract(ref.id, "(?<=_).*(?=_[^_]*$)"),
           samp.epitope = str_extract(samp.id, "(?<=_).*(?=_[^_]*$)"))

# Remove negative v negative (two decoy) and positive v positive (two non-decoy) comparisons (negative data only)
negative.data <- negative.data %>%
    filter(!(str_detect(ref.id, "decoy") & str_detect(samp.id, "decoy"))) %>%
    filter(str_detect(ref.id, "decoy") | str_detect(samp.id, "decoy"))


# If removing subsets of data, remove here
remove.substring <- "YVLDHLIVV"
negative.data <- negative.data %>%
    filter(!(str_detect(ref.epitope, remove.substring) | str_detect(samp.epitope, remove.substring)))
positive.data <- positive.data %>%
    filter(!(str_detect(ref.epitope, remove.substring) | str_detect(samp.epitope, remove.substring)))

# remove.substring <- "NLVPMVATV"
# negative.data <- negative.data %>%
#     filter(!(str_detect(ref.epitope, remove.substring) | str_detect(samp.epitope, remove.substring)))
# positive.data <- positive.data %>%
#     filter(!(str_detect(ref.epitope, remove.substring) | str_detect(samp.epitope, remove.substring)))



# Combine data for various results files
# Separate receptor:receptor pairs into true positive matches (same epitope), decoy receptor pairs (different epitopes), and unlikely receptor pairs (non-matches)
# Separate positive data into true positive matches (same epitope) and unlikely matches (different epitopes)
# Define a function to check for complex matches
isMatch <- function(ref_epitope, samp_epitope) {
    # Check if ref_epitope is a substring of samp_epitope or vice versa
    return(grepl(ref_epitope, samp_epitope) || grepl(samp_epitope, ref_epitope))
}

# Separate positive data into true positive matches (complex matches) and unlikely matches (non-matches)
positive.data <- positive.data %>%
    mutate(is.match = mapply(isMatch, ref.epitope, samp.epitope))

mismatch.data <- positive.data %>%
    filter(!is.match) %>% 
    select(-is.match)

positive.data <- positive.data %>%
    filter(is.match) %>% 
    select(-is.match)

# Create a combined data frame 
all.data <- bind_rows(
    positive.data %>% mutate(type = "True receptor paired kernels"),
    negative.data %>% mutate(type = "Decoy receptor paired kernels"),
    mismatch.data %>% mutate(type = "Unlikely receptor paired kernels")
)

# Remove erroneous YVL-GIL decoys
all.data <- all.data %>%
    filter(!(ref.epitope == "YVLDHLIVV" & samp.epitope == "GILGFVFTL" & type == "Decoy receptor pairs"))

# Make subgraph, receptor, and master information tables
# Data table for receptor:receptor subgraphs
all.subgraph.data <- all.data %>%
    group_by(id, sg_group) %>%
    mutate(avg.distance = mean(distance), avg.kdist = mean(kdist), avg.corr = mean(corr), count = n()) %>%
    filter(kdist == max(kdist)) %>%
    ungroup()
all.subgraph.data <- all.subgraph.data %>%
    rename_with(~ paste0("top.", .), .cols = 1:5)

# Data table for receptor:receptor matches
all.receptor.data <- all.data %>%
    group_by(id) %>%
    mutate(avg.distance = mean(distance), avg.kdist = mean(kdist), avg.corr = mean(corr), count = n()) %>%
    filter(kdist == max(kdist)) %>%
    ungroup()

all.receptor.data <- all.receptor.data %>%
    rename_with(~ paste0("top.", .), .cols = 1:9) %>%
    mutate(type = str_replace(type, "paired kernels", "pairs"))

# Join all data with initial data (for sequence similarity calculation)
# Import initial data
# positive.reference.data <- read.csv("./data/piema-benchmark-positive-sequences.csv") 
# negative.reference.data <- read.csv("./data/piema-benchmark-negative-sequences.csv")
# positive.reference.data <- read.csv("./data/piema-benchmark-positive-sequences-full.csv") 
# negative.reference.data <- read.csv("./data/piema-benchmark-negative-sequences-full.csv")
positive.reference.data <- read.csv("./data/piema-benchmark-positive-sequences-expanded.csv") 
negative.reference.data <- read.csv("./data/piema-benchmark-negative-sequences-expanded.csv")
reference.data <- positive.reference.data %>%
    select(all_of(names(negative.reference.data %>% 
        select(-True.epitope, -True.epitope.gene, -True.epitope.species)))) %>%
    rbind(negative.reference.data %>% 
        select(-True.epitope, -True.epitope.gene, -True.epitope.species)) %>% 
    select(clone.id, CDR1a, CDR2a, CDR2.5a, CDR3a, CDR1b, CDR2b, CDR2.5b, CDR3b, full.seq)

# Join to results
all.receptor.data <- all.receptor.data %>%
    left_join(reference.data, by = c("ref.id" = "clone.id")) %>%
    left_join(reference.data, by = c("samp.id" = "clone.id"), suffix = c(".ref", ".samp"))

# Calculate CDR3 sequence similarity for each receptor pair (i.e., row)
all.receptor.data <- all.receptor.data %>%
    mutate(CDR3combined.ref = paste0(CDR3a.ref, CDR3b.ref),
           CDR3combined.samp = paste0(CDR3a.samp, CDR3b.samp),
           CDR3.similarity = 1 - stringdist::stringdist(CDR3combined.ref, CDR3combined.samp, method = "lv") / nchar(CDR3combined.ref))

# Calculate combined CDR sequence similarity for each receptor pair
all.receptor.data <- all.receptor.data %>%
    mutate(CDRcombined.ref = paste0(CDR1a.ref, CDR2a.ref, CDR2.5a.ref, CDR3combined.ref, CDR1b.ref, CDR2b.ref, CDR2.5b.ref, CDR3combined.ref),
           CDRcombined.samp = paste0(CDR1a.samp, CDR2a.samp, CDR2.5a.samp, CDR3combined.samp, CDR1b.samp, CDR2b.samp, CDR2.5b.samp, CDR3combined.samp),
           CDR.similarity = 1 - stringdist::stringdist(CDRcombined.ref, CDRcombined.samp, method = "lv") / nchar(CDRcombined.ref))

# Calculate full sequence similarity for each receptor pair
all.receptor.data <- all.receptor.data %>%
    mutate(full.similarity = 1 - stringdist::stringdist(full.seq.ref, full.seq.samp, method = "lv") / nchar(full.seq.ref))

# Remove unnecessary columns, retaining CDR3a and CDR3b sequences in addition to new computed columns
all.receptor.data <- all.receptor.data %>% 
    select(-CDR1a.ref, -CDR2a.ref, -CDR2.5a.ref, -CDR1b.ref, -CDR2b.ref, -CDR2.5b.ref,
        -CDR1a.samp, -CDR2a.samp, -CDR2.5a.samp, -CDR1b.samp, -CDR2b.samp, -CDR2.5b.samp,
        -full.seq.ref, -full.seq.samp, -CDRcombined.ref, -CDRcombined.samp, -CDR3combined.ref, -CDR3combined.samp)

# Perform a final filtering step to exclude receptor pairs that slipped by the initial filtering
threshold <- 0.9
threshold.feature <- "CDR3.similarity"
# threshold.feature <- "full.similarity"
all.receptor.data <- all.receptor.data %>%
    filter(!!sym(threshold.feature) <= threshold)

# Filter kernel matches based on receptor pair matches
all.data <- all.data %>%
    filter(id %in% all.receptor.data$id)


# Write final data to file
write.csv(all.data, file = paste0(out.path, "final_all_kernels_data.csv"), row.names = FALSE)
write.csv(all.receptor.data, file = paste0(out.path, "final_receptor_data.csv"), row.names = FALSE)

# Count the number of receptors and receptor data sources for each epitope
# Combine unique values of ref.id and ref.epitope with samp.id and samp.epitope
unique.receptor.data <- bind_rows(
    all.receptor.data %>% select(id = ref.id, epitope = ref.epitope),
    all.receptor.data %>% select(id = samp.id, epitope = samp.epitope)
) %>% distinct()

# Add decoy substring to epitope where decoy is present in id
unique.receptor.data <- unique.receptor.data %>%
    mutate(epitope = case_when(
        str_detect(id, "decoy") ~ paste0(epitope, " (decoy)"),
        TRUE ~ epitope
    ))

# Count the number of receptors for each epitope
count.receptor.data <- unique.receptor.data %>%
    group_by(epitope) %>%
    summarise(count = n_distinct(id))

# Reshape so that counts for epitopes with substring "(decoy)"" are added as a separate column to the original epitope
count.receptor.data <- count.receptor.data %>%
    separate(epitope, into = c("epitope", "decoy"), sep = " \\(") %>%
    pivot_wider(names_from = decoy, values_from = count, values_fill = list(count = 0)) %>%
    dplyr::rename_with(~ c("positive", "decoy"), .cols = 2:3)

# Then, add counts where epitope == "GILGFVFTL_YVLDHLIVV" to counts where epitope == "YVLDHLIVV" and remove original row
# count.receptor.data <- count.receptor.data %>%
#     mutate(positive = ifelse(epitope == "YVLDHLIVV", positive + count.receptor.data$positive[count.receptor.data$epitope == "GILGFVFTL_YVLDHLIVV"], positive)) %>%
#     filter(epitope != "GILGFVFTL_YVLDHLIVV")

# Bind to Epitope, Epitope gene, and Epitope species
epitope.data <- positive.reference.data %>% 
    distinct(Epitope, Epitope.species, Epitope.gene)

count.receptor.data <- count.receptor.data %>%
    left_join(epitope.data, by = c("epitope" = "Epitope"))

# Additionally, for positive data, add data counts for a new table
count.method.data <- unique.receptor.data %>%
    mutate(method = str_extract(id, "^[^_]+")) %>%
    filter(method != "decoy") %>%
    group_by(epitope, method) %>%
    summarise(count = n_distinct(id)) %>%
    pivot_wider(names_from = method, values_from = count, values_fill = list(count = 0))

# Write data 
write.csv(unique.receptor.data, file = paste0(out.path, "final_unique_receptors.csv"), row.names = FALSE)
write.csv(count.receptor.data, file = paste0(out.path, "final_receptor_counts.csv"), row.names = FALSE)
write.csv(count.method.data, file = paste0(out.path, "final_method_counts.csv"), row.names = FALSE)

# Master data table
# Includes number of kernel matches, number of receptor matches, average kernel matches per receptor, average KAS per receptor for each ref.epitope:samp.epitope pairing
all.data.master <- all.receptor.data %>%
    group_by(ref.epitope, samp.epitope, type) %>%
    summarise(total.kernels = sum(count),
              total.receptor.pairs = n(),
              total.ref.receptors = n_distinct(ref.id),
              total.samp.receptors = n_distinct(samp.id),
              avg.kernels.per.receptor.pair = total.kernels / total.receptor.pairs,
              avg.distance.per.receptor.pair = mean(avg.distance),
              avg.top.kdist.per.receptor.pair = mean(top.kdist),
              avg.kdist.per.receptor.pair = mean(avg.kdist),
              avg.corr.per.receptor.pair = mean(avg.corr),
              avg.top.corr.per.receptor.pair = mean(top.corr),
              avg.subgraph.distance = mean(top.patDist),
              avg.subgraph.shape.sim = mean(top.shapSim),
              avg.median.kdist.per.receptor.pair = mean(top.med_scr),
              avg.CDR3.sequence.sim = mean(CDR3.similarity),
              avg.CDR.sequence.sim = mean(CDR.similarity),
              avg.full.sequence.sim = mean(full.similarity)) %>%
    ungroup() 

all.data.master <- bind_rows(all.data.master, 
    all.data.master %>% 
        filter(type == "Unlikely receptor pairs") %>%
        group_by(ref.epitope) %>%
        summarise(
            total.kernels = sum(total.kernels),
            total.receptor.pairs = sum(total.receptor.pairs),
            total.ref.receptors = max(total.ref.receptors),
            total.samp.receptors = sum(total.samp.receptors),
            avg.kernels.per.receptor.pair = total.kernels / total.receptor.pairs,
            avg.distance.per.receptor.pair = mean(avg.distance.per.receptor.pair),
            avg.top.kdist.per.receptor.pair = mean(avg.top.kdist.per.receptor.pair),
            avg.kdist.per.receptor.pair = mean(avg.kdist.per.receptor.pair),
            avg.corr.per.receptor.pair = mean(avg.corr.per.receptor.pair),
            avg.top.corr.per.receptor.pair = mean(avg.top.corr.per.receptor.pair),
            avg.subgraph.distance = mean(avg.subgraph.distance),
            avg.subgraph.shape.sim = mean(avg.subgraph.shape.sim),
            avg.median.kdist.per.receptor.pair = mean(avg.median.kdist.per.receptor.pair),
            avg.CDR3.sequence.sim = mean(avg.CDR3.sequence.sim),
            avg.CDR.sequence.sim = mean(avg.CDR.sequence.sim),
            avg.full.sequence.sim = mean(avg.full.sequence.sim)
        ) %>% 
        mutate(samp.epitope = "Combined", type = "Unlikely receptor pairs") %>% 
        ungroup() 
    ) %>%
    mutate(match.epitope = case_when(
        type == "Decoy receptor pairs" ~ "Decoy",
        type == "True receptor pairs" ~ "True positive",
        TRUE ~ samp.epitope
    ))


# Export data
write.csv(all.data.master, file = paste0(out.path, "results_master.csv"), row.names = FALSE)


# Select plots of particular interest from plot-results.R