# ============================================================
# File: process-results.R
# Description: Script for processing and analyzing PIEMA benchmark results
# Author: Maia Bennett-Boehm
# University of Nebraska at Omaha
# ============================================================


# Load necessary libraries
library(tidyverse)

# Import data
negative.data.path <- "./results/negative/apbs"
positive.data.path <- "./results/positive/apbs"
out.path <- "./analysis/apbs"

negative.data <- data.frame()
positive.data <- read.csv(paste0(positive.data.path, "positive_final_results.csv"))

# Import all .csv files in the folder "./results/negative/apbs" and append them to negative.data
files <- list.files(path = negative.data.path, pattern = "*.csv", full.names = TRUE)
for (file in files) {
    temp.data <- read.csv(file)
    negative.data <- rbind(negative.data, temp.data)
}

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

# Create basic density plot
ggplot(all.data, aes(x = kdist, color = type)) +
    geom_density() +
    labs(title = "Density plot of KAS values",
                x = "KAS",
                y = "Density",
                color = "Receptor:receptor pair type") +
    theme_minimal()


# Make per-epitope subtables
# Create a combined data frame with duplicated rows for 'Unlikely binding pairs' to allow facet wrapping
all.data.faceted <- all.data %>%
    mutate(epitope = ref.epitope) %>%
    bind_rows(
        all.data %>%
            filter(type == "Unlikely receptor paired kernels") %>%
            mutate(epitope = samp.epitope)
    )

all.data.minus.cross <- all.data.faceted %>%
    filter(!str_detect(epitope, "_"))

# Create density plot with facet wrapping by epitope
ggplot(all.data.minus.cross, aes(x = kdist, color = type)) +
    geom_density() +
    labs(title = "Density plot of KAS values by epitope, not including cross-reactive receptors",
         x = "KAS",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

# Make information tables
# Data table for receptor:receptor matches
all.receptor.data <- all.data %>%
    group_by(id) %>%
    mutate(avg.kdist = mean(kdist), count = n()) %>%
    filter(kdist == max(kdist)) %>%
    ungroup()

all.receptor.data <- all.receptor.data %>%
    rename_with(~ paste0("top.", .), .cols = 1:9) %>%
    mutate(type = str_replace(type, "paired kernels", "pairs"))

# Create density plot of top KAS values 
ggplot(all.receptor.data, aes(x = top.kdist, color = type)) +
    geom_density() +
    labs(title = "Density plot of top KAS values for each receptor:receptor pair",
         x = "KAS",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

# Facet wrap by epitope
all.receptor.data.faceted <- all.receptor.data %>%
    mutate(epitope = ref.epitope) %>%
    bind_rows(
        all.receptor.data %>%
            filter(type == "Unlikely receptor pairs") %>%
            mutate(epitope = samp.epitope)
    )

# Create density plot of top KAS values for each receptor:receptor pair, faceted by epitope
ggplot(all.receptor.data.faceted, aes(x = top.kdist, color = type)) +
    geom_density() +
    labs(title = "Density plot of top KAS values for each receptor:receptor pair, faceted by epitope",
         x = "KAS",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

# Density plot of kernel match counts
ggplot(all.receptor.data, aes(x = count, color = type)) +
    geom_density() +
    labs(title = "Kernel match counts per receptor:receptor pair",
         x = "Count",
         y = "Density",
         fill = "Receptor:receptor pair type") +
    theme_minimal()

# Facet wrap by epitope
ggplot(all.receptor.data.faceted, aes(x = count, color = type)) +
    geom_density() +
    labs(title = "Kernel match counts per receptor:receptor pair, faceted by epitope",
         x = "Count",
         y = "Density",
         fill = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)




#per.epitope.tables <- all.data %>%
#    group_by(ref.epitope, samp.epitope, type) %>%
#    summarize(count = n()) %>%
#    spread(key = type, value = count, fill = 0)

# Create bar plot with counts by epitope
