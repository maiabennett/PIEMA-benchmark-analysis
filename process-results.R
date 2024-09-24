# ============================================================
# File: process-results.R
# Description: Script for processing and analyzing PIEMA benchmark results
# Author: Maia Bennett-Boehm
# University of Nebraska at Omaha
# ============================================================


# Load necessary libraries
library(tidyverse)

# Specify paths
negative.data.path <- "./results/negative/apbs"
positive.data.path <- "./results/positive/apbs"
out.path <- "./analysis/apbs"

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


# Make information tables
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

all.receptor.data.faceted <- all.receptor.data %>%
    mutate(epitope = ref.epitope) %>%
    bind_rows(
        all.receptor.data %>%
            filter(type == "Unlikely receptor pairs") %>%
            mutate(epitope = samp.epitope)
    )

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
              avg.subgraph.distance = mean(top.patDist),
              avg.subgraph.shape.sim = mean(top.shapSim)) %>%
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
            avg.subgraph.distance = mean(avg.subgraph.distance),
            avg.subgraph.shape.sim = mean(avg.subgraph.shape.sim)
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
write.csv(all.data.master, file = paste0(out.path, "/results_master.csv"), row.names = FALSE)

# Plots
# Density plots
# KAS scores for all kernel matches
ggplot(all.data, aes(x = kdist, color = type)) +
    geom_density() +
    labs(title = "Density plot of KAS values",
                x = "KAS",
                y = "Density",
                color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/density_plot_kas_values.png"), width = 6, height = 4)

# KAS scores for all kernel matches, faceted by epitope
ggplot(all.data.faceted, aes(x = kdist, color = type)) +
    geom_density() +
    labs(title = "Density plot of KAS values by epitope, not including cross-reactive receptors",
         x = "KAS",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "/density_plot_kas_values_per_epitope.png"), width = 6, height = 4)

# Spearman correlation for all kernel matches
ggplot(all.data, aes(x = corr, color = type)) +
    geom_density() +
    labs(title = "Spearman correlation for all kernel matches",
         x = "Spearman correlation",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    xlim(0, 1)

ggsave(paste0(out.path, "/density_plot_spearman_correlation.png"), width = 6, height = 4)

# Spearman correlation for all kernel matches, faceted by epitope
ggplot(all.data.faceted, aes(x = corr, color = type)) +
    geom_density() +
    labs(title = "Spearman correlation for all kernel matches by epitope, not including cross-reactive receptors",
         x = "Spearman correlation",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    xlim(0, 1) + 
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "/density_plot_spearman_correlation_per_epitope.png"), width = 6, height = 4)

# Euclidean distance for all kernel matches
ggplot(all.data, aes(x = distance, color = type)) +
    geom_density() +
    labs(title = "Euclidean distance for all kernel matches",
         x = "Euclidean distance",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/density_plot_euclidean_distance.png"), width = 6, height = 4)

# Euclidean distance for all kernel matches, faceted by epitope
ggplot(all.data.faceted, aes(x = distance, color = type)) +
    geom_density() +
    labs(title = "Euclidean distance for all kernel matches by epitope, not including cross-reactive receptors",
         x = "Euclidean distance",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "/density_plot_euclidean_distance_per_epitope.png"), width = 6, height = 4)

# Shape similarity (cosine distance between subgraphs)
ggplot(all.data, aes(x = shapSim, color = type)) +
    geom_density() +
    labs(title = "Shape similarity for all kernel matches",
         x = "Shape similarity",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/density_plot_shape_similarity.png"), width = 6, height = 4)

# Shape similarity (cosine distance between subgraphs), faceted by epitope
ggplot(all.data.faceted, aes(x = shapSim, color = type)) +
    geom_density() +
    labs(title = "Shape similarity for all kernel matches by epitope, not including cross-reactive receptors",
         x = "Shape similarity",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "/density_plot_shape_similarity_per_epitope.png"), width = 6, height = 4)

# Subgraph Euclidean distance
ggplot(all.data, aes(x = patDist, color = type)) +
    geom_density() +
    labs(title = "Subgraph Euclidean distance for all kernel matches",
         x = "Subgraph Euclidean distance",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/density_plot_subgraph_euclidean_distance.png"), width = 6, height = 4)

# Subgraph Euclidean distance, faceted by epitope
ggplot(all.data.faceted, aes(x = patDist, color = type)) +
    geom_density() +
    labs(title = "Subgraph Euclidean distance for all kernel matches by epitope, not including cross-reactive receptors",
         x = "Subgraph Euclidean distance",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "/density_plot_subgraph_euclidean_distance_per_epitope.png"), width = 6, height = 4)

# Top KAS values per receptor:receptor pair
ggplot(all.receptor.data, aes(x = top.kdist, color = type)) +
    geom_density() +
    labs(title = "Density plot of top KAS values for each receptor:receptor pair",
         x = "KAS",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/density_plot_top_kas_values.png"), width = 6, height = 4)

# Top KAS values per receptor:receptor pair, faceted by epitope
ggplot(all.receptor.data.faceted, aes(x = top.kdist, color = type)) +
    geom_density() +
    labs(title = "Density plot of top KAS values for each receptor:receptor pair, faceted by epitope",
         x = "KAS",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "/density_plot_top_kas_values_per_epitope.png"), width = 6, height = 4)

# Number of kernel matches per receptor:receptor pair
ggplot(all.receptor.data, aes(x = count, color = type)) +
    geom_density() +
    labs(title = "Kernel match counts per receptor:receptor pair",
         x = "Count",
         y = "Density",
         fill = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/density_plot_kernel_match_counts.png"), width = 6, height = 4)

# Number of kernel matches per receptor:receptor pair, faceted by epitope
ggplot(all.receptor.data.faceted, aes(x = count, color = type)) +
    geom_density() +
    labs(title = "Kernel match counts per receptor:receptor pair, faceted by epitope",
         x = "Count",
         y = "Density",
         fill = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "/density_plot_kernel_match_counts_per_epitope.png"), width = 6, height = 4)


# Point plots
# Top KAS value per reference receptor epitope (x axis) colored by sample receptor epitope 
ggplot(all.data.master, aes(x = ref.epitope, y = avg.top.kdist.per.receptor.pair, color = match.epitope)) +
    geom_point(aes(size = ifelse(match.epitope %in% c("True positive", "Decoy", "Combined"), 2, 1))) +
    labs(title = "Top KAS value per reference receptor epitope, colored by sample receptor epitope",
         x = "Reference epitope",
         y = "Top KAS value",
         color = "Paired receptor epitope") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_size_continuous(range = c(2, 4), guide = "none")

ggsave(paste0(out.path, "/point_plot_top_kas_values.png"), width = 6, height = 4)

# Average KAS value per reference receptor epitope (x axis) colored by sample receptor epitope
ggplot(all.data.master, aes(x = ref.epitope, y = avg.kdist.per.receptor.pair, color = match.epitope)) +
    geom_point(aes(size = ifelse(match.epitope %in% c("True positive", "Decoy", "Combined"), 2, 1))) +
    labs(title = "Average KAS value per reference receptor epitope, colored by sample receptor epitope",
         x = "Reference epitope",
         y = "Average KAS value",
         color = "Paired receptor epitope") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_size_continuous(range = c(2, 4), guide = "none")

ggsave(paste0(out.path, "/point_plot_average_kas_values.png"), width = 6, height = 4)

# Average Spearman correlation per reference receptor epitope (x axis) colored by sample receptor epitope
ggplot(all.data.master, aes(x = ref.epitope, y = avg.corr.per.receptor.pair, color = match.epitope)) +
    geom_point(aes(size = ifelse(match.epitope %in% c("True positive", "Decoy", "Combined"), 2, 1))) +
    labs(title = "Average Spearman correlation per reference receptor epitope, colored by sample receptor epitope",
         x = "Reference epitope",
         y = "Average Spearman correlation",
         color = "Paired receptor epitope") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_size_continuous(range = c(2, 4), guide = "none")

ggsave(paste0(out.path, "/point_plot_average_spearman_correlation.png"), width = 6, height = 4)

# Average Euclidean distance per reference receptor epitope (x axis) colored by sample receptor epitope
ggplot(all.data.master, aes(x = ref.epitope, y = avg.distance.per.receptor.pair, color = match.epitope)) +
    geom_point(aes(size = ifelse(match.epitope %in% c("True positive", "Decoy", "Combined"), 2, 1))) +
    labs(title = "Average Euclidean distance per reference receptor epitope, colored by sample receptor epitope",
         x = "Reference epitope",
         y = "Average Euclidean distance",
         color = "Paired receptor epitope") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_size_continuous(range = c(2, 4), guide = "none")

ggsave(paste0(out.path, "/point_plot_average_euclidean_distance.png"), width = 6, height = 4)

# Average subgraph Euclidean distance per reference receptor epitope (x axis) colored by sample receptor epitope
ggplot(all.data.master, aes(x = ref.epitope, y = avg.subgraph.distance, color = match.epitope)) +
    geom_point(aes(size = ifelse(match.epitope %in% c("True positive", "Decoy", "Combined"), 2, 1))) +
    labs(title = "Average subgraph Euclidean distance per reference receptor epitope, colored by sample receptor epitope",
         x = "Reference epitope",
         y = "Average subgraph Euclidean distance",
         color = "Paired receptor epitope") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_size_continuous(range = c(2, 4), guide = "none")

ggsave(paste0(out.path, "/point_plot_average_subgraph_euclidean_distance.png"), width = 6, height = 4)

# Average subgraph shape similarity per reference receptor epitope (x axis) colored by sample receptor epitope
ggplot(all.data.master, aes(x = ref.epitope, y = avg.subgraph.shape.sim, color = match.epitope)) +
    geom_point(aes(size = ifelse(match.epitope %in% c("True positive", "Decoy", "Combined"), 2, 1))) +
    labs(title = "Average subgraph shape similarity per reference receptor epitope, colored by sample receptor epitope",
         x = "Reference epitope",
         y = "Average subgraph shape similarity",
         color = "Paired receptor epitope") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_size_continuous(range = c(2, 4), guide = "none")

ggsave(paste0(out.path, "/point_plot_average_subgraph_shape_similarity.png"), width = 6, height = 4)

# Average number of kernel matches per receptor pair per reference receptor epitope (x axis) colored by sample receptor epitope
ggplot(all.data.master, aes(x = ref.epitope, y = avg.kernels.per.receptor.pair, color = match.epitope)) +
    geom_point(aes(size = ifelse(match.epitope %in% c("True positive", "Decoy", "Combined"), 2, 1))) +
    labs(title = "Average number of kernel matches per receptor pair per reference receptor epitope, colored by sample receptor epitope",
         x = "Reference epitope",
         y = "Average number of kernel matches per receptor pair",
         color = "Paired receptor epitope") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_size_continuous(range = c(2, 4), guide = "none")

ggsave(paste0(out.path, "/point_plot_average_kernel_matches_per_receptor_pair.png"), width = 6, height = 4)


# Dot plots
# Subgraph shape similarity vs subgraph Euclidean distance, colored by receptor pair type
ggplot(all.data.master, aes(x = avg.subgraph.distance, y = avg.subgraph.shape.sim, color = type)) +
    geom_point() +
    labs(title = "Subgraph shape similarity vs subgraph Euclidean distance",
         x = "Average subgraph Euclidean distance",
         y = "Average subgraph shape similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_subgraph_shape_similarity_vs_euclidean_distance.png"), width = 6, height = 4)

# Top KAS vs average KAS, colored by receptor pair type
ggplot(all.data.master, aes(x = avg.kdist.per.receptor.pair, y = avg.top.kdist.per.receptor.pair, color = type)) +
    geom_point() +
    labs(title = "Top KAS vs average KAS",
         x = "Average KAS",
         y = "Top KAS",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_top_kas_vs_average_kas.png"), width = 6, height = 4)

# Average Spearman vs average KAS, colored by receptor pair type
ggplot(all.data.master, aes(x = avg.kdist.per.receptor.pair, y = avg.corr.per.receptor.pair, color = type)) +
    geom_point() +
    labs(title = "Average Spearman correlation vs average KAS",
         x = "Average KAS",
         y = "Average Spearman correlation",
            color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_average_spearman_correlation_vs_average_kas.png"), width = 6, height = 4)

# Average Euclidean distance vs average KAS, colored by receptor pair type
ggplot(all.data.master, aes(x = avg.kdist.per.receptor.pair, y = avg.distance.per.receptor.pair, color = type)) +
    geom_point() +
    labs(title = "Average Euclidean distance vs average KAS",
         x = "Average KAS",
         y = "Average Euclidean distance",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_average_euclidean_distance_vs_average_kas.png"), width = 6, height = 4)

# Average subgraph shape similarity vs top KAS, colored by receptor pair type
ggplot(all.data.master, aes(x = avg.top.kdist.per.receptor.pair, y = avg.subgraph.shape.sim, color = type)) +
    geom_point() +
    labs(title = "Top KAS vs average subgraph shape similarity",
         x = "Top KAS",
         y = "Average subgraph shape similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_top_kas_vs_average_subgraph_shape_similarity.png"), width = 6, height = 4)

# Average subgraph shape similarity vs average KAS, colored by receptor pair type
ggplot(all.data.master, aes(x = avg.kdist.per.receptor.pair, y = avg.subgraph.shape.sim, color = type)) +
    geom_point() +
    labs(title = "Average KAS vs average subgraph shape similarity",
         x = "Average KAS",
         y = "Average subgraph shape similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_average_kas_vs_average_subgraph_shape_similarity.png"), width = 6, height = 4)

# Average shape similarity vs average number of kernel matches, colored by receptor pair type
ggplot(all.data.master, aes(x = avg.kernels.per.receptor.pair, y = avg.subgraph.shape.sim, color = type)) +
    geom_point() +
    labs(title = "Average number of kernel matches vs average subgraph shape similarity",
         x = "Average number of kernel matches",
         y = "Average subgraph shape similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_average_kernel_matches_vs_average_subgraph_shape_similarity.png"), width = 6, height = 4)

# Average number of kernel matches versus average Euclidean distance
ggplot(all.data.master, aes(x = avg.distance.per.receptor.pair, y = avg.kernels.per.receptor.pair, color = type)) +
    geom_point() +
    labs(title = "Average number of kernel matches vs average Euclidean distance",
         x = "Average Euclidean distance",
         y = "Average number of kernel matches",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_average_kernel_matches_vs_average_euclidean_distance.png"), width = 6, height = 4)

# Average number of kernel matches versus Spearman correlation coefficient
ggplot(all.data.master, aes(x = avg.corr.per.receptor.pair, 
    y = avg.kernels.per.receptor.pair, color = type)) +
    geom_point() +
    labs(title = "Average number of kernel matches vs Spearman correlation coefficient",
         x = "Average Spearman correlation",
         y = "Average number of kernel matches",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_average_kernel_matches_vs_average_spearman_correlation.png"), width = 6, height = 4)


# Sankey diagram 
library(networkD3)
library(htmlwidgets)

# Sankey diagram for ref.epitope to sample.epitope, with magnitude of kernel matches
# This does not include decoy receptor pairs
kernel.counts <- all.data.master %>% 
    filter(type != "Decoy receptor pairs") %>%
    filter(samp.epitope != "Combined") %>%
    select(ref.epitope, samp.epitope, total.kernels)

nodes <- data.frame(name = c(unique(kernel.counts$ref.epitope), unique(kernel.counts$samp.epitope)))

links <- data.frame(target = match(kernel.counts$samp.epitope, nodes$name) + length(unique(receptor.counts$ref.epitope)) - 1,
                    source = match(kernel.counts$ref.epitope, nodes$name) - 1,
                    value = kernel.counts$total.kernels,
                    stringsAsFactors = FALSE)

nodes <- nodes %>%
    mutate(source = 0:(nrow(nodes) - 1))

links <- links %>%
    left_join(nodes, join_by(source))

sankeyplot <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source", 
    Target = "target", Value = "value", NodeID = "name", sinksRight = FALSE, LinkGroup = "name")

(sankeyplot <- prependContent(sankeyplot, htmltools::tags$h3("Total kernel matches per epitope pair")))

saveWidget(sankeyplot, file = paste0(out.path, "/sankey_plot_kernel_matches.html"), selfcontained = TRUE)

# Sankey diagram for ref.epitope to sample.epitope, with magnitude of receptor pairs
receptor.counts <- all.data.master %>% 
    filter(type != "Decoy receptor pairs") %>%
    filter(samp.epitope != "Combined") %>%
    select(ref.epitope, samp.epitope, total.receptor.pairs)

nodes <- data.frame(name = c(unique(receptor.counts$ref.epitope), unique(receptor.counts$samp.epitope)))

links <- data.frame(target = match(receptor.counts$samp.epitope, nodes$name) + length(unique(receptor.counts$ref.epitope)) - 1,
                    source = match(receptor.counts$ref.epitope, nodes$name) - 1,
                    value = receptor.counts$total.receptor.pairs,
                    stringsAsFactors = FALSE)

nodes <- nodes %>%
    mutate(source = 0:(nrow(nodes) - 1))

links <- links %>%
    left_join(nodes, join_by(source))

sankeyplot <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source",
    Target = "target", Value = "value", NodeID = "name", sinksRight = FALSE, LinkGroup = "name")

(sankeyplot <- prependContent(sankeyplot, htmltools::tags$h3("Total receptor pairs per epitope pair")))

saveWidget(sankeyplot, file = paste0(out.path, "/sankey_plot_receptor_pairs.html"), selfcontained = TRUE)


# PCA of all features (top KAS, average KAS, Spearman correlation, Euclidean distance, subgraph shape similarity, average number of kernels per pairing)

library(ggfortify)
library(factoextra)

# Select data
pca.data <- all.data.master %>%
    filter(samp.epitope != "Combined") %>%
    select(ref.epitope, samp.epitope, type, avg.top.kdist.per.receptor.pair, avg.kdist.per.receptor.pair, avg.corr.per.receptor.pair, avg.distance.per.receptor.pair, avg.subgraph.shape.sim, avg.kernels.per.receptor.pair)

# Perform PCA
pca.result <- prcomp(pca.data %>% select(-ref.epitope, -samp.epitope, -type), scale. = TRUE)

# Plot PCA
autoplot(pca.result, data = pca.data, colour = 'type') +
    labs(title = "PCA of all features",
            x = "Principal Component 1",
            y = "Principal Component 2",
            color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/pca_plot.png"), width = 6, height = 4)

autoplot(pca.result, data = pca.data, colour = 'type', loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3) +
    labs(title = "PCA of all features, with eigenvectors",
            x = "Principal Component 1",
            y = "Principal Component 2",
            color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/pca_plot_with_eigenvectors.png"), width = 6, height = 4)

fviz_pca_ind(pca.result, col.ind = pca.data$type, geom = "point", addEllipses=TRUE, ellipse.level=0.95) +
    labs(title = "PCA of all features, individual factor map",
            x = "Principal Component 1",
            y = "Principal Component 2",
            color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/pca_plot_individual_factor_map.png"), width = 6, height = 4)
