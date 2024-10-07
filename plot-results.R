# ============================================================
# File: process-results.R
# Description: Script for processing and analyzing PIEMA benchmark results
# Author: Maia Bennett-Boehm
# University of Nebraska at Omaha
# ============================================================


# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(networkD3)
library(htmlwidgets)
library(ggfortify)
library(factoextra)

# Set working directory
# out.path <- "./analysis/apbs/run1"
# out.path <- "./analysis/easymifs/run1/OP"
out.path <- "./analysis/apbs/run2"

# Load data
all.data <- read.csv(paste0(out.path, "/final_all_kernel_data.csv"))

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

all.data.minus.nlv <- all.data.faceted %>%
    filter(ref.epitope != "NLVPMVATV") %>%
    filter(samp.epitope != "NLVPMVATV")

all.receptor.data <- read.csv(paste0(out.path, "/final_receptor_data.csv"))

# Faceted data for epitope plots
all.receptor.data.faceted <- all.receptor.data %>%
    mutate(epitope = ref.epitope) %>%
    bind_rows(
        all.receptor.data %>%
            filter(type == "Unlikely receptor pairs") %>%
            mutate(epitope = samp.epitope)
    )

all.receptor.data.minus.nlv <- all.receptor.data.faceted %>%
    filter(ref.epitope != "NLVPMVATV") %>%
    filter(samp.epitope != "NLVPMVATV")

all.data.master <- read.csv(paste0(out.path, "/results_master.csv"))

all.data.master.minus.nlv <- all.data.master %>%
    filter(ref.epitope != "NLVPMVATV") %>%
    filter(samp.epitope != "NLVPMVATV")

# All plots
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

# All minus NLV
ggplot(all.data.minus.nlv, aes(x = kdist, color = type)) +
    geom_density() +
    labs(title = "Density plot of KAS values, not including NLVPMVATV epitope",
                x = "KAS",
                y = "Density",
                color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/density_plot_kas_values_minus_nlv.png"), width = 6, height = 4)

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

# Minus NLV
ggplot(all.data.minus.nlv, aes(x = kdist, color = type)) +
    geom_density() +
    labs(title = "Density plot of KAS values by epitope",
         x = "KAS",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/density_plot_kas_values_per_epitope_minus_nlv.png"), width = 6, height = 4)

# Spearman correlation for all kernel matches
ggplot(all.data, aes(x = corr, color = type)) +
    geom_density() +
    labs(title = "Spearman correlation for all kernel matches",
         x = "Spearman correlation",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    xlim(0.5, 1)

ggsave(paste0(out.path, "/density_plot_spearman_correlation.png"), width = 6, height = 4)

# Spearman correlation for all kernel matches, faceted by epitope
ggplot(all.data.faceted, aes(x = corr, color = type)) +
    geom_density() +
    labs(title = "Spearman correlation for all kernel matches by epitope, not including cross-reactive receptors",
         x = "Spearman correlation",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    xlim(0.5, 1) + 
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

# Minus NLV
ggplot(all.data.minus.nlv, aes(x = patDist, color = type)) +
    geom_density() +
    labs(title = "Subgraph Euclidean distance for all kernel matches",
         x = "Subgraph Euclidean distance",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/density_plot_subgraph_euclidean_distance_minus_nlv.png"), width = 6, height = 4)

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

# Minus NLV
ggplot(all.receptor.data.minus.nlv, aes(x = top.kdist, color = type)) +
    geom_density() +
    labs(title = "Density plot of top KAS values for each receptor:receptor pair, not including NLVPMVATV epitope",
         x = "KAS",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/density_plot_top_kas_values_minus_nlv.png"), width = 6, height = 4)

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

# Minus NLV
ggplot(all.receptor.data.minus.nlv, aes(x = count, color = type)) +
    geom_density() +
    labs(title = "Kernel match counts per receptor:receptor pair",
         x = "Count",
         y = "Density",
         fill = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/density_plot_kernel_match_counts_minus_nlv.png"), width = 6, height = 4)

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

# CDR3 sequence similarity
ggplot(all.receptor.data, aes(x = CDR3.similarity, color = type)) +
    geom_density() +
    labs(title = "CDR3 sequence similarity for each receptor:receptor pair",
         x = "CDR3 sequence similarity",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/density_plot_cdr3_sequence_similarity.png"), width = 6, height = 4)

# Minus NLV
ggplot(all.receptor.data.minus.nlv, aes(x = CDR3.similarity, color = type)) +
    geom_density() +
    labs(title = "CDR3 sequence similarity for each receptor:receptor pair",
         x = "CDR3 sequence similarity",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/density_plot_cdr3_sequence_similarity_minus_nlv.png"), width = 6, height = 4)

# CDR3 sequence similarity, faceted by epitope- shows difference in per-epitope performance is not necessarily just a result of CDR3 similarity, but is definitely related (same trends but not as strong)
ggplot(all.receptor.data.faceted, aes(x = CDR3.similarity, color = type)) +
    geom_density() +
    labs(title = "CDR3 sequence similarity for each receptor:receptor pair, faceted by epitope",
         x = "CDR3 sequence similarity",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "/density_plot_cdr3_sequence_similarity_per_epitope.png"), width = 6, height = 4)

# CDR sequence similarity
ggplot(all.receptor.data, aes(x = CDR.similarity, color = type)) +
    geom_density() +
    labs(title = "CDR sequence similarity for each receptor:receptor pair",
         x = "CDR sequence similarity",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/density_plot_cdr_sequence_similarity.png"), width = 6, height = 4)

# Minus NLV
ggplot(all.receptor.data.minus.nlv, aes(x = CDR.similarity, color = type)) +
    geom_density() +
    labs(title = "CDR sequence similarity for each receptor:receptor pair",
         x = "CDR sequence similarity",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/density_plot_cdr_sequence_similarity_minus_nlv.png"), width = 6, height = 4)

# CDR sequence similarity, faceted by epitope
ggplot(all.receptor.data.faceted, aes(x = CDR.similarity, color = type)) +
    geom_density() +
    labs(title = "CDR sequence similarity for each receptor:receptor pair, faceted by epitope",
         x = "CDR sequence similarity",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "/density_plot_cdr_sequence_similarity_per_epitope.png"), width = 6, height = 4)

# Full sequence similarity
ggplot(all.receptor.data, aes(x = full.similarity, color = type)) +
    geom_density() +
    labs(title = "Full sequence similarity for each receptor:receptor pair",
         x = "Full sequence similarity",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/density_plot_full_sequence_similarity.png"), width = 6, height = 4)

# Minus NLV
ggplot(all.receptor.data.minus.nlv, aes(x = full.similarity, color = type)) +
    geom_density() +
    labs(title = "Full sequence similarity for each receptor:receptor pair",
         x = "Full sequence similarity",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/density_plot_full_sequence_similarity_minus_nlv.png"), width = 6, height = 4)

# Full sequence similarity, faceted by epitope
ggplot(all.receptor.data.faceted, aes(x = full.similarity, color = type)) +
    geom_density() +
    labs(title = "Full sequence similarity for each receptor:receptor pair, faceted by epitope",
         x = "Full sequence similarity",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "/density_plot_full_sequence_similarity_per_epitope.png"), width = 6, height = 4)

# Average KAS per receptor
ggplot(all.receptor.data, aes(x = avg.kdist, color = type)) +
    geom_density() +
    labs(title = "Average KAS per receptor:receptor pair",
         x = "Average KAS",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/density_plot_average_kas_values.png"), width = 6, height = 4)

# Minus NLV
ggplot(all.receptor.data.minus.nlv, aes(x = avg.kdist, color = type)) +
    geom_density() +
    labs(title = "Average KAS per receptor:receptor pair",
         x = "Average KAS",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/density_plot_average_kas_values_minus_nlv.png"), width = 6, height = 4)

# Average KAS per receptor, faceted by epitope
ggplot(all.receptor.data.faceted, aes(x = avg.kdist, color = type)) +
    geom_density() +
    labs(title = "Average KAS per receptor:receptor pair, faceted by epitope",
         x = "Average KAS",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "/density_plot_average_kas_values_per_epitope.png"), width = 6, height = 4)

# Median subgraph KAS
ggplot(all.receptor.data, aes(x = top.med_scr, color = type)) +
    geom_density() +
    labs(title = "Median subgraph KAS per receptor:receptor pair",
         x = "Median subgraph KAS",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/density_plot_median_subgraph_kas_values.png"), width = 6, height = 4)

# Minus NLV
ggplot(all.receptor.data.minus.nlv, aes(x = top.med_scr, color = type)) +
    geom_density() +
    labs(title = "Median subgraph KAS per receptor:receptor pair",
         x = "Median subgraph KAS",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/density_plot_median_subgraph_kas_values_minus_nlv.png"), width = 6, height = 4)

# Median subgraph KAS, faceted by epitope
ggplot(all.receptor.data.faceted, aes(x = top.med_scr, color = type)) +
    geom_density() +
    labs(title = "Median subgraph KAS per receptor:receptor pair, faceted by epitope",
         x = "Median subgraph KAS",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "/density_plot_median_subgraph_kas_values_per_epitope.png"), width = 6, height = 4)

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

ggsave(paste0(out.path, "/point_plot_top_kas_values.png"), width = 12, height = 8)

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

ggsave(paste0(out.path, "/point_plot_average_kas_values.png"), width = 12, height = 8)

# Median subgraph KAS value per reference receptor epitope (x axis) colored by sample receptor epitope
ggplot(all.data.master, aes(x = ref.epitope, y = avg.median.kdist.per.receptor.pair, color = match.epitope)) +
    geom_point(aes(size = ifelse(match.epitope %in% c("True positive", "Decoy", "Combined"), 2, 1))) +
    labs(title = "Median subgraph KAS value per reference receptor epitope, colored by sample receptor epitope",
         x = "Reference epitope",
         y = "Median subgraph KAS value",
         color = "Paired receptor epitope") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_size_continuous(range = c(2, 4), guide = "none")

ggsave(paste0(out.path, "/point_plot_median_subgraph_kas_values.png"), width = 12, height = 8)

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

ggsave(paste0(out.path, "/point_plot_average_spearman_correlation.png"), width = 12, height = 8)

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

ggsave(paste0(out.path, "/point_plot_average_euclidean_distance.png"), width = 12, height = 8)

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

ggsave(paste0(out.path, "/point_plot_average_subgraph_euclidean_distance.png"), width = 12, height = 8)

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

ggsave(paste0(out.path, "/point_plot_average_subgraph_shape_similarity.png"), width = 12, height = 8)

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

ggsave(paste0(out.path, "/point_plot_average_kernel_matches_per_receptor_pair.png"), width = 12, height = 8)

# Average CDR3 sequence similarity per reference receptor epitope (x axis) colored by sample receptor epitope
# So interesting, shows that sequence similarity might work well for some epitopes (YLQ, RPI) but not others (GLC), supporting need for combination of approaches
ggplot(all.data.master, aes(x = ref.epitope, y = avg.CDR3.sequence.sim, color = match.epitope)) +
    geom_point(aes(size = ifelse(match.epitope %in% c("True positive", "Decoy", "Combined"), 2, 1))) +
    labs(title = "Average CDR3 sequence similarity per reference receptor epitope, colored by sample receptor epitope",
         x = "Reference epitope",
         y = "Average CDR3 sequence similarity",
         color = "Paired receptor epitope") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_size_continuous(range = c(2, 4), guide = "none")

ggsave(paste0(out.path, "/point_plot_average_cdr3_sequence_similarity.png"), width = 12, height = 8)

# Average CDR sequence similarity per reference receptor epitope (x axis) colored by sample receptor epitope
# Combined CDR sequence similarity works better, but not for all epitopes and could be an artifact of V/J gene preservation as well
ggplot(all.data.master, aes(x = ref.epitope, y = avg.CDR.sequence.sim, color = match.epitope)) +
    geom_point(aes(size = ifelse(match.epitope %in% c("True positive", "Decoy", "Combined"), 2, 1))) +
    labs(title = "Average CDR sequence similarity per reference receptor epitope, colored by sample receptor epitope",
         x = "Reference epitope",
         y = "Average CDR sequence similarity",
         color = "Paired receptor epitope") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_size_continuous(range = c(2, 4), guide = "none")

ggsave(paste0(out.path, "/point_plot_average_cdr_sequence_similarity.png"), width = 12, height = 8)

# Average full sequence similarity per reference receptor epitope (x axis) colored by sample receptor epitope
# Full sequence similarity works for all but YLQ
ggplot(all.data.master, aes(x = ref.epitope, y = avg.full.sequence.sim, color = match.epitope)) +
    geom_point(aes(size = ifelse(match.epitope %in% c("True positive", "Decoy", "Combined"), 2, 1))) +
    labs(title = "Average full sequence similarity per reference receptor epitope, colored by sample receptor epitope",
         x = "Reference epitope",
         y = "Average full sequence similarity",
         color = "Paired receptor epitope") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_size_continuous(range = c(2, 4), guide = "none")

ggsave(paste0(out.path, "/point_plot_average_full_sequence_similarity.png"), width = 12, height = 8)


# Dot plots
# Subgraph shape similarity vs subgraph Euclidean distance, colored by receptor pair type
ggplot(all.data.master, aes(x = avg.subgraph.distance, y = avg.subgraph.shape.sim, color = type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.subgraph.shape.sim > 0.99 & avg.subgraph.distance < 3.75 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.subgraph.shape.sim > 0.99 & avg.subgraph.distance < 3.75 & type == "Decoy receptor pairs")), 
              aes(label = paste0(ref.epitope)),
              vjust = -1, 
              size = 3) +
    labs(title = "Subgraph shape similarity vs subgraph Euclidean distance",
         x = "Average subgraph Euclidean distance",
         y = "Average subgraph shape similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_subgraph_shape_similarity_vs_euclidean_distance.png"), width = 12, height = 8)

# Top KAS vs average KAS, colored by receptor pair type
ggplot(all.data.master, aes(x = avg.kdist.per.receptor.pair, y = avg.top.kdist.per.receptor.pair, color = type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.top.kdist.per.receptor.pair > 0.2 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.top.kdist.per.receptor.pair > 0.2 & type == "Decoy receptor pairs")), 
              aes(label = paste0(ref.epitope)),
              vjust = -1, 
              size = 3) +
    labs(title = "Top KAS vs average KAS",
         x = "Average KAS",
         y = "Top KAS",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_top_kas_vs_average_kas.png"), width = 12, height = 8)

# Average Spearman vs average KAS, colored by receptor pair type
ggplot(all.data.master, aes(x = avg.kdist.per.receptor.pair, y = avg.corr.per.receptor.pair, color = type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.corr.per.receptor.pair > 0.7325 & avg.kdist.per.receptor.pair > -0.075 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.corr.per.receptor.pair > 0.7325 & avg.kdist.per.receptor.pair > -0.075 & type == "Decoy receptor pairs")), 
              aes(label = paste0(ref.epitope)),
              vjust = -1, 
              size = 3) +
    labs(title = "Average Spearman correlation vs average KAS",
         x = "Average KAS",
         y = "Average Spearman correlation",
            color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_average_spearman_correlation_vs_average_kas.png"), width = 12, height = 8)

# Average Euclidean distance vs average KAS, colored by receptor pair type
ggplot(all.data.master, aes(x = avg.kdist.per.receptor.pair, y = avg.distance.per.receptor.pair, color = type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.distance.per.receptor.pair < 6.25 & type == "Unlikely receptor pairs")),
                aes(label = paste0(ref.epitope, "+", samp.epitope)),
                vjust = -1, 
                size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.distance.per.receptor.pair < 6.25 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Average Euclidean distance vs average KAS",
         x = "Average KAS",
         y = "Average Euclidean distance",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_average_euclidean_distance_vs_average_kas.png"), width = 12, height = 8)

# Average subgraph shape similarity vs top KAS, colored by receptor pair type
ggplot(all.data.master, aes(x = avg.top.kdist.per.receptor.pair, y = avg.subgraph.shape.sim, color = type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.subgraph.shape.sim > 0.99 & avg.top.kdist.per.receptor.pair > 0.175 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.subgraph.shape.sim > 0.99 & avg.top.kdist.per.receptor.pair > 0.175 & type == "Decoy receptor pairs")), 
              aes(label = paste0(ref.epitope)),
              vjust = -1, 
              size = 3) +
    labs(title = "Top KAS vs average subgraph shape similarity",
         x = "Top KAS",
         y = "Average subgraph shape similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_top_kas_vs_average_subgraph_shape_similarity.png"), width = 12, height = 8)

# Average subgraph shape similarity vs average KAS, colored by receptor pair type
ggplot(all.data.master, aes(x = avg.kdist.per.receptor.pair, y = avg.subgraph.shape.sim, color = type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.subgraph.shape.sim > 0.99 & avg.kdist.per.receptor.pair > -0.075 & type == "Unlikely receptor pairs")),
                aes(label = paste0(ref.epitope, "+", samp.epitope)),
                vjust = -1, 
                size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.subgraph.shape.sim > 0.99 & avg.kdist.per.receptor.pair > -0.08 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Average KAS vs average subgraph shape similarity",
         x = "Average KAS",
         y = "Average subgraph shape similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_average_kas_vs_average_subgraph_shape_similarity.png"), width = 12, height = 8)

# Average shape similarity vs average number of kernel matches, colored by receptor pair type
ggplot(all.data.master, aes(x = avg.kernels.per.receptor.pair, y = avg.subgraph.shape.sim, color = type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.subgraph.shape.sim > 0.99 & avg.kernels.per.receptor.pair > 6.5 & type == "Unlikely receptor pairs")),
                aes(label = paste0(ref.epitope, "+", samp.epitope)),
                vjust = -1, 
                size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.subgraph.shape.sim > 0.99 & avg.kernels.per.receptor.pair > 6.5 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Average number of kernel matches vs average subgraph shape similarity",
         x = "Average number of kernel matches",
         y = "Average subgraph shape similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_average_kernel_matches_vs_average_subgraph_shape_similarity.png"), width = 12, height = 8)

# Average number of kernel matches versus average Euclidean distance
ggplot(all.data.master, aes(x = avg.distance.per.receptor.pair, y = avg.kernels.per.receptor.pair, color = type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.kernels.per.receptor.pair > 7 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.kernels.per.receptor.pair > 7 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Average number of kernel matches vs average Euclidean distance",
         x = "Average Euclidean distance",
         y = "Average number of kernel matches",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_average_kernel_matches_vs_average_euclidean_distance.png"), width = 12, height = 8)

# Average number of kernel matches versus Spearman correlation coefficient
ggplot(all.data.master, aes(x = avg.corr.per.receptor.pair, 
    y = avg.kernels.per.receptor.pair, color = type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.kernels.per.receptor.pair > 7 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.kernels.per.receptor.pair > 7 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Average number of kernel matches vs Spearman correlation coefficient",
         x = "Average Spearman correlation",
         y = "Average number of kernel matches",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_average_kernel_matches_vs_average_spearman_correlation.png"), width = 12, height = 8)

# Average number of kernel matches versus top KAS
ggplot(all.data.master, aes(x = avg.top.kdist.per.receptor.pair, 
    y = avg.kernels.per.receptor.pair, color = type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.kernels.per.receptor.pair > 7 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.kernels.per.receptor.pair > 7 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Average number of kernel matches vs top KAS",
         x = "Top KAS",
         y = "Average number of kernel matches",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_average_kernel_matches_vs_top_kas.png"), width = 12, height = 8)

# Top KAS versus Spearman correlation
ggplot(all.data.master, aes(x = avg.corr.per.receptor.pair, 
    y = avg.top.kdist.per.receptor.pair, color = type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.top.kdist.per.receptor.pair > 0.2 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.top.kdist.per.receptor.pair > 0.2 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Top KAS vs Spearman correlation coefficient",
         x = "Average Spearman correlation",
         y = "Top KAS",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_top_kas_vs_spearman_correlation.png"), width = 12, height = 8)

# CDR3 sequence similarity versus Spearman correlation
ggplot(all.data.master, aes(x = avg.corr.per.receptor.pair, 
    y = avg.CDR3.sequence.sim, color = type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.CDR3.sequence.sim > 0.375 & avg.corr.per.receptor.pair > 0.73 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.CDR3.sequence.sim > 0.375 & avg.corr.per.receptor.pair > 0.73 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "CDR3 sequence similarity vs Spearman correlation coefficient",
         x = "Average Spearman correlation",
         y = "Average CDR3 sequence similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_cdr3_sequence_similarity_vs_spearman_correlation.png"), width = 12, height = 8)

# CDR3 sequence similarity versus top KAS
ggplot(all.data.master, aes(x = avg.top.kdist.per.receptor.pair, 
    y = avg.CDR3.sequence.sim, color = type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.CDR3.sequence.sim > 0.375 & avg.top.kdist.per.receptor.pair > 0.2 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.CDR3.sequence.sim > 0.375 & avg.top.kdist.per.receptor.pair > 0.2 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "CDR3 sequence similarity vs top KAS",
         x = "Top KAS",
         y = "Average CDR3 sequence similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_cdr3_sequence_similarity_vs_top_kas.png"), width = 12, height = 8)

# CDR3 sequence similarity versus average KAS
ggplot(all.data.master, aes(x = avg.kdist.per.receptor.pair, 
    y = avg.CDR3.sequence.sim, color = type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.CDR3.sequence.sim > 0.375 & avg.kdist.per.receptor.pair > -0.075 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.CDR3.sequence.sim > 0.375 & avg.kdist.per.receptor.pair > -0.075 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "CDR3 sequence similarity vs average KAS",
         x = "Average KAS",
         y = "Average CDR3 sequence similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_cdr3_sequence_similarity_vs_average_kas.png"), width = 12, height = 8)

# CDR3 sequence similarity versus Euclidean distance
ggplot(all.data.master, aes(x = avg.distance.per.receptor.pair, 
    y = avg.CDR3.sequence.sim, color = type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.CDR3.sequence.sim > 0.375 & avg.distance.per.receptor.pair < 6.25 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.CDR3.sequence.sim > 0.375 & avg.distance.per.receptor.pair < 6.25 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "CDR3 sequence similarity vs average Euclidean distance",
         x = "Average Euclidean distance",
         y = "Average CDR3 sequence similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_cdr3_sequence_similarity_vs_average_euclidean_distance.png"), width = 12, height = 8)

# CDR3 sequence similarity versus average number of kernel matches
ggplot(all.data.master, aes(x = avg.kernels.per.receptor.pair, 
    y = avg.CDR3.sequence.sim, color = type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.CDR3.sequence.sim > 0.375 & avg.kernels.per.receptor.pair > 7 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.CDR3.sequence.sim > 0.375 & avg.kernels.per.receptor.pair > 7 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "CDR3 sequence similarity vs average number of kernel matches",
         x = "Average number of kernel matches",
         y = "Average CDR3 sequence similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_cdr3_sequence_similarity_vs_average_kernel_matches.png"), width = 12, height = 8)

# Full sequence similarity versus Spearman correlation
ggplot(all.data.master, aes(x = avg.corr.per.receptor.pair, 
    y = avg.full.sequence.sim, color = type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.full.sequence.sim > 0.425 & avg.corr.per.receptor.pair > 0.73 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.full.sequence.sim > 0.425 & avg.corr.per.receptor.pair > 0.73 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Full sequence similarity vs Spearman correlation coefficient",
         x = "Average Spearman correlation",
         y = "Average full sequence similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_full_sequence_similarity_vs_spearman_correlation.png"), width = 12, height = 8)

# Full sequence similarity versus top KAS
ggplot(all.data.master, aes(x = avg.top.kdist.per.receptor.pair, 
    y = avg.full.sequence.sim, color = type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.full.sequence.sim > 0.425 & avg.top.kdist.per.receptor.pair > 0.2 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.full.sequence.sim > 0.425 & avg.top.kdist.per.receptor.pair > 0.2 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Full sequence similarity vs top KAS",
         x = "Top KAS",
         y = "Average full sequence similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_full_sequence_similarity_vs_top_kas.png"), width = 12, height = 8)

# Full sequence similarity versus average KAS
ggplot(all.data.master, aes(x = avg.kdist.per.receptor.pair, 
    y = avg.full.sequence.sim, color = type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.full.sequence.sim > 0.425 & avg.kdist.per.receptor.pair > -0.075 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.full.sequence.sim > 0.425 & avg.kdist.per.receptor.pair > -0.075 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Full sequence similarity vs average KAS",
         x = "Average KAS",
         y = "Average full sequence similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_full_sequence_similarity_vs_average_kas.png"), width = 12, height = 8)

# Full sequence similarity versus Euclidean distance
ggplot(all.data.master, aes(x = avg.distance.per.receptor.pair, 
    y = avg.full.sequence.sim, color = type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.full.sequence.sim > 0.425 & avg.distance.per.receptor.pair < 6.25 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.full.sequence.sim > 0.425 & avg.distance.per.receptor.pair < 6.25 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Full sequence similarity vs average Euclidean distance",
         x = "Average Euclidean distance",
         y = "Average full sequence similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/dot_plot_full_sequence_similarity_vs_average_euclidean_distance.png"), width = 12, height = 8)

# Dot plots of all data (not averages)
# Top KAS vs Spearman
ggplot(all.receptor.data, aes(x = top.kdist, y = top.corr, color = type)) +
    geom_point(alpha = 0.5) +
    labs(title = "Top KAS vs Spearman correlation coefficient",
         x = "Top KAS",
         y = "Spearman correlation",
         color = "Receptor:receptor pair type") +
    ylim(0.6, 1) +
    theme_minimal()

ggsave(paste0(out.path, "/all_data_dot_plot_top_kas_vs_spearman_correlation.png"), width = 12, height = 8)

# Top KAS vs Euclidean distance
ggplot(all.receptor.data, aes(x = top.kdist, y = top.distance, color = type)) +
    geom_point(alpha = 0.5) +
    labs(title = "Top KAS vs Euclidean distance",
         x = "Top KAS",
         y = "Euclidean distance",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/all_data_dot_plot_top_kas_vs_euclidean_distance.png"), width = 12, height = 8)

# Top KAS vs subgraph shape similarity
ggplot(all.receptor.data, aes(x = top.kdist, y = top.shapSim, color = type)) +
    geom_point(alpha = 0.5) +
    labs(title = "Top KAS vs subgraph shape similarity",
         x = "Top KAS",
         y = "Subgraph shape similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/all_data_dot_plot_top_kas_vs_subgraph_shape_similarity.png"), width = 12, height = 8)

# Top KAS vs number of kernels
ggplot(all.receptor.data, aes(x = top.kdist, y = count, color = type)) +
    geom_point(alpha = 0.5) +
    labs(title = "Top KAS vs number of kernels",
         x = "Top KAS",
         y = "Number of kernels",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/all_data_dot_plot_top_kas_vs_number_of_kernels.png"), width = 12, height = 8)

# CDR3 sequence similarity vs Spearman correlation
ggplot(all.receptor.data, aes(x = CDR3.similarity, y = top.corr, color = type)) +
    geom_point(alpha = 0.5) +
    labs(title = "CDR3 sequence similarity vs Spearman correlation coefficient",
         x = "CDR3 sequence similarity",
         y = "Spearman correlation",
         color = "Receptor:receptor pair type") +
    ylim(0.6, 1) +
    theme_minimal()

ggsave(paste0(out.path, "/all_data_dot_plot_cdr3_sequence_similarity_vs_spearman_correlation.png"), width = 12, height = 8)

# CDR3 sequence similarity vs top KAS
ggplot(all.receptor.data, aes(x = CDR3.similarity, y = top.kdist, color = type)) +
    geom_point(alpha = 0.5) +
    labs(title = "CDR3 sequence similarity vs top KAS",
         x = "CDR3 sequence similarity",
         y = "Top KAS",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/all_data_dot_plot_cdr3_sequence_similarity_vs_top_kas.png"), width = 12, height = 8)

# CDR3 sequence similarity vs Euclidean distance
ggplot(all.receptor.data, aes(x = CDR3.similarity, y = top.distance, color = type)) +
    geom_point(alpha = 0.5) +
    labs(title = "CDR3 sequence similarity vs Euclidean distance",
         x = "CDR3 sequence similarity",
         y = "Euclidean distance",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/all_data_dot_plot_cdr3_sequence_similarity_vs_euclidean_distance.png"), width = 12, height = 8)

# CDR3 sequence similarity vs number of kernels
ggplot(all.receptor.data, aes(x = CDR3.similarity, y = count, color = type)) +
    geom_point(alpha = 0.5) +
    labs(title = "CDR3 sequence similarity vs number of kernels",
         x = "CDR3 sequence similarity",
         y = "Number of kernels",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/all_data_dot_plot_cdr3_sequence_similarity_vs_number_of_kernels.png"), width = 12, height = 8)

# Full sequence similarity vs Spearman correlation 
ggplot(all.receptor.data, aes(x = full.similarity, y = top.corr, color = type)) +
    geom_point(alpha = 0.5) +
    labs(title = "Full sequence similarity vs Spearman correlation coefficient",
         x = "Full sequence similarity",
         y = "Spearman correlation",
         color = "Receptor:receptor pair type") +
    ylim(0.6, 1) +
    theme_minimal()

ggsave(paste0(out.path, "/all_data_dot_plot_full_sequence_similarity_vs_spearman_correlation.png"), width = 12, height = 8)

# Full sequence similarity vs top KAS
ggplot(all.receptor.data, aes(x = full.similarity, y = top.kdist, color = type)) +
    geom_point(alpha = 0.5) +
    labs(title = "Full sequence similarity vs top KAS",
         x = "Full sequence similarity",
         y = "Top KAS",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/all_data_dot_plot_full_sequence_similarity_vs_top_kas.png"), width = 12, height = 8)

# Full sequence similarity vs Euclidean distance
ggplot(all.receptor.data, aes(x = full.similarity, y = top.distance, color = type)) +
    geom_point(alpha = 0.5) +
    labs(title = "Full sequence similarity vs Euclidean distance",
         x = "Full sequence similarity",
         y = "Euclidean distance",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/all_data_dot_plot_full_sequence_similarity_vs_euclidean_distance.png"), width = 12, height = 8)

# Full sequence similarity vs number of kernels
ggplot(all.receptor.data, aes(x = full.similarity, y = count, color = type)) +
    geom_point(alpha = 0.5) +
    labs(title = "Full sequence similarity vs number of kernels",
         x = "Full sequence similarity",
         y = "Number of kernels",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/all_data_dot_plot_full_sequence_similarity_vs_number_of_kernels.png"), width = 12, height = 8)


# Sankey diagrams
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
# Select data
pca.data <- all.data.master %>%
    filter(samp.epitope != "Combined") %>%
    select(ref.epitope, samp.epitope, type, avg.top.kdist.per.receptor.pair, avg.kdist.per.receptor.pair, avg.corr.per.receptor.pair, avg.distance.per.receptor.pair, avg.subgraph.shape.sim, avg.kernels.per.receptor.pair, avg.CDR3.sequence.sim, avg.full.sequence.sim)

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

# PCA of all features (top KAS, average KAS, Spearman correlation, Euclidean distance, subgraph shape similarity, average number of kernels per pairing) without NLV data
# Select data
pca.data <- all.data.master.minus.nlv %>%
    filter(samp.epitope != "Combined") %>%
    select(ref.epitope, samp.epitope, type, avg.top.kdist.per.receptor.pair, avg.kdist.per.receptor.pair, avg.corr.per.receptor.pair, avg.distance.per.receptor.pair, avg.subgraph.shape.sim, avg.kernels.per.receptor.pair, avg.CDR3.sequence.sim, avg.full.sequence.sim)

# Perform PCA
pca.result <- prcomp(pca.data %>% select(-ref.epitope, -samp.epitope, -type), scale. = TRUE)

# Plot PCA
autoplot(pca.result, data = pca.data, colour = 'type') +
    labs(title = "PCA of all features",
            x = "Principal Component 1",
            y = "Principal Component 2",
            color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/pca_plot_minus_NLV.png"), width = 6, height = 4)

autoplot(pca.result, data = pca.data, colour = 'type', loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3) +
    labs(title = "PCA of all features, with eigenvectors",
            x = "Principal Component 1",
            y = "Principal Component 2",
            color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/pca_plot_with_eigenvectors_minus_NLV.png"), width = 6, height = 4)

fviz_pca_ind(pca.result, col.ind = pca.data$type, geom = "point", addEllipses=TRUE, ellipse.level=0.95) +
    labs(title = "PCA of all features, individual factor map",
            x = "Principal Component 1",
            y = "Principal Component 2",
            color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/pca_plot_individual_factor_map_minus_NLV.png"), width = 6, height = 4)


# PCA on all data and features without NLV data
# Select data
pca.data <- all.receptor.data.minus.nlv %>% 
    select(ref.epitope, samp.epitope, type, top.kdist, avg.kdist, avg.distance, avg.corr, top.shapSim, count, CDR3.similarity, full.similarity) 

# Perform PCA
pca.result <- prcomp(pca.data %>% select(-ref.epitope, -samp.epitope, -type), scale. = TRUE)

# Plot PCA
autoplot(pca.result, data = pca.data, colour = 'type') +
    labs(title = "PCA of all features",
            x = "Principal Component 1",
            y = "Principal Component 2",
            color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/pca_plot_all_data.png"), width = 6, height = 4)

autoplot(pca.result, data = pca.data, colour = 'type', loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3) +
    labs(title = "PCA of all features, with eigenvectors",
            x = "Principal Component 1",
            y = "Principal Component 2",
            color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/pca_plot_all_data_with_eigenvectors.png"), width = 6, height = 4)

fviz_pca_ind(pca.result, col.ind = pca.data$type, geom = "point", addEllipses=TRUE, ellipse.level=0.95) +
    labs(title = "PCA of all features, individual factor map",
            x = "Principal Component 1",
            y = "Principal Component 2",
            color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "/pca_plot_all_data_individual_factor_map.png"), width = 6, height = 4)