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
library(DataExplorer)

# Set working directory
# out.path <- "./analysis/apbs/without-cross-reactives/CDR3-similarity/all/"
out.path <- "./analysis/apbs/CDR3-similarity/all/"


# Load data
all.data <- read.csv(paste0(out.path, "final_all_kernels_data.csv")) %>%
    dplyr::rename(
        # Metrics averaged across all kernel matches per receptor pair
        euc.dist = distance, kas = kdist, spearman = corr, sg.median.kas = med_scr, 
        # Descriptors
        binding.pair.type = type, pair.id = id, sg.id = sg_group)

# Create a combined data frame with duplicated rows for 'Unlikely binding pairs' to allow facet wrapping
all.data.faceted <- all.data %>%
    mutate(epitope = ref.epitope) %>%
    bind_rows(
        all.data %>%
            filter(binding.pair.type == "Unlikely receptor paired kernels") %>%
            mutate(epitope = samp.epitope)
    )

all.receptor.data <- read.csv(paste0(out.path, "final_receptor_data.csv")) %>%
    dplyr::rename(
        # Metrics specific to top-scoring (KAS) kernel match
        top.kas = top.kdist, top.spearman = top.corr, top.euc.dist = top.distance, 
        # Metrics specific to subgraph with top-scoring (KAS) kernel match
        top.sg.euc.dist = top.patDist, top.sg.cosine.sim = top.shapSim, top.sg.median.kas = top.med_scr, 
        # Metrics averaged across all kernel matches per receptor pair
        avg.euc.dist = avg.distance, avg.kas = avg.kdist, avg.spearman = avg.corr,
        # Kernel count
        kernel.count = count, 
        # Descriptors
        binding.pair.type = type, pair.id = id, top.sg.id = top.sg_group)

# Faceted data for epitope plots
all.receptor.data.faceted <- all.receptor.data %>%
    mutate(epitope = ref.epitope) %>%
    bind_rows(
        all.receptor.data %>%
            filter(binding.pair.type == "Unlikely receptor pairs") %>%
            mutate(epitope = samp.epitope)
    )

all.data.master <- read.csv(paste0(out.path, "results_master.csv"))


# All plots
# Density plots
# KAS scores for all kernel matches
ggplot(all.data, aes(x = kas, color = binding.pair.type)) +
    geom_density() +
    labs(title = "Density plot of KAS values",
                x = "KAS",
                y = "Density",
                color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "density_plot_kas_values.png"), width = 6, height = 4)

# KAS scores for all kernel matches, faceted by epitope
ggplot(all.data.faceted, aes(x = kas, color = binding.pair.type)) +
    geom_density() +
    labs(title = "Density plot of KAS values by epitope, not including cross-reactive receptors",
         x = "KAS",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "density_plot_kas_values_per_epitope.png"), width = 6, height = 4)

# Spearman correlation for all kernel matches
ggplot(all.data, aes(x = spearman, color = binding.pair.type)) +
    geom_density() +
    labs(title = "Spearman correlation for all kernel matches",
         x = "Spearman correlation",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    xlim(0.5, 1)

ggsave(paste0(out.path, "density_plot_spearman_correlation.png"), width = 6, height = 4)

# Spearman correlation for all kernel matches, faceted by epitope
ggplot(all.data.faceted, aes(x = spearman, color = binding.pair.type)) +
    geom_density() +
    labs(title = "Spearman correlation for all kernel matches by epitope, not including cross-reactive receptors",
         x = "Spearman correlation",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    xlim(0.5, 1) + 
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "density_plot_spearman_correlation_per_epitope.png"), width = 6, height = 4)

# Euclidean distance for all kernel matches
ggplot(all.data, aes(x = euc.dist, color = binding.pair.type)) +
    geom_density() +
    labs(title = "Euclidean distance for all kernel matches",
         x = "Euclidean distance",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "density_plot_euclidean_distance.png"), width = 6, height = 4)

# Euclidean distance for all kernel matches, faceted by epitope
ggplot(all.data.faceted, aes(x = euc.dist, color = binding.pair.type)) +
    geom_density() +
    labs(title = "Euclidean distance for all kernel matches by epitope, not including cross-reactive receptors",
         x = "Euclidean distance",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "density_plot_euclidean_distance_per_epitope.png"), width = 6, height = 4)

# Shape similarity (cosine distance between subgraphs)
ggplot(all.data, aes(x = shapSim, color = binding.pair.type)) +
    geom_density() +
    labs(title = "Shape similarity for all kernel matches",
         x = "Shape similarity",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "density_plot_shape_similarity.png"), width = 6, height = 4)

# Shape similarity (cosine distance between subgraphs), faceted by epitope
ggplot(all.data.faceted, aes(x = shapSim, color = binding.pair.type)) +
    geom_density() +
    labs(title = "Shape similarity for all kernel matches by epitope, not including cross-reactive receptors",
         x = "Shape similarity",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "density_plot_shape_similarity_per_epitope.png"), width = 6, height = 4)

# Subgraph Euclidean distance
ggplot(all.data, aes(x = patDist, color = binding.pair.type)) +
    geom_density() +
    labs(title = "Subgraph Euclidean distance for all kernel matches",
         x = "Subgraph Euclidean distance",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "density_plot_subgraph_euclidean_distance.png"), width = 6, height = 4)


# Subgraph Euclidean distance, faceted by epitope
ggplot(all.data.faceted, aes(x = patDist, color = binding.pair.type)) +
    geom_density() +
    labs(title = "Subgraph Euclidean distance for all kernel matches by epitope, not including cross-reactive receptors",
         x = "Subgraph Euclidean distance",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "density_plot_subgraph_euclidean_distance_per_epitope.png"), width = 6, height = 4)

# Top KAS values per receptor:receptor pair
ggplot(all.receptor.data, aes(x = top.kas, color = binding.pair.type)) +
    geom_density() +
    labs(title = "Density plot of top KAS values for each receptor:receptor pair",
         x = "KAS",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "density_plot_top_kas_values.png"), width = 6, height = 4)


# Top KAS values per receptor:receptor pair, faceted by epitope
ggplot(all.receptor.data.faceted, aes(x = top.kas, color = binding.pair.type)) +
    geom_density() +
    labs(title = "Density plot of top KAS values for each receptor:receptor pair, faceted by epitope",
         x = "KAS",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "density_plot_top_kas_values_per_epitope.png"), width = 6, height = 4)

# Number of kernel matches per receptor:receptor pair
ggplot(all.receptor.data, aes(x = kernel.count, color = binding.pair.type)) +
    geom_density() +
    labs(title = "Kernel match counts per receptor:receptor pair",
         x = "Count",
         y = "Density",
         fill = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "density_plot_kernel_match_counts.png"), width = 6, height = 4)



# Number of kernel matches per receptor:receptor pair, faceted by epitope
ggplot(all.receptor.data.faceted, aes(x = kernel.count, color = binding.pair.type)) +
    geom_density() +
    labs(title = "Kernel match counts per receptor:receptor pair, faceted by epitope",
         x = "Count",
         y = "Density",
         fill = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "density_plot_kernel_match_counts_per_epitope.png"), width = 6, height = 4)

# CDR3 sequence similarity
ggplot(all.receptor.data, aes(x = CDR3.similarity, color = binding.pair.type)) +
    geom_density() +
    labs(title = "CDR3 sequence similarity for each receptor:receptor pair",
         x = "CDR3 sequence similarity",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "density_plot_cdr3_sequence_similarity.png"), width = 6, height = 4)


# CDR3 sequence similarity, faceted by epitope- shows difference in per-epitope performance is not necessarily just a result of CDR3 similarity, but is definitely related (same trends but not as strong)
ggplot(all.receptor.data.faceted, aes(x = CDR3.similarity, color = binding.pair.type)) +
    geom_density() +
    labs(title = "CDR3 sequence similarity for each receptor:receptor pair, faceted by epitope",
         x = "CDR3 sequence similarity",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "density_plot_cdr3_sequence_similarity_per_epitope.png"), width = 6, height = 4)

# CDR sequence similarity
ggplot(all.receptor.data, aes(x = CDR.similarity, color = binding.pair.type)) +
    geom_density() +
    labs(title = "CDR sequence similarity for each receptor:receptor pair",
         x = "CDR sequence similarity",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "density_plot_cdr_sequence_similarity.png"), width = 6, height = 4)


# CDR sequence similarity, faceted by epitope
ggplot(all.receptor.data.faceted, aes(x = CDR.similarity, color = binding.pair.type)) +
    geom_density() +
    labs(title = "CDR sequence similarity for each receptor:receptor pair, faceted by epitope",
         x = "CDR sequence similarity",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "density_plot_cdr_sequence_similarity_per_epitope.png"), width = 6, height = 4)

# Full sequence similarity
ggplot(all.receptor.data, aes(x = full.similarity, color = binding.pair.type)) +
    geom_density() +
    labs(title = "Full sequence similarity for each receptor:receptor pair",
         x = "Full sequence similarity",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "density_plot_full_sequence_similarity.png"), width = 6, height = 4)


# Full sequence similarity, faceted by epitope
ggplot(all.receptor.data.faceted, aes(x = full.similarity, color = binding.pair.type)) +
    geom_density() +
    labs(title = "Full sequence similarity for each receptor:receptor pair, faceted by epitope",
         x = "Full sequence similarity",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "density_plot_full_sequence_similarity_per_epitope.png"), width = 6, height = 4)

# Average KAS per receptor
ggplot(all.receptor.data, aes(x = avg.kas, color = binding.pair.type)) +
    geom_density() +
    labs(title = "Average KAS per receptor:receptor pair",
         x = "Average KAS",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "density_plot_average_kas_values.png"), width = 6, height = 4)



# Average KAS per receptor, faceted by epitope
ggplot(all.receptor.data.faceted, aes(x = avg.kas, color = binding.pair.type)) +
    geom_density() +
    labs(title = "Average KAS per receptor:receptor pair, faceted by epitope",
         x = "Average KAS",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "density_plot_average_kas_values_per_epitope.png"), width = 6, height = 4)

# Median subgraph KAS
ggplot(all.receptor.data, aes(x = top.sg.median.kas, color = binding.pair.type)) +
    geom_density() +
    labs(title = "Median subgraph KAS per receptor:receptor pair",
         x = "Median subgraph KAS",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "density_plot_median_subgraph_kas_values.png"), width = 6, height = 4)

# Median subgraph KAS, faceted by epitope
ggplot(all.receptor.data.faceted, aes(x = top.sg.median.kas, color = binding.pair.type)) +
    geom_density() +
    labs(title = "Median subgraph KAS per receptor:receptor pair, faceted by epitope",
         x = "Median subgraph KAS",
         y = "Density",
         color = "Receptor:receptor pair type") +
    theme_minimal() +
    facet_wrap(~ epitope)

ggsave(paste0(out.path, "density_plot_median_subgraph_kas_values_per_epitope.png"), width = 6, height = 4)


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

ggsave(paste0(out.path, "point_plot_top_kas_values.png"), width = 12, height = 8)

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

ggsave(paste0(out.path, "point_plot_average_kas_values.png"), width = 12, height = 8)

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

ggsave(paste0(out.path, "point_plot_median_subgraph_kas_values.png"), width = 12, height = 8)

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

ggsave(paste0(out.path, "point_plot_average_spearman_correlation.png"), width = 12, height = 8)

# Top Spearman correlation
ggplot(all.data.master, aes(x = ref.epitope, y = avg.top.corr.per.receptor.pair, color = match.epitope)) +
    geom_point(aes(size = ifelse(match.epitope %in% c("True positive", "Decoy", "Combined"), 2, 1))) +
    labs(title = "Top Spearman correlation per reference receptor epitope, colored by sample receptor epitope",
         x = "Reference epitope",
         y = "Top Spearman correlation",
            color = "Paired receptor epitope") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_size_continuous(range = c(2, 4), guide = "none")

ggsave(paste0(out.path, "point_plot_top_spearman_correlation.png"), width = 12, height = 8)

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

ggsave(paste0(out.path, "point_plot_average_euclidean_distance.png"), width = 12, height = 8)

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

ggsave(paste0(out.path, "point_plot_average_subgraph_euclidean_distance.png"), width = 12, height = 8)

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

ggsave(paste0(out.path, "point_plot_average_subgraph_shape_similarity.png"), width = 12, height = 8)

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

ggsave(paste0(out.path, "point_plot_average_kernel_matches_per_receptor_pair.png"), width = 12, height = 8)

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

ggsave(paste0(out.path, "point_plot_average_cdr3_sequence_similarity.png"), width = 12, height = 8)

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

ggsave(paste0(out.path, "point_plot_average_cdr_sequence_similarity.png"), width = 12, height = 8)

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

ggsave(paste0(out.path, "point_plot_average_full_sequence_similarity.png"), width = 12, height = 8)


# Box plots
# All receptor data (i.e., the data going into the logistic classifier), per feature, per epitope, colored by receptor pair type
# Make significance tables: Three-way comparisons (i.e., true receptor pairs vs unlikely receptor pairs and true vs decoy receptor pairs)
# All data
features.of.interest <- c("top.euc.dist", "top.kas", "top.spearman", "top.sg.cosine.sim", "top.sg.median.kas", "top.sg.euc.dist", "avg.euc.dist", "avg.kas", "avg.spearman", "kernel.count", "CDR3.similarity", "CDR.similarity", "full.similarity")
comparisons <- list(
    c("True receptor pairs", "Unlikely receptor pairs"),
    c("True receptor pairs", "Decoy receptor pairs")
)
significance.results.three.groups <- data.frame()

for (feature in features.of.interest) {
    significance.results <- all.receptor.data %>% 
        do({
            data <- .
            results <- lapply(comparisons, function(comp) {
                test <- wilcox.test(as.formula(paste(feature, "~ binding.pair.type")), data = data %>% filter(binding.pair.type %in% comp))
                data.frame(
                    Feature = feature,
                    Comparison = paste(comp, collapse = " vs "),
                    P.Value = test$p.value
                )
            })
            bind_rows(results)
        }) %>% 
        ungroup()

    significance.results.three.groups <- bind_rows(significance.results.three.groups, significance.results)
}

write.csv(significance.results.three.groups, paste0(out.path, "significance_results_three_groups.csv"))

# Make significance tables: Two-way comparisons (i.e., true receptor pairs vs unlikely and decoy receptor pairs)
comparisons <- c("True receptor pairs", "Negative receptor pairs")
significance.results.two.groups <- data.frame()

for (feature in features.of.interest) {
    two.way.data <- all.receptor.data %>% 
        mutate(binding.pair.type = ifelse(binding.pair.type %in% c("Unlikely receptor pairs", "Decoy receptor pairs"), "Negative receptor pairs", binding.pair.type)) 
    significance.results <- two.way.data %>% 
        do({
            data <- .
            test <- wilcox.test(as.formula(paste(feature, "~ binding.pair.type")), data = data)
            data.frame(
                Feature = feature,
                Comparison = paste(comparisons, collapse = " vs "),
                P.Value = test$p.value
            )
        }) %>% 
        ungroup()

    significance.results.two.groups <- bind_rows(significance.results.two.groups, significance.results)
}

write.csv(significance.results.two.groups, paste0(out.path, "significance_results_two_groups.csv"))


# Per epitope
comparisons <- list(
    c("True receptor pairs", "Unlikely receptor pairs"),
    c("True receptor pairs", "Decoy receptor pairs")
)
significance.results.three.groups <- data.frame()

for (feature in features.of.interest) {
    significance.results <- all.receptor.data.faceted %>%
        group_by(epitope) %>%
        do({
            data <- .
            results <- lapply(comparisons, function(comp) {
                test <- wilcox.test(as.formula(paste(feature, "~ binding.pair.type")), data = data %>% filter(binding.pair.type %in% comp))
                data.frame(
                    Feature = feature,
                    Epitope = unique(data$epitope),
                    Comparison = paste(comp, collapse = " vs "),
                    P.Value = test$p.value
                )
            })
            bind_rows(results)
        }) %>%
        ungroup()

    significance.results.three.groups <- bind_rows(significance.results.three.groups, significance.results)
}

write.csv(significance.results.three.groups, paste0(out.path, "significance_results_per_epitope_three_groups.csv"))

# Make significance tables: Two-way comparisons (i.e., true receptor pairs vs unlikely and decoy receptor pairs)
comparisons <- c("True receptor pairs", "Negative receptor pairs")
significance.results.two.groups <- data.frame()

for (feature in features.of.interest) {
    two.way.data <- all.receptor.data.faceted %>%
        mutate(binding.pair.type = ifelse(binding.pair.type %in% c("Unlikely receptor pairs", "Decoy receptor pairs"), "Negative receptor pairs", binding.pair.type)) 
    significance.results <- two.way.data %>%
        group_by(epitope) %>%
        do({
            data <- .
            test <- wilcox.test(as.formula(paste(feature, "~ binding.pair.type")), data = data)
            data.frame(
                Feature = feature,
                Epitope = unique(data$epitope),
                Comparison = paste(comparisons, collapse = " vs "),
                P.Value = test$p.value
            )
        }) %>%
        ungroup()

    significance.results.two.groups <- bind_rows(significance.results.two.groups, significance.results)
}

write.csv(significance.results.two.groups, paste0(out.path, "significance_results_per_epitope_two_groups.csv"))

# Box plots for 2 and 3-way comparisons
# Top KAS
library(ggsignif)
all.receptor.data.two.way <- all.receptor.data %>%
        mutate(binding.pair.type = ifelse(binding.pair.type %in% c("Unlikely receptor pairs", "Decoy receptor pairs"), "Negative receptor pairs", binding.pair.type))


ggplot(all.receptor.data, aes(x = binding.pair.type, y = top.kas, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Top KAS values per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Top KAS") +
    geom_signif(comparisons = list(
        c("True receptor pairs", "Unlikely receptor pairs"),
        c("True receptor pairs", "Decoy receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3) 

ggsave(paste0(out.path, "box_plot_top_kas_values_per_epitope_three_group.png"), width = 12, height = 8)

ggplot(all.receptor.data.two.way, 
    aes(x = binding.pair.type, y = top.kas, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("True receptor pairs" = "#00BA38", "Negative receptor pairs" = "#F8766D")) +
    labs(title = "Top KAS values per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Top KAS") +
    geom_signif(comparisons = list(c("True receptor pairs", "Negative receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_top_kas_values_per_epitope_two_group.png"), width = 12, height = 8)

# Average KAS
ggplot(all.receptor.data, aes(x = binding.pair.type, y = avg.kas, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Average KAS values per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Average KAS") +
    geom_signif(comparisons = list(
        c("True receptor pairs", "Unlikely receptor pairs"),
        c("True receptor pairs", "Decoy receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_average_kas_values_per_epitope_three_group.png"), width = 12, height = 8)

ggplot(all.receptor.data.two.way, 
    aes(x = binding.pair.type, y = avg.kas, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("True receptor pairs" = "#00BA38", "Negative receptor pairs" = "#F8766D")) +
    labs(title = "Average KAS values per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Average KAS") +
    geom_signif(comparisons = list(c("True receptor pairs", "Negative receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_average_kas_values_per_epitope_two_group.png"), width = 12, height = 8)


# Median subgraph KAS
ggplot(all.receptor.data, aes(x = binding.pair.type, y = top.sg.median.kas, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Median subgraph KAS values per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Median subgraph KAS") +
    geom_signif(comparisons = list(
        c("True receptor pairs", "Unlikely receptor pairs"),
        c("True receptor pairs", "Decoy receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_median_subgraph_kas_values_per_epitope_three_group.png"), width = 12, height = 8)

ggplot(all.receptor.data.two.way, 
    aes(x = binding.pair.type, y = top.sg.median.kas, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    scale_fill_manual(values = c("True receptor pairs" = "#00BA38", "Negative receptor pairs" = "#F8766D")) +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Median subgraph KAS values per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Median subgraph KAS") +
    geom_signif(comparisons = list(c("True receptor pairs", "Negative receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_median_subgraph_kas_values_per_epitope_two_group.png"), width = 12, height = 8)

# Top Spearman correlation
ggplot(all.receptor.data, aes(x = binding.pair.type, y = top.spearman, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Top Spearman correlation values per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Top Spearman correlation") +
    # ylim(0.5,1) +
    geom_signif(comparisons = list(
        c("True receptor pairs", "Unlikely receptor pairs"),
        c("True receptor pairs", "Decoy receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_top_spearman_correlation_per_epitope_three_group.png"), width = 12, height = 8)

ggplot(all.receptor.data.two.way, 
    aes(x = binding.pair.type, y = top.spearman, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    scale_fill_manual(values = c("True receptor pairs" = "#00BA38", "Negative receptor pairs" = "#F8766D")) +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Top Spearman correlation values per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Top Spearman correlation") +
    # ylim(0.5,1) +
    geom_signif(comparisons = list(c("True receptor pairs", "Negative receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_top_spearman_correlation_per_epitope_two_group.png"), width = 12, height = 8)

# Average Spearman correlation
ggplot(all.receptor.data, aes(x = binding.pair.type, y = avg.spearman, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Average Spearman correlation values per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Average Spearman correlation") +
    # ylim(0,1) +
    geom_signif(comparisons = list(
        c("True receptor pairs", "Unlikely receptor pairs"),
        c("True receptor pairs", "Decoy receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_average_spearman_correlation_per_epitope_three_group.png"), width = 12, height = 8)

ggplot(all.receptor.data.two.way, 
    aes(x = binding.pair.type, y = avg.spearman, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    scale_fill_manual(values = c("True receptor pairs" = "#00BA38", "Negative receptor pairs" = "#F8766D")) +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Average Spearman correlation values per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Average Spearman correlation") +
    # ylim(0,1) +
    geom_signif(comparisons = list(c("True receptor pairs", "Negative receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_average_spearman_correlation_per_epitope_two_group.png"), width = 12, height = 8)

# Average Euclidean distance
ggplot(all.receptor.data, aes(x = binding.pair.type, y = avg.euc.dist, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Average Euclidean distance values per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Average Euclidean distance") +
    geom_signif(comparisons = list(
        c("True receptor pairs", "Unlikely receptor pairs"),
        c("True receptor pairs", "Decoy receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_average_euclidean_distance_per_epitope_three_group.png"), width = 12, height = 8)

ggplot(all.receptor.data.two.way, 
    aes(x = binding.pair.type, y = avg.euc.dist, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    scale_fill_manual(values = c("True receptor pairs" = "#00BA38", "Negative receptor pairs" = "#F8766D")) +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Average Euclidean distance values per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Average Euclidean distance") +
    geom_signif(comparisons = list(c("True receptor pairs", "Negative receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_average_euclidean_distance_per_epitope_two_group.png"), width = 12, height = 8)

# Average subgraph Euclidean distance
ggplot(all.receptor.data, aes(x = binding.pair.type, y = top.sg.median.kas, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Average subgraph Euclidean distance values per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Average subgraph Euclidean distance") +
    geom_signif(comparisons = list(
        c("True receptor pairs", "Unlikely receptor pairs"),
        c("True receptor pairs", "Decoy receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_average_subgraph_euclidean_distance_per_epitope_three_group.png"), width = 12, height = 8)

ggplot(all.receptor.data.two.way, 
    aes(x = binding.pair.type, y = avg.euc.dist, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    scale_fill_manual(values = c("True receptor pairs" = "#00BA38", "Negative receptor pairs" = "#F8766D")) +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Average subgraph Euclidean distance values per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Average subgraph Euclidean distance") +
    geom_signif(comparisons = list(c("True receptor pairs", "Negative receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_average_subgraph_euclidean_distance_per_epitope_two_group.png"), width = 12, height = 8)


# Average subgraph shape similarity
ggplot(all.receptor.data, aes(x = binding.pair.type, y = top.sg.cosine.sim, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Average subgraph shape similarity values per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Average subgraph shape similarity") +
    geom_signif(comparisons = list(
        c("True receptor pairs", "Unlikely receptor pairs"),
        c("True receptor pairs", "Decoy receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_average_subgraph_shape_similarity_per_epitope_three_group.png"), width = 12, height = 8)

ggplot(all.receptor.data.two.way, 
    aes(x = binding.pair.type, y = top.sg.cosine.sim, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    scale_fill_manual(values = c("True receptor pairs" = "#00BA38", "Negative receptor pairs" = "#F8766D")) +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Average subgraph shape similarity values per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Average subgraph shape similarity") +
    geom_signif(comparisons = list(c("True receptor pairs", "Negative receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_average_subgraph_shape_similarity_per_epitope_two_group.png"), width = 12, height = 8)

# Top subgraph euclidean distance
ggplot(all.receptor.data, aes(x = binding.pair.type, y = top.sg.euc.dist, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Top subgraph Euclidean distance values per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Top subgraph Euclidean distance") +
    geom_signif(comparisons = list(
        c("True receptor pairs", "Unlikely receptor pairs"),
        c("True receptor pairs", "Decoy receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_top_subgraph_euclidean_distance_per_epitope_three_group.png"), width = 12, height = 8)

ggplot(all.receptor.data.two.way, 
    aes(x = binding.pair.type, y = top.sg.euc.dist, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    scale_fill_manual(values = c("True receptor pairs" = "#00BA38", "Negative receptor pairs" = "#F8766D")) +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Top subgraph Euclidean distance values per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Top subgraph Euclidean distance") +
    geom_signif(comparisons = list(c("True receptor pairs", "Negative receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_top_subgraph_euclidean_distance_per_epitope_two_group.png"), width = 12, height = 8)

# Average number of kernel matches
ggplot(all.receptor.data, aes(x = binding.pair.type, y = kernel.count, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Average number of kernel matches per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Average number of kernel matches") +
    geom_signif(comparisons = list(
        c("True receptor pairs", "Unlikely receptor pairs"),
        c("True receptor pairs", "Decoy receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_average_number_of_kernel_matches_per_epitope_three_group.png"), width = 12, height = 8)

ggplot(all.receptor.data.two.way, 
    aes(x = binding.pair.type, y = kernel.count, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    scale_fill_manual(values = c("True receptor pairs" = "#00BA38", "Negative receptor pairs" = "#F8766D")) +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Average number of kernel matches per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Average number of kernel matches") +
    geom_signif(comparisons = list(c("True receptor pairs", "Negative receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_average_number_of_kernel_matches_per_epitope_two_group.png"), width = 12, height = 8)

# Average CDR3 sequence similarity
ggplot(all.receptor.data, aes(x = binding.pair.type, y = CDR3.similarity, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Average CDR3 sequence similarity per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Average CDR3 sequence similarity") +
    geom_signif(comparisons = list(
        c("True receptor pairs", "Unlikely receptor pairs"),
        c("True receptor pairs", "Decoy receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_average_cdr3_sequence_similarity_per_epitope_three_group.png"), width = 12, height = 8)

ggplot(all.receptor.data.two.way, 
    aes(x = binding.pair.type, y = CDR3.similarity, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    scale_fill_manual(values = c("True receptor pairs" = "#00BA38", "Negative receptor pairs" = "#F8766D")) +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Average CDR3 sequence similarity per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Average CDR3 sequence similarity") +
    geom_signif(comparisons = list(c("True receptor pairs", "Negative receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_average_cdr3_sequence_similarity_per_epitope_two_group.png"), width = 12, height = 8)

# Average CDR sequence similarity
ggplot(all.receptor.data, aes(x = binding.pair.type, y = CDR.similarity, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Average CDR sequence similarity per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Average CDR sequence similarity") +
    geom_signif(comparisons = list(
        c("True receptor pairs", "Unlikely receptor pairs"),
        c("True receptor pairs", "Decoy receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_average_cdr_sequence_similarity_per_epitope_three_group.png"), width = 12, height = 8)

ggplot(all.receptor.data.two.way, 
    aes(x = binding.pair.type, y = CDR.similarity, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    scale_fill_manual(values = c("True receptor pairs" = "#00BA38", "Negative receptor pairs" = "#F8766D")) +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Average CDR sequence similarity per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Average CDR sequence similarity") +
    geom_signif(comparisons = list(c("True receptor pairs", "Negative receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_average_cdr_sequence_similarity_per_epitope_two_group.png"), width = 12, height = 8)

# Average full sequence similarity
ggplot(all.receptor.data, aes(x = binding.pair.type, y = full.similarity, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Average full sequence similarity per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Average full sequence similarity") +
    geom_signif(comparisons = list(
        c("True receptor pairs", "Unlikely receptor pairs"),
        c("True receptor pairs", "Decoy receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_average_full_sequence_similarity_per_epitope_three_group.png"), width = 12, height = 8)

ggplot(all.receptor.data.two.way, 
    aes(x = binding.pair.type, y = full.similarity, fill = binding.pair.type)) +
    geom_boxplot() +
    facet_wrap(~ref.epitope, nrow = 2) +
    theme_minimal() +
    scale_fill_manual(values = c("True receptor pairs" = "#00BA38", "Negative receptor pairs" = "#F8766D")) +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Average full sequence similarity per receptor:receptor pair, colored by receptor pair type",
         x = "Receptor:receptor pair type",
         y = "Average full sequence similarity") +
    geom_signif(comparisons = list(c("True receptor pairs", "Negative receptor pairs")),
        map_signif_level = TRUE, 
        textsize = 3)

ggsave(paste0(out.path, "box_plot_average_full_sequence_similarity_per_epitope_two_group.png"), width = 12, height = 8)

# All features, all data (2 and 3-way comparisons), not faceted by epitope
ggplot((all.receptor.data %>% 
        select(binding.pair.type, features.of.interest) %>%
        pivot_longer(cols = features.of.interest, names_to = "feature", values_to = "value")),
    aes(x = binding.pair.type, y = value, fill = binding.pair.type)) +
    geom_boxplot() +
    theme_minimal() +
    facet_wrap(~feature, scales = "free", nrow = 2) +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
    labs(title = "All features",
         x = "Receptor:receptor pair type",
         y = "Value") 

ggsave(paste0(out.path, "box_plot_all_features_all_data.png"), width = 12, height = 8)

ggplot((all.receptor.data.two.way %>% 
        select(binding.pair.type, features.of.interest) %>%
        pivot_longer(cols = features.of.interest, names_to = "feature", values_to = "value")),
    aes(x = binding.pair.type, y = value, fill = feature)) +
    geom_boxplot() +
    theme_minimal() +
    facet_wrap(~feature, scales = "free", nrow = 2) +
    theme(strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
    labs(title = "All features",
         x = "Receptor:receptor pair type",
         y = "Value")

ggsave(paste0(out.path, "box_plot_all_features_all_data_two_group.png"), width = 12, height = 8)



# Dot plots
# Subgraph shape similarity vs subgraph Euclidean distance, colored by receptor pair type
ggplot(all.data.master, aes(x = avg.subgraph.distance, y = avg.subgraph.shape.sim, color = binding.pair.type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.subgraph.shape.sim > 0.991 & avg.subgraph.distance < 3.5 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.subgraph.shape.sim > 0.991 & avg.subgraph.distance < 3.5 & type == "Decoy receptor pairs")), 
              aes(label = paste0(ref.epitope)),
              vjust = -1, 
              size = 3) +
    labs(title = "Subgraph shape similarity vs subgraph Euclidean distance",
         x = "Average subgraph Euclidean distance",
         y = "Average subgraph shape similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "dot_plot_subgraph_shape_similarity_vs_euclidean_distance.png"), width = 12, height = 8)

# Top KAS vs average KAS, colored by receptor pair type
ggplot(all.data.master, aes(x = avg.kdist.per.receptor.pair, y = avg.top.kdist.per.receptor.pair, color = binding.pair.type)) +
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

ggsave(paste0(out.path, "dot_plot_top_kas_vs_average_kas.png"), width = 12, height = 8)

# Average Spearman vs average KAS, colored by receptor pair type
ggplot(all.data.master, aes(x = avg.kdist.per.receptor.pair, y = avg.corr.per.receptor.pair, color = binding.pair.type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.corr.per.receptor.pair > 0.74 & avg.kdist.per.receptor.pair > -0.05 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.corr.per.receptor.pair > 0.74 & avg.kdist.per.receptor.pair > -0.05 & type == "Decoy receptor pairs")), 
              aes(label = paste0(ref.epitope)),
              vjust = -1, 
              size = 3) +
    labs(title = "Average Spearman correlation vs average KAS",
         x = "Average KAS",
         y = "Average Spearman correlation",
            color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "dot_plot_average_spearman_correlation_vs_average_kas.png"), width = 12, height = 8)

# Average Euclidean distance vs average KAS, colored by receptor pair type
ggplot(all.data.master, aes(x = avg.kdist.per.receptor.pair, y = avg.distance.per.receptor.pair, color = binding.pair.type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.distance.per.receptor.pair < 6 & avg.kdist.per.receptor.pair > -0.05 & type == "Unlikely receptor pairs")),
                aes(label = paste0(ref.epitope, "+", samp.epitope)),
                vjust = -1, 
                size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.distance.per.receptor.pair < 6 & avg.kdist.per.receptor.pair > -0.05 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Average Euclidean distance vs average KAS",
         x = "Average KAS",
         y = "Average Euclidean distance",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "dot_plot_average_euclidean_distance_vs_average_kas.png"), width = 12, height = 8)

# Average subgraph shape similarity vs top KAS, colored by receptor pair type
ggplot(all.data.master, aes(x = avg.top.kdist.per.receptor.pair, y = avg.subgraph.shape.sim, color = binding.pair.type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.subgraph.shape.sim > 0.991 & avg.top.kdist.per.receptor.pair > 0.2 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.subgraph.shape.sim > 0.991 & avg.top.kdist.per.receptor.pair > 0.2 & type == "Decoy receptor pairs")), 
              aes(label = paste0(ref.epitope)),
              vjust = -1, 
              size = 3) +
    labs(title = "Top KAS vs average subgraph shape similarity",
         x = "Top KAS",
         y = "Average subgraph shape similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "dot_plot_top_kas_vs_average_subgraph_shape_similarity.png"), width = 12, height = 8)

# Average subgraph shape similarity vs average KAS, colored by receptor pair type
ggplot(all.data.master, aes(x = avg.kdist.per.receptor.pair, y = avg.subgraph.shape.sim, color = binding.pair.type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.subgraph.shape.sim > 0.991 & avg.kdist.per.receptor.pair > -0.05 & type == "Unlikely receptor pairs")),
                aes(label = paste0(ref.epitope, "+", samp.epitope)),
                vjust = -1, 
                size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.subgraph.shape.sim > 0.991 & avg.kdist.per.receptor.pair > -0.08 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Average KAS vs average subgraph shape similarity",
         x = "Average KAS",
         y = "Average subgraph shape similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "dot_plot_average_kas_vs_average_subgraph_shape_similarity.png"), width = 12, height = 8)

# Average shape similarity vs average number of kernel matches, colored by receptor pair type
ggplot(all.data.master, aes(x = avg.kernels.per.receptor.pair, y = avg.subgraph.shape.sim, color = binding.pair.type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.subgraph.shape.sim > 0.991 & avg.kernels.per.receptor.pair > 8 & type == "Unlikely receptor pairs")),
                aes(label = paste0(ref.epitope, "+", samp.epitope)),
                vjust = -1, 
                size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.subgraph.shape.sim > 0.991 & avg.kernels.per.receptor.pair > 8 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Average number of kernel matches vs average subgraph shape similarity",
         x = "Average number of kernel matches",
         y = "Average subgraph shape similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "dot_plot_average_kernel_matches_vs_average_subgraph_shape_similarity.png"), width = 12, height = 8)

# Average number of kernel matches versus average Euclidean distance
ggplot(all.data.master, aes(x = avg.distance.per.receptor.pair, y = avg.kernels.per.receptor.pair, color = binding.pair.type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.kernels.per.receptor.pair > 8 & avg.distance.per.receptor.pair < 6 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.kernels.per.receptor.pair > 8 & avg.distance.per.receptor.pair < 6 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Average number of kernel matches vs average Euclidean distance",
         x = "Average Euclidean distance",
         y = "Average number of kernel matches",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "dot_plot_average_kernel_matches_vs_average_euclidean_distance.png"), width = 12, height = 8)

# Average number of kernel matches versus Spearman correlation coefficient
ggplot(all.data.master, aes(x = avg.corr.per.receptor.pair, 
    y = avg.kernels.per.receptor.pair, color = binding.pair.type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.kernels.per.receptor.pair > 8 & avg.corr.per.receptor.pair > 0.74 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.kernels.per.receptor.pair > 8 & avg.corr.per.receptor.pair > 0.74 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Average number of kernel matches vs Spearman correlation coefficient",
         x = "Average Spearman correlation",
         y = "Average number of kernel matches",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "dot_plot_average_kernel_matches_vs_average_spearman_correlation.png"), width = 12, height = 8)

# Average number of kernel matches versus top KAS
ggplot(all.data.master, aes(x = avg.top.kdist.per.receptor.pair, 
    y = avg.kernels.per.receptor.pair, color = binding.pair.type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.kernels.per.receptor.pair > 8 & avg.top.kdist.per.receptor.pair > 0.2 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.kernels.per.receptor.pair > 8 & avg.top.kdist.per.receptor.pair > 0.2 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Average number of kernel matches vs top KAS",
         x = "Top KAS",
         y = "Average number of kernel matches",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "dot_plot_average_kernel_matches_vs_top_kas.png"), width = 12, height = 8)

# Top KAS versus Spearman correlation
ggplot(all.data.master, aes(x = avg.corr.per.receptor.pair, 
    y = avg.top.kdist.per.receptor.pair, color = binding.pair.type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.top.kdist.per.receptor.pair > 0.2 & avg.corr.per.receptor.pair > 0.74 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.top.kdist.per.receptor.pair > 0.2 & avg.corr.per.receptor.pair > 0.74 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Top KAS vs Spearman correlation coefficient",
         x = "Average Spearman correlation",
         y = "Top KAS",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "dot_plot_top_kas_vs_spearman_correlation.png"), width = 12, height = 8)

# CDR3 sequence similarity versus Spearman correlation
ggplot(all.data.master, aes(x = avg.corr.per.receptor.pair, 
    y = avg.CDR3.sequence.sim, color = binding.pair.type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.CDR3.sequence.sim > 0.4 & avg.corr.per.receptor.pair > 0.74 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.CDR3.sequence.sim > 0.4 & avg.corr.per.receptor.pair > 0.74 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "CDR3 sequence similarity vs Spearman correlation coefficient",
         x = "Average Spearman correlation",
         y = "Average CDR3 sequence similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "dot_plot_cdr3_sequence_similarity_vs_spearman_correlation.png"), width = 12, height = 8)

# CDR3 sequence similarity versus top KAS
ggplot(all.data.master, aes(x = avg.top.kdist.per.receptor.pair, 
    y = avg.CDR3.sequence.sim, color = binding.pair.type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.CDR3.sequence.sim > 0.4 & avg.top.kdist.per.receptor.pair > 0.2 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.CDR3.sequence.sim > 0.4 & avg.top.kdist.per.receptor.pair > 0.2 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "CDR3 sequence similarity vs top KAS",
         x = "Top KAS",
         y = "Average CDR3 sequence similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "dot_plot_cdr3_sequence_similarity_vs_top_kas.png"), width = 12, height = 8)

# CDR3 sequence similarity versus average KAS
ggplot(all.data.master, aes(x = avg.kdist.per.receptor.pair, 
    y = avg.CDR3.sequence.sim, color = binding.pair.type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.CDR3.sequence.sim > 0.4 & avg.kdist.per.receptor.pair > -0.05 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.CDR3.sequence.sim > 0.4 & avg.kdist.per.receptor.pair > -0.05 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "CDR3 sequence similarity vs average KAS",
         x = "Average KAS",
         y = "Average CDR3 sequence similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "dot_plot_cdr3_sequence_similarity_vs_average_kas.png"), width = 12, height = 8)

# CDR3 sequence similarity versus Euclidean distance
ggplot(all.data.master, aes(x = avg.distance.per.receptor.pair, 
    y = avg.CDR3.sequence.sim, color = binding.pair.type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.CDR3.sequence.sim > 0.4 & avg.distance.per.receptor.pair < 6 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.CDR3.sequence.sim > 0.4 & avg.distance.per.receptor.pair < 6 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "CDR3 sequence similarity vs average Euclidean distance",
         x = "Average Euclidean distance",
         y = "Average CDR3 sequence similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "dot_plot_cdr3_sequence_similarity_vs_average_euclidean_distance.png"), width = 12, height = 8)

# CDR3 sequence similarity versus average number of kernel matches
ggplot(all.data.master, aes(x = avg.kernels.per.receptor.pair, 
    y = avg.CDR3.sequence.sim, color = binding.pair.type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.CDR3.sequence.sim > 0.4 & avg.kernels.per.receptor.pair > 8 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.CDR3.sequence.sim > 0.4 & avg.kernels.per.receptor.pair > 8 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "CDR3 sequence similarity vs average number of kernel matches",
         x = "Average number of kernel matches",
         y = "Average CDR3 sequence similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "dot_plot_cdr3_sequence_similarity_vs_average_kernel_matches.png"), width = 12, height = 8)

# Full sequence similarity versus Spearman correlation
ggplot(all.data.master, aes(x = avg.corr.per.receptor.pair, 
    y = avg.full.sequence.sim, color = binding.pair.type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.full.sequence.sim > 0.475 & avg.corr.per.receptor.pair > 0.74 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.full.sequence.sim > 0.475 & avg.corr.per.receptor.pair > 0.74 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Full sequence similarity vs Spearman correlation coefficient",
         x = "Average Spearman correlation",
         y = "Average full sequence similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "dot_plot_full_sequence_similarity_vs_spearman_correlation.png"), width = 12, height = 8)

# Full sequence similarity versus top KAS
ggplot(all.data.master, aes(x = avg.top.kdist.per.receptor.pair, 
    y = avg.full.sequence.sim, color = binding.pair.type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.full.sequence.sim > 0.475 & avg.top.kdist.per.receptor.pair > 0.2 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.full.sequence.sim > 0.475 & avg.top.kdist.per.receptor.pair > 0.2 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Full sequence similarity vs top KAS",
         x = "Top KAS",
         y = "Average full sequence similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "dot_plot_full_sequence_similarity_vs_top_kas.png"), width = 12, height = 8)

# Full sequence similarity versus average KAS
ggplot(all.data.master, aes(x = avg.kdist.per.receptor.pair, 
    y = avg.full.sequence.sim, color = binding.pair.type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.full.sequence.sim > 0.475 & avg.kdist.per.receptor.pair > -0.05 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.full.sequence.sim > 0.475 & avg.kdist.per.receptor.pair > -0.05 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Full sequence similarity vs average KAS",
         x = "Average KAS",
         y = "Average full sequence similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "dot_plot_full_sequence_similarity_vs_average_kas.png"), width = 12, height = 8)

# Full sequence similarity versus Euclidean distance
ggplot(all.data.master, aes(x = avg.distance.per.receptor.pair, 
    y = avg.full.sequence.sim, color = binding.pair.type)) +
    geom_point() +
    geom_text_repel(data = subset(all.data.master, type == "True receptor pairs"), 
              aes(label = ifelse(ref.epitope != samp.epitope, paste0(ref.epitope, "+", samp.epitope), ref.epitope)), 
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.full.sequence.sim > 0.475 & avg.distance.per.receptor.pair < 6 & type == "Unlikely receptor pairs")), 
              aes(label = paste0(ref.epitope, "+", samp.epitope)),
              vjust = -1, 
              size = 3) +
    geom_text_repel(data = subset(all.data.master, (avg.full.sequence.sim > 0.475 & avg.distance.per.receptor.pair < 6 & type == "Decoy receptor pairs")),
                aes(label = paste0(ref.epitope)),
                vjust = -1, 
                size = 3) +
    labs(title = "Full sequence similarity vs average Euclidean distance",
         x = "Average Euclidean distance",
         y = "Average full sequence similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "dot_plot_full_sequence_similarity_vs_average_euclidean_distance.png"), width = 12, height = 8)

# Dot plots of all data (not averages)
# Top KAS vs Spearman
ggplot(all.receptor.data, aes(x = top.kas, y = top.spearman, color = binding.pair.type)) +
    geom_point(alpha = 0.5) +
    labs(title = "Top KAS vs Spearman correlation coefficient",
         x = "Top KAS",
         y = "Spearman correlation",
         color = "Receptor:receptor pair type") +
    ylim(0.6, 1) +
    theme_minimal()

ggsave(paste0(out.path, "all_data_dot_plot_top_kas_vs_spearman_correlation.png"), width = 12, height = 8)

# Top KAS vs Euclidean distance
ggplot(all.receptor.data, aes(x = top.kas, y = top.euc.dist, color = binding.pair.type)) +
    geom_point(alpha = 0.5) +
    labs(title = "Top KAS vs Euclidean distance",
         x = "Top KAS",
         y = "Euclidean distance",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "all_data_dot_plot_top_kas_vs_euclidean_distance.png"), width = 12, height = 8)

# Top KAS vs subgraph shape similarity
ggplot(all.receptor.data, aes(x = top.kas, y = top.sg.cosine.sim, color = binding.pair.type)) +
    geom_point(alpha = 0.5) +
    labs(title = "Top KAS vs subgraph shape similarity",
         x = "Top KAS",
         y = "Subgraph shape similarity",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "all_data_dot_plot_top_kas_vs_subgraph_shape_similarity.png"), width = 12, height = 8)

# Top KAS vs number of kernels
ggplot(all.receptor.data, aes(x = top.kas, y = kernel.count, color = binding.pair.type)) +
    geom_point(alpha = 0.5) +
    labs(title = "Top KAS vs number of kernels",
         x = "Top KAS",
         y = "Number of kernels",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "all_data_dot_plot_top_kas_vs_number_of_kernels.png"), width = 12, height = 8)

# CDR3 sequence similarity vs Spearman correlation
ggplot(all.receptor.data, aes(x = CDR3.similarity, y = top.spearman, color = binding.pair.type)) +
    geom_point(alpha = 0.5) +
    labs(title = "CDR3 sequence similarity vs Spearman correlation coefficient",
         x = "CDR3 sequence similarity",
         y = "Spearman correlation",
         color = "Receptor:receptor pair type") +
    ylim(0.6, 1) +
    theme_minimal()

ggsave(paste0(out.path, "all_data_dot_plot_cdr3_sequence_similarity_vs_spearman_correlation.png"), width = 12, height = 8)

# CDR3 sequence similarity vs top KAS
ggplot(all.receptor.data, aes(x = CDR3.similarity, y = top.kas, color = binding.pair.type)) +
    geom_point(alpha = 0.5) +
    labs(title = "CDR3 sequence similarity vs top KAS",
         x = "CDR3 sequence similarity",
         y = "Top KAS",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "all_data_dot_plot_cdr3_sequence_similarity_vs_top_kas.png"), width = 12, height = 8)

# CDR3 sequence similarity vs Euclidean distance
ggplot(all.receptor.data, aes(x = CDR3.similarity, y = top.euc.dist, color = binding.pair.type)) +
    geom_point(alpha = 0.5) +
    labs(title = "CDR3 sequence similarity vs Euclidean distance",
         x = "CDR3 sequence similarity",
         y = "Euclidean distance",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "all_data_dot_plot_cdr3_sequence_similarity_vs_euclidean_distance.png"), width = 12, height = 8)

# CDR3 sequence similarity vs number of kernels
ggplot(all.receptor.data, aes(x = CDR3.similarity, y = kernel.count, color = binding.pair.type)) +
    geom_point(alpha = 0.5) +
    labs(title = "CDR3 sequence similarity vs number of kernels",
         x = "CDR3 sequence similarity",
         y = "Number of kernels",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "all_data_dot_plot_cdr3_sequence_similarity_vs_number_of_kernels.png"), width = 12, height = 8)

# Full sequence similarity vs Spearman correlation 
ggplot(all.receptor.data, aes(x = full.similarity, y = top.spearman, color = binding.pair.type)) +
    geom_point(alpha = 0.5) +
    labs(title = "Full sequence similarity vs Spearman correlation coefficient",
         x = "Full sequence similarity",
         y = "Spearman correlation",
         color = "Receptor:receptor pair type") +
    ylim(0.6, 1) +
    theme_minimal()

ggsave(paste0(out.path, "all_data_dot_plot_full_sequence_similarity_vs_spearman_correlation.png"), width = 12, height = 8)

# Full sequence similarity vs top KAS
ggplot(all.receptor.data, aes(x = full.similarity, y = top.kas, color = binding.pair.type)) +
    geom_point(alpha = 0.5) +
    labs(title = "Full sequence similarity vs top KAS",
         x = "Full sequence similarity",
         y = "Top KAS",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "all_data_dot_plot_full_sequence_similarity_vs_top_kas.png"), width = 12, height = 8)

# Full sequence similarity vs Euclidean distance
ggplot(all.receptor.data, aes(x = full.similarity, y = top.euc.dist, color = binding.pair.type)) +
    geom_point(alpha = 0.5) +
    labs(title = "Full sequence similarity vs Euclidean distance",
         x = "Full sequence similarity",
         y = "Euclidean distance",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "all_data_dot_plot_full_sequence_similarity_vs_euclidean_distance.png"), width = 12, height = 8)

# Full sequence similarity vs number of kernels
ggplot(all.receptor.data, aes(x = full.similarity, y = kernel.count, color = binding.pair.type)) +
    geom_point(alpha = 0.5) +
    labs(title = "Full sequence similarity vs number of kernels",
         x = "Full sequence similarity",
         y = "Number of kernels",
         color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "all_data_dot_plot_full_sequence_similarity_vs_number_of_kernels.png"), width = 12, height = 8)


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

saveWidget(sankeyplot, file = paste0(out.path, "sankey_plot_kernel_matches.html"), selfcontained = TRUE)

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

saveWidget(sankeyplot, file = paste0(out.path, "sankey_plot_receptor_pairs.html"), selfcontained = TRUE)


# PCA of all features (top KAS, average KAS, Spearman correlation, Euclidean distance, subgraph shape similarity, average number of kernels per pairing) from master data list
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

ggsave(paste0(out.path, "pca_plot.png"), width = 6, height = 4)

autoplot(pca.result, data = pca.data, colour = 'type', loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3) +
    labs(title = "PCA of all features, with eigenvectors",
            x = "Principal Component 1",
            y = "Principal Component 2",
            color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "pca_plot_with_eigenvectors.png"), width = 6, height = 4)

fviz_pca_ind(pca.result, col.ind = pca.data$type, geom = "point", addEllipses=TRUE, ellipse.level=0.95) +
    labs(title = "PCA of all features, individual factor map",
            x = "Principal Component 1",
            y = "Principal Component 2",
            color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "pca_plot_individual_factor_map.png"), width = 6, height = 4)

# PCA of all features (top KAS, average KAS, Spearman correlation, Euclidean distance, subgraph shape similarity, average number of kernels per pairing) on all receptors
# Select data
pca.data <- all.receptor.data %>%
    select(features.of.interest, type, ref.epitope, samp.epitope)

# Perform PCA
pca.result <- prcomp(pca.data %>% select(-ref.epitope, -samp.epitope, -type), scale. = TRUE)

# Plot PCA
autoplot(pca.result, data = pca.data, colour = 'type') +
    labs(title = "PCA of all features",
            x = "Principal Component 1",
            y = "Principal Component 2",
            color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "pca_plot_all_receptors.png"), width = 6, height = 4)

autoplot(pca.result, data = pca.data, colour = 'type', loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3) +
    labs(title = "PCA of all features, with eigenvectors",
            x = "Principal Component 1",
            y = "Principal Component 2",
            color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "pca_plot_with_eigenvectors_all_receptors.png"), width = 6, height = 4)

fviz_pca_ind(pca.result, col.ind = pca.data$type, geom = "point", addEllipses=TRUE, ellipse.level=0.95) +
    labs(title = "PCA of all features, individual factor map",
            x = "Principal Component 1",
            y = "Principal Component 2",
            color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "pca_plot_individual_factor_map_all_receptors.png"), width = 6, height = 4)


# PCA on all data and features (all kernels)
# Select data
pca.data <- all.data %>% 
    select(c(3:5, 8:9, 13:15))

# Perform PCA
pca.result <- prcomp(pca.data %>% select(-ref.epitope, -samp.epitope, -type), scale. = TRUE)

# Plot PCA
autoplot(pca.result, data = pca.data, colour = 'type') +
    labs(title = "PCA of all features",
            x = "Principal Component 1",
            y = "Principal Component 2",
            color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "pca_plot_all_kernels.png"), width = 6, height = 4)

autoplot(pca.result, data = pca.data, colour = 'type', loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3) +
    labs(title = "PCA of all features, with eigenvectors",
            x = "Principal Component 1",
            y = "Principal Component 2",
            color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "pca_plot_with_eigenvectors_all_kernels.png"), width = 6, height = 4)

fviz_pca_ind(pca.result, col.ind = pca.data$type, geom = "point", addEllipses=TRUE, ellipse.level=0.95) +
    labs(title = "PCA of all features, individual factor map",
            x = "Principal Component 1",
            y = "Principal Component 2",
            color = "Receptor:receptor pair type") +
    theme_minimal()

ggsave(paste0(out.path, "pca_plot_individual_factor_map_all_kernels.png"), width = 6, height = 4)


# Look at correlation between sequence and structural features
plot_correlation(all.receptor.data %>% select(-ref.epitope, -samp.epitope, -binding.pair.type, -top.sg.id))

ggsave(paste0(out.path, "feature_correlation_plot.png"), width = 10, height = 10)
