# ============================================================
# File: case-study-selection.R
# Description: Script for selecting receptors for further investigation based on the results of PIEMA and analyzing them in the context of the logistic regression classifiers
# Author: Maia Bennett-Boehm
# University of Nebraska at Omaha
# ============================================================

# Load the necessary libraries
library(tidyverse)

# Import the data
input.path <- "./analysis/classifier/"
out.path <- "./analysis/case-study/"
model1.preds <- read_csv(paste0(input.path, "all_model1_preds.csv"))
model2.preds <- read_csv(paste0(input.path, "all_model2_preds.csv"))

dir.create(out.path, showWarnings = FALSE)

# Subset to unique receptor pairs for investigation via PIEMA step 7
unique.pairs <- model2.preds %>% 
    distinct(pair.id, .keep_all = TRUE)

unique.true.pairs <- unique.pairs %>%
    filter(shared.specificity == "Yes")

unique.unlikely.pairs <- unique.pairs %>% 
    filter(ref.epitope != samp.epitope) %>%
    filter(!str_detect("decoy", pair.id))

# Subset further to true pairs with high KAS scores (top 10%), true pairs with low KAS scores (top 10%)
dissim.seq.pairs <- unique.true.pairs %>%
    arrange(full.similarity) %>%
    head(150)

sim.seq.pairs <- unique.true.pairs %>%
    arrange(desc(full.similarity)) %>%
    head(150)

# Subset further to high sequence similarity and low sequence similarity in each previous group (15 per, for ease of analysis)
sim.seq.high.scoring.pairs <- sim.seq.pairs %>%
    arrange(desc(top.kas)) %>%
    head(75)

dissim.seq.high.scoring.pairs <- dissim.seq.pairs %>%
    arrange(desc(top.kas)) %>%
    head(75)

sim.seq.low.scoring.pairs <- sim.seq.pairs %>%
    arrange(top.kas) %>%
    head(75)

dissim.seq.low.scoring.pairs <- dissim.seq.pairs %>%
    arrange(top.kas) %>%
    head(75)

# Write data along with the pair IDs of the top 10% of each subset in comma-sep format to be used in PIEMA
write.csv(sim.seq.high.scoring.pairs, paste0(out.path, "sim_seq_high_scoring_pairs.csv"), row.names = FALSE)
write.csv(dissim.seq.high.scoring.pairs, paste0(out.path, "dissim_seq_high_scoring_pairs.csv"), row.names = FALSE)
write.csv(sim.seq.low.scoring.pairs, paste0(out.path, "sim_seq_low_scoring_pairs.csv"), row.names = FALSE)
write.csv(dissim.seq.low.scoring.pairs, paste0(out.path, "dissim_seq_low_scoring_pairs.csv"), row.names = FALSE)

# Top 10% pair IDs for PIEMA
pair.ids <- sim.seq.high.scoring.pairs %>% 
    head(15) %>% 
    pull(pair.id)
writeLines(paste(pair.ids, collapse = ","), con = paste0(out.path, "similar-seq-high-scores.txt"))

pair.ids <- dissim.seq.high.scoring.pairs %>% 
    head(15) %>% 
    pull(pair.id)
writeLines(paste(pair.ids, collapse = ","), con = paste0(out.path, "dissimilar-seq-high-scores.txt"))

pair.ids <- sim.seq.low.scoring.pairs %>% 
    head(15) %>% 
    pull(pair.id)
writeLines(paste(pair.ids, collapse = ","), con = paste0(out.path, "similar-seq-low-scores.txt"))

pair.ids <- dissim.seq.low.scoring.pairs %>% 
    head(15) %>% 
    pull(pair.id)
writeLines(paste(pair.ids, collapse = ","), con = paste0(out.path, "dissimilar-seq-low-scores.txt"))


# Subset unlikely pairs with high sequence similarity (top 10%) and high KAS scores (top 10%) for further investigation in PIEMA
sim.seq.unlikely.pairs <- unique.unlikely.pairs %>%
    arrange(desc(full.similarity)) %>%
    head(150)

high.scoring.unlikely.pairs <- unique.unlikely.pairs %>%
    arrange(desc(top.kas)) %>%
    head(150)

write.csv(sim.seq.unlikely.pairs, paste0(out.path, "sim_seq_unlikely_pairs.csv"), row.names = FALSE)
write.csv(high.scoring.unlikely.pairs, paste0(out.path, "high_scoring_unlikely_pairs.csv"), row.names = FALSE)

# Write top 10% pair IDs for PIEMA
pair.ids <- sim.seq.unlikely.pairs %>% 
    head(15) %>% 
    pull(pair.id)
writeLines(paste(pair.ids, collapse = ","), con = paste0(out.path, "similar-seq-unlikely.txt"))

pair.ids <- high.scoring.unlikely.pairs %>% 
    head(15) %>% 
    pull(pair.id)
writeLines(paste(pair.ids, collapse = ","), con = paste0(out.path, "high-scoring-unlikely.txt"))

# Now, investigate the performance of the logistic regression classifiers on these pairs
# First, subset the data to only the pairs of interest including all model approaches (combined, sequence, structure, for both model tasks)
sim.seq.high.scoring.preds <- model1.preds %>% 
    filter(pair.id %in% sim.seq.high.scoring.pairs$pair.id) %>%
    bind_rows(model2.preds %>% filter(pair.id %in% sim.seq.high.scoring.pairs$pair.id))

dissim.seq.high.scoring.preds <- model1.preds %>% 
    filter(pair.id %in% dissim.seq.high.scoring.pairs$pair.id) %>%
    bind_rows(model2.preds %>% filter(pair.id %in% dissim.seq.high.scoring.pairs$pair.id)) 

sim.seq.low.scoring.preds <- model1.preds %>%
    filter(pair.id %in% sim.seq.low.scoring.pairs$pair.id) %>%
    bind_rows(model2.preds %>% filter(pair.id %in% sim.seq.low.scoring.pairs$pair.id))

dissim.seq.low.scoring.preds <- model1.preds %>%
    filter(pair.id %in% dissim.seq.low.scoring.pairs$pair.id) %>%
    bind_rows(model2.preds %>% filter(pair.id %in% dissim.seq.low.scoring.pairs$pair.id))

all.highlighted.preds <- sim.seq.high.scoring.preds %>%
    mutate(group = "sim.seq.high.scoring") %>%
    bind_rows(dissim.seq.high.scoring.preds %>%
        mutate(group = "dissim.seq.high.scoring")) %>%
    bind_rows(sim.seq.low.scoring.preds%>%
        mutate(group = "sim.seq.low.scoring")) %>%
    bind_rows(dissim.seq.low.scoring.preds%>%
        mutate(group = "dissim.seq.low.scoring"))

sim.seq.unlikely.preds <- model1.preds %>%
    filter(pair.id %in% sim.seq.unlikely.pairs$pair.id) %>%
    bind_rows(model2.preds %>% filter(pair.id %in% sim.seq.unlikely.pairs$pair.id))

high.scoring.unlikely.preds <- model1.preds %>%
    filter(pair.id %in% high.scoring.unlikely.pairs$pair.id) %>%
    bind_rows(model2.preds %>% filter(pair.id %in% high.scoring.unlikely.pairs$pair.id))

all.unlikely.preds <- sim.seq.unlikely.preds %>%
    mutate(group = "sim.seq.unlikely") %>%
    bind_rows(high.scoring.unlikely.preds %>%
        mutate(group = "high.scoring.unlikely"))

write.csv(all.highlighted.preds, paste0(out.path, "highlighted_preds.csv"), row.names = FALSE)
write.csv(all.unlikely.preds, paste0(out.path, "unlikely_preds.csv"), row.names = FALSE)

# Plot scores classifier predictions per model/approach/etc. for each group
# All structural and sequence scores per group
piema.scores <- c("avg.kas", "avg.euc.dist", "avg.spearman", "top.sg.euc.dist", "top.sg.cosine.sim", "kernel.count", "top.kas", "CDR3.similarity", "full.similarity")
ggplot(all.highlighted.preds %>%
    distinct(pair.id, .keep_all = TRUE) %>%
    pivot_longer(cols = c(piema.scores), names_to = "Metric", values_to = "Score"), 
    aes(x = Metric, y = Score, fill = group)) +
    geom_boxplot() +
    facet_wrap(~Metric, scales = "free") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PIEMA scores for highlighted receptor pair groups", x = "Metric", y = "Score") +
    scale_fill_viridis_d(option = "magma", begin = 0.2, end = 0.9)

ggsave(paste0(out.path, "highlighted_piema_scores.png"), width = 12, height = 12)

# Proportion of Yes and No
# By model and approach
ggplot(all.highlighted.preds %>%
    filter(Dataset == "CDR3") %>%
    mutate(`LogReg Approach` = paste0(Model, ", ", Approach)) %>%
    group_by(`LogReg Approach`, group, .pred_class) %>%
    summarise(count = n()) %>%
    group_by(`LogReg Approach`, group) %>%
    mutate(proportion = count / sum(count) * 100), 
    aes(x = .pred_class, y = proportion, fill = group)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ `LogReg Approach`) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Predictions for highlighted receptor pair groups", x = "Prediction of shared specificity", y = "% of observations") +
    scale_fill_viridis_d(option = "magma", begin = 0.2, end = 0.9)

ggsave(paste0(out.path, "highlighted_pred_classes.png"), width = 15, height = 10)

# By model, approach, and features
ggplot(all.highlighted.preds %>%
    filter(Dataset == "CDR3") %>%
    mutate(`LogReg Approach` = paste0(Model, ", ", Approach, " (", Feature, ")")) %>%
    group_by(`LogReg Approach`, group, .pred_class) %>%
    summarise(count = n()) %>%
    group_by(`LogReg Approach`, group) %>%
    mutate(proportion = count / sum(count) * 100), 
    aes(x = .pred_class, y = proportion, fill = group)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ `LogReg Approach`, ncol = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Predictions for highlighted receptor pair groups", x = "Prediction of shared specificity", y = "% of observations") +
    scale_fill_viridis_d(option = "magma", begin = 0.2, end = 0.9)

ggsave(paste0(out.path, "highlighted_pred_classes_by_features.png"), width = 10, height = 15)

# Distribution of .pred_yes and .pred_no
# By model and approach
ggplot(all.highlighted.preds %>%
    filter(Dataset == "CDR3") %>%
    mutate(`LogReg Approach` = paste0(Model, ", ", Approach)) %>%
    pivot_longer(cols = c(".pred_Yes", ".pred_No"), names_to = "Predicted class", values_to = "Probability"),
    aes(x = `Predicted class`, y = Probability, fill = group)) +
    geom_boxplot() +
    facet_wrap( ~ `LogReg Approach`) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Predicted probabilities for highlighted receptor pair groups", x = "Predicted class", y = "Probability") +
    scale_fill_viridis_d(option = "magma", begin = 0.2, end = 0.9)

ggsave(paste0(out.path, "highlighted_pred_probs.png"), width = 15, height = 10)

# By model, approach, and features
ggplot(all.highlighted.preds %>%
    filter(Dataset == "CDR3") %>%
    mutate(`LogReg Approach` = paste0(Model, ", ", Approach, " (", Feature, ")")) %>%
    pivot_longer(cols = c(".pred_Yes", ".pred_No"), names_to = "Predicted class", values_to = "Probability"),
    aes(x = `Predicted class`, y = Probability, fill = group)) +
    geom_boxplot() +
    facet_wrap( ~ `LogReg Approach`, ncol = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Predicted probabilities for highlighted receptor pair groups", x = "Predicted class", y = "Probability") +
    scale_fill_viridis_d(option = "magma", begin = 0.2, end = 0.9)

ggsave(paste0(out.path, "highlighted_pred_probs_by_features.png"), width = 15, height = 25)

# Density plots of class probabilities 
# By model and approach
ggplot(all.highlighted.preds %>%
    filter(Dataset == "CDR3") %>%
    mutate(`LogReg Approach` = paste0(Model, ", ", Approach)) %>%
    pivot_longer(cols = c(".pred_Yes", ".pred_No"), names_to = "Predicted class", values_to = "Probability"),
    aes(x = Probability, fill = group, linetype = `Predicted class`)) +
    geom_density(alpha = 0.5) +
    facet_wrap( ~ `LogReg Approach`) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Density of predicted probabilities for highlighted receptor pair groups", x = "Probability", y = "Density") +
    scale_fill_viridis_d(option = "magma", begin = 0.2, end = 0.9)

ggsave(paste0(out.path, "highlighted_pred_probs_density.png"), width = 15, height = 10)

# By model, approach, and features
ggplot(all.highlighted.preds %>%
    filter(Dataset == "CDR3") %>%
    mutate(`LogReg Approach` = paste0(Model, ", ", Approach, " (", Feature, ")")) %>%
    pivot_longer(cols = c(".pred_Yes", ".pred_No"), names_to = "Predicted class", values_to = "Probability"),
    aes(x = Probability, fill = group, linetype = `Predicted class`)) +
    geom_density(alpha = 0.5) +
    facet_wrap( ~ `LogReg Approach`, ncol = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Density of predicted probabilities for highlighted receptor pair groups", x = "Probability", y = "Density") +
    scale_fill_viridis_d(option = "magma", begin = 0.2, end = 0.9)

ggsave(paste0(out.path, "highlighted_pred_probs_density_by_features.png"), width = 15, height = 25)
