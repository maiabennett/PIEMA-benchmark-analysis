# ============================================================
# File: logistic-classifier-analysis.R
# Description: Script for analyzing the performance of all trained logistic regression classifier approaches
# Author: Maia Bennett-Boehm
# University of Nebraska at Omaha
# ============================================================

library(tidyverse)
library(tidymodels)
library(themis)
library(broom)
library(viridis)
library(ggrepel)
library(scales)

# Combine metrics
# One column with metrics
# One column with original dataset (CDR3 similarity (i.e., CDR3) or full sequence similarity (i.e., Full))
# One column with model approach (Model 1 = positive vs decoy, Model 2 = positive vs unlikely + decoy)
# One column with the feature type (Combined, Structure, or Sequence)
# One column with the feature selection approach (all, limited, or tune for combined and structure, both, CDR3, or full for sequence)
# One column with data inclusion (NLV, Yes or No)
# One column with the metric value
run.paths <- list(CDR3 = "./analysis/classifier/CDR3-similarity/", Full = "./analysis/classifier/full-similarity/")
# nlv.paths <- list(Yes = "with-NLV/", No = "without-NLV/")
epitope.paths <- list(All = "all/")
feature.paths <- list(Combined = "combined_features_", Structure = "structure_features_", Sequence = "sequence_features_")
model.paths <- list(`Model 1` = "model1_", `Model 2` = "model2_")

all.metrics <- data.frame()
limited.epitope.metrics <- data.frame()

for (run in names(run.paths)) {
    # for (nlv in names(nlv.paths)) {
    for (epitope in names(epitope.paths)) {
        for (feature in names(feature.paths)) {
            for (model in names(model.paths)) {
                metrics <- read.csv(paste0(run.paths[[run]], epitope.paths[[epitope]],
                    # nlv.paths[[nlv]], 
                    feature.paths[[feature]], model.paths[[model]], "test_metrics.csv")) %>%
                    mutate(Metrics = .metric, Dataset = run, Approach = feature, Model = model, Epitopes = epitope
                        # NLV = nlv
                        ) %>% 
                    select(-.metric) %>%
                    pivot_longer(cols = -c(Metrics, Dataset, Approach, Model, Epitopes
                        # NLV
                        ), names_to = "Feature", values_to = "Value") %>%
                    pivot_wider(names_from = Metrics, values_from = Value) 

                all.metrics <- rbind(all.metrics, metrics)
            }
        }
    }
}


        

# Combine per-epitope metrics
all.epitope.metrics <- data.frame()
limited.epitope.metrics <- data.frame()

for (run in names(run.paths)) {
    for (epitope in names(epitope.paths)) {
        for (feature in names(feature.paths)) {
            for (model in names(model.paths)) {
                metrics <- read.csv(paste0(run.paths[[run]], epitope.paths[[epitope]], feature.paths[[feature]], model.paths[[model]], "epitope_metrics.csv")) %>%
                    mutate(Metrics = .metric, Dataset = run, Approach = feature, Model = model, Epitopes = epitope) %>% 
                    select(-.metric) 
                epitope_col <- if ("epitope" %in% colnames(metrics)) "epitope" else "ref.epitope"
                metrics <- metrics %>%
                    pivot_longer(cols = -c(Metrics, Dataset, Approach, Model, !!sym(epitope_col), Epitopes), names_to = "Feature", values_to = "Value") %>%
                    pivot_wider(names_from = Metrics, values_from = Value) %>% 
                    dplyr::rename(Epitope = !!sym(epitope_col))

                all.epitope.metrics <- rbind(all.epitope.metrics, metrics)
            }
        }
    }
}


# Combine coefficients
all.coefficients.features <- data.frame()
all.coefficients.pca <- data.frame()

for (run in names(run.paths)) {
    for (epitope in names(epitope.paths)) {
        for (feature in names(feature.paths)) {
            for (model in names(model.paths)) {
                coefficients <- read.csv(paste0(run.paths[[run]], epitope.paths[[epitope]], feature.paths[[feature]], model.paths[[model]], "coefficients.csv")) %>%
                    mutate(Coefficient = term, Dataset = run, Approach = feature, Model = model, Epitopes = epitope) %>% 
                    select(-term) 
                p.val.cols <- grep("p.value|penalty", colnames(coefficients), value = TRUE)
                coefficients <- coefficients %>%
                    select(-all_of(p.val.cols)) %>%
                    pivot_longer(cols = -c(Coefficient, Dataset, Approach, Model, Epitopes), names_to = "Feature", values_to = "Value") %>%
                    pivot_wider(names_from = Coefficient, values_from = Value, names_sep = "_") 

                all.coefficients.features <- bind_rows(all.coefficients.features, coefficients)
            }
        }
    }
}

# Save results
out.path <- "./analysis/classifier/"
write.csv(all.metrics, paste0(out.path, "all_metrics.csv"), row.names = FALSE)
write.csv(all.epitope.metrics, paste0(out.path, "all_epitope_metrics.csv"), row.names = FALSE)
# write.csv(limited.epitope.metrics, paste0(out.path, "limited_epitope_metrics.csv"), row.names = FALSE)
write.csv(all.coefficients.features, paste0(out.path, "all_coefficients_features.csv"), row.names = FALSE)
# write.csv(all.coefficients.pca, paste0(out.path, "all_coefficients_pca.csv"), row.names = FALSE)


# Plot results to identify best approaches, feature influence, etc
# Plotting accuracy (y axis), with bars colored by the values of approach + number of features, with a x axis being a combination of dataset + inclusion of NLV
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        #mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        mutate(Epitopes = ifelse(Epitopes == "All", "All Epitopes", "Limited Epitopes")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            #Approach == "Combined" & Feature == "limited.features" ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, Epitopes)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset and approach", subtitle = "Model 1: Positive vs. Decoy pairs", y = "Accuracy") +
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.45,0.75),oob = rescale_none)

ggsave(paste0(out.path, "accuracy_model1.png"), width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        #mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        mutate(Epitopes = ifelse(Epitopes == "All", "All Epitopes", "Limited Epitopes")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            #Approach == "Combined" & Feature == "limited.features" ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, Epitopes)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset and approach", subtitle = "Model 2: Positive vs. Unlikely and Decoy pairs", y = "Accuracy") +
    scale_fill_viridis_d(option = "viridis", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.45,0.75),oob = rescale_none)

ggsave(paste0(out.path, "accuracy_model2.png"), width = 20, height = 15, dpi = 300)

# Plotting ROC_AUC
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        #mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        mutate(Epitopes = ifelse(Epitopes == "All", "All Epitopes", "Limited Epitopes")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            #Approach == "Combined" & Feature == "limited.features" ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, Epitopes)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC_AUC by dataset and approach", subtitle = "Model 1: Positive vs. Decoy pairs", y = "ROC_AUC") +
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.45,0.7),oob = rescale_none)

ggsave(paste0(out.path, "roc_auc_model1.png"), width = 20, height = 10, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        #mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        mutate(Epitopes = ifelse(Epitopes == "All", "All Epitopes", "Limited Epitopes")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            #Approach == "Combined" & Feature == "limited.features" ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, Epitopes)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC_AUC by dataset and approach", subtitle = "Model 2: Positive vs. Unlikely and Decoy pairs", y = "ROC_AUC") +
    scale_fill_viridis_d(option = "viridis", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.45,0.7),oob = rescale_none)

ggsave(paste0(out.path, "roc_auc_model2.png"), width = 20, height = 10, dpi = 300)

# Plotting PR_AUC
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        #mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        mutate(Epitopes = ifelse(Epitopes == "All", "All Epitopes", "Limited Epitopes")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            #Approach == "Combined" & Feature == "limited.features" ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, Epitopes)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR_AUC by dataset and approach", subtitle = "Model 1: Positive vs. Decoy pairs", y = "PR_AUC") +
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.1,0.6),oob = rescale_none)

ggsave(paste0(out.path, "pr_auc_model1.png"), width = 20, height = 10, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        #mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        mutate(Epitopes = ifelse(Epitopes == "All", "All Epitopes", "Limited Epitopes")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            #Approach == "Combined" & Feature == "limited.features" ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, Epitopes)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR_AUC by dataset and approach", subtitle = "Model 2: Positive vs. Unlikely and Decoy pairs", y = "PR_AUC") +
    scale_fill_viridis_d(option = "viridis", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.1,0.6),oob = rescale_none)

ggsave(paste0(out.path, "pr_auc_model2.png"), width = 20, height = 10, dpi = 300)

# Plotting recall
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        #mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        mutate(Epitopes = ifelse(Epitopes == "All", "All Epitopes", "Limited Epitopes")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            #Approach == "Combined" & Feature == "limited.features" ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, Epitopes)),
    aes(x = Approach, y = recall, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Recall by dataset and approach", subtitle = "Model 1: Positive vs. Decoy pairs", y = "Recall") +
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.35,0.6),oob = rescale_none)

ggsave(paste0(out.path, "recall_model1.png"), width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        #mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        mutate(Epitopes = ifelse(Epitopes == "All", "All Epitopes", "Limited Epitopes")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            #Approach == "Combined" & Feature == "limited.features" ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, Epitopes)),
    aes(x = Approach, y = recall, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Recall by dataset and approach", subtitle = "Model 2: Positive vs. Unlikely and Decoy pairs", y = "Recall") +
    scale_fill_viridis_d(option = "viridis", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.35,0.6),oob = rescale_none)

ggsave(paste0(out.path, "recall_model2.png"), width = 20, height = 15, dpi = 300)

# Precision
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        #mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        mutate(Epitopes = ifelse(Epitopes == "All", "All Epitopes", "Limited Epitopes")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            #Approach == "Combined" & Feature == "limited.features" ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, Epitopes)),
    aes(x = Approach, y = precision, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Precision by dataset and approach", subtitle = "Model 1: Positive vs. Decoy pairs", y = "Precision") +
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.05,0.5),oob = rescale_none)

ggsave(paste0(out.path, "precision_model1.png"), width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        #mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        mutate(Epitopes = ifelse(Epitopes == "All", "All Epitopes", "Limited Epitopes")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            #Approach == "Combined" & Feature == "limited.features" ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, Epitopes)),
    aes(x = Approach, y = precision, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Precision by dataset and approach", subtitle = "Model 2: Positive vs. Unlikely and Decoy pairs", y = "Precision") +
    scale_fill_viridis_d(option = "viridis", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.05,0.5),oob = rescale_none)

ggsave(paste0(out.path, "precision_model2.png"), width = 20, height = 15, dpi = 300)


# Feature coefficients
# Convert to long format for plotting, and plot features on X, Coefficients on Y, colored by Approach + Feature
# Combined feature sets
# These still need more work, but not necessary rn
# library(ggpattern)

# Model 1
ggplot(
    all.coefficients.features %>% 
        filter(Model == "Model 1") %>%
        pivot_longer(cols = -c(Dataset, Approach, Model, Epitopes, Feature), names_to = "Features", values_to = "Coefficient") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(Epitopes = ifelse(Epitopes == "All", "All Epitopes", "Limited Epitopes")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, Epitopes)),
    aes(x = Features, y = Coefficient, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Feature coefficients by dataset and approach", subtitle = "Model 1: Positive vs Decoy pairs", y = "Coefficient") +
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = .9)

ggsave(paste0(out.path, "feature_coefficients_model1.png"), width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.coefficients.features %>% 
        filter(Model == "Model 2") %>%
        pivot_longer(cols = -c(Dataset, Approach, Model, Epitopes, Feature), names_to = "Features", values_to = "Coefficient") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(Epitopes = ifelse(Epitopes == "All", "All Epitopes", "Limited Epitopes")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, Epitopes)),
    aes(x = Features, y = Coefficient, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Feature coefficients by dataset and approach", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "Coefficient") +
    scale_fill_viridis_d(option = "viridis", begin = 0.1, end = .9)

ggsave(paste0(out.path, "feature_coefficients_model2.png"), width = 20, height = 15, dpi = 300)


# Per-epitope performance metrics
# Accuracy
# Model 1
ggplot(
    all.epitope.metrics %>%
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        # filter(Approach == "Structure" | 
        #     (Approach == "Sequence" & Feature == "cdr3") | 
        #     (Approach == "Combined" & Feature == "limited.features")) %>%
        # #mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        mutate(Epitopes = ifelse(Epitopes == "All", "All Epitopes", "Limited Epitopes")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, Epitopes)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    #facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    facet_wrap( ~ Epitope, nrow = 2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset and epitope", subtitle = "Model 1: Positive vs Decoy pairs", y = "Accuracy") +
    #scale_fill_manual(values = model1.approach.colors) +
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.45,0.8),oob = rescale_none)

ggsave(paste0(out.path, "accuracy_by_epitope_model1.png"), width = 20, height = 10, dpi = 300)


# Model 2
ggplot(
    all.epitope.metrics %>%
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        # filter(Approach == "Structure" | 
        #     (Approach == "Sequence" & Feature == "cdr3") | 
        #     (Approach == "Combined" & Feature == "limited.features")) %>%
        # #mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        mutate(Epitopes = ifelse(Epitopes == "All", "All Epitopes", "Limited Epitopes")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, Epitopes)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    facet_wrap( ~ Epitope, nrow = 2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset and epitope", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "Accuracy") +
    #scale_fill_manual(values = model2.approach.colors) +
    scale_fill_viridis_d(option = "viridis", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.45,0.8),oob = rescale_none)

ggsave(paste0(out.path, "accuracy_by_epitope_model2.png"), width = 20, height = 10, dpi = 300)


# ROC AUC
# Model 1
# With NLV
ggplot(
    all.epitope.metrics %>%
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        # filter(Approach == "Structure" | 
        #     (Approach == "Sequence" & Feature == "cdr3") | 
        #     (Approach == "Combined" & Feature == "limited.features")) %>%
        # #mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        mutate(Epitopes = ifelse(Epitopes == "All", "All Epitopes", "Limited Epitopes")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, Epitopes)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    facet_wrap( ~ Epitope, nrow = 2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC AUC by dataset and epitope", subtitle = "Model 1: Positive vs Decoy pairs", y = "ROC AUC") +
    #scale_fill_manual(values = model1.approach.colors) +
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.45,0.9),oob = rescale_none)

ggsave(paste0(out.path, "roc_auc_by_epitope_model1.png"), width = 20, height = 10, dpi = 300)


# Model 2
ggplot(
    all.epitope.metrics %>%
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        # filter(Approach == "Structure" | 
        #     (Approach == "Sequence" & Feature == "cdr3") | 
        #     (Approach == "Combined" & Feature == "limited.features")) %>%
        # #mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        mutate(Epitopes = ifelse(Epitopes == "All", "All Epitopes", "Limited Epitopes")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, Epitopes)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    facet_wrap( ~ Epitope, nrow = 2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC AUC by dataset and epitope", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "ROC AUC") +
    #scale_fill_manual(values = model2.approach.colors) +
    scale_fill_viridis_d(option = "viridis", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.45,0.9),oob = rescale_none)

ggsave(paste0(out.path, "roc_auc_by_epitope_model2.png"), width = 20, height = 10, dpi = 300)


# PR AUC
# Model 1 
ggplot(
    all.epitope.metrics %>%
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        # filter(Approach == "Structure" | 
        #     (Approach == "Sequence" & Feature == "cdr3") | 
        #     (Approach == "Combined" & Feature == "limited.features")) %>%
        # #mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        mutate(Epitopes = ifelse(Epitopes == "All", "All Epitopes", "Limited Epitopes")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, Epitopes)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    facet_wrap( ~ Epitope, nrow = 2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR AUC by dataset and epitope", subtitle = "Model 1: Positive vs Decoy pairs", y = "PR AUC") +
    #scale_fill_manual(values = model1.approach.colors) +
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.05,0.8),oob = rescale_none)

ggsave(paste0(out.path, "pr_auc_by_epitope_model1.png"), width = 20, height = 10, dpi = 300)


# Model 2
ggplot(
    all.epitope.metrics %>%
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        # filter(Approach == "Structure" | 
        #     (Approach == "Sequence" & Feature == "cdr3") | 
        #     (Approach == "Combined" & Feature == "limited.features")) %>%
        # #mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        mutate(Epitopes = ifelse(Epitopes == "All", "All Epitopes", "Limited Epitopes")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, Epitopes)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    facet_wrap( ~ Epitope, nrow = 2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR AUC by dataset and epitope", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "PR AUC") +
    #scale_fill_manual(values = model2.approach.colors) +
    scale_fill_viridis_d(option = "viridis", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.05,0.8),oob = rescale_none)

ggsave(paste0(out.path, "pr_auc_by_epitope_model2.png"), width = 20, height = 10, dpi = 300)


# Import data for ROC AUC and PR AUC curves
model1.preds <- data.frame()
model2.preds <- data.frame()

sequence.features <- list(`CDR3 only` = "CDR3", `Full seq. only` = "full_seq", `CDR3 + full seq.` = "both_features")
struct.combined.features <- list(`All features` = "all_features", `PCA features` = "pca_features", `Tuned features` = "best_tune")

for (run in names(run.paths)) {
    for (epitope in names(epitope.paths)) {
        for (feature in names(feature.paths)) {
            if (feature == "Sequence") {
                for (sequence.feature in names(sequence.features)) {
                    preds1 <- read.csv(paste0(run.paths[[run]], epitope.paths[[epitope]], feature.paths[[feature]], "model1_test_pred_", sequence.features[[sequence.feature]], ".csv")) %>% 
                        mutate(pca1 = NA, pca2 = NA, Dataset = run, Approach = feature, Model = "Model 1", Epitopes = epitope, Feature = sequence.feature)
                    preds2 <- read.csv(paste0(run.paths[[run]], epitope.paths[[epitope]], feature.paths[[feature]], "model2_test_pred_", sequence.features[[sequence.feature]], ".csv")) %>% 
                        mutate(pca1 = NA, pca2 = NA, Dataset = run, Approach = feature, Model = "Model 2", Epitopes = epitope, Feature = sequence.feature)
                    model1.preds <- rbind(model1.preds, preds1)
                    model2.preds <- rbind(model2.preds, preds2) 
                }
            }
            else {
                for (struct.combined.feature in names(struct.combined.features)) {
                    preds1 <- read.csv(paste0(run.paths[[run]], epitope.paths[[epitope]], feature.paths[[feature]], "model1_test_pred_", struct.combined.features[[struct.combined.feature]], ".csv")) %>% 
                        mutate(Dataset = run, Approach = feature, Model = "Model 1", Epitopes = epitope, Feature = struct.combined.feature)
                    preds2 <- read.csv(paste0(run.paths[[run]], epitope.paths[[epitope]], feature.paths[[feature]], "model2_test_pred_", struct.combined.features[[struct.combined.feature]], ".csv")) %>% 
                        mutate(Dataset = run, Approach = feature, Model = "Model 2", Epitopes = epitope, Feature = struct.combined.feature)
                    model1.preds <- rbind(model1.preds, preds1) 
                    model2.preds <- rbind(model2.preds, preds2)
                }
            }   
        }
    }
}

model1.preds <- model1.preds %>% 
    mutate(shared.specificity = factor(shared.specificity, levels = c("No", "Yes")))
model2.preds <- model2.preds %>%
    mutate(shared.specificity = factor(shared.specificity, levels = c("No", "Yes")))

write.csv(model1.preds, paste0(out.path, "all_model1_preds.csv"), row.names = FALSE)
write.csv(model2.preds, paste0(out.path, "all_model2_preds.csv"), row.names = FALSE)

# Model 1
# ROC AUC curves
# Color by LogReg approach, linetype by Feature
model1.preds %>% 
    mutate(`LogReg Approach` = paste(Approach, "approach,", Feature)) %>% 
    group_by(`LogReg Approach`, Approach) %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>% 
    ggplot(aes(x = 1 - specificity, y = sensitivity, color = `LogReg Approach`, linetype = Approach)) +
    geom_path(lwd = 1.5, alpha = 0.6) +
    geom_abline(lty = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = 0.9) +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted")) +
    labs(title = "ROC curve, Model 1: Positive vs Decoy pairs", x = "1 - Specificity", y = "Sensitivity")

ggsave(paste0(out.path, "roc_auc_curves_model1.png"), width = 20, height = 10, dpi = 300)

# ROC AUC curves per epitope
model1.preds %>% 
    mutate(`LogReg Approach` = paste(Approach, "approach,", Feature)) %>% 
    group_by(`LogReg Approach`, Approach, ref.epitope) %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, color = `LogReg Approach`, linetype = Approach)) +
    geom_path(lwd = 1.5, alpha = 0.6) +
    geom_abline(lty = 3) +
    coord_equal() +
    facet_wrap(~ ref.epitope, ncol = 4) +
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = 0.9) +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted")) +
    labs(title = "ROC curve, Model 1: Positive vs Decoy pairs", x = "1 - Specificity", y = "Sensitivity")

ggsave(paste0(out.path, "roc_auc_curves_model1_per_epitope.png"), width = 20, height = 15, dpi = 300)


# PR AUC curves
model1.preds %>%
    mutate(`LogReg Approach` = paste(Approach, "approach,", Feature)) %>% 
    group_by(`LogReg Approach`, Approach) %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>% 
    ggplot(aes(x = recall, y = precision, color = `LogReg Approach`, linetype = Approach)) +
    geom_path(lwd = 1.5, alpha = 0.6) +
    coord_equal() +
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = 0.9) +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted")) +
    labs(title = "PR curve, Model 1: Positive vs Decoy pairs", x = "Recall", y = "Precision")

ggsave(paste0(out.path, "pr_auc_curves_model1.png"), width = 20, height = 10, dpi = 300)

# PR AUC curves per epitope
model1.preds %>%
    mutate(`LogReg Approach` = paste(Approach, "approach,", Feature)) %>% 
    group_by(`LogReg Approach`, Approach, ref.epitope) %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    ggplot(aes(x = recall, y = precision, color = `LogReg Approach`, linetype = Approach)) +
    geom_path(lwd = 1.5, alpha = 0.6) +
    coord_equal() +
    facet_wrap(~ ref.epitope, ncol = 4) +
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = 0.9) +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted")) +
    labs(title = "PR curve, Model 1: Positive vs Decoy pairs", x = "Recall", y = "Precision")

ggsave(paste0(out.path, "pr_auc_curves_model1_per_epitope.png"), width = 20, height = 15, dpi = 300)

# Model 2
# ROC AUC curves
model2.preds %>% 
    mutate(`LogReg Approach` = paste(Approach, "approach,", Feature)) %>% 
    group_by(`LogReg Approach`, Approach) %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>% 
    ggplot(aes(x = 1 - specificity, y = sensitivity, color = `LogReg Approach`, linetype = Approach)) +
    geom_path(lwd = 1.5, alpha = 0.6) +
    geom_abline(lty = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "viridis", begin = 0.1, end = 0.9) +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted")) +
    labs(title = "ROC curve, Model 2: Positve vs Unlikely and Decoy pairs", x = "1 - Specificity", y = "Sensitivity")

ggsave(paste0(out.path, "roc_auc_curves_model2.png"), width = 20, height = 10, dpi = 300)

# ROC AUC curves per epitope
model2.preds %>% 
    mutate(`LogReg Approach` = paste(Approach, "approach,", Feature)) %>%
    # mutate(Epitope = ref.epitope) %>%
    #     bind_rows(model2.preds %>%
    #         filter(ref.epitope != samp.epitope) %>%
    #         mutate(Epitope = samp.epitope)) %>%
    group_by(`LogReg Approach`, Approach, ref.epitope) %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, color = `LogReg Approach`, linetype = Approach)) +
    geom_path(lwd = 1.5, alpha = 0.6) +
    geom_abline(lty = 3) +
    coord_equal() +
    facet_wrap(~ ref.epitope, ncol = 4) +
    scale_color_viridis_d(option = "viridis", begin = 0.1, end = 0.9) +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted")) +
    labs(title = "ROC curve, Model 2: Positive vs Unlikely and Decoy pairs", x = "1 - Specificity", y = "Sensitivity")

ggsave(paste0(out.path, "roc_auc_curves_model2_per_epitope.png"), width = 20, height = 15, dpi = 300)

# PR AUC curves
model2.preds %>%
    mutate(`LogReg Approach` = paste(Approach, "approach,", Feature)) %>% 
    group_by(`LogReg Approach`, Approach) %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>% 
    ggplot(aes(x = recall, y = precision, color = `LogReg Approach`, linetype = Approach)) +
    geom_path(lwd = 1.5, alpha = 0.6) +
    coord_equal() +
    scale_color_viridis_d(option = "viridis", begin = 0.1, end = 0.9) +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted")) +
    labs(title = "PR curve, Model 2: Positive vs Unlikely and Decoy pairs", x = "Recall", y = "Precision")

ggsave(paste0(out.path, "pr_auc_curves_model2.png"), width = 20, height = 10, dpi = 300)

# PR AUC curves per epitope
model2.preds %>%
    mutate(`LogReg Approach` = paste(Approach, "approach,", Feature)) %>% 
    group_by(`LogReg Approach`, Approach, ref.epitope) %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    ggplot(aes(x = recall, y = precision, color = `LogReg Approach`, linetype = Approach)) +
    geom_path(lwd = 1.5, alpha = 0.6) +
    coord_equal() +
    facet_wrap(~ ref.epitope, ncol = 4) +
    scale_color_viridis_d(option = "viridis", begin = 0.1, end = 0.9) +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted")) +
    labs(title = "PR curve, Model 2: Positive vs Unlikely and Decoy pairs", x = "Recall", y = "Precision")

ggsave(paste0(out.path, "pr_auc_curves_model2_per_epitope.png"), width = 20, height = 15, dpi = 300)
