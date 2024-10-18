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
run.paths <- list(CDR3 = "./analysis/classifier/without-cross-reactives/run1/", Full = "./analysis/classifier/without-cross-reactives/run2/")
nlv.paths <- list(Yes = "with-NLV/", No = "without-NLV/")
feature.paths <- list(Combined = "combined_features_", Structure = "structure_features_", Sequence = "sequence_features_")
model.paths <- list(`Model 1` = "model1_", `Model 2` = "model2_")

all.metrics <- data.frame()

for (run in names(run.paths)) {
    for (nlv in names(nlv.paths)) {
        for (feature in names(feature.paths)) {
            for (model in names(model.paths)) {
                metrics <- read.csv(paste0(run.paths[[run]], nlv.paths[[nlv]], feature.paths[[feature]], model.paths[[model]], "test_metrics.csv")) %>%
                    mutate(Metrics = .metric, Dataset = run, Approach = feature, Model = model, NLV = nlv) %>% 
                    select(-.metric) %>%
                    pivot_longer(cols = -c(Metrics, Dataset, Approach, Model, NLV), names_to = "Feature", values_to = "Value") %>%
                    pivot_wider(names_from = Metrics, values_from = Value) 

                all.metrics <- rbind(all.metrics, metrics)
            }
        }
    }
}


        

# Combine per-epitope metrics- Need to figure out how to represent per-epitope metrics first
all.epitope.metrics.with.nlv <- data.frame()
all.epitope.metrics.without.nlv <- data.frame()

for (run in names(run.paths)) {
    for (feature in names(feature.paths)) {
        for (model in names(model.paths)) {
            metrics <- read.csv(paste0(run.paths[[run]], nlv.paths[["Yes"]], feature.paths[[feature]], model.paths[[model]], "epitope_metrics.csv")) %>%
                mutate(Metrics = .metric, Dataset = run, Approach = feature, Model = model, NLV = "Yes") %>% 
                select(-.metric) 
            epitope_col <- if ("epitope" %in% colnames(metrics)) "epitope" else "ref.epitope"
            metrics <- metrics %>%
                pivot_longer(cols = -c(Metrics, Dataset, Approach, Model, !!sym(epitope_col), NLV), names_to = "Feature", values_to = "Value") %>%
                pivot_wider(names_from = Metrics, values_from = Value) %>% 
                dplyr::rename(Epitope = !!sym(epitope_col))

            all.epitope.metrics.with.nlv <- rbind(all.epitope.metrics.with.nlv, metrics)

            metrics <- read.csv(paste0(run.paths[[run]], nlv.paths[["No"]], feature.paths[[feature]], model.paths[[model]], "epitope_metrics.csv")) %>%
                mutate(Metrics = .metric, Dataset = run, Approach = feature, Model = model, NLV = "No") %>% 
                select(-.metric) 
            epitope_col <- if ("epitope" %in% colnames(metrics)) "epitope" else "ref.epitope"
            metrics <- 
            metrics %>%
                pivot_longer(cols = -c(Metrics, Dataset, Approach, Model, !!sym(epitope_col), NLV), names_to = "Feature", values_to = "Value") %>%
                pivot_wider(names_from = Metrics, values_from = Value) %>% 
                dplyr::rename(Epitope = !!sym(epitope_col))

            all.epitope.metrics.without.nlv <- rbind(all.epitope.metrics.without.nlv, metrics)
        }
    }
}



# Combine coefficients
all.coefficients.with.nlv <- data.frame()
all.coefficients.without.nlv <- data.frame()

for (run in names(run.paths)) {
    for (feature in names(feature.paths)) {
        for (model in names(model.paths)) {
            coefficients <- read.csv(paste0(run.paths[[run]], nlv.paths[["Yes"]], feature.paths[[feature]], model.paths[[model]], "coefficients.csv")) %>%
                mutate(Coefficient = term, Dataset = run, Approach = feature, Model = model, NLV = "Yes") %>% 
                select(-term) 
            p.val.cols <- grep("p.value|penalty", colnames(coefficients), value = TRUE)
            coefficients <- coefficients %>%
                select(-all_of(p.val.cols)) %>%
                pivot_longer(cols = -c(Coefficient, Dataset, Approach, Model, NLV), names_to = "Feature", values_to = "Value") %>%
                pivot_wider(names_from = Coefficient, values_from = Value, names_sep = "_") 

            all.coefficients.with.nlv <- bind_rows(all.coefficients.with.nlv, coefficients)

            coefficients <- read.csv(paste0(run.paths[[run]], nlv.paths[["No"]], feature.paths[[feature]], model.paths[[model]], "coefficients.csv")) %>%
                mutate(Coefficient = term, Dataset = run, Approach = feature, Model = model, NLV = "No") %>% 
                select(-term)
            p.val.cols <- grep("p.value|penalty", colnames(coefficients), value = TRUE)
            coefficients <- coefficients %>%
                select(-all_of(p.val.cols)) %>%
                pivot_longer(cols = -c(Coefficient, Dataset, Approach, Model, NLV), names_to = "Feature", values_to = "Value") %>%
                pivot_wider(names_from = Coefficient, values_from = Value, names_sep = "_") 

            all.coefficients.without.nlv <- bind_rows(all.coefficients.without.nlv, coefficients)
        }
    }
}

# Save results
out.path <- "./analysis/classifier/without-cross-reactives/"
write.csv(all.metrics, paste0(out.path, "all_metrics.csv"), row.names = FALSE)
write.csv(all.epitope.metrics.with.nlv, paste0(out.path, "all_epitope_metrics_with_nlv.csv"), row.names = FALSE)
write.csv(all.epitope.metrics.without.nlv, paste0(out.path, "all_epitope_metrics_without_nlv.csv"), row.names = FALSE)
write.csv(all.coefficients.with.nlv, paste0(out.path, "all_coefficients_with_nlv.csv"), row.names = FALSE)
write.csv(all.coefficients.without.nlv, paste0(out.path, "all_coefficients_without_nlv.csv"), row.names = FALSE)

#Try only plotting non-NLV performance

# Plot results to identify best approaches, feature influence, etc
# Plotting accuracy (y axis), with bars colored by the values of approach + number of features, with a x axis being a combination of dataset + inclusion of NLV
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            Approach == "Combined" & Feature == "limited.features" ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset and approach", subtitle = "Model 1: Positive vs. Decoy pairs", y = "Accuracy") +
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "accuracy_by_dataset_and_approach_model1.png"), width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            Approach == "Combined" & Feature == "limited.features" ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset and approach", subtitle = "Model 2: Positive vs. Unlikely and Decoy pairs", y = "Accuracy") +
    scale_fill_viridis_d(option = "viridis", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "accuracy_by_dataset_and_approach_model2.png"), width = 20, height = 15, dpi = 300)

# Plotting ROC_AUC
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            Approach == "Combined" & Feature == "limited.features" ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC_AUC by dataset and approach", subtitle = "Model 1: Positive vs. Decoy pairs", y = "ROC_AUC") +
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "roc_auc_by_dataset_and_approach_model1.png"), width = 20, height = 10, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            Approach == "Combined" & Feature == "limited.features" ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC_AUC by dataset and approach", subtitle = "Model 2: Positive vs. Unlikely and Decoy pairs", y = "ROC_AUC") +
    scale_fill_viridis_d(option = "viridis", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "roc_auc_by_dataset_and_approach_model2.png"), width = 20, height = 10, dpi = 300)

# Plotting PR_AUC
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            Approach == "Combined" & Feature == "limited.features" ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR_AUC by dataset and approach", subtitle = "Model 1: Positive vs. Decoy pairs", y = "PR_AUC") +
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.15,0.6),oob = rescale_none)

ggsave(paste0(out.path, "pr_auc_by_dataset_and_approach_model1.png"), width = 20, height = 10, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            Approach == "Combined" & Feature == "limited.features" ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR_AUC by dataset and approach", subtitle = "Model 2: Positive vs. Unlikely and Decoy pairs", y = "PR_AUC") +
    scale_fill_viridis_d(option = "viridis", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.15,0.6),oob = rescale_none)

ggsave(paste0(out.path, "pr_auc_by_dataset_and_approach_model2.png"), width = 20, height = 10, dpi = 300)

# Plotting recall
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            Approach == "Combined" & Feature == "limited.features" ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = recall, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Recall by dataset and approach", subtitle = "Model 1: Positive vs. Decoy pairs", y = "Recall") +
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.35,0.6),oob = rescale_none)

ggsave(paste0(out.path, "recall_by_dataset_and_approach_model1.png"), width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            Approach == "Combined" & Feature == "limited.features" ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = recall, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Recall by dataset and approach", subtitle = "Model 2: Positive vs. Unlikely and Decoy pairs", y = "Recall") +
    scale_fill_viridis_d(option = "viridis", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.35,0.6),oob = rescale_none)

ggsave(paste0(out.path, "recall_by_dataset_and_approach_model2.png"), width = 20, height = 15, dpi = 300)

# Precision
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            Approach == "Combined" & Feature == "limited.features" ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = precision, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Precision by dataset and approach", subtitle = "Model 1: Positive vs. Decoy pairs", y = "Precision") +
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.15,0.5),oob = rescale_none)

ggsave(paste0(out.path, "precision_by_dataset_and_approach_model1.png"), width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            Approach == "Combined" & Feature == "limited.features" ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = precision, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Precision by dataset and approach", subtitle = "Model 2: Positive vs. Unlikely and Decoy pairs", y = "Precision") +
    scale_fill_viridis_d(option = "viridis", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.15,0.5),oob = rescale_none)

ggsave(paste0(out.path, "precision_by_dataset_and_approach_model2.png"), width = 20, height = 15, dpi = 300)

# Plotting the same metrics based on the inclusion of full sequence similarity in the model
# Includes all structural methods (for reference), the relevant sequence method, and the relevant combined method(s)
# Set colors for consistency in plots excluding certain methods
model.approaches <- unique(all.metrics %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        arrange(`LogReg Approach`) %>%
        pull(`LogReg Approach`)) 
model1.approach.colors <- viridis(length(model.approaches), option = "plasma", begin = 0.1, end = .9)
model1.approach.colors <- setNames(model1.approach.colors, model.approaches)
model2.approach.colors <- viridis(length(model.approaches), option = "viridis", begin = 0.1, end = .9)
model2.approach.colors <- setNames(model2.approach.colors, model.approaches)

# Accuracy 
# With full sequence
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & (Feature == "full" | Feature == "both")) | 
            (Approach == "Combined" & (Feature == "all.features" | Feature == "tuned.features"))) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset for approaches including full sequence similarity", subtitle = "Model 1: Positive vs Decoy pairs", y = "Accuracy") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "accuracy_by_dataset_full_sequence_model1.png"), width = 20, height = 15, dpi = 300)

# Without NLV 
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        filter(NLV == "No") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & (Feature == "full" | Feature == "both")) | 
            (Approach == "Combined" & (Feature == "all.features" | Feature == "tuned.features"))) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset for approaches including full sequence similarity, without NLV", subtitle = "Model 1: Positive vs Decoy pairs", y = "Accuracy") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "accuracy_by_dataset_full_sequence_model1_without_NLV.png"), width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & (Feature == "full" | Feature == "both")) | 
            (Approach == "Combined" & (Feature == "all.features" | Feature == "tuned.features"))) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset for approaches including full sequence similarity", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "Accuracy") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "accuracy_by_dataset_full_sequence_model2.png"), width = 20, height = 15, dpi = 300)

# Without NLV
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        filter(NLV == "No") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & (Feature == "full" | Feature == "both")) | 
            (Approach == "Combined" & (Feature == "all.features" | Feature == "tuned.features"))) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset for approaches including full sequence similarity, without NLV", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "Accuracy") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "accuracy_by_dataset_full_sequence_model2_without_NLV.png"), width = 20, height = 15, dpi = 300)

# Without full sequence
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Feature == "cdr3") | 
            (Approach == "Combined" & Feature == "limited.features")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset for approaches excluding full sequence similarity", subtitle = "Model 1: Positive vs Decoy pairs", y = "Accuracy") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "accuracy_by_dataset_cdr3_sequence_model1.png"), width = 20, height = 15, dpi = 300)

# Without NLV  
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        filter(NLV == "No") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Feature == "cdr3") | 
            (Approach == "Combined" & Feature == "limited.features")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset for approaches excluding full sequence similarity, without NLV", subtitle = "Model 1: Positive vs Decoy pairs", y = "Accuracy") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "accuracy_by_dataset_cdr3_sequence_model1_without_NLV.png"), width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Feature == "cdr3") | 
            (Approach == "Combined" & Feature == "limited.features")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset for approaches excluding full sequence similarity", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "Accuracy") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "accuracy_by_dataset_cdr3_sequence_model2.png"), width = 20, height = 15, dpi = 300)

# Without NLV
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        filter(NLV == "No") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Feature == "cdr3") | 
            (Approach == "Combined" & Feature == "limited.features")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset for approaches excluding full sequence similarity, without NLV", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "Accuracy") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "accuracy_by_dataset_cdr3_sequence_model2_without_NLV.png"), width = 20, height = 15, dpi = 300)

# ROC_AUC
# With full sequence
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & (Feature == "full" | Feature == "both")) | 
            (Approach == "Combined" & (Feature == "all.features" | Feature == "tuned.features"))) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC_AUC by dataset for approaches including full sequence similarity", subtitle = "Model 1: Positive vs Decoy pairs", y = "ROC_AUC") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "roc_auc_by_dataset_full_sequence_model1.png"), width = 20, height = 15, dpi = 300)

# Without NLV
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        filter(NLV == "No") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & (Feature == "full" | Feature == "both")) | 
            (Approach == "Combined" & (Feature == "all.features" | Feature == "tuned.features"))) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC_AUC by dataset for approaches including full sequence similarity, without NLV", subtitle = "Model 1: Positive vs Decoy pairs", y = "ROC_AUC") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "roc_auc_by_dataset_full_sequence_model1_without_NLV.png"), width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & (Feature == "full" | Feature == "both")) | 
            (Approach == "Combined" & (Feature == "all.features" | Feature == "tuned.features"))) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC_AUC by dataset for approaches including full sequence similarity", subtitle = "Model 1: Positive vs Decoy pairs", y = "ROC_AUC") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "roc_auc_by_dataset_full_sequence_model2.png"), width = 20, height = 15, dpi = 300)

# Without NLV
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        filter(NLV == "No") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & (Feature == "full" | Feature == "both")) | 
            (Approach == "Combined" & (Feature == "all.features" | Feature == "tuned.features"))) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC_AUC by dataset for approaches including full sequence similarity, without NLV", subtitle = "Model 1: Positive vs Decoy pairs", y = "ROC_AUC") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "roc_auc_by_dataset_full_sequence_model2_without_NLV.png"), width = 20, height = 15, dpi = 300)

# Without full sequence
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Feature == "cdr3") | 
            (Approach == "Combined" & Feature == "limited.features")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC_AUC by dataset for approaches excluding full sequence similarity", subtitle = "Model 1: Positive vs Decoy pairs", y = "ROC_AUC") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "roc_auc_by_dataset_cdr3_sequence_model1.png"), width = 20, height = 15, dpi = 300)

# Without NLV
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        filter(NLV == "No") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Feature == "cdr3") | 
            (Approach == "Combined" & Feature == "limited.features")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC_AUC by dataset for approaches excluding full sequence similarity, without NLV", subtitle = "Model 1: Positive vs Decoy pairs", y = "ROC_AUC") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "roc_auc_by_dataset_cdr3_sequence_model1_without_NLV.png"), width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Feature == "cdr3") | 
            (Approach == "Combined" & Feature == "limited.features")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC_AUC by dataset for approaches excluding full sequence similarity", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "ROC_AUC") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "roc_auc_by_dataset_cdr3_sequence_model2.png"), width = 20, height = 15, dpi = 300)

# Without NLV
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        filter(NLV == "No") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Feature == "cdr3") | 
            (Approach == "Combined" & Feature == "limited.features")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC_AUC by dataset for approaches excluding full sequence similarity, without NLV", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "ROC_AUC") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "roc_auc_by_dataset_cdr3_sequence_model2_without_NLV.png"), width = 20, height = 15, dpi = 300)

# PR_AUC
# With full sequence
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & (Feature == "full" | Feature == "both")) | 
            (Approach == "Combined" & (Feature == "all.features" | Feature == "tuned.features"))) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR_AUC by dataset for approaches including full sequence similarity", subtitle = "Model 1: Positive vs Decoy pairs", y = "PR_AUC") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.15,0.6),oob = rescale_none)

ggsave(paste0(out.path, "pr_auc_by_dataset_full_sequence_model1.png"), width = 20, height = 15, dpi = 300)

# Without NLV
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        filter(NLV == "No") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & (Feature == "full" | Feature == "both")) | 
            (Approach == "Combined" & (Feature == "all.features" | Feature == "tuned.features"))) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR_AUC by dataset for approaches including full sequence similarity, without NLV", subtitle = "Model 1: Positive vs Decoy pairs", y = "PR_AUC") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.15,0.6),oob = rescale_none)

ggsave(paste0(out.path, "pr_auc_by_dataset_full_sequence_model1_without_NLV.png"), width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & (Feature == "full" | Feature == "both")) | 
            (Approach == "Combined" & (Feature == "all.features" | Feature == "tuned.features"))) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR_AUC by dataset for approaches including full sequence similarity", subtitle = "Model 1: Positive vs Decoy pairs", y = "PR_AUC") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.15,0.6),oob = rescale_none)

ggsave(paste0(out.path, "pr_auc_by_dataset_full_sequence_model2.png"), width = 20, height = 15, dpi = 300)

# Without NLV
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        filter(NLV == "No") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & (Feature == "full" | Feature == "both")) | 
            (Approach == "Combined" & (Feature == "all.features" | Feature == "tuned.features"))) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR_AUC by dataset for approaches including full sequence similarity, without NLV", subtitle = "Model 1: Positive vs Decoy pairs", y = "PR_AUC") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.15,0.6),oob = rescale_none)

ggsave(paste0(out.path, "pr_auc_by_dataset_full_sequence_model2_without_NLV.png"), width = 20, height = 15, dpi = 300)

# Without full sequence
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Feature == "cdr3") | 
            (Approach == "Combined" & Feature == "limited.features")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR_AUC by dataset for approaches excluding full sequence similarity", subtitle = "Model 1: Positive vs Decoy pairs", y = "PR_AUC") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.15,0.6),oob = rescale_none)

ggsave(paste0(out.path, "pr_auc_by_dataset_cdr3_sequence_model1.png"), width = 20, height = 15, dpi = 300)

# Without NLV
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        filter(NLV == "No") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Feature == "cdr3") | 
            (Approach == "Combined" & Feature == "limited.features")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR_AUC by dataset for approaches excluding full sequence similarity, without NLV", subtitle = "Model 1: Positive vs Decoy pairs", y = "PR_AUC") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.15,0.6),oob = rescale_none)

ggsave(paste0(out.path, "pr_auc_by_dataset_cdr3_sequence_model1_without_NLV.png"), width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Feature == "cdr3") | 
            (Approach == "Combined" & Feature == "limited.features")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR_AUC by dataset for approaches excluding full sequence similarity", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "PR_AUC") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.15,0.6),oob = rescale_none)

ggsave(paste0(out.path, "pr_auc_by_dataset_cdr3_sequence_model2.png"), width = 20, height = 15, dpi = 300)

# Without NLV
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        filter(NLV == "No") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Feature == "cdr3") | 
            (Approach == "Combined" & Feature == "limited.features")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR_AUC by dataset for approaches excluding full sequence similarity, without NLV", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "PR_AUC") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.15,0.6),oob = rescale_none)

ggsave(paste0(out.path, "pr_auc_by_dataset_cdr3_sequence_model2_without_NLV.png"), width = 20, height = 15, dpi = 300)

# Feature coefficients
# Convert to long format for plotting, and plot features on X, Coefficients on Y, colored by Approach + Feature
# Start HERE
# Combined feature sets
# These still need more work, but not necessary rn
library(ggpattern)
ggplot(
    all.coefficients.with.nlv %>%
        filter(Approach == "Combined") %>%
        filter(Model == "Model 1") %>%
        pivot_longer(cols = -c(Dataset, Approach, Model, NLV, Feature), names_to = "Variable", values_to = "Coefficient") %>%
        #rename_with(~ gsub("^Value_", "", .), starts_with("Value_")) %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        #mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Variable, y = Coefficient, fill = `LogReg Approach`, alpha = Data)) +
    geom_bar(stat = "identity", position = "dodge") +
    #facet_wrap( ~ Data) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Feature coefficients by dataset and approach", subtitle = "Model 1: Positive vs Decoy pairs", y = "Coefficient") +
    scale_fill_manual(values = model1.approach.colors)

ggsave(paste0(out.path, "combined_feature_coefficients_model1.png"), width = 20, height = 15, dpi = 300)



# Per-epitope performance metrics
# Accuracy
# Model 1
# With NLV
ggplot(
    all.epitope.metrics.with.nlv %>%
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        # filter(Approach == "Structure" | 
        #     (Approach == "Sequence" & Feature == "cdr3") | 
        #     (Approach == "Combined" & Feature == "limited.features")) %>%
        # mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    #facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    facet_wrap( ~ Epitope, nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset and epitope", subtitle = "Model 1: Positive vs Decoy pairs", y = "Accuracy") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "accuracy_by_epitope_model1_with_NLV.png"), width = 20, height = 10, dpi = 300)

# Without NLV
ggplot(
    all.epitope.metrics.without.nlv %>%
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        # filter(Approach == "Structure" | 
        #     (Approach == "Sequence" & Feature == "cdr3") | 
        #     (Approach == "Combined" & Feature == "limited.features")) %>%
        # mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    facet_wrap( ~ Epitope, nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset and epitope", subtitle = "Model 1: Positive vs Decoy pairs", y = "Accuracy") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave(paste0(out.path, "accuracy_by_epitope_model1_without_NLV.png"), width = 20, height = 10, dpi = 300)

# Model 2
# With NLV
ggplot(
    all.epitope.metrics.with.nlv %>%
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        # filter(Approach == "Structure" | 
        #     (Approach == "Sequence" & Feature == "cdr3") | 
        #     (Approach == "Combined" & Feature == "limited.features")) %>%
        # mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    facet_wrap( ~ Epitope, nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset and epitope", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "Accuracy") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.8),oob = rescale_none)

ggsave(paste0(out.path, "accuracy_by_epitope_model2_with_NLV.png"), width = 20, height = 10, dpi = 300)

# Without NLV
ggplot( 
    all.epitope.metrics.without.nlv %>%
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        # filter(Approach == "Structure" | 
        #     (Approach == "Sequence" & Feature == "cdr3") | 
        #     (Approach == "Combined" & Feature == "limited.features")) %>%
        # mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    facet_wrap( ~ Epitope, nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset and epitope", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "Accuracy") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.8),oob = rescale_none)

ggsave(paste0(out.path, "accuracy_by_epitope_model2_without_NLV.png"), width = 20, height = 10, dpi = 300)

# ROC AUC
# Model 1
# With NLV
ggplot(
    all.epitope.metrics.with.nlv %>%
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        # filter(Approach == "Structure" | 
        #     (Approach == "Sequence" & Feature == "cdr3") | 
        #     (Approach == "Combined" & Feature == "limited.features")) %>%
        # mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    facet_wrap( ~ Epitope, nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC AUC by dataset and epitope", subtitle = "Model 1: Positive vs Decoy pairs", y = "ROC AUC") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.8),oob = rescale_none)

ggsave(paste0(out.path, "roc_auc_by_epitope_model1_with_NLV.png"), width = 20, height = 10, dpi = 300)

# With NLV, only considering CDR3 sequence data
ggplot(
    all.epitope.metrics.with.nlv %>%
        filter(Model == "Model 1") %>%
        filter(Dataset == "CDR3") %>%
        # filter(Approach == "Structure" | 
        #     (Approach == "Sequence" & Feature == "cdr3") | 
        #     (Approach == "Combined" & Feature == "limited.features")) %>%
        # mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    facet_wrap( ~ Epitope, nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC AUC by dataset and epitope", subtitle = "Model 1: Positive vs Decoy pairs", y = "ROC AUC") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.8),oob = rescale_none)

ggsave(paste0(out.path, "roc_auc_by_epitope_model1_with_NLV_CDR3_sim.png"), width = 20, height = 10, dpi = 300)

# Without NLV
ggplot(
    all.epitope.metrics.without.nlv %>%
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        # filter(Approach == "Structure" | 
        #     (Approach == "Sequence" & Feature == "cdr3") | 
        #     (Approach == "Combined" & Feature == "limited.features")) %>%
        # mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    facet_wrap( ~ Epitope, nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC AUC by dataset and epitope", subtitle = "Model 1: Positive vs Decoy pairs", y = "ROC AUC") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.8),oob = rescale_none)

ggsave(paste0(out.path, "roc_auc_by_epitope_model1_without_NLV.png"), width = 20, height = 10, dpi = 300)

# Model 2
# With NLV
ggplot(
    all.epitope.metrics.with.nlv %>%
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        # filter(Approach == "Structure" | 
        #     (Approach == "Sequence" & Feature == "cdr3") | 
        #     (Approach == "Combined" & Feature == "limited.features")) %>%
        # mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    facet_wrap( ~ Epitope, nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC AUC by dataset and epitope", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "ROC AUC") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.8),oob = rescale_none)

ggsave(paste0(out.path, "roc_auc_by_epitope_model2_with_NLV.png"), width = 20, height = 10, dpi = 300)

# Without NLV
ggplot(
    all.epitope.metrics.without.nlv %>%
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        # filter(Approach == "Structure" | 
        #     (Approach == "Sequence" & Feature == "cdr3") | 
        #     (Approach == "Combined" & Feature == "limited.features")) %>%
        # mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    facet_wrap( ~ Epitope, nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC AUC by dataset and epitope", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "ROC AUC") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.8),oob = rescale_none)

ggsave(paste0(out.path, "roc_auc_by_epitope_model2_without_NLV.png"), width = 20, height = 10, dpi = 300)

# PR AUC
# Model 1 
# With NLV
ggplot(
    all.epitope.metrics.with.nlv %>%
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        # filter(Approach == "Structure" | 
        #     (Approach == "Sequence" & Feature == "cdr3") | 
        #     (Approach == "Combined" & Feature == "limited.features")) %>%
        # mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    facet_wrap( ~ Epitope, nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR AUC by dataset and epitope", subtitle = "Model 1: Positive vs Decoy pairs", y = "PR AUC") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.2,0.7),oob = rescale_none)

ggsave(paste0(out.path, "pr_auc_by_epitope_model1_with_NLV.png"), width = 20, height = 10, dpi = 300)

# Without NLV
ggplot(
    all.epitope.metrics.without.nlv %>%
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        # filter(Approach == "Structure" | 
        #     (Approach == "Sequence" & Feature == "cdr3") | 
        #     (Approach == "Combined" & Feature == "limited.features")) %>%
        # mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    facet_wrap( ~ Epitope, nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR AUC by dataset and epitope", subtitle = "Model 1: Positive vs Decoy pairs", y = "PR AUC") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.2,0.7),oob = rescale_none)

ggsave(paste0(out.path, "pr_auc_by_epitope_model1_without_NLV.png"), width = 20, height = 10, dpi = 300)

# Model 2
# With NLV
ggplot(
    all.epitope.metrics.with.nlv %>%
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        # filter(Approach == "Structure" | 
        #     (Approach == "Sequence" & Feature == "cdr3") | 
        #     (Approach == "Combined" & Feature == "limited.features")) %>%
        # mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    facet_wrap( ~ Epitope, nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR AUC by dataset and epitope", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "PR AUC") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.05,0.55),oob = rescale_none)

ggsave(paste0(out.path, "pr_auc_by_epitope_model2_with_NLV.png"), width = 20, height = 10, dpi = 300)

# Without NLV
ggplot(
    all.epitope.metrics.without.nlv %>%
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        # filter(Approach == "Structure" | 
        #     (Approach == "Sequence" & Feature == "cdr3") | 
        #     (Approach == "Combined" & Feature == "limited.features")) %>%
        # mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Feature, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Feature == "cdr3" ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Feature == "both" ~ paste(Approach, "approach (CDR3 + full seq.)"),
            Approach == "Sequence" & Feature == "full" ~ paste(Approach, "approach (full seq. only)"),
            .default = paste(Approach, "approach,", Feature))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    facet_wrap( ~ Epitope, nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR AUC by dataset and epitope", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "PR AUC") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.05,0.55),oob = rescale_none)

ggsave(paste0(out.path, "pr_auc_by_epitope_model2_without_NLV.png"), width = 20, height = 10, dpi = 300)

# Selection of best models: one each per approach and model type
# Model 1
# ROC AUC curves
# ROC AUC curves per epitope
# PR AUC curves
# PR AUC curves per epitope

# Model 2
# ROC AUC curves
# ROC AUC curves per epitope
# PR AUC curves
# PR AUC curves per epitope

