# ============================================================
# File: process-results.R
# Description: Script for training and analyzing the PIEMA logistic regression classifier
# Author: Maia Bennett-Boehm
# University of Nebraska at Omaha
# ============================================================

# Load the necessary libraries
library(GGally)
library(skimr)
library(DataExplorer)
library(caret) # for nicer confusion matrix
library(discrim) # for discriminant analysis
library(tidyverse)
library(tidymodels)
library(themis)
library(broom)
library(viridis)
library(ggrepel)

# out.path <- "./analysis/classifier/run1/9-feature-set/with-NLV/"
# out.path <- "./analysis/classifier/run1/9-feature-set/without-NLV/"
# out.path <- "./analysis/classifier/run1/7-feature-set/with-NLV/"
# out.path <- "./analysis/classifier/run1/7-feature-set/without-NLV/"
# out.path <- "./analysis/classifier/run1/5-feature-set/with-NLV/"
# out.path <- "./analysis/classifier/run1/5-feature-set/without-NLV/"
# out.path <- "./analysis/classifier/run2/9-feature-set/with-NLV/"
out.path <- "./analysis/classifier/run2/9-feature-set/without-NLV/"
# out.path <- "./analysis/classifier/run2/7-feature-set/with-NLV/"
# out.path <- "./analysis/classifier/run2/7-feature-set/without-NLV/"
# out.path <- "./analysis/classifier/run2/5-feature-set/with-NLV/"
# out.path <- "./analysis/classifier/run2/5-feature-set/without-NLV/"
# input.path <- "./analysis/apbs/run1/"
input.path <- "./analysis/apbs/run2/"

# Save and load the workspace
# save.image(file = "workspace.RData")
# load(file = "workspace.RData")


# This script uses tidymodels to train a logistic regression classifier on the PIEMA data
# For this purpose, two models will be trained: 
# 1. A model that uses only "True binding pairs" (positive data) and "Decoy binding pairs" (negative data)
# 2. A model that uses all data, including "True binding pairs" (positive data) and "Decoy binding pairs" and "Unlikely binding pairs" (negative data)

# Load the data
piema.data <- read.csv(paste0(input.path, "final_receptor_data.csv")) %>%
    dplyr::rename(top.kas = top.kdist, top.spearman = top.corr, top.euc.dist = top.distance, 
        top.sg.euc.dist = top.patDist, top.cosine.sim = top.shapSim, top.sg.median.kas = top.med_scr, 
        binding.pair.type = type, avg.euc.dist = avg.distance, avg.kas = avg.kdist, avg.spearman = avg.corr,
        kernel.count = count, pair.id = id, top.sg.id = top.sg_group)

# If running without NLV
piema.data <- piema.data %>%
    filter(ref.epitope != "NLVPMVATV") %>%
    filter(samp.epitope != "NLVPMVATV")

# Prepare the data, adding a column to denote whether the receptor pairs share binding specificity (Yes) or not (No)
# First,  we only want to use the "True binding pairs" and "Decoy binding pairs"
model1.data <- piema.data %>% 
    filter(binding.pair.type == "True receptor pairs" | binding.pair.type == "Decoy receptor pairs") %>%
    mutate(shared.specificity = ifelse(binding.pair.type == "True receptor pairs", "Yes", "No")) %>%
    mutate(shared.specificity = as.factor(shared.specificity))

# Next, we want to use all data
model2.data <- piema.data %>%
    mutate(shared.specificity = ifelse(binding.pair.type == "True receptor pairs", "Yes", "No"))%>%
    mutate(shared.specificity = as.factor(shared.specificity))

model1.long <- model1.data %>% 
    select(shared.specificity, top.kas, avg.euc.dist, top.sg.euc.dist, avg.kas, avg.spearman, top.cosine.sim, top.sg.median.kas, kernel.count, CDR3.similarity, full.similarity) %>%
    pivot_longer(cols = c(top.kas, avg.euc.dist, top.sg.euc.dist, avg.kas, avg.spearman, top.cosine.sim, top.sg.median.kas, kernel.count, CDR3.similarity, full.similarity), names_to = "feature", values_to = "value")

ggplot(model1.long, aes(x = shared.specificity, y = value, fill = feature)) +
    geom_boxplot() +
    facet_wrap(~feature, scales = "free", ncol = 5) +
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = 0.9) +
    theme_minimal() +
    theme(legend.position = "none", axis.title.x = element_blank(), 
        strip.background = element_rect(fill = "lightgrey"), 
        strip.text = element_text(color = "black", face = "bold"))

ggsave(paste0(out.path, "model1_feature_boxplot.png"), width = 20, height = 10)

model2.long <- model2.data %>%
    select(shared.specificity, top.kas, avg.euc.dist, top.sg.euc.dist, avg.kas, avg.spearman, top.cosine.sim, top.sg.median.kas, kernel.count, CDR3.similarity, full.similarity) %>%
    pivot_longer(cols = c(top.kas, avg.euc.dist, top.sg.euc.dist, avg.kas, avg.spearman, top.cosine.sim, top.sg.median.kas, kernel.count, CDR3.similarity, full.similarity), names_to = "feature", values_to = "value")

ggplot(model2.long, aes(x = shared.specificity, y = value, fill = feature)) +
    geom_boxplot() +
    facet_wrap(~feature, scales = "free", ncol = 5) +
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = 0.9) +
    theme_minimal() +
    theme(legend.position = "none", axis.title.x = element_blank(), 
        strip.background = element_rect(fill = "grey60"), 
        strip.text = element_text(color = "white", face = "bold"))

ggsave(paste0(out.path, "model2_feature_boxplot.png"), width = 20, height = 10)


# Set seed for reproducibility
set.seed(77482951)

# Split the data into training and testing sets
model1.split <- initial_split(model1.data, prop = 0.75, strata = shared.specificity)
model1.train <- training(model1.split)
model1.test <- testing(model1.split)

model2.split <- initial_split(model2.data, prop = 0.75, strata = shared.specificity)
model2.train <- training(model2.split)
model2.test <- testing(model2.split)

# Define the variables of interest and recipe for the logistic regression model
target.var <- "shared.specificity"
positive.class <- "Yes"
negative.class <- "No"

# Structure + sequence model formula
# Nine feature set
combined.model.form <- shared.specificity ~ avg.euc.dist + avg.kas + top.sg.median.kas + top.sg.euc.dist + avg.spearman + top.cosine.sim + kernel.count + CDR3.similarity + full.similarity
# Seven feature set
# combined.model.form <- shared.specificity ~ top.sg.euc.dist + top.sg.median.kas + avg.spearman + top.cosine.sim + kernel.count + CDR3.similarity + full.similarity
# Five feature set
# combined.model.form <- shared.specificity ~ top.sg.euc.dist + top.sg.median.kas + avg.spearman + kernel.count + CDR3.similarity
# Structure-only model formula
# 7 features (in nine feature set)
structure.model.form <- shared.specificity ~ avg.euc.dist + avg.kas + top.sg.median.kas + top.sg.euc.dist + avg.spearman + top.cosine.sim + kernel.count
# 5 features (in seven feature set)
# structure.model.form <- shared.specificity ~ top.sg.euc.dist + top.sg.median.kas + avg.spearman + top.cosine.sim + kernel.count 
# 4 features (in five feature set)
# structure.model.form <- shared.specificity ~  top.sg.euc.dist + top.sg.median.kas + avg.spearman + kernel.count
# Sequence-only model formula
# 2 features (in nine and seven feature sets)
sequence.model.form <- shared.specificity ~ CDR3.similarity + full.similarity
# 1 feature (in five feature set)
# sequence.model.form <- shared.specificity ~ CDR3.similarity

# Define F beta metric
# f_meas <- function(data, truth, estimate, na_rm = TRUE, ...) {
#     yardstick::f_meas(
#         data = data,
#         truth = !!rlang::enquo(truth),
#         estimate = !!rlang::enquo(estimate),
#         beta = 0.5,
#         na_rm = na_rm,
#         ...
#     )
# }
# f_meas <- new_class_metric(f_meas, direction = "maximize")

classification.metrics <- metric_set(
    yardstick::accuracy, 
    yardstick::kap, 
    yardstick::roc_auc, 
    yardstick::pr_auc,
    yardstick::mn_log_loss, 
    yardstick::sens, 
    yardstick::spec, 
    yardstick::f_meas,
    yardstick::precision,
    yardstick::recall,
    yardstick::bal_accuracy,
    # f_meas
    ) 

# Define the recipe for the logistic regression models
combined.model1.recipe <- recipe(combined.model.form, data = model1.train) %>%
    step_dummy(all_nominal(), -all_outcomes()) %>%
    step_normalize(all_predictors())
combined.model2.recipe <- recipe(combined.model.form, data = model2.train) %>%
    step_dummy(all_nominal(), -all_outcomes()) %>%
    step_normalize(all_predictors())

structure.model1.recipe <- recipe(structure.model.form, data = model1.train) %>%
    step_dummy(all_nominal(), -all_outcomes()) %>%
    step_normalize(all_predictors())
structure.model2.recipe <- recipe(structure.model.form, data = model2.train) %>%
    step_dummy(all_nominal(), -all_outcomes()) %>%
    step_normalize(all_predictors())

sequence.model1.recipe <- recipe(sequence.model.form, data = model1.train) %>%
    step_dummy(all_nominal(), -all_outcomes()) %>%
    step_normalize(all_predictors())
sequence.model2.recipe <- recipe(sequence.model.form, data = model2.train) %>%
    step_dummy(all_nominal(), -all_outcomes()) %>%
    step_normalize(all_predictors())

# Recipes with upsampling
combined.model1.upsampling.recipe <- combined.model1.recipe %>%
    step_smotenc(target.var)
combined.model2.upsampling.recipe <- combined.model2.recipe %>%
    step_smotenc(target.var)

structure.model1.upsampling.recipe <- structure.model1.recipe %>%
    step_smotenc(target.var)
structure.model2.upsampling.recipe <- structure.model2.recipe %>%
    step_smotenc(target.var)

sequence.model1.upsampling.recipe <- sequence.model1.recipe %>%
    step_smotenc(target.var)
sequence.model2.upsampling.recipe <- sequence.model2.recipe %>%
    step_smotenc(target.var)

# Recipes with downsampling
combined.model1.downsampling.recipe <- combined.model1.recipe %>%
    step_downsample(target.var)
combined.model2.downsampling.recipe <- combined.model2.recipe %>%
    step_downsample(target.var)

structure.model1.downsampling.recipe <- structure.model1.recipe %>%
    step_downsample(target.var)
structure.model2.downsampling.recipe <- structure.model2.recipe %>%
    step_downsample(target.var)

sequence.model1.downsampling.recipe <- sequence.model1.recipe %>%
    step_downsample(target.var)
sequence.model2.downsampling.recipe <- sequence.model2.recipe %>%
    step_downsample(target.var)

# Define recipe set
combined.model1.recipes <- list(basic = combined.model1.recipe, upsampled = combined.model1.upsampling.recipe, downsampled = combined.model1.downsampling.recipe)
combined.model2.recipes <- list(basic = combined.model2.recipe, upsampled = combined.model2.upsampling.recipe, downsampled = combined.model2.downsampling.recipe)

structure.model1.recipes <- list(basic = structure.model1.recipe, upsampled = structure.model1.upsampling.recipe, downsampled = structure.model1.downsampling.recipe)
structure.model2.recipes <- list(basic = structure.model2.recipe, upsampled = structure.model2.upsampling.recipe, downsampled = structure.model2.downsampling.recipe)

sequence.model1.recipes <- list(basic = sequence.model1.recipe, upsampled = sequence.model1.upsampling.recipe, downsampled = sequence.model1.downsampling.recipe)
sequence.model2.recipes <- list(basic = sequence.model2.recipe, upsampled = sequence.model2.upsampling.recipe, downsampled = sequence.model2.downsampling.recipe)

# Define the logistic regression model for both models
logreg.model <- 
  logistic_reg() %>%
  set_engine("glm") %>%
  set_mode("classification") 


# Define multi-recipe workflow sets for each model
# Combined, structure-only, and sequence-only models are trained separately so their results can be more easily compared
combined.model1.workflow <- workflow_set(combined.model1.recipes,
    models = list(logreg = logreg.model),
    cross = FALSE)
combined.model2.workflow <- workflow_set(combined.model2.recipes,
    models = list(logreg = logreg.model),
    cross = FALSE)

structure.model1.workflow <- workflow_set(structure.model1.recipes,
    models = list(logreg = logreg.model),
    cross = FALSE)
structure.model2.workflow <- workflow_set(structure.model2.recipes,
    models = list(logreg = logreg.model),
    cross = FALSE)

sequence.model1.workflow <- workflow_set(sequence.model1.recipes,
    models = list(logreg = logreg.model),
    cross = FALSE)
sequence.model2.workflow <- workflow_set(sequence.model2.recipes,
    models = list(logreg = logreg.model),
    cross = FALSE)

# Specify resamples and CV folds
set.seed(77482951)
control <- control_resamples(save_workflow = TRUE,
                                 save_pred = TRUE,
                                 event_level = "second")

model1.cv.folds <- vfold_cv(model1.train, v = 10)
model2.cv.folds <- vfold_cv(model2.train, v = 10)


# Fit the model sets
combined.model1.fit <- combined.model1.workflow %>%
    workflow_map("fit_resamples", 
        seed = 77482951,
        verbose = TRUE,
        resamples = model1.cv.folds,
        metrics = classification.metrics,
        control = control)
combined.model2.fit <- combined.model2.workflow %>%
    workflow_map("fit_resamples",
        seed = 77482951,
        verbose = TRUE,
        resamples = model2.cv.folds,
        metrics = classification.metrics,
        control = control)

structure.model1.fit <- structure.model1.workflow %>%
    workflow_map("fit_resamples",
        seed = 77482951,
        verbose = TRUE,
        resamples = model1.cv.folds,
        metrics = classification.metrics,
        control = control)
structure.model2.fit <- structure.model2.workflow %>%
    workflow_map("fit_resamples",
        seed = 77482951,
        verbose = TRUE,
        resamples = model2.cv.folds,
        metrics = classification.metrics,
        control = control)

sequence.model1.fit <- sequence.model1.workflow %>%
    workflow_map("fit_resamples",
        seed = 77482951,
        verbose = TRUE,
        resamples = model1.cv.folds,
        metrics = classification.metrics,
        control = control)
sequence.model2.fit <- sequence.model2.workflow %>%
    workflow_map("fit_resamples",
        seed = 77482951,
        verbose = TRUE,
        resamples = model2.cv.folds,
        metrics = classification.metrics,
        control = control)


# Store workflow set predictions
combined.model1.train.predictions <- combined.model1.fit %>%
    collect_predictions()
combined.model2.train.predictions <- combined.model2.fit %>%
    collect_predictions()

structure.model1.train.predictions <- structure.model1.fit %>%
    collect_predictions()
structure.model2.train.predictions <- structure.model2.fit %>%
    collect_predictions()

sequence.model1.train.predictions <- sequence.model1.fit %>%
    collect_predictions()
sequence.model2.train.predictions <- sequence.model2.fit %>%
    collect_predictions()

# Calculate workflow set metrics
combined.model1.train.metrics <- combined.model1.fit %>%
    collect_metrics() 
combined.model2.train.metrics <- combined.model2.fit %>%
    collect_metrics()

structure.model1.train.metrics <- structure.model1.fit %>%
    collect_metrics()
structure.model2.train.metrics <- structure.model2.fit %>%
    collect_metrics()

sequence.model1.train.metrics <- sequence.model1.fit %>%
    collect_metrics()
sequence.model2.train.metrics <- sequence.model2.fit %>%
    collect_metrics()

# Fetch the best model by the specified metric 
# metric.for.selection <- "f_meas"
metric.for.selection <- "recall"

combined.model1.best.fit <- combined.model1.fit %>%
    fit_best(metric = metric.for.selection)
combined.model2.best.fit <- combined.model2.fit %>%
    fit_best(metric = metric.for.selection)

structure.model1.best.fit <- structure.model1.fit %>%
    fit_best(metric = metric.for.selection)
structure.model2.best.fit <- structure.model2.fit %>%
    fit_best(metric = metric.for.selection)

sequence.model1.best.fit <- sequence.model1.fit %>%
    fit_best(metric = metric.for.selection)
sequence.model2.best.fit <- sequence.model2.fit %>%
    fit_best(metric = metric.for.selection)

# Add column to training metrics to identify the selected model 
combined.model1.selected.method <- combined.model1.train.metrics %>%
    filter(.metric == metric.for.selection) %>%
    filter(mean == max(mean)) %>%
    pull(wflow_id)
combined.model1.train.metrics <- combined.model1.train.metrics %>%
    mutate(selected = ifelse(wflow_id == combined.model1.selected.method, "Yes", "No"))

combined.model2.selected.method <- combined.model2.train.metrics %>%
    filter(.metric == metric.for.selection) %>%
    filter(mean == max(mean)) %>%
    pull(wflow_id)
combined.model2.train.metrics <- combined.model2.train.metrics %>%
    mutate(selected = ifelse(wflow_id == combined.model2.selected.method, "Yes", "No"))

structure.model1.selected.method <- structure.model1.train.metrics %>%
    filter(.metric == metric.for.selection) %>%
    filter(mean == max(mean)) %>%
    pull(wflow_id)
structure.model1.train.metrics <- structure.model1.train.metrics %>%
    mutate(selected = ifelse(wflow_id == structure.model1.selected.method, "Yes", "No"))

structure.model2.selected.method <- structure.model2.train.metrics %>%
    filter(.metric == metric.for.selection) %>%
    filter(mean == max(mean)) %>%
    pull(wflow_id)
structure.model2.train.metrics <- structure.model2.train.metrics %>%
    mutate(selected = ifelse(wflow_id == structure.model2.selected.method, "Yes", "No"))

sequence.model1.selected.method <- sequence.model1.train.metrics %>%
    filter(.metric == metric.for.selection) %>%
    filter(mean == max(mean)) %>%
    pull(wflow_id)
sequence.model1.train.metrics <- sequence.model1.train.metrics %>%
    mutate(selected = ifelse(wflow_id == sequence.model1.selected.method, "Yes", "No"))

sequence.model2.selected.method <- sequence.model2.train.metrics %>%
    filter(.metric == metric.for.selection) %>%
    filter(mean == max(mean)) %>%
    pull(wflow_id)
sequence.model2.train.metrics <- sequence.model2.train.metrics %>%
    mutate(selected = ifelse(wflow_id == sequence.model2.selected.method, "Yes", "No"))

# Fit the model sets to the test data and store test predictions
combined.model1.final.fit <- combined.model1.best.fit %>%
    augment(new_data = model1.test)
combined.model2.final.fit <- combined.model2.best.fit %>%
    augment(new_data = model2.test)

structure.model1.final.fit <- structure.model1.best.fit %>%
    augment(new_data = model1.test)
structure.model2.final.fit <- structure.model2.best.fit %>%
    augment(new_data = model2.test)

sequence.model1.final.fit <- sequence.model1.best.fit %>%
    augment(new_data = model1.test)
sequence.model2.final.fit <- sequence.model2.best.fit %>%
    augment(new_data = model2.test)

# Calculate test metrics
combined.model1.test.metrics <- classification.metrics(data = combined.model1.final.fit,
    truth = shared.specificity,
    estimate = .pred_class,
    .pred_Yes,
    event_level = 'second')
combined.model2.test.metrics <- classification.metrics(data = combined.model2.final.fit,
    truth = shared.specificity,
    estimate = .pred_class,
    .pred_Yes,
    event_level = 'second')

structure.model1.test.metrics <- classification.metrics(data = structure.model1.final.fit,
    truth = shared.specificity,
    estimate = .pred_class,
    .pred_Yes,
    event_level = 'second')
structure.model2.test.metrics <- classification.metrics(data = structure.model2.final.fit,
    truth = shared.specificity,
    estimate = .pred_class,
    .pred_Yes,
    event_level = 'second')

sequence.model1.test.metrics <- classification.metrics(data = sequence.model1.final.fit,
    truth = shared.specificity,
    estimate = .pred_class,
    .pred_Yes,
    event_level = 'second')
sequence.model2.test.metrics <- classification.metrics(data = sequence.model2.final.fit,
    truth = shared.specificity,
    estimate = .pred_class,
    .pred_Yes,
    event_level = 'second')

# Combine metrics for easy comparison
# Combined and sequence-only are very similar, but structure-only has much higher sensitivity- interesting
model1.metrics <- combined.model1.test.metrics %>% 
    select(.metric, .estimate) %>%
    dplyr::rename(combined.estimate = .estimate) %>%
    cbind(structure.model1.test.metrics %>% select(.estimate) %>% 
        dplyr::rename(structure.estimate = .estimate)) %>%
    cbind(sequence.model1.test.metrics %>% select(.estimate) %>%
        dplyr::rename(sequence.estimate = .estimate))

model2.metrics <- combined.model2.test.metrics %>% 
    select(.metric, .estimate) %>%
    dplyr::rename(combined.estimate = .estimate) %>%
    cbind(structure.model2.test.metrics %>% select(.estimate) %>% 
        dplyr::rename(structure.estimate = .estimate)) %>%
    cbind(sequence.model2.test.metrics %>% select(.estimate) %>%
        dplyr::rename(sequence.estimate = .estimate))

# Save the results
write.csv(combined.model1.train.metrics, paste0(out.path, "combined_model1_train_metrics.csv"), row.names = FALSE)
write.csv(combined.model2.train.metrics, paste0(out.path, "combined_model2_train_metrics.csv"), row.names = FALSE)
write.csv(combined.model1.final.fit, paste0(out.path, "combined_model1_predictions.csv"), row.names = FALSE)
write.csv(combined.model2.final.fit, paste0(out.path, "combined_model2_predictions.csv"), row.names = FALSE)
write.csv(structure.model1.final.fit, paste0(out.path, "structure_model1_predictions.csv"), row.names = FALSE)
write.csv(structure.model2.final.fit, paste0(out.path, "structure_model2_predictions.csv"), row.names = FALSE)
write.csv(sequence.model1.final.fit, paste0(out.path, "sequence_model1_predictions.csv"), row.names = FALSE)
write.csv(sequence.model2.final.fit, paste0(out.path, "sequence_model2_predictions.csv"), row.names = FALSE)
write.csv(model1.metrics, paste0(out.path, "model1_metrics.csv"), row.names = FALSE)
write.csv(model2.metrics, paste0(out.path, "model2_metrics.csv"), row.names = FALSE)

# Plot whole model ROC + PR
# Model 1
combined.model1.rocauc <- combined.model1.final.fit %>% 
    roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    pull(.estimate)
combined.model1.roc <- combined.model1.final.fit %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>% 
    mutate(model = paste0("Best model combining structure and sequence features\n", 
        "ROC AUC: ", round(combined.model1.rocauc, 3)))
structure.model1.rocauc <- structure.model1.final.fit %>%
    roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    pull(.estimate)
structure.model1.roc <- structure.model1.final.fit %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>% 
    mutate(model = paste0("Best model using structure features\n", 
        "ROC AUC: ", round(structure.model1.rocauc, 3)))
sequence.model1.rocauc <- sequence.model1.final.fit %>%
    roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    pull(.estimate)
sequence.model1.roc <- sequence.model1.final.fit %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>% 
    mutate(model = paste0("Best model using sequence features\n", 
        "ROC AUC: ", round(sequence.model1.rocauc, 3)))
bind_rows(combined.model1.roc, structure.model1.roc, sequence.model1.roc) %>% 
    ggplot(aes(x = 1 - specificity, y = sensitivity, col = model, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.6) +
    geom_abline(lty = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "mako", begin = 0.2, end = 0.8) 

ggsave(paste0(out.path, "model1_roc.png"), width = 10, height = 10)

# Model 2
combined.model2.roc <- combined.model2.final.fit %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(model = paste0("Best model combining structure and sequence features\n"),
        "ROC AUC: ", combined.model1.final.fit %>% 
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate))
structure.model2.roc <- structure.model2.final.fit %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(model = paste0("Best model using structure features\n"),
        "ROC AUC: ", round(structure.model1.final.fit %>% 
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3))
sequence.model2.roc <- sequence.model2.final.fit %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(model = paste0("Best model using sequence features\n"),
        "ROC AUC: ", round(sequence.model1.final.fit %>% 
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3))
bind_rows(combined.model2.roc, structure.model2.roc, sequence.model2.roc) %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, col = model, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.8) +
    geom_abline(lty = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "mako", begin = 0.2, end = 0.8) 

ggsave(paste0(out.path, "model2_roc.png"), width = 10, height = 10)

# Model 1
combined.model1.pr <- combined.model1.final.fit %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(model = paste0("Best model combining structure and sequence features\n", 
        "PR AUC: ", round(combined.model1.final.fit %>% 
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model1.pr <- structure.model1.final.fit %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(model = paste0("Best model using structure features\n", 
        "PR AUC: ", round(structure.model1.final.fit %>% 
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model1.pr <- sequence.model1.final.fit %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(model = paste0("Best model using sequence features\n", 
        "PR AUC: ", round(sequence.model1.final.fit %>% 
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
bind_rows(combined.model1.pr, structure.model1.pr, sequence.model1.pr) %>%
    ggplot(aes(x = recall, y = precision, col = model, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.8) +
    coord_equal() + 
    scale_color_viridis_d(option = "mako", begin = 0.2, end = 0.8) 
ggsave(paste0(out.path, "model1_pr.png"), width = 10, height = 10)

# Model 2
combined.model2.pr <- combined.model2.final.fit %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(model = paste0("Best model combining structure and sequence features\n",
        "PR AUC: ", round(combined.model2.final.fit %>% 
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model2.pr <- structure.model2.final.fit %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(model = paste0("Best model using structure features\n",
        "PR AUC: ", round(structure.model2.final.fit %>% 
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model2.pr <- sequence.model2.final.fit %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(model = paste0("Best model using sequence features\n",
        "PR AUC: ", round(sequence.model2.final.fit %>% 
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))

bind_rows(combined.model2.pr, structure.model2.pr, sequence.model2.pr) %>%
    ggplot(aes(x = recall, y = precision, col = model, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.8) +
    coord_equal() + 
    scale_color_viridis_d(option = "mako", begin = 0.2, end = 0.8)

ggsave(paste0(out.path, "model2_pr.png"), width = 10, height = 10)

# Extract and store final model formulas, other information
combined.model1.best.fit$fit[2]
combined.model1.coefficients <- combined.model1.best.fit %>%
    extract_fit_engine() %>%
    broom::tidy()
combined.model2.best.fit$fit[2]
combined.model2.coefficients <- combined.model2.best.fit %>%
    extract_fit_engine() %>%
    broom::tidy()

structure.model1.best.fit$fit[2]
structure.model1.coefficients <- structure.model1.best.fit %>%
    extract_fit_engine() %>%
    broom::tidy()   
structure.model2.best.fit$fit[2]
structure.model2.coefficients <- structure.model2.best.fit %>%
    extract_fit_engine() %>%
    broom::tidy()

sequence.model1.best.fit$fit[2]
sequence.model1.coefficients <- sequence.model1.best.fit %>%
    extract_fit_engine() %>%
    broom::tidy()
sequence.model2.best.fit$fit[2]
sequence.model2.coefficients <- sequence.model2.best.fit %>%
    extract_fit_engine() %>%
    broom::tidy()

# Combine coefficients for easy comparison
model1.coefficients <- combined.model1.coefficients %>% 
    select(term, estimate, p.value) %>%
    dplyr::rename(combined.estimate = estimate, combined.p.value = p.value) %>%
    full_join(structure.model1.coefficients %>% 
            select(term, estimate, p.value) %>% 
            dplyr::rename(structure.estimate = estimate, structure.p.value = p.value), 
            by = "term") %>%
    full_join(sequence.model1.coefficients %>% 
            select(term, estimate, p.value) %>% 
            dplyr::rename(sequence.estimate = estimate, sequence.p.value = p.value), 
            by = "term") 

model2.coefficients <- combined.model2.coefficients %>%
    select(term, estimate, p.value) %>%
    dplyr::rename(combined.estimate = estimate, combined.p.value = p.value) %>%
    full_join(structure.model2.coefficients %>% 
            select(term, estimate, p.value) %>% 
            dplyr::rename(structure.estimate = estimate, structure.p.value = p.value), 
            by = "term") %>%
    full_join(sequence.model2.coefficients %>%
            select(term, estimate, p.value) %>% 
            dplyr::rename(sequence.estimate = estimate, sequence.p.value = p.value), 
            by = "term")

# Save coefficients
write.csv(model1.coefficients, paste0(out.path, "model1_coefficients.csv"), row.names = FALSE)
write.csv(model2.coefficients, paste0(out.path, "model2_coefficients.csv"), row.names = FALSE)

# Check performance of models on specific epitopes
# Fetch predictions from the best models (i.e., final fit objects), group by epitope, and calculate metrics
# For model 1, this is simple, as ref.epitope and samp.epitope are (generally) the same 
# However, we do need to duplicate rows where the ref and samp epitopes are different AND shared.specificity is 'Yes' (i.e., not missing the cross-reactive positive cases)
# Metrics sets- without crossreactive accounting
combined.model1.epitope.metrics <- combined.model1.final.fit %>%
    group_by(ref.epitope) %>%
    classification.metrics(truth = shared.specificity,
        estimate = .pred_class,
        .pred_Yes,
        event_level = 'second')
structure.model1.epitope.metrics <- structure.model1.final.fit %>%
    group_by(ref.epitope) %>%
    classification.metrics(truth = shared.specificity,
        estimate = .pred_class,
        .pred_Yes,
        event_level = 'second')
sequence.model1.epitope.metrics <- sequence.model1.final.fit %>%
    group_by(ref.epitope) %>%
    classification.metrics(truth = shared.specificity,
        estimate = .pred_class,
        .pred_Yes,
        event_level = 'second')

model1.epitope.metrics <- combined.model1.epitope.metrics %>% 
    select(ref.epitope, .metric, .estimate) %>% 
    dplyr::rename(combined.estimate = .estimate) %>% 
    cbind(structure.model1.epitope.metrics %>% select(.estimate) %>% 
        dplyr::rename(structure.estimate = .estimate)) %>% 
    cbind(sequence.model1.epitope.metrics %>% select(.estimate) %>% 
        dplyr::rename(sequence.estimate = .estimate)) %>% 
    arrange(ref.epitope)

# Metrics sets- with crossreactive accounting
combined.model1.final.fit.with.crossreactives <- combined.model1.final.fit %>%
    mutate(epitope = ref.epitope) %>%
    bind_rows(combined.model1.final.fit %>% 
        filter(shared.specificity == "Yes" & ref.epitope != samp.epitope) %>%
        mutate(epitope = samp.epitope))
combined.model1.epitope.metrics.with.crossreactives <- combined.model1.final.fit.with.crossreactives %>%
    group_by(epitope) %>%
    classification.metrics(truth = shared.specificity,
        estimate = .pred_class,
        .pred_Yes,
        event_level = 'second')
structure.model1.final.fit.with.crossreactives <- structure.model1.final.fit %>%
    mutate(epitope = ref.epitope) %>%
    bind_rows(combined.model1.final.fit %>% 
        filter(shared.specificity == "Yes" & ref.epitope != samp.epitope) %>%
        mutate(epitope = samp.epitope))
structure.model1.epitope.metrics.with.crossreactives <- structure.model1.final.fit.with.crossreactives %>%
    group_by(epitope) %>%
    classification.metrics(truth = shared.specificity,
        estimate = .pred_class,
        .pred_Yes,
        event_level = 'second')
sequence.model1.final.fit.with.crossreactives <- sequence.model1.final.fit %>%
    mutate(epitope = ref.epitope) %>%
    bind_rows(combined.model1.final.fit %>% 
        filter(shared.specificity == "Yes" & ref.epitope != samp.epitope) %>%
        mutate(epitope = samp.epitope))
sequence.model1.epitope.metrics.with.crossreactives <- sequence.model1.final.fit.with.crossreactives %>%
    group_by(epitope) %>%
    classification.metrics(truth = shared.specificity,
        estimate = .pred_class,
        .pred_Yes,
        event_level = 'second')

model1.epitope.metrics.with.crossreactives <- combined.model1.epitope.metrics.with.crossreactives %>% 
    select(epitope, .metric, .estimate) %>%
    dplyr::rename(combined.estimate = .estimate) %>%
    cbind(structure.model1.epitope.metrics.with.crossreactives %>% select(.estimate) %>% 
        dplyr::rename(structure.estimate = .estimate)) %>%
    cbind(sequence.model1.epitope.metrics.with.crossreactives %>% select(.estimate) %>%
        dplyr::rename(sequence.estimate = .estimate)) %>%
    arrange(epitope)


# For model 2, we need to ensure that all epitopes are represented where they are ref.epitopes AND samp.epitopes
# The logreg model doesn't necessarily care about which is a 'ref' and which is a 'samp', it cares about the pairwise relationships between receptors
# So, simply grouping by ref.epitope misses negative data associated with each receptor where it is a samp.epitope
# As such, we need to duplicate data where the ref and samp epitopes are different AND shared.specificity is 'No' (i.e., not accidentally counting cross-reactives as negative) or where the epitopes are different and shared.specificity is 'Yes' (i.e., not missing cross-reactives), and then group by epitope
combined.model2.final.fit.expanded <- combined.model2.final.fit %>%
    mutate(epitope = ref.epitope) %>%
    bind_rows(combined.model1.final.fit %>% 
        filter(ref.epitope != samp.epitope) %>%
        mutate(epitope = samp.epitope))
combined.model2.epitope.metrics <- combined.model2.final.fit.expanded %>%
    group_by(epitope) %>%
    classification.metrics(truth = shared.specificity,
        estimate = .pred_class,
        .pred_Yes,
        event_level = 'second')
structure.model2.final.fit.expanded <- structure.model2.final.fit %>%
    mutate(epitope = ref.epitope) %>%
    bind_rows(combined.model1.final.fit %>% 
        filter(ref.epitope != samp.epitope) %>%
        mutate(epitope = samp.epitope))
structure.model2.epitope.metrics <- structure.model2.final.fit.expanded %>%
    group_by(epitope) %>%
    classification.metrics(truth = shared.specificity,
        estimate = .pred_class,
        .pred_Yes,
        event_level = 'second')
sequence.model2.final.fit.expanded <- sequence.model2.final.fit %>%
    mutate(epitope = ref.epitope) %>%
    bind_rows(combined.model1.final.fit %>% 
        filter(ref.epitope != samp.epitope) %>%
        mutate(epitope = samp.epitope))
sequence.model2.epitope.metrics <- sequence.model2.final.fit.expanded %>%
    group_by(epitope) %>%
    classification.metrics(truth = shared.specificity,
        estimate = .pred_class,
        .pred_Yes,
        event_level = 'second')


model2.epitope.metrics <- combined.model2.epitope.metrics %>%
    select(epitope, .metric, .estimate) %>%
    dplyr::rename(combined.estimate = .estimate) %>%
    cbind(structure.model2.epitope.metrics %>% select(.estimate) %>% 
        dplyr::rename(structure.estimate = .estimate)) %>%
    cbind(sequence.model2.epitope.metrics %>% select(.estimate) %>%
        dplyr::rename(sequence.estimate = .estimate)) %>%
    arrange(epitope)

# Save epitope metrics
write.csv(model1.epitope.metrics, paste0(out.path, "model1_epitope_metrics.csv"), row.names = FALSE)
write.csv(model1.epitope.metrics.with.crossreactives, paste0(out.path, "model1_epitope_metrics_with_crossreactives.csv"), row.names = FALSE)
write.csv(model2.epitope.metrics, paste0(out.path, "model2_epitope_metrics.csv"), row.names = FALSE)

# Per-model, per-epitope ROCs
# Define function to label points at best fbeta
beta <- 0.5
getfbeta <- function(data, beta = 1.0) {
    if ("specificity" %in% names(data)) {
        data %>%
            dplyr::mutate(precision = sensitivity / (sensitivity + (1 - specificity)),
                fbeta = (1 + beta^2) * (precision * sensitivity) / ((beta^2 * precision) + sensitivity))
    }
    else if ("recall" %in% names(data)) {
        data %>%
            dplyr::mutate(fbeta = (1 + beta^2) * (precision * recall) / ((beta^2 * precision) + recall))
    }
}

combined.model1.gil.roc <- combined.model1.final.fit %>%
    filter(ref.epitope == "GILGFVFTL") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL\n",
        "ROC AUC: ", round(combined.model1.final.fit %>% 
            filter(ref.epitope == "GILGFVFTL") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model1.gil.yvl.roc <- combined.model1.final.fit %>%
    filter(ref.epitope == "GILGFVFTL_YVLDHLIVV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL_YVLDHLIVV\n",
        "ROC AUC: ", round(combined.model1.final.fit %>% 
            filter(ref.epitope == "GILGFVFTL_YVLDHLIVV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model1.glc.roc <- combined.model1.final.fit %>%
    filter(ref.epitope == "GLCTLVAML") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GLCTLVAML\n",
        "ROC AUC: ", round(combined.model1.final.fit %>% 
            filter(ref.epitope == "GLCTLVAML") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model1.nlv.roc <- combined.model1.final.fit %>%
    filter(ref.epitope == "NLVPMVATV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("NLVPMVATV\n",
        "ROC AUC: ", round(combined.model1.final.fit %>% 
            filter(ref.epitope == "NLVPMVATV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model1.rpi.roc <- combined.model1.final.fit %>%
    filter(ref.epitope == "RPIIRPATL") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("RPIIRPATL\n",
        "ROC AUC: ", round(combined.model1.final.fit %>% 
            filter(ref.epitope == "RPIIRPATL") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model1.ylq.roc <- combined.model1.final.fit %>%
    filter(ref.epitope == "YLQPRTFLL") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YLQPRTFLL\n",
        "ROC AUC: ", round(combined.model1.final.fit %>% 
            filter(ref.epitope == "YLQPRTFLL") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model1.yvl.roc <- combined.model1.final.fit %>%
    filter(ref.epitope == "YVLDHLIVV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YVLDHLIVV\n",
        "ROC AUC: ", round(combined.model1.final.fit %>% 
            filter(ref.epitope == "YVLDHLIVV") %>% 
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))

combined.model1.epitope.roc <- bind_rows(combined.model1.gil.roc, 
    combined.model1.gil.yvl.roc, 
    combined.model1.glc.roc, 
    # combined.model1.nlv.roc, 
    combined.model1.rpi.roc, 
    combined.model1.ylq.roc, 
    combined.model1.yvl.roc)  %>%
    getfbeta(beta)
label.data <- combined.model1.epitope.roc %>%
    group_by(epitope) %>%
   filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
combined.model1.epitope.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, col = epitope)) + 
    geom_path(lwd = 1.5, alpha = 0.8) +
    geom_abline(lty = 3) + 
    geom_text_repel(data = label.data, 
            aes(label = epitope), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9)

ggsave(paste0(out.path, "combined_model1_roc_per_epitope.png"), width = 10, height = 10)

# Same, but accounting for cross-reactive epitopes
combined.model1.gil.roc.wc <- combined.model1.final.fit.with.crossreactives %>%
    filter(epitope == "GILGFVFTL") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL\n",
        "ROC AUC: ", round(combined.model1.final.fit.with.crossreactives %>% 
            filter(epitope == "GILGFVFTL") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model1.gil.yvl.roc.wc <- combined.model1.final.fit.with.crossreactives %>%
    filter(epitope == "GILGFVFTL_YVLDHLIVV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL_YVLDHLIVV\n",
        "ROC AUC: ", round(combined.model1.final.fit.with.crossreactives %>% 
            filter(epitope == "GILGFVFTL_YVLDHLIVV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model1.yvl.roc.wc <- combined.model1.final.fit.with.crossreactives %>%
    filter(epitope == "YVLDHLIVV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YVLDHLIVV\n",
        "ROC AUC: ", round(combined.model1.final.fit.with.crossreactives %>% 
            filter(epitope == "YVLDHLIVV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))

combined.model1.epitope.roc.wc <- bind_rows(combined.model1.gil.roc.wc, 
    combined.model1.gil.yvl.roc.wc, 
    combined.model1.glc.roc, 
    # combined.model1.nlv.roc, 
    combined.model1.rpi.roc, 
    combined.model1.ylq.roc, 
    combined.model1.yvl.roc.wc)  %>%
    getfbeta(beta)
label.data <- combined.model1.epitope.roc.wc %>%
    group_by(epitope) %>%
   filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
combined.model1.epitope.roc.wc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, col = epitope)) + 
    geom_path(lwd = 1.5, alpha = 0.8) +
    geom_abline(lty = 3) + 
    geom_text_repel(data = label.data, 
            aes(label = epitope), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9)

ggsave(paste0(out.path, "combined_model1_roc_per_epitope_with_crossreactives.png"), width = 10, height = 10)

# Structure-only model
structure.model1.gil.roc <- structure.model1.final.fit %>%
    filter(ref.epitope == "GILGFVFTL") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL\n",
        "ROC AUC: ", round(structure.model1.final.fit %>% 
            filter(ref.epitope == "GILGFVFTL") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model1.gil.yvl.roc <- structure.model1.final.fit %>%
    filter(ref.epitope == "GILGFVFTL_YVLDHLIVV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL_YVLDHLIVV\n",
        "ROC AUC: ", round(structure.model1.final.fit %>% 
            filter(ref.epitope == "GILGFVFTL_YVLDHLIVV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model1.glc.roc <- structure.model1.final.fit %>%
    filter(ref.epitope == "GLCTLVAML") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GLCTLVAML\n",
        "ROC AUC: ", round(structure.model1.final.fit %>% 
            filter(ref.epitope == "GLCTLVAML") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model1.nlv.roc <- structure.model1.final.fit %>%
    filter(ref.epitope == "NLVPMVATV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("NLVPMVATV\n",
        "ROC AUC: ", round(structure.model1.final.fit %>% 
            filter(ref.epitope == "NLVPMVATV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model1.rpi.roc <- structure.model1.final.fit %>%
    filter(ref.epitope == "RPIIRPATL") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("RPIIRPATL\n",
        "ROC AUC: ", round(structure.model1.final.fit %>% 
            filter(ref.epitope == "RPIIRPATL") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model1.ylq.roc <- structure.model1.final.fit %>%
    filter(ref.epitope == "YLQPRTFLL") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YLQPRTFLL\n",
        "ROC AUC: ", round(structure.model1.final.fit %>% 
            filter(ref.epitope == "YLQPRTFLL") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model1.yvl.roc <- structure.model1.final.fit %>%
    filter(ref.epitope == "YVLDHLIVV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YVLDHLIVV\n",
        "ROC AUC: ", round(structure.model1.final.fit %>% 
            filter(ref.epitope == "YVLDHLIVV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))

structure.model1.epitope.roc <- bind_rows(structure.model1.gil.roc, 
    structure.model1.gil.yvl.roc, 
    structure.model1.glc.roc, 
    # structure.model1.nlv.roc, 
    structure.model1.rpi.roc,
    structure.model1.ylq.roc, 
    structure.model1.yvl.roc)  %>%
    getfbeta(beta)
label.data <- structure.model1.epitope.roc %>%
    group_by(epitope) %>%
   filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
structure.model1.epitope.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, col = epitope)) + 
    geom_path(lwd = 1.5, alpha = 0.8) +
    geom_abline(lty = 3) + 
    geom_text_repel(data = label.data, 
            aes(label = epitope), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9)

ggsave(paste0(out.path, "structure_model1_roc_per_epitope.png"), width = 10, height = 10)

# Same, but accounting for cross-reactive epitopes
structure.model1.gil.roc.wc <- structure.model1.final.fit.with.crossreactives %>%
    filter(epitope == "GILGFVFTL") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL\n",
        "ROC AUC: ", round(structure.model1.final.fit.with.crossreactives %>% 
            filter(epitope == "GILGFVFTL") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model1.gil.yvl.roc.wc <- structure.model1.final.fit.with.crossreactives %>%
    filter(epitope == "GILGFVFTL_YVLDHLIVV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL_YVLDHLIVV\n",
        "ROC AUC: ", round(structure.model1.final.fit.with.crossreactives %>% 
            filter(epitope == "GILGFVFTL_YVLDHLIVV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model1.yvl.roc.wc <- structure.model1.final.fit.with.crossreactives %>%
    filter(epitope == "YVLDHLIVV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YVLDHLIVV\n",
        "ROC AUC: ", round(structure.model1.final.fit.with.crossreactives %>% 
            filter(epitope == "YVLDHLIVV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))

structure.model1.epitope.roc.wc <- bind_rows(structure.model1.gil.roc.wc, 
    structure.model1.gil.yvl.roc.wc, 
    structure.model1.glc.roc, 
    # structure.model1.nlv.roc, 
    structure.model1.rpi.roc, 
    structure.model1.ylq.roc, 
    structure.model1.yvl.roc.wc)  %>%
    getfbeta(beta)
label.data <- structure.model1.epitope.roc.wc %>%
    group_by(epitope) %>%
   filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
structure.model1.epitope.roc.wc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, col = epitope)) + 
    geom_path(lwd = 1.5, alpha = 0.8) +
    geom_abline(lty = 3) + 
    geom_text_repel(data = label.data, 
            aes(label = epitope), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9)

ggsave(paste0(out.path, "structure_model1_roc_per_epitope_with_crossreactives.png"), width = 10, height = 10)

# Sequence-only model
sequence.model1.gil.roc <- sequence.model1.final.fit %>%
    filter(ref.epitope == "GILGFVFTL") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL\n",
        "ROC AUC: ", round(sequence.model1.final.fit %>% 
            filter(ref.epitope == "GILGFVFTL") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model1.gil.yvl.roc <- sequence.model1.final.fit %>%
    filter(ref.epitope == "GILGFVFTL_YVLDHLIVV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL_YVLDHLIVV\n",
        "ROC AUC: ", round(sequence.model1.final.fit %>% 
            filter(ref.epitope == "GILGFVFTL_YVLDHLIVV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model1.glc.roc <- sequence.model1.final.fit %>%
    filter(ref.epitope == "GLCTLVAML") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GLCTLVAML\n",
        "ROC AUC: ", round(sequence.model1.final.fit %>% 
            filter(ref.epitope == "GLCTLVAML") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model1.nlv.roc <- sequence.model1.final.fit %>%
    filter(ref.epitope == "NLVPMVATV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("NLVPMVATV\n",
        "ROC AUC: ", round(sequence.model1.final.fit %>% 
            filter(ref.epitope == "NLVPMVATV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model1.rpi.roc <- sequence.model1.final.fit %>%
    filter(ref.epitope == "RPIIRPATL") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("RPIIRPATL\n",
        "ROC AUC: ", round(sequence.model1.final.fit %>% 
            filter(ref.epitope == "RPIIRPATL") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model1.ylq.roc <- sequence.model1.final.fit %>%
    filter(ref.epitope == "YLQPRTFLL") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YLQPRTFLL\n",
        "ROC AUC: ", round(sequence.model1.final.fit %>% 
            filter(ref.epitope == "YLQPRTFLL") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model1.yvl.roc <- sequence.model1.final.fit %>%
    filter(ref.epitope == "YVLDHLIVV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YVLDHLIVV\n",
        "ROC AUC: ", round(sequence.model1.final.fit %>% 
            filter(ref.epitope == "YVLDHLIVV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))

sequence.model1.epitope.roc <- bind_rows(sequence.model1.gil.roc, 
    sequence.model1.gil.yvl.roc, 
    sequence.model1.glc.roc, 
    # sequence.model1.nlv.roc, 
    sequence.model1.rpi.roc, 
    sequence.model1.ylq.roc, 
    sequence.model1.yvl.roc)  %>%
    getfbeta(beta)
label.data <- sequence.model1.epitope.roc %>%
    group_by(epitope) %>%
   filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
sequence.model1.epitope.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, col = epitope)) + 
    geom_path(lwd = 1.5, alpha = 0.8) +
    geom_abline(lty = 3) + 
    geom_text_repel(data = label.data, 
            aes(label = epitope), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9)

ggsave(paste0(out.path, "sequence_model1_roc_per_epitope.png"), width = 10, height = 10)

# Same, but accounting for cross-reactive epitopes
sequence.model1.gil.roc.wc <- sequence.model1.final.fit.with.crossreactives %>%
    filter(epitope == "GILGFVFTL") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL\n",
        "ROC AUC: ", round(sequence.model1.final.fit.with.crossreactives %>% 
            filter(epitope == "GILGFVFTL") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model1.gil.yvl.roc.wc <- sequence.model1.final.fit.with.crossreactives %>%
    filter(epitope == "GILGFVFTL_YVLDHLIVV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL_YVLDHLIVV\n",
        "ROC AUC: ", round(sequence.model1.final.fit.with.crossreactives %>% 
            filter(epitope == "GILGFVFTL_YVLDHLIVV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model1.yvl.roc.wc <- sequence.model1.final.fit.with.crossreactives %>%
    filter(epitope == "YVLDHLIVV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YVLDHLIVV\n",
        "ROC AUC: ", round(sequence.model1.final.fit.with.crossreactives %>% 
            filter(epitope == "YVLDHLIVV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))

sequence.model1.epitope.roc.wc <- bind_rows(sequence.model1.gil.roc.wc, 
    sequence.model1.gil.yvl.roc.wc, 
    sequence.model1.glc.roc, 
    # sequence.model1.nlv.roc, 
    sequence.model1.rpi.roc, 
    sequence.model1.ylq.roc, 
    sequence.model1.yvl.roc.wc)  %>%
    getfbeta(beta)
label.data <- sequence.model1.epitope.roc.wc %>%
    group_by(epitope) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
sequence.model1.epitope.roc.wc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, col = epitope)) + 
    geom_path(lwd = 1.5, alpha = 0.8) +
    geom_abline(lty = 3) + 
    geom_text_repel(data = label.data, 
            aes(label = epitope), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9)

ggsave(paste0(out.path, "sequence_model1_roc_per_epitope_with_crossreactives.png"), width = 10, height = 10)


# Per-epitope, per-model ROC curves
# Store epitope colors from above plots
epitopes <- piema.data %>% distinct(ref.epitope) %>% pull(ref.epitope)
epitope.colors <- viridis(length(epitopes), option = "plasma", begin = 0.1, end = 0.9)
epitope.colors <- setNames(epitope.colors, epitopes)

# Bind data
model1.gil.roc <- bind_rows(combined.model1.gil.roc %>% mutate(model = "Combined"), structure.model1.gil.roc %>% mutate(model = "Structure"), sequence.model1.gil.roc %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model1.gil.roc %>% 
    group_by(epitope, model) %>% 
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>% 
    slice(1)
model1.gil.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.5, color = epitope.colors[["GILGFVFTL"]]) +
    geom_abline(lty = 3) + 
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["GILGFVFTL"]]) + 
    coord_equal() + 
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_GIL_roc.png"), width = 20, height = 20)

# With cross-reactives
model1.gil.roc.wc <- bind_rows(combined.model1.gil.roc.wc %>% mutate(model = "Combined"), structure.model1.gil.roc.wc %>% mutate(model = "Structure"), sequence.model1.gil.roc.wc %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model1.gil.roc.wc %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model1.gil.roc.wc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.5, color = epitope.colors[["GILGFVFTL"]]) +
    geom_abline(lty = 3) + 
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["GILGFVFTL"]]) + 
    coord_equal() + 
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_GIL_roc_with_crossreactives.png"), width = 20, height = 20)


model1.gil.yvl.roc <- bind_rows(combined.model1.gil.yvl.roc %>% mutate(model = "Combined"), structure.model1.gil.yvl.roc %>% mutate(model = "Structure"), sequence.model1.gil.yvl.roc %>% mutate(model = "Sequence")) %>%
    getfbeta(beta)
label.data <- model1.gil.yvl.roc %>% 
    group_by(epitope, model) %>% 
   filter(fbeta == max(fbeta, na.rm = TRUE)) %>% 
    slice(1)
model1.gil.yvl.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.5, color = epitope.colors[["GILGFVFTL_YVLDHLIVV"]]) +
    geom_abline(lty = 3) + 
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["GILGFVFTL_YVLDHLIVV"]]) + 
    coord_equal() + 
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_GIL-YVL_roc.png"), width = 20, height = 20)

# With cross-reactives
model1.gil.yvl.roc.wc <- bind_rows(combined.model1.gil.yvl.roc.wc %>% mutate(model = "Combined"), structure.model1.gil.yvl.roc.wc %>% mutate(model = "Structure"), sequence.model1.gil.yvl.roc.wc %>% mutate(model = "Sequence")) %>%
    getfbeta(beta)
label.data <- model1.gil.yvl.roc.wc %>% 
    group_by(epitope, model) %>% 
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>% 
    slice(1)
model1.gil.yvl.roc.wc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.5, color = epitope.colors[["GILGFVFTL_YVLDHLIVV"]]) +
    geom_abline(lty = 3) + 
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["GILGFVFTL_YVLDHLIVV"]]) + 
    coord_equal() + 
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_GIL-YVL_roc_with_crossreactives.png"), width = 20, height = 20)

model1.glc.roc <- bind_rows(combined.model1.glc.roc %>% mutate(model = "Combined"), structure.model1.glc.roc %>% mutate(model = "Structure"), sequence.model1.glc.roc %>% mutate(model = "Sequence")) %>%
    getfbeta(beta)
label.data <- model1.glc.roc %>% 
    group_by(epitope, model) %>% 
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%   
    slice(1)
model1.glc.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.5, color = epitope.colors[["GLCTLVAML"]]) +
    geom_abline(lty = 3) + 
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["GLCTLVAML"]]) + 
    coord_equal() + 
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_GLC_roc.png"), width = 20, height = 20)

model1.nlv.roc <- bind_rows(combined.model1.nlv.roc %>% mutate(model = "Combined"), structure.model1.nlv.roc %>% mutate(model = "Structure"), sequence.model1.nlv.roc %>% mutate(model = "Sequence")) %>%
    getfbeta(beta)
label.data <- model1.nlv.roc %>% 
    group_by(epitope, model) %>% 
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model1.nlv.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, linetype = model)) +
    geom_path(lwd = 1.5, alpha = 0.5, color = epitope.colors[["NLVPMVATV"]]) +
    geom_abline(lty = 3) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["NLVPMVATV"]]) +
    coord_equal() +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_NLV_roc.png"), width = 20, height = 20)

model1.rpi.roc <- bind_rows(combined.model1.rpi.roc %>% mutate(model = "Combined"), structure.model1.rpi.roc %>% mutate(model = "Structure"), sequence.model1.rpi.roc %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model1.rpi.roc %>% 
    group_by(epitope, model) %>% 
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model1.rpi.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.5, color = epitope.colors[["RPIIRPATL"]]) +
    geom_abline(lty = 3) + 
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["RPIIRPATL"]]) + 
    coord_equal() + 
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_RPI_roc.png"), width = 20, height = 20)

model1.ylq.roc <- bind_rows(combined.model1.ylq.roc %>% mutate(model = "Combined"), structure.model1.ylq.roc %>% mutate(model = "Structure"), sequence.model1.ylq.roc %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model1.ylq.roc %>% 
    group_by(epitope, model) %>% 
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>% 
    slice(1)
model1.ylq.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.5, color = epitope.colors[["YLQPRTFLL"]]) +
    geom_abline(lty = 3) + 
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["YLQPRTFLL"]]) + 
    coord_equal() + 
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_YLQ_roc.png"), width = 20, height = 20)

model1.yvl.roc <- bind_rows(combined.model1.yvl.roc %>% mutate(model = "Combined"), structure.model1.yvl.roc %>% mutate(model = "Structure"), sequence.model1.yvl.roc %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model1.yvl.roc %>% 
    group_by(epitope, model) %>% 
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model1.yvl.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.5, color = epitope.colors[["YVLDHLIVV"]]) +
    geom_abline(lty = 3) + 
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["YVLDHLIVV"]]) + 
    coord_equal() + 
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_YVL_roc.png"), width = 20, height = 20)

# With cross-reactives
model1.yvl.roc.wc <- bind_rows(combined.model1.yvl.roc.wc %>% mutate(model = "Combined"), structure.model1.yvl.roc.wc %>% mutate(model = "Structure"), sequence.model1.yvl.roc.wc %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model1.yvl.roc.wc %>% 
    group_by(epitope, model) %>% 
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model1.yvl.roc.wc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.5, color = epitope.colors[["YVLDHLIVV"]]) +
    geom_abline(lty = 3) + 
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["YVLDHLIVV"]]) + 
    coord_equal() + 
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_YVL_roc_with_crossreactives.png"), width = 20, height = 20)


# All strategies and epitopes combined (messy but informative, for now)
# Combined model is solid line
# Structure model is dashed line
# Sequence model is dotted line
# MAKE PER EPITOPE NOT PER MODEL STRATEGY
model1.epitope.roc <- bind_rows(combined.model1.epitope.roc %>% mutate(model = "Combined"), structure.model1.epitope.roc %>% mutate(model = "Structure"), sequence.model1.epitope.roc %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta) %>%
    mutate(epitope = gsub("\nPR AUC: [0-9.]+", "", epitope))
label.data <- model1.epitope.roc %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model1.epitope.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, col = epitope, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.5) +
    geom_abline(lty = 3) + 
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_roc_per_epitope.png"), width = 20, height = 20)

# With cross-reactives (and remove the the ROC AUCs)
model1.epitope.roc.wc <- bind_rows(combined.model1.epitope.roc.wc %>% mutate(model = "Combined"), structure.model1.epitope.roc.wc %>% mutate(model = "Structure"), sequence.model1.epitope.roc.wc %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta) %>%
    mutate(epitope = gsub("\nPR AUC: [0-9.]+", "", epitope))
label.data <- model1.epitope.roc.wc %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model1.epitope.roc.wc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, col = epitope, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.5) +
    geom_abline(lty = 3) + 
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_roc_per_epitope_with_crossreactives.png"), width = 20, height = 20)


# More PR curves
# Per model, per epitope
combined.model1.gil.pr <- combined.model1.final.fit %>%
    filter(ref.epitope == "GILGFVFTL") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL\n",
        "PR AUC: ", round(combined.model1.final.fit %>% 
            filter(ref.epitope == "GILGFVFTL") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model1.gil.yvl.pr <- combined.model1.final.fit %>%
    filter(ref.epitope == "GILGFVFTL_YVLDHLIVV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL_YVLDHLIVV\n",
        "PR AUC: ", round(combined.model1.final.fit %>% 
            filter(ref.epitope == "GILGFVFTL_YVLDHLIVV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model1.glc.pr <- combined.model1.final.fit %>%
    filter(ref.epitope == "GLCTLVAML") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GLCTLVAML\n",
        "PR AUC: ", round(combined.model1.final.fit %>% 
            filter(ref.epitope == "GLCTLVAML") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model1.nlv.pr <- combined.model1.final.fit %>%
    filter(ref.epitope == "NLVPMVATV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("NLVPMVATV\n",
        "PR AUC: ", round(combined.model1.final.fit %>% 
            filter(ref.epitope == "NLVPMVATV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model1.rpi.pr <- combined.model1.final.fit %>%
    filter(ref.epitope == "RPIIRPATL") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("RPIIRPATL\n",
        "PR AUC: ", round(combined.model1.final.fit %>% 
            filter(ref.epitope == "RPIIRPATL") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model1.ylq.pr <- combined.model1.final.fit %>%
    filter(ref.epitope == "YLQPRTFLL") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YLQPRTFLL\n",
        "PR AUC: ", round(combined.model1.final.fit %>% 
            filter(ref.epitope == "YLQPRTFLL") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model1.yvl.pr <- combined.model1.final.fit %>%
    filter(ref.epitope == "YVLDHLIVV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YVLDHLIVV\n",
        "PR AUC: ", round(combined.model1.final.fit %>% 
            filter(ref.epitope == "YVLDHLIVV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))

combined.model1.epitope.pr <- bind_rows(combined.model1.gil.pr, 
    combined.model1.gil.yvl.pr, 
    combined.model1.glc.pr, 
    # combined.model1.nlv.pr, 
    combined.model1.rpi.pr, 
    combined.model1.ylq.pr, 
    combined.model1.yvl.pr)  %>%
    getfbeta(beta)
label.data <- combined.model1.epitope.pr %>%
    group_by(epitope) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
combined.model1.epitope.pr %>%
    ggplot(aes(x = recall, y = precision, color = epitope)) + 
    geom_path(lwd = 1.5, alpha = 0.8) +
    geom_text_repel(data = label.data, 
            aes(label = epitope), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9)

ggsave(paste0(out.path, "combined_model1_pr_per_epitope.png"), width = 10, height = 10)

# With cross-reactives
combined.model1.gil.pr.wc <- combined.model1.final.fit.with.crossreactives %>%
    filter(epitope == "GILGFVFTL") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL\n",
        "PR AUC: ", round(combined.model1.final.fit.with.crossreactives %>% 
            filter(epitope == "GILGFVFTL") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model1.gil.yvl.pr.wc <- combined.model1.final.fit.with.crossreactives %>%
    filter(epitope == "GILGFVFTL_YVLDHLIVV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL_YVLDHLIVV\n",
        "PR AUC: ", round(combined.model1.final.fit.with.crossreactives %>% 
            filter(epitope == "GILGFVFTL_YVLDHLIVV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model1.yvl.pr.wc <- combined.model1.final.fit.with.crossreactives %>%
    filter(epitope == "YVLDHLIVV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YVLDHLIVV\n",
        "PR AUC: ", round(combined.model1.final.fit.with.crossreactives %>% 
            filter(epitope == "YVLDHLIVV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))

combined.model1.epitope.pr.wc <- bind_rows(combined.model1.gil.pr.wc, 
    combined.model1.gil.yvl.pr.wc, 
    combined.model1.glc.pr, 
    # combined.model1.nlv.pr, 
    combined.model1.rpi.pr, 
    combined.model1.ylq.pr, 
    combined.model1.yvl.pr.wc)  %>%
    getfbeta(beta)
label.data <- combined.model1.epitope.pr.wc %>%
    group_by(epitope) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
combined.model1.epitope.pr.wc %>%
    ggplot(aes(x = recall, y = precision, color = epitope)) + 
    geom_path(lwd = 1.5, alpha = 0.8) +
    geom_text_repel(data = label.data, 
            aes(label = epitope), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9)

ggsave(paste0(out.path, "combined_model1_pr_per_epitope_with_crossreactives.png"), width = 10, height = 10)


# Structure model
structure.model1.gil.pr <- structure.model1.final.fit %>%
    filter(ref.epitope == "GILGFVFTL") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL\n",
        "PR AUC: ", round(structure.model1.final.fit %>% 
            filter(ref.epitope == "GILGFVFTL") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model1.gil.yvl.pr <- structure.model1.final.fit %>%
    filter(ref.epitope == "GILGFVFTL_YVLDHLIVV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL_YVLDHLIVV\n",
        "PR AUC: ", round(structure.model1.final.fit %>% 
            filter(ref.epitope == "GILGFVFTL_YVLDHLIVV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model1.glc.pr <- structure.model1.final.fit %>%
    filter(ref.epitope == "GLCTLVAML") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GLCTLVAML\n",
        "PR AUC: ", round(structure.model1.final.fit %>% 
            filter(ref.epitope == "GLCTLVAML") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model1.nlv.pr <- structure.model1.final.fit %>%
    filter(ref.epitope == "NLVPMVATV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("NLVPMVATV\n",
        "PR AUC: ", round(structure.model1.final.fit %>% 
            filter(ref.epitope == "NLVPMVATV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model1.rpi.pr <- structure.model1.final.fit %>%
    filter(ref.epitope == "RPIIRPATL") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("RPIIRPATL\n",
        "PR AUC: ", round(structure.model1.final.fit %>%
            filter(ref.epitope == "RPIIRPATL") %>%   
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model1.ylq.pr <- structure.model1.final.fit %>%
    filter(ref.epitope == "YLQPRTFLL") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YLQPRTFLL\n",
        "PR AUC: ", round(structure.model1.final.fit %>% 
            filter(ref.epitope == "YLQPRTFLL") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model1.yvl.pr <- structure.model1.final.fit %>%
    filter(ref.epitope == "YVLDHLIVV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YVLDHLIVV\n",
        "PR AUC: ", round(structure.model1.final.fit %>% 
            filter(ref.epitope == "YVLDHLIVV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))

structure.model1.epitope.pr <- bind_rows(structure.model1.gil.pr, 
    structure.model1.gil.yvl.pr, 
    structure.model1.glc.pr, 
    # structure.model1.nlv.pr, 
    structure.model1.rpi.pr, 
    structure.model1.ylq.pr, 
    structure.model1.yvl.pr)  %>%
    getfbeta(beta)
label.data <- structure.model1.epitope.pr %>%
    group_by(epitope) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
structure.model1.epitope.pr %>%
    ggplot(aes(x = recall, y = precision, col = epitope)) + 
    geom_path(lwd = 1.5, alpha = 0.8) +
    geom_text_repel(data = label.data, 
            aes(label = epitope), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9)

ggsave(paste0(out.path, "structure_model1_pr_per_epitope.png"), width = 10, height = 10)

# With cross-reactives
structure.model1.gil.pr.wc <- structure.model1.final.fit.with.crossreactives %>%
    filter(epitope == "GILGFVFTL") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL\n",
        "PR AUC: ", round(structure.model1.final.fit.with.crossreactives %>% 
            filter(epitope == "GILGFVFTL") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model1.gil.yvl.pr.wc <- structure.model1.final.fit.with.crossreactives %>%
    filter(epitope == "GILGFVFTL_YVLDHLIVV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL_YVLDHLIVV\n",
        "PR AUC: ", round(structure.model1.final.fit.with.crossreactives %>% 
            filter(epitope == "GILGFVFTL_YVLDHLIVV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model1.yvl.pr.wc <- structure.model1.final.fit.with.crossreactives %>%
    filter(epitope == "YVLDHLIVV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YVLDHLIVV\n",
        "PR AUC: ", round(structure.model1.final.fit.with.crossreactives %>% 
            filter(epitope == "YVLDHLIVV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))

structure.model1.epitope.pr.wc <- bind_rows(structure.model1.gil.pr.wc, 
    structure.model1.gil.yvl.pr.wc, 
    structure.model1.glc.pr, 
    # structure.model1.nlv.pr, 
    structure.model1.rpi.pr, 
    structure.model1.ylq.pr, 
    structure.model1.yvl.pr.wc)  %>%
    getfbeta(beta)
label.data <- structure.model1.epitope.pr.wc %>%
    group_by(epitope) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
structure.model1.epitope.pr.wc %>%
    ggplot(aes(x = recall, y = precision, col = epitope)) + 
    geom_path(lwd = 1.5, alpha = 0.8) +
    geom_text_repel(data = label.data, 
            aes(label = epitope), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9)

ggsave(paste0(out.path, "structure_model1_pr_per_epitope_with_crossreactives.png"), width = 10, height = 10)


# Sequence model
sequence.model1.gil.pr <- sequence.model1.final.fit %>%
    filter(ref.epitope == "GILGFVFTL") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL\n",
        "PR AUC: ", round(sequence.model1.final.fit %>% 
            filter(ref.epitope == "GILGFVFTL") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model1.gil.yvl.pr <- sequence.model1.final.fit %>%
    filter(ref.epitope == "GILGFVFTL_YVLDHLIVV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL_YVLDHLIVV\n",
        "PR AUC: ", round(sequence.model1.final.fit %>% 
            filter(ref.epitope == "GILGFVFTL_YVLDHLIVV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model1.glc.pr <- sequence.model1.final.fit %>%
    filter(ref.epitope == "GLCTLVAML") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GLCTLVAML\n",
        "PR AUC: ", round(sequence.model1.final.fit %>% 
            filter(ref.epitope == "GLCTLVAML") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model1.nlv.pr <- sequence.model1.final.fit %>%
    filter(ref.epitope == "NLVPMVATV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("NLVPMVATV\n",
        "PR AUC: ", round(sequence.model1.final.fit %>% 
            filter(ref.epitope == "NLVPMVATV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model1.rpi.pr <- sequence.model1.final.fit %>%
    filter(ref.epitope == "RPIIRPATL") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("RPIIRPATL\n",
        "PR AUC: ", round(sequence.model1.final.fit %>% 
            filter(ref.epitope == "RPIIRPATL") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model1.ylq.pr <- sequence.model1.final.fit %>%
    filter(ref.epitope == "YLQPRTFLL") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YLQPRTFLL\n",
        "PR AUC: ", round(sequence.model1.final.fit %>% 
            filter(ref.epitope == "YLQPRTFLL") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model1.yvl.pr <- sequence.model1.final.fit %>%
    filter(ref.epitope == "YVLDHLIVV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YVLDHLIVV\n",
        "PR AUC: ", round(sequence.model1.final.fit %>% 
            filter(ref.epitope == "YVLDHLIVV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))

sequence.model1.epitope.pr <- bind_rows(sequence.model1.gil.pr, 
    sequence.model1.gil.yvl.pr, 
    sequence.model1.glc.pr, 
    # sequence.model1.nlv.pr, 
    sequence.model1.rpi.pr, 
    sequence.model1.ylq.pr, 
    sequence.model1.yvl.pr)  %>%
    getfbeta(beta)
label.data <- sequence.model1.epitope.pr %>%
    group_by(epitope) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
sequence.model1.epitope.pr %>%
    ggplot(aes(x = recall, y = precision, col = epitope)) + 
    geom_path(lwd = 1.5, alpha = 0.8) +
    geom_text_repel(data = label.data, 
            aes(label = epitope), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9)

ggsave(paste0(out.path, "sequence_model1_pr_per_epitope.png"), width = 10, height = 10)

# With cross-reactives
sequence.model1.gil.pr.wc <- sequence.model1.final.fit.with.crossreactives %>%
    filter(epitope == "GILGFVFTL") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL\n",
        "PR AUC: ", round(sequence.model1.final.fit.with.crossreactives %>% 
            filter(epitope == "GILGFVFTL") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model1.gil.yvl.pr.wc <- sequence.model1.final.fit.with.crossreactives %>%
    filter(epitope == "GILGFVFTL_YVLDHLIVV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL_YVLDHLIVV\n",
        "PR AUC: ", round(sequence.model1.final.fit.with.crossreactives %>% 
            filter(epitope == "GILGFVFTL_YVLDHLIVV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model1.yvl.pr.wc <- sequence.model1.final.fit.with.crossreactives %>%
    filter(epitope == "YVLDHLIVV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YVLDHLIVV\n",
        "PR AUC: ", round(sequence.model1.final.fit.with.crossreactives %>% 
            filter(epitope == "YVLDHLIVV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))

sequence.model1.epitope.pr.wc <- bind_rows(sequence.model1.gil.pr.wc, 
    sequence.model1.gil.yvl.pr.wc, 
    sequence.model1.glc.pr, 
    # sequence.model1.nlv.pr, 
    sequence.model1.rpi.pr, 
    sequence.model1.ylq.pr, 
    sequence.model1.yvl.pr.wc)  %>%
    getfbeta(beta)
label.data <- sequence.model1.epitope.pr.wc %>%
    group_by(epitope) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
sequence.model1.epitope.pr.wc %>%
    ggplot(aes(x = recall, y = precision, col = epitope)) + 
    geom_path(lwd = 1.5, alpha = 0.8) +
    geom_text_repel(data = label.data, 
            aes(label = epitope), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9)

ggsave(paste0(out.path, "sequence_model1_pr_per_epitope_with_crossreactives.png"), width = 10, height = 10)

# Per epitope, per model
model1.gil.pr <- bind_rows(combined.model1.gil.pr %>% mutate(model = "Combined"), structure.model1.gil.pr %>% mutate(model = "Structure"), sequence.model1.gil.pr %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model1.gil.pr %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model1.gil.pr %>%
    ggplot(aes(x = recall, y = precision, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["GILGFVFTL"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["GILGFVFTL"]]) + 
    coord_equal() + 
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_GIL_pr.png"), width = 10, height = 10)

# With cross-reactives
model1.gil.pr.wc <- bind_rows(combined.model1.gil.pr.wc %>% mutate(model = "Combined"), structure.model1.gil.pr.wc %>% mutate(model = "Structure"), sequence.model1.gil.pr.wc %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model1.gil.pr.wc %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model1.gil.pr.wc %>%
    ggplot(aes(x = recall, y = precision, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["GILGFVFTL"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["GILGFVFTL"]]) + 
    coord_equal() + 
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_GIL_pr_with_crossreactives.png"), width = 10, height = 10)

model1.gil.yvl.pr <- bind_rows(combined.model1.gil.yvl.pr %>% mutate(model = "Combined"), structure.model1.gil.yvl.pr %>% mutate(model = "Structure"), sequence.model1.gil.yvl.pr %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model1.gil.yvl.pr %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model1.gil.yvl.pr %>%
    ggplot(aes(x = recall, y = precision, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["GILGFVFTL_YVLDHLIVV"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["GILGFVFTL_YVLDHLIVV"]]) +
    coord_equal() +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_GIL-YVL_pr.png"), width = 10, height = 10)

# With cross-reactives
model1.gil.yvl.pr.wc <- bind_rows(combined.model1.gil.yvl.pr.wc %>% mutate(model = "Combined"), structure.model1.gil.yvl.pr.wc %>% mutate(model = "Structure"), sequence.model1.gil.yvl.pr.wc %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model1.gil.yvl.pr.wc %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model1.gil.yvl.pr.wc %>%
    ggplot(aes(x = recall, y = precision, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["GILGFVFTL_YVLDHLIVV"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["GILGFVFTL_YVLDHLIVV"]]) +
    coord_equal() +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_GIL-YVL_pr_with_crossreactives.png"), width = 10, height = 10)


model1.glc.pr <- bind_rows(combined.model1.glc.pr %>% mutate(model = "Combined"), structure.model1.glc.pr %>% mutate(model = "Structure"), sequence.model1.glc.pr %>% mutate(model = "Sequence"))    %>%
    getfbeta(beta)
label.data <- model1.glc.pr %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model1.glc.pr %>%
    ggplot(aes(x = recall, y = precision, linetype = model)) +
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["GLCTLVAML"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["GLCTLVAML"]]) +
    coord_equal() +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_GLC_pr.png"), width = 10, height = 10)

model1.nlv.pr <- bind_rows(combined.model1.nlv.pr %>% mutate(model = "Combined"), structure.model1.nlv.pr %>% mutate(model = "Structure"), sequence.model1.nlv.pr %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model1.nlv.pr %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model1.nlv.pr %>%
    ggplot(aes(x = recall, y = precision, linetype = model)) +
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["NLVPMVATV"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["NLVPMVATV"]]) +
    coord_equal() +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_NLV_pr.png"), width = 10, height = 10)

model1.rpi.pr <- bind_rows(combined.model1.rpi.pr %>% mutate(model = "Combined"), structure.model1.rpi.pr %>% mutate(model = "Structure"), sequence.model1.rpi.pr %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model1.rpi.pr %>% 
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model1.rpi.pr %>%
    ggplot(aes(x = recall, y = precision, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["RPIIRPATL"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["RPIIRPATL"]]) +
    coord_equal() +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_RPI_pr.png"), width = 10, height = 10)

model1.ylq.pr <- bind_rows(combined.model1.ylq.pr %>% mutate(model = "Combined"), structure.model1.ylq.pr %>% mutate(model = "Structure"), sequence.model1.ylq.pr %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model1.ylq.pr %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model1.ylq.pr %>%
    ggplot(aes(x = recall, y = precision, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["YLQPRTFLL"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["YLQPRTFLL"]]) +
    coord_equal() +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_YLQ_pr.png"), width = 10, height = 10)

model1.yvl.pr <- bind_rows(combined.model1.yvl.pr %>% mutate(model = "Combined"), structure.model1.yvl.pr %>% mutate(model = "Structure"), sequence.model1.yvl.pr %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model1.yvl.pr %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model1.yvl.pr %>%
    ggplot(aes(x = recall, y = precision, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["YVLDHLIVV"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["YVLDHLIVV"]]) +
    coord_equal() +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_YVL_pr.png"), width = 10, height = 10)

# With cross-reactives
model1.yvl.pr.wc <- bind_rows(combined.model1.yvl.pr.wc %>% mutate(model = "Combined"), structure.model1.yvl.pr.wc %>% mutate(model = "Structure"), sequence.model1.yvl.pr.wc %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model1.yvl.pr.wc %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model1.yvl.pr.wc %>%
    ggplot(aes(x = recall, y = precision, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["YVLDHLIVV"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["YVLDHLIVV"]]) +
    coord_equal() +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_YVL_pr_with_crossreactives.png"), width = 10, height = 10)


# All strategies and epitopes combined (messy but informative, for now)
# Combined model is solid line
# Structure model is dashed line
# Sequence model is dotted line
model1.epitope.pr <- bind_rows(combined.model1.epitope.pr %>% mutate(model = "Combined"), structure.model1.epitope.pr %>% mutate(model = "Structure"), sequence.model1.epitope.pr %>% mutate(model = "Sequence")) %>%
    getfbeta(beta) %>%
    mutate(epitope = gsub("\nROC AUC: [0-9.]+", "", epitope))
label.data <- model1.epitope.pr %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model1.epitope.pr %>%
    ggplot(aes(x = recall, y = precision, col = epitope, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.5) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_pr_per_epitope.png"), width = 20, height = 20)

# With cross-reactives
model1.epitope.pr.wc <- bind_rows(combined.model1.epitope.pr.wc %>% mutate(model = "Combined"), structure.model1.epitope.pr.wc %>% mutate(model = "Structure"), sequence.model1.epitope.pr.wc %>% mutate(model = "Sequence")) %>%
    getfbeta(beta) %>%
    mutate(epitope = gsub("\nROC AUC: [0-9.]+", "", epitope))
label.data <- model1.epitope.pr.wc %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model1.epitope.pr.wc %>%
    ggplot(aes(x = recall, y = precision, col = epitope, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.5) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model1_pr_per_epitope_with_crossreactives.png"), width = 20, height = 20)


# Model 2 ROC
# Per model, per epitope
combined.model2.gil.roc <- combined.model2.final.fit.expanded %>%
    filter(epitope == "GILGFVFTL") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL\n",
        "ROC AUC: ", round(combined.model2.final.fit.expanded %>%
            filter(epitope == "GILGFVFTL") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model2.gil.yvl.roc <- combined.model2.final.fit.expanded %>%
    filter(epitope == "GILGFVFTL_YVLDHLIVV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL_YVLDHLIVV\n",
        "ROC AUC: ", round(combined.model2.final.fit.expanded %>% 
            filter(epitope == "GILGFVFTL_YVLDHLIVV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model2.glc.roc <- combined.model2.final.fit.expanded %>%
    filter(epitope == "GLCTLVAML") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GLCTLVAML\n",
        "ROC AUC: ", round(combined.model2.final.fit.expanded %>% 
            filter(epitope == "GLCTLVAML") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model2.nlv.roc <- combined.model2.final.fit.expanded %>%
    filter(epitope == "NLVPMVATV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("NLVPMVATV\n",
        "ROC AUC: ", round(combined.model2.final.fit.expanded %>% 
            filter(epitope == "NLVPMVATV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model2.rpi.roc <- combined.model2.final.fit.expanded %>%
    filter(epitope == "RPIIRPATL") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("RPIIRPATL\n",
        "ROC AUC: ", round(combined.model2.final.fit.expanded %>% 
            filter(epitope == "RPIIRPATL") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model2.ylq.roc <- combined.model2.final.fit.expanded %>%
    filter(epitope == "YLQPRTFLL") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YLQPRTFLL\n",
        "ROC AUC: ", round(combined.model2.final.fit.expanded %>% 
            filter(epitope == "YLQPRTFLL") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model2.yvl.roc <- combined.model2.final.fit.expanded %>%
    filter(epitope == "YVLDHLIVV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YVLDHLIVV\n",
        "ROC AUC: ", round(combined.model2.final.fit.expanded %>% 
            filter(epitope == "YVLDHLIVV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))

combined.model2.epitope.roc <- bind_rows(combined.model2.gil.roc, 
    combined.model2.gil.yvl.roc, 
    combined.model2.glc.roc, 
    # combined.model2.nlv.roc, 
    combined.model2.rpi.roc, 
    combined.model2.ylq.roc, 
    combined.model2.yvl.roc) %>%
    getfbeta(beta)
label.data <- combined.model2.epitope.roc %>%
    group_by(epitope) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
combined.model2.epitope.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, col = epitope)) + 
    geom_path(lwd = 1.5, alpha = 0.8) +
    geom_text_repel(data = label.data, 
            aes(label = epitope), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9)

ggsave(paste0(out.path, "combined_model2_roc_per_epitope.png"), width = 10, height = 10)

# Structure model
structure.model2.gil.roc <- structure.model2.final.fit %>%
    filter(ref.epitope == "GILGFVFTL") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL\n",
        "ROC AUC: ", round(structure.model2.final.fit %>% 
            filter(ref.epitope == "GILGFVFTL") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model2.gil.yvl.roc <- structure.model2.final.fit %>%
    filter(ref.epitope == "GILGFVFTL_YVLDHLIVV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL_YVLDHLIVV\n",
        "ROC AUC: ", round(structure.model2.final.fit %>% 
            filter(ref.epitope == "GILGFVFTL_YVLDHLIVV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model2.glc.roc <- structure.model2.final.fit %>%
    filter(ref.epitope == "GLCTLVAML") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GLCTLVAML\n",
        "ROC AUC: ", round(structure.model2.final.fit %>% 
            filter(ref.epitope == "GLCTLVAML") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model2.nlv.roc <- structure.model2.final.fit %>%
    filter(ref.epitope == "NLVPMVATV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("NLVPMVATV\n",
        "ROC AUC: ", round(structure.model2.final.fit %>% 
            filter(ref.epitope == "NLVPMVATV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model2.rpi.roc <- structure.model2.final.fit %>%
    filter(ref.epitope == "RPIIRPATL") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("RPIIRPATL\n",
        "ROC AUC: ", round(structure.model2.final.fit %>% 
            filter(ref.epitope == "RPIIRPATL") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model2.ylq.roc <- structure.model2.final.fit %>%
    filter(ref.epitope == "YLQPRTFLL") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YLQPRTFLL\n",
        "ROC AUC: ", round(structure.model2.final.fit %>% 
            filter(ref.epitope == "YLQPRTFLL") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model2.yvl.roc <- structure.model2.final.fit %>%
    filter(ref.epitope == "YVLDHLIVV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YVLDHLIVV\n",
        "ROC AUC: ", round(structure.model2.final.fit %>% 
            filter(ref.epitope == "YVLDHLIVV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))

structure.model2.epitope.roc <- bind_rows(structure.model2.gil.roc, 
    structure.model2.gil.yvl.roc, 
    structure.model2.glc.roc, 
    # structure.model2.nlv.roc, 
    structure.model2.rpi.roc, 
    structure.model2.ylq.roc, 
    structure.model2.yvl.roc)  %>%
    getfbeta(beta)
label.data <- structure.model2.epitope.roc %>%
    group_by(epitope) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
structure.model2.epitope.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, col = epitope)) + 
    geom_path(lwd = 1.5, alpha = 0.8) +
    geom_text_repel(data = label.data, 
            aes(label = epitope), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9)

ggsave(paste0(out.path, "structure_model2_roc_per_epitope.png"), width = 10, height = 10)

# Sequence model
sequence.model2.gil.roc <- sequence.model2.final.fit %>%
    filter(ref.epitope == "GILGFVFTL") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL\n",
        "ROC AUC: ", round(sequence.model2.final.fit %>% 
            filter(ref.epitope == "GILGFVFTL") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model2.gil.yvl.roc <- sequence.model2.final.fit %>%
    filter(ref.epitope == "GILGFVFTL_YVLDHLIVV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL_YVLDHLIVV\n",
        "ROC AUC: ", round(sequence.model2.final.fit %>% 
            filter(ref.epitope == "GILGFVFTL_YVLDHLIVV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model2.glc.roc <- sequence.model2.final.fit %>%
    filter(ref.epitope == "GLCTLVAML") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GLCTLVAML\n",
        "ROC AUC: ", round(sequence.model2.final.fit %>% 
            filter(ref.epitope == "GLCTLVAML") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model2.nlv.roc <- sequence.model2.final.fit %>%
    filter(ref.epitope == "NLVPMVATV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("NLVPMVATV\n",
        "ROC AUC: ", round(sequence.model2.final.fit %>% 
            filter(ref.epitope == "NLVPMVATV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model2.rpi.roc <- sequence.model2.final.fit %>%
    filter(ref.epitope == "RPIIRPATL") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("RPIIRPATL\n",
        "ROC AUC: ", round(sequence.model2.final.fit %>% 
            filter(ref.epitope == "RPIIRPATL") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model2.ylq.roc <- sequence.model2.final.fit %>%
    filter(ref.epitope == "YLQPRTFLL") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YLQPRTFLL\n",
        "ROC AUC: ", round(sequence.model2.final.fit %>% 
            filter(ref.epitope == "YLQPRTFLL") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model2.yvl.roc <- sequence.model2.final.fit %>%
    filter(ref.epitope == "YVLDHLIVV") %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YVLDHLIVV\n",
        "ROC AUC: ", round(sequence.model2.final.fit %>% 
            filter(ref.epitope == "YVLDHLIVV") %>%
            roc_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))

sequence.model2.epitope.roc <- bind_rows(sequence.model2.gil.roc, 
    sequence.model2.gil.yvl.roc, 
    sequence.model2.glc.roc, 
    # sequence.model2.nlv.roc, 
    sequence.model2.rpi.roc, 
    sequence.model2.ylq.roc, 
    sequence.model2.yvl.roc)  %>%
    getfbeta(beta)
label.data <- sequence.model2.epitope.roc %>%
    group_by(epitope) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
sequence.model2.epitope.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, col = epitope)) + 
    geom_path(lwd = 1.5, alpha = 0.8) +
    geom_text_repel(data = label.data, 
            aes(label = epitope), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9)

ggsave(paste0(out.path, "sequence_model2_roc_per_epitope.png"), width = 10, height = 10)

# Per epitope, per model
model2.gil.roc <- bind_rows(combined.model2.gil.roc %>% mutate(model = "Combined"), structure.model2.gil.roc %>% mutate(model = "Structure"), sequence.model2.gil.roc %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model2.gil.roc %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model2.gil.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["GILGFVFTL"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["GILGFVFTL"]]) + 
    coord_equal() + 
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model2_GIL_roc.png"), width = 10, height = 10)

model2.gil.yvl.roc <- bind_rows(combined.model2.gil.yvl.roc %>% mutate(model = "Combined"), structure.model2.gil.yvl.roc %>% mutate(model = "Structure"), sequence.model2.gil.yvl.roc %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model2.gil.yvl.roc %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model2.gil.yvl.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["GILGFVFTL_YVLDHLIVV"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["GILGFVFTL_YVLDHLIVV"]]) +
    coord_equal() +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model2_GIL-YVL_roc.png"), width = 10, height = 10)

model2.glc.roc <- bind_rows(combined.model2.glc.roc %>% mutate(model = "Combined"), structure.model2.glc.roc %>% mutate(model = "Structure"), sequence.model2.glc.roc %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model2.glc.roc %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model2.glc.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, linetype = model)) +
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["GLCTLVAML"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["GLCTLVAML"]]) +
    coord_equal() +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model2_GLC_roc.png"), width = 10, height = 10)

model2.nlv.roc <- bind_rows(combined.model2.nlv.roc %>% mutate(model = "Combined"), structure.model2.nlv.roc %>% mutate(model = "Structure"), sequence.model2.nlv.roc %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model2.nlv.roc %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model2.nlv.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, linetype = model)) +
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["NLVPMVATV"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["NLVPMVATV"]]) +
    coord_equal() +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model2_NLV_roc.png"), width = 10, height = 10)

model2.rpi.roc <- bind_rows(combined.model2.rpi.roc %>% mutate(model = "Combined"), structure.model2.rpi.roc %>% mutate(model = "Structure"), sequence.model2.rpi.roc %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model2.rpi.roc %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model2.rpi.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["RPIIRPATL"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["RPIIRPATL"]]) +
    coord_equal() +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model2_RPI_roc.png"), width = 10, height = 10)

model2.ylq.roc <- bind_rows(combined.model2.ylq.roc %>% mutate(model = "Combined"), structure.model2.ylq.roc %>% mutate(model = "Structure"), sequence.model2.ylq.roc %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model2.ylq.roc %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model2.ylq.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["YLQPRTFLL"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["YLQPRTFLL"]]) +
    coord_equal() +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model2_YLQ_roc.png"), width = 10, height = 10)

model2.yvl.roc <- bind_rows(combined.model2.yvl.roc %>% mutate(model = "Combined"), structure.model2.yvl.roc %>% mutate(model = "Structure"), sequence.model2.yvl.roc %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model2.yvl.roc %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model2.yvl.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["YVLDHLIVV"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["YVLDHLIVV"]]) +
    coord_equal() +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model2_YVL_roc.png"), width = 10, height = 10)

# All models combined
model2.epitope.roc <- bind_rows(combined.model2.epitope.roc %>% mutate(model = "Combined"), structure.model2.epitope.roc %>% mutate(model = "Structure"), sequence.model2.epitope.roc %>% mutate(model = "Sequence"))
    getfbeta(beta) %>%
    mutate(epitope = gsub("\nROC AUC: [0-9.]+", "", epitope))
label.data <- model2.epitope.roc %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model2.epitope.roc %>%
    ggplot(aes(x = 1 - specificity, y = sensitivity, col = epitope, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.5) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

# Model 2 PR
# Per model, per epitope
combined.model2.gil.pr <- combined.model2.final.fit.expanded %>%
    filter(epitope == "GILGFVFTL") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL\n",
        "PR AUC: ", round(combined.model2.final.fit.expanded %>% 
            filter(epitope == "GILGFVFTL") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model2.gil.yvl.pr <- combined.model2.final.fit.expanded %>%
    filter(epitope == "GILGFVFTL_YVLDHLIVV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL_YVLDHLIVV\n",
        "PR AUC: ", round(combined.model2.final.fit.expanded %>% 
            filter(epitope == "GILGFVFTL_YVLDHLIVV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model2.glc.pr <- combined.model2.final.fit.expanded %>%
    filter(epitope == "GLCTLVAML") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GLCTLVAML\n",
        "PR AUC: ", round(combined.model2.final.fit.expanded %>% 
            filter(epitope == "GLCTLVAML") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model2.nlv.pr <- combined.model2.final.fit.expanded %>%
    filter(epitope == "NLVPMVATV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("NLVPMVATV\n",
        "PR AUC: ", round(combined.model2.final.fit.expanded %>% 
            filter(epitope == "NLVPMVATV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model2.rpi.pr <- combined.model2.final.fit.expanded %>%
    filter(epitope == "RPIIRPATL") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("RPIIRPATL\n",
        "PR AUC: ", round(combined.model2.final.fit.expanded %>% 
            filter(epitope == "RPIIRPATL") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model2.ylq.pr <- combined.model2.final.fit.expanded %>%
    filter(epitope == "YLQPRTFLL") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YLQPRTFLL\n",
        "PR AUC: ", round(combined.model2.final.fit.expanded %>% 
            filter(epitope == "YLQPRTFLL") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
combined.model2.yvl.pr <- combined.model2.final.fit.expanded %>%
    filter(epitope == "YVLDHLIVV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YVLDHLIVV\n",
        "PR AUC: ", round(combined.model2.final.fit.expanded %>% 
            filter(epitope == "YVLDHLIVV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
            
combined.model2.epitope.pr <- bind_rows(combined.model2.gil.pr, 
    combined.model2.gil.yvl.pr, 
    combined.model2.glc.pr, 
    # combined.model2.nlv.pr, 
    combined.model2.rpi.pr, 
    combined.model2.ylq.pr, 
    combined.model2.yvl.pr) %>%
    getfbeta(beta)
label.data <- combined.model2.epitope.pr %>%
    group_by(epitope) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
combined.model2.epitope.pr %>%
    ggplot(aes(x = recall, y = precision, col = epitope)) + 
    geom_path(lwd = 1.5, alpha = 0.8) +
    geom_text_repel(data = label.data, 
            aes(label = epitope), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9)

ggsave(paste0(out.path, "combined_model2_pr_per_epitope.png"), width = 10, height = 10)

structure.model2.gil.pr <- structure.model2.final.fit %>%
    filter(ref.epitope == "GILGFVFTL") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL\n",
        "PR AUC: ", round(structure.model2.final.fit %>% 
            filter(ref.epitope == "GILGFVFTL") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model2.gil.yvl.pr <- structure.model2.final.fit %>%
    filter(ref.epitope == "GILGFVFTL_YVLDHLIVV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL_YVLDHLIVV\n",
        "PR AUC: ", round(structure.model2.final.fit %>% 
            filter(ref.epitope == "GILGFVFTL_YVLDHLIVV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model2.glc.pr <- structure.model2.final.fit %>%
    filter(ref.epitope == "GLCTLVAML") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GLCTLVAML\n",
        "PR AUC: ", round(structure.model2.final.fit %>% 
            filter(ref.epitope == "GLCTLVAML") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model2.nlv.pr <- structure.model2.final.fit %>%
    filter(ref.epitope == "NLVPMVATV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("NLVPMVATV\n",
        "PR AUC: ", round(structure.model2.final.fit %>% 
            filter(ref.epitope == "NLVPMVATV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model2.rpi.pr <- structure.model2.final.fit %>%
    filter(ref.epitope == "RPIIRPATL") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("RPIIRPATL\n",
        "PR AUC: ", round(structure.model2.final.fit %>% 
            filter(ref.epitope == "RPIIRPATL") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model2.ylq.pr <- structure.model2.final.fit %>%
    filter(ref.epitope == "YLQPRTFLL") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YLQPRTFLL\n",
        "PR AUC: ", round(structure.model2.final.fit %>% 
            filter(ref.epitope == "YLQPRTFLL") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
structure.model2.yvl.pr <- structure.model2.final.fit %>%
    filter(ref.epitope == "YVLDHLIVV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YVLDHLIVV\n",
        "PR AUC: ", round(structure.model2.final.fit %>% 
            filter(ref.epitope == "YVLDHLIVV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))

structure.model2.epitope.pr <- bind_rows(structure.model2.gil.pr, 
    structure.model2.gil.yvl.pr, 
    structure.model2.glc.pr, 
    # structure.model2.nlv.pr, 
    structure.model2.rpi.pr, 
    structure.model2.ylq.pr, 
    structure.model2.yvl.pr) %>%
    getfbeta(beta)
label.data <- structure.model2.epitope.pr %>%
    group_by(epitope) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
structure.model2.epitope.pr %>%
    ggplot(aes(x = recall, y = precision, col = epitope)) + 
    geom_path(lwd = 1.5, alpha = 0.8) +
    geom_text_repel(data = label.data, 
            aes(label = epitope), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9)

ggsave(paste0(out.path, "structure_model2_pr_per_epitope.png"), width = 10, height = 10)

# Sequence model
sequence.model2.gil.pr <- sequence.model2.final.fit %>%
    filter(ref.epitope == "GILGFVFTL") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL\n",
        "PR AUC: ", round(sequence.model2.final.fit %>% 
            filter(ref.epitope == "GILGFVFTL") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model2.gil.yvl.pr <- sequence.model2.final.fit %>%
    filter(ref.epitope == "GILGFVFTL_YVLDHLIVV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GILGFVFTL_YVLDHLIVV\n",
        "PR AUC: ", round(sequence.model2.final.fit %>% 
            filter(ref.epitope == "GILGFVFTL_YVLDHLIVV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model2.glc.pr <- sequence.model2.final.fit %>%
    filter(ref.epitope == "GLCTLVAML") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("GLCTLVAML\n",
        "PR AUC: ", round(sequence.model2.final.fit %>% 
            filter(ref.epitope == "GLCTLVAML") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model2.nlv.pr <- sequence.model2.final.fit %>%
    filter(ref.epitope == "NLVPMVATV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("NLVPMVATV\n",
        "PR AUC: ", round(sequence.model2.final.fit %>% 
            filter(ref.epitope == "NLVPMVATV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model2.rpi.pr <- sequence.model2.final.fit %>%
    filter(ref.epitope == "RPIIRPATL") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("RPIIRPATL\n",
        "PR AUC: ", round(sequence.model2.final.fit %>% 
            filter(ref.epitope == "RPIIRPATL") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model2.ylq.pr <- sequence.model2.final.fit %>%
    filter(ref.epitope == "YLQPRTFLL") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YLQPRTFLL\n",
        "PR AUC: ", round(sequence.model2.final.fit %>% 
            filter(ref.epitope == "YLQPRTFLL") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))
sequence.model2.yvl.pr <- sequence.model2.final.fit %>%
    filter(ref.epitope == "YVLDHLIVV") %>%
    pr_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    mutate(epitope = paste0("YVLDHLIVV\n",
        "PR AUC: ", round(sequence.model2.final.fit %>% 
            filter(ref.epitope == "YVLDHLIVV") %>%
            pr_auc(truth = shared.specificity, .pred_Yes, event_level = "second") %>% pull(.estimate), 3)))

sequence.model2.epitope.pr <- bind_rows(sequence.model2.gil.pr, 
    sequence.model2.gil.yvl.pr, 
    sequence.model2.glc.pr, 
    # sequence.model2.nlv.pr, 
    sequence.model2.rpi.pr, 
    sequence.model2.ylq.pr, 
    sequence.model2.yvl.pr) %>%
    getfbeta(beta)
label.data <- sequence.model2.epitope.pr %>%
    group_by(epitope) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
sequence.model2.epitope.pr %>%
    ggplot(aes(x = recall, y = precision, col = epitope)) + 
    geom_path(lwd = 1.5, alpha = 0.8) +
    geom_text_repel(data = label.data, 
            aes(label = epitope), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9)

ggsave(paste0(out.path, "sequence_model2_pr_per_epitope.png"), width = 10, height = 10)

# Per epitope, per model
model2.gil.pr <- bind_rows(combined.model2.gil.pr %>% mutate(model = "Combined"), structure.model2.gil.pr %>% mutate(model = "Structure"), sequence.model2.gil.pr %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model2.gil.pr %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model2.gil.pr %>%
    ggplot(aes(x = recall, y = precision, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["GILGFVFTL"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["GILGFVFTL"]]) + 
    coord_equal() + 
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model2_GIL_pr.png"), width = 10, height = 10)

model2.gil.yvl.pr <- bind_rows(combined.model2.gil.yvl.pr %>% mutate(model = "Combined"), structure.model2.gil.yvl.pr %>% mutate(model = "Structure"), sequence.model2.gil.yvl.pr %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model2.gil.yvl.pr %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model2.gil.yvl.pr %>%
    ggplot(aes(x = recall, y = precision, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["GILGFVFTL_YVLDHLIVV"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["GILGFVFTL_YVLDHLIVV"]]) +
    coord_equal() +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model2_GIL-YVL_pr.png"), width = 10, height = 10)

model2.glc.pr <- bind_rows(combined.model2.glc.pr %>% mutate(model = "Combined"), structure.model2.glc.pr %>% mutate(model = "Structure"), sequence.model2.glc.pr %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model2.glc.pr %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model2.glc.pr %>%
    ggplot(aes(x = recall, y = precision, linetype = model)) +
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["GLCTLVAML"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["GLCTLVAML"]]) +
    coord_equal() +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model2_GLC_pr.png"), width = 10, height = 10)

model2.nlv.pr <- bind_rows(combined.model2.nlv.pr %>% mutate(model = "Combined"), structure.model2.nlv.pr %>% mutate(model = "Structure"), sequence.model2.nlv.pr %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model2.nlv.pr %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model2.nlv.pr %>%
    ggplot(aes(x = recall, y = precision, linetype = model)) +
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["NLVPMVATV"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["NLVPMVATV"]]) +
    coord_equal() +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model2_NLV_pr.png"), width = 10, height = 10)

model2.rpi.pr <- bind_rows(combined.model2.rpi.pr %>% mutate(model = "Combined"), structure.model2.rpi.pr %>% mutate(model = "Structure"), sequence.model2.rpi.pr %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model2.rpi.pr %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model2.rpi.pr %>%
    ggplot(aes(x = recall, y = precision, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["RPIIRPATL"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["RPIIRPATL"]]) +
    coord_equal() +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model2_RPI_pr.png"), width = 10, height = 10)

model2.ylq.pr <- bind_rows(combined.model2.ylq.pr %>% mutate(model = "Combined"), structure.model2.ylq.pr %>% mutate(model = "Structure"), sequence.model2.ylq.pr %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model2.ylq.pr %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model2.ylq.pr %>%
    ggplot(aes(x = recall, y = precision, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["YLQPRTFLL"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["YLQPRTFLL"]]) +
    coord_equal() +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model2_YLQ_pr.png"), width = 10, height = 10)

model2.yvl.pr <- bind_rows(combined.model2.yvl.pr %>% mutate(model = "Combined"), structure.model2.yvl.pr %>% mutate(model = "Structure"), sequence.model2.yvl.pr %>% mutate(model = "Sequence"))  %>%
    getfbeta(beta)
label.data <- model2.yvl.pr %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model2.yvl.pr %>%
    ggplot(aes(x = recall, y = precision, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.8, color = epitope.colors[["YVLDHLIVV"]]) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3,
            color = epitope.colors[["YVLDHLIVV"]]) +
    coord_equal() +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model2_YVL_pr.png"), width = 10, height = 10)


# All models combined
model2.epitope.pr <- bind_rows(combined.model2.epitope.pr %>% mutate(model = "Combined"), structure.model2.epitope.pr %>% mutate(model = "Structure"), sequence.model2.epitope.pr %>% mutate(model = "Sequence")) %>%
    getfbeta(beta) %>%
    mutate(epitope = gsub("\nPR AUC: [0-9.]+", "", epitope))
label.data <- model2.epitope.pr %>%
    group_by(epitope, model) %>%
    filter(fbeta == max(fbeta, na.rm = TRUE)) %>%
    slice(1)
model2.epitope.pr %>%
    ggplot(aes(x = recall, y = precision, col = epitope, linetype = model)) + 
    geom_path(lwd = 1.5, alpha = 0.5) +
    geom_text_repel(data = label.data, 
            aes(label = paste0(epitope, "\n", "(", model, ")")), 
            vjust = -3, 
            size = 3) + 
    coord_equal() + 
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_linetype_manual(values = c("Combined" = "solid", "Structure" = "dashed", "Sequence" = "dotted"))

ggsave(paste0(out.path, "model2_pr_per_epitope.png"), width = 10, height = 10)


# Combine metrics


# Combine coefficients
