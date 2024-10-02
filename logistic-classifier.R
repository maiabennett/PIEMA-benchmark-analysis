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

out.path <- "./analysis/classifier/"

# This script uses tidymodels to train a logistic regression classifier on the PIEMA data
# For this purpose, two models will be trained: 
# 1. A model that uses only "True binding pairs" (positive data) and "Decoy binding pairs" (negative data)
# 2. A model that uses all data, including "True binding pairs" (positive data) and "Decoy binding pairs" and "Unlikely binding pairs" (negative data)

# Load the data
piema.data <- read.csv("./analysis/apbs/final_receptor_data.csv") %>%
    dplyr::rename(top.kas = top.kdist, top.spearman = top.corr, top.euc.dist = top.distance, 
        top.sg.euc.dist = top.patDist, top.cosine.sim = top.shapSim, top.sg.median.kas = top.med_scr, 
        binding.pair.type = type, avg.euc.dist = avg.distance, avg.kas = avg.kdist, avg.spearman = avg.corr,
        kernel.count = count, pair.id = id, top.sg.id = top.sg_group)

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
    select(shared.specificity, top.kas, avg.euc.dist, avg.kas, avg.spearman, top.cosine.sim, top.sg.median.kas, kernel.count, CDR3.similarity, full.similarity) %>%
    pivot_longer(cols = c(top.kas, avg.euc.dist, avg.kas, avg.spearman, top.cosine.sim, top.sg.median.kas, kernel.count, CDR3.similarity, full.similarity), names_to = "feature", values_to = "value")

ggplot(model1.long, aes(x = shared.specificity, y = value, fill = feature)) +
    geom_boxplot() +
    facet_wrap(~feature, scales = "free", ncol = 4) +
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = 0.9) +
    theme_minimal()

ggsave(paste0(out.path, "model1_feature_boxplot.png"), width = 10, height = 15)

model2.long <- model2.data %>%
    select(shared.specificity, top.kas, avg.euc.dist, avg.kas, avg.spearman, top.cosine.sim, top.sg.median.kas, kernel.count, CDR3.similarity, full.similarity) %>%
    pivot_longer(cols = c(top.kas, avg.euc.dist, avg.kas, avg.spearman, top.cosine.sim, top.sg.median.kas, kernel.count, CDR3.similarity, full.similarity), names_to = "feature", values_to = "value")

ggplot(model2.long, aes(x = shared.specificity, y = value, fill = feature)) +
    geom_boxplot() +
    facet_wrap(~feature, scales = "free", ncol = 4) +
    scale_color_viridis_d(option = "plasma", begin = 0.1, end = 0.9) +
    theme_minimal()

ggsave(paste0(out.path, "model2_feature_boxplot.png"), width = 10, height = 15)


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
# model.form <- shared.specificity ~ top.kas + avg.euc.dist + avg.kas + avg.spearman + top.cosine.sim + kernel.count + CDR3.similarity + full.similarity
## START HERE
model.form <- shared.specificity ~ avg.euc.dist + avg.kas + avg.spearman + top.cosine.sim + kernel.count + CDR3.similarity + full.similarity
positive.class <- "Yes"
negative.class <- "No"

# Define F beta metric
f_meas <- function(data, truth, estimate, na_rm = TRUE, ...) {
    yardstick::f_meas(
        data = data,
        truth = !!rlang::enquo(truth),
        estimate = !!rlang::enquo(estimate),
        beta = 0.5,
        na_rm = na_rm,
        ...
    )
}
f_meas <- new_class_metric(f_meas, direction = "maximize")

classification.metrics <- metric_set(yardstick::accuracy, 
                                     yardstick::kap, 
                                     yardstick::roc_auc, 
                                     yardstick::mn_log_loss, 
                                     yardstick::sens, 
                                     yardstick::spec, 
                                     #yardstick::f_meas,
                                     yardstick::precision,
                                     yardstick::recall,
                                     yardstick::bal_accuracy,
                                     f_meas) 

# Define the recipe for the logistic regression models
model1.recipe <- recipe(model.form, data = model1.train) %>%
    step_dummy(all_nominal(), -all_outcomes()) %>%
    step_normalize(all_predictors())

model2.recipe <- recipe(model.form, data = model2.train) %>%
    step_dummy(all_nominal(), -all_outcomes()) %>%
    step_normalize(all_predictors())

# Recipes with upsampling
model1.upsampling.recipe <- model1.recipe %>%
    step_smotenc(target.var)

model2.upsampling.recipe <- model2.recipe %>%
    step_smotenc(target.var)

# Recipes with downsampling
model1.downsampling.recipe <- model1.recipe %>%
    step_downsample(target.var)

model2.downsampling.recipe <- model2.recipe %>%
    step_downsample(target.var)

# Define recipe set
model1.recipes <- list(basic = model1.recipe, upsampled = model1.upsampling.recipe, downsampled = model1.downsampling.recipe)
model2.recipes <- list(basic = model2.recipe, upsampled = model2.upsampling.recipe, downsampled = model2.downsampling.recipe)

# Define the logistic regression model for both models
logreg.model <- 
  logistic_reg() %>%
  set_engine("glm") %>%
  set_mode("classification") 


# Define multi-recipe workflow sets for both models
model1.workflow <- workflow_set(model1.recipes,
    models = list(logreg = logreg.model),
    cross = FALSE)
model2.workflow <- workflow_set(model2.recipes,
    models = list(logreg = logreg.model),
    cross = FALSE)

# Specify resamples and CV folds
control <- control_resamples(save_workflow = TRUE,
                                 save_pred = TRUE,
                                 event_level = "second")

model1.cv.folds <- vfold_cv(model1.train, v = 5)
model2.cv.folds <- vfold_cv(model2.train, v = 5)


# Fit the model sets
model1.fit <- model1.workflow %>%
    workflow_map("fit_resamples", 
        seed = 77482951,
        verbose = TRUE,
        resamples = model1.cv.folds,
        metrics = classification.metrics,
        control = control)

model2.fit <- model2.workflow %>%
    workflow_map("fit_resamples",
        seed = 77482951,
        verbose = TRUE,
        resamples = model2.cv.folds,
        metrics = classification.metrics,
        control = control)


# Store workflow set predictions
model1.train.predictions <- model1.fit %>%
    collect_predictions()
model2.train.predictions <- model2.fit %>%
    collect_predictions()

# Calculate workflow set metrics
model1.train.metrics <- model1.fit %>%
    collect_metrics() 
model2.train.metrics <- model2.fit %>%
    collect_metrics()

# Fetch the best model by F beta metric 
metric.for.selection <- "f_meas"

model1.best.fit <- model1.fit %>%
    fit_best(metric = metric.for.selection)
model2.best.fit <- model2.fit %>%
    fit_best(metric = metric.for.selection)


# Fit the model sets to the test data and store test predictions
model1.final.fit <- model1.best.fit %>%
    augment(new_data = model1.test)
model2.final.fit <- model2.best.fit %>%
    augment(new_data = model2.test)


# Calculate test metrics
model1.test.metrics <- classification.metrics(data = model1.final.fit,
    truth = shared.specificity,
    estimate = .pred_class,
    .pred_Yes,
    event_level = 'second')
                                         
model2.test.metrics <- classification.metrics(data = model2.final.fit,
    truth = shared.specificity,
    estimate = .pred_class,
    .pred_Yes,
    event_level = 'second')

# Plot ROC 
# Look into storing models based on various metrics and plotting them together
(model1.roc <- model1.final.fit %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    autoplot())

(model2.roc <- model2.final.fit %>%
    roc_curve(truth = shared.specificity, .pred_Yes, event_level = "second") %>%
    autoplot())

# Extract and store final model formulas, other information
model1.best.fit$fit[2]
model1.coefficients <- model1.best.fit %>%
    extract_fit_engine() %>%
    broom::tidy()

model2.best.fit$fit[2]
model2.coefficients <- model2.best.fit %>%
    extract_fit_engine() %>%
    broom::tidy()

# Looks like I should remove top KAS