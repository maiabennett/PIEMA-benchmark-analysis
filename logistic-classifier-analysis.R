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
# One column with the feature selection approach (Combined, Structure, or Sequence)
# One column with number of features (9, 7, or 5 for Combined, 7, 5 or 4 for Structure, 2 or 1 for Sequence)
# One column with data inclusion (NLV, Yes or No)
# One column with the metric value
all.metrics <- read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/with-NLV/model1_metrics.csv") %>% 
    mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 1", Approach = "Combined", 
        Features = 9, NLV = "Yes", Value = combined.estimate) %>%
    select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
    pivot_wider(names_from = Metrics, values_from = Value) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/with-NLV/model1_metrics.csv") %>% 
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 1", Approach = "Structure", 
            Features = 7, NLV = "Yes", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/with-NLV/model1_metrics.csv") %>% 
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 1", Approach = "Sequence", 
            Features = 2, NLV = "Yes", Value = sequence.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/without-NLV/model1_metrics.csv") %>% 
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 1", Approach = "Combined", 
            Features = 9, NLV = "No", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/without-NLV/model1_metrics.csv") %>% 
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 1", Approach = "Structure", 
            Features = 7, NLV = "No", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%  
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/without-NLV/model1_metrics.csv") %>% 
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 1", Approach = "Sequence", 
            Features = 2, NLV = "No", Value = sequence.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/with-NLV/model2_metrics.csv") %>% 
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 2", Approach = "Combined", 
            Features = 9, NLV = "Yes", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/with-NLV/model2_metrics.csv") %>% 
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 2", Approach = "Structure", 
            Features = 7, NLV = "Yes", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/with-NLV/model2_metrics.csv") %>% 
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 2", Approach = "Sequence", 
            Features = 2, NLV = "Yes", Value = sequence.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/without-NLV/model2_metrics.csv") %>% 
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 2", Approach = "Combined", 
            Features = 9, NLV = "No", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/without-NLV/model2_metrics.csv") %>% 
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 2", Approach = "Structure", 
            Features = 7, NLV = "No", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/without-NLV/model2_metrics.csv") %>% 
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 2", Approach = "Sequence", 
            Features = 2, NLV = "No", Value = sequence.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/7-feature-set/with-NLV/model1_metrics.csv") %>% 
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 1", Approach = "Combined", 
            Features = 7, NLV = "Yes", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/7-feature-set/with-NLV/model1_metrics.csv") %>% 
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 1", Approach = "Structure", 
            Features = 5, NLV = "Yes", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/7-feature-set/without-NLV/model1_metrics.csv") %>% 
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 1", Approach = "Combined", 
            Features = 7, NLV = "No", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/7-feature-set/without-NLV/model1_metrics.csv") %>% 
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 1", Approach = "Structure", 
            Features = 5, NLV = "No", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/7-feature-set/with-NLV/model2_metrics.csv") %>% 
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 2", Approach = "Combined", 
            Features = 7, NLV = "Yes", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/7-feature-set/with-NLV/model2_metrics.csv") %>% 
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 2", Approach = "Structure", 
            Features = 5, NLV = "Yes", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/7-feature-set/without-NLV/model2_metrics.csv") %>% 
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 2", Approach = "Combined", 
            Features = 7, NLV = "No", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/7-feature-set/without-NLV/model2_metrics.csv") %>% 
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 2", Approach = "Structure", 
            Features = 5, NLV = "No", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/with-NLV/model1_metrics.csv") %>% 
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 1", Approach = "Combined", 
            Features = 5, NLV = "Yes", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/with-NLV/model1_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 1", Approach = "Structure", 
            Features = 4, NLV = "Yes", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/with-NLV/model1_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 1", Approach = "Sequence", 
            Features = 1, NLV = "Yes", Value = sequence.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/without-NLV/model1_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 1", Approach = "Combined", 
            Features = 5, NLV = "No", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/without-NLV/model1_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 1", Approach = "Structure", 
            Features = 4, NLV = "No", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/without-NLV/model1_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 1", Approach = "Sequence", 
            Features = 1, NLV = "No", Value = sequence.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/with-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 2", Approach = "Combined", 
            Features = 5, NLV = "Yes", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/with-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 2", Approach = "Structure", 
            Features = 4, NLV = "Yes", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/with-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 2", Approach = "Sequence", 
            Features = 1, NLV = "Yes", Value = sequence.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/without-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 2", Approach = "Combined", 
            Features = 5, NLV = "No", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/without-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 2", Approach = "Structure", 
            Features = 4, NLV = "No", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/without-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "CDR3", Model = "Model 2", Approach = "Sequence", 
            Features = 1, NLV = "No", Value = sequence.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/with-NLV/model1_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 1", Approach = "Combined", 
            Features = 9, NLV = "Yes", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/with-NLV/model1_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 1", Approach = "Structure", 
            Features = 7, NLV = "Yes", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/with-NLV/model1_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 1", Approach = "Sequence", 
            Features = 2, NLV = "Yes", Value = sequence.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/without-NLV/model1_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 1", Approach = "Combined", 
            Features = 9, NLV = "No", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/without-NLV/model1_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 1", Approach = "Structure", 
            Features = 7, NLV = "No", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/without-NLV/model1_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 1", Approach = "Sequence", 
            Features = 2, NLV = "No", Value = sequence.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/with-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 2", Approach = "Combined", 
            Features = 9, NLV = "Yes", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/with-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 2", Approach = "Structure", 
            Features = 7, NLV = "Yes", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/with-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 2", Approach = "Sequence", 
            Features = 2, NLV = "Yes", Value = sequence.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/without-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 2", Approach = "Combined", 
            Features = 9, NLV = "No", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/without-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 2", Approach = "Structure", 
            Features = 7, NLV = "No", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/without-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 2", Approach = "Sequence", 
            Features = 2, NLV = "No", Value = sequence.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/7-feature-set/with-NLV/model1_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 1", Approach = "Combined", 
            Features = 7, NLV = "Yes", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/7-feature-set/with-NLV/model1_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 1", Approach = "Structure", 
            Features = 5, NLV = "Yes", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/7-feature-set/without-NLV/model1_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 1", Approach = "Combined", 
            Features = 7, NLV = "No", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/7-feature-set/without-NLV/model1_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 1", Approach = "Structure", 
            Features = 5, NLV = "No", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/7-feature-set/with-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 2", Approach = "Combined", 
            Features = 7, NLV = "Yes", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/7-feature-set/with-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 2", Approach = "Structure", 
            Features = 5, NLV = "Yes", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/7-feature-set/without-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 2", Approach = "Combined", 
            Features = 7, NLV = "No", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/7-feature-set/without-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 2", Approach = "Structure", 
            Features = 5, NLV = "No", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/with-NLV/model1_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 1", Approach = "Combined", 
            Features = 5, NLV = "Yes", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/with-NLV/model1_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 1", Approach = "Structure", 
            Features = 4, NLV = "Yes", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/with-NLV/model1_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 1", Approach = "Sequence", 
            Features = 1, NLV = "Yes", Value = sequence.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/without-NLV/model1_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 1", Approach = "Combined", 
            Features = 5, NLV = "No", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/without-NLV/model1_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 1", Approach = "Structure", 
            Features = 4, NLV = "No", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/without-NLV/model1_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 1", Approach = "Sequence", 
            Features = 1, NLV = "No", Value = sequence.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/with-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 2", Approach = "Combined", 
            Features = 5, NLV = "Yes", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/with-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 2", Approach = "Structure", 
            Features = 4, NLV = "Yes", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/with-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 2", Approach = "Sequence", 
            Features = 1, NLV = "Yes", Value = sequence.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/without-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 2", Approach = "Combined", 
            Features = 5, NLV = "No", Value = combined.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/without-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 2", Approach = "Structure", 
            Features = 4, NLV = "No", Value = structure.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value)) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/without-NLV/model2_metrics.csv") %>%
        mutate(Metrics = .metric, Dataset = "Full", Model = "Model 2", Approach = "Sequence", 
            Features = 1, NLV = "No", Value = sequence.estimate) %>%
        select(-.metric, -combined.estimate, -structure.estimate, -sequence.estimate) %>%
        pivot_wider(names_from = Metrics, values_from = Value))

        

# Combine per-epitope metrics- Need to figure out how to represent per-epitope metrics first
all.epitope.metrics

# Combine coefficients
all.coef.nine.features <- read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/with-NLV/model1_coefficients.csv") %>% 
    mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 1", Approach = "Combined", 
        NLV = "Yes", Value = combined.estimate, Significance = combined.p.value) %>%
    select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
        -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
    pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_") %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/with-NLV/model1_coefficients.csv") %>% 
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 1", Approach = "Structure", 
            NLV = "Yes", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/with-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 1", Approach = "Sequence", 
            NLV = "Yes", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/without-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 1", Approach = "Combined", 
            NLV = "No", Value = combined.estimate, Significance = combined.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/without-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 1", Approach = "Structure", 
            NLV = "No", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/without-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 1", Approach = "Sequence", 
            NLV = "No", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/with-NLV/model2_coefficients.csv") %>% 
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 2", Approach = "Combined", 
            NLV = "Yes", Value = combined.estimate, Significance = combined.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/with-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 2", Approach = "Structure", 
            NLV = "Yes", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/with-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 2", Approach = "Sequence", 
            NLV = "Yes", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/without-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 2", Approach = "Combined", 
            NLV = "No", Value = combined.estimate, Significance = combined.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/without-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 2", Approach = "Structure", 
            NLV = "No", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/9-feature-set/without-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 2", Approach = "Sequence", 
            NLV = "No", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/with-NLV/model1_coefficients.csv") %>% 
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 1", Approach = "Combined", 
            NLV = "Yes", Value = combined.estimate, Significance = combined.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/with-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 1", Approach = "Structure", 
            NLV = "Yes", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/with-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 1", Approach = "Sequence", 
            NLV = "Yes", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/without-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 1", Approach = "Combined", 
            NLV = "No", Value = combined.estimate, Significance = combined.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/without-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 1", Approach = "Structure", 
            NLV = "No", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/without-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 1", Approach = "Sequence", 
            NLV = "No", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/with-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 2", Approach = "Combined", 
            NLV = "Yes", Value = combined.estimate, Significance = combined.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/with-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 2", Approach = "Structure", 
            NLV = "Yes", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/with-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 2", Approach = "Sequence", 
            NLV = "Yes", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/without-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 2", Approach = "Combined", 
            NLV = "No", Value = combined.estimate, Significance = combined.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/without-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 2", Approach = "Structure", 
            NLV = "No", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/9-feature-set/without-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 2", Approach = "Sequence", 
            NLV = "No", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate, 
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_"))

all.coef.seven.features <- read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/7-feature-set/with-NLV/model1_coefficients.csv") %>%
    mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 1", Approach = "Combined", 
        NLV = "Yes", Value = combined.estimate, Significance = combined.p.value) %>%
    select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
        -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
    pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_") %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/7-feature-set/with-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 1", Approach = "Structure", 
            NLV = "Yes", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/7-feature-set/with-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 1", Approach = "Sequence", 
            NLV = "Yes", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/7-feature-set/without-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 1", Approach = "Combined", 
            NLV = "No", Value = combined.estimate, Significance = combined.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/7-feature-set/without-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 1", Approach = "Structure", 
            NLV = "No", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/7-feature-set/without-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 1", Approach = "Sequence", 
            NLV = "No", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/7-feature-set/with-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 2", Approach = "Combined", 
            NLV = "Yes", Value = combined.estimate, Significance = combined.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/7-feature-set/with-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 2", Approach = "Structure", 
            NLV = "Yes", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/7-feature-set/with-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 2", Approach = "Sequence", 
            NLV = "Yes", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/7-feature-set/without-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 2", Approach = "Combined", 
            NLV = "No", Value = combined.estimate, Significance = combined.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/7-feature-set/without-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 2", Approach = "Structure", 
            NLV = "No", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/7-feature-set/without-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 2", Approach = "Sequence", 
            NLV = "No", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/7-feature-set/with-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 1", Approach = "Combined", 
            NLV = "Yes", Value = combined.estimate, Significance = combined.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/7-feature-set/with-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 1", Approach = "Structure", 
            NLV = "Yes", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/7-feature-set/with-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 1", Approach = "Sequence", 
            NLV = "Yes", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/7-feature-set/without-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 1", Approach = "Combined", 
            NLV = "No", Value = combined.estimate, Significance = combined.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/7-feature-set/without-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 1", Approach = "Structure", 
            NLV = "No", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/7-feature-set/without-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 1", Approach = "Sequence", 
            NLV = "No", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/7-feature-set/with-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 2", Approach = "Combined", 
            NLV = "Yes", Value = combined.estimate, Significance = combined.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/7-feature-set/with-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 2", Approach = "Structure", 
            NLV = "Yes", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/7-feature-set/with-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 2", Approach = "Sequence", 
            NLV = "Yes", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/7-feature-set/without-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 2", Approach = "Combined", 
            NLV = "No", Value = combined.estimate, Significance = combined.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/7-feature-set/without-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 2", Approach = "Structure", 
            NLV = "No", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/7-feature-set/without-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 2", Approach = "Sequence", 
            NLV = "No", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_"))

all.coef.five.features <- read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/with-NLV/model1_coefficients.csv") %>%
    mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 1", Approach = "Combined", 
        NLV = "Yes", Value = combined.estimate, Significance = combined.p.value) %>%
    select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
        -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
    pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_") %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/with-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 1", Approach = "Structure", 
            NLV = "Yes", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/with-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 1", Approach = "Sequence", 
            NLV = "Yes", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/without-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 1", Approach = "Combined", 
            NLV = "No", Value = combined.estimate, Significance = combined.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/without-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 1", Approach = "Structure", 
            NLV = "No", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/without-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 1", Approach = "Sequence", 
            NLV = "No", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/with-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 2", Approach = "Combined", 
            NLV = "Yes", Value = combined.estimate, Significance = combined.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/with-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 2", Approach = "Structure", 
            NLV = "Yes", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/with-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 2", Approach = "Sequence", 
            NLV = "Yes", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/without-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 2", Approach = "Combined", 
            NLV = "No", Value = combined.estimate, Significance = combined.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/without-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 2", Approach = "Structure", 
            NLV = "No", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run1/5-feature-set/without-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "CDR3", Model = "Model 2", Approach = "Sequence", 
            NLV = "No", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/with-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 1", Approach = "Combined", 
            NLV = "Yes", Value = combined.estimate, Significance = combined.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/with-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 1", Approach = "Structure", 
            NLV = "Yes", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/with-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 1", Approach = "Sequence", 
            NLV = "Yes", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/without-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 1", Approach = "Combined", 
            NLV = "No", Value = combined.estimate, Significance = combined.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/without-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 1", Approach = "Structure", 
            NLV = "No", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/without-NLV/model1_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 1", Approach = "Sequence", 
            NLV = "No", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/with-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 2", Approach = "Combined", 
            NLV = "Yes", Value = combined.estimate, Significance = combined.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/with-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 2", Approach = "Structure", 
            NLV = "Yes", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/with-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 2", Approach = "Sequence", 
            NLV = "Yes", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/without-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 2", Approach = "Combined", 
            NLV = "No", Value = combined.estimate, Significance = combined.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/without-NLV/model2_coefficients.csv") %>%
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 2", Approach = "Structure", 
            NLV = "No", Value = structure.estimate, Significance = structure.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_")) %>%
    rbind(read.csv("./analysis/classifier/with-cross-reactives/with-avg-eucdist/run2/5-feature-set/without-NLV/model2_coefficients.csv") %>%  
        mutate(Coefficient = term, Dataset = "Full", Model = "Model 2", Approach = "Sequence", 
            NLV = "No", Value = sequence.estimate, Significance = sequence.p.value) %>%
        select(-term, -combined.estimate, -combined.p.value, -structure.estimate,
            -structure.p.value, -sequence.estimate, -sequence.p.value) %>%
        pivot_wider(names_from = Coefficient, values_from = c(Value, Significance), names_sep = "_"))

all.coef.combined.features <- all.coef.nine.features %>%
    filter(Approach == "Combined") %>%
    mutate(Features = 9) %>%
    bind_rows(all.coef.seven.features %>%
        filter(Approach == "Combined") %>%
        mutate(Features = 7)) %>%
    bind_rows(all.coef.five.features %>%
        filter(Approach == "Combined") %>%
        mutate(Features = 5))

all.coef.structure.features <- all.coef.nine.features %>%
    filter(Approach == "Structure") %>%
    mutate(Features = 7) %>%
    bind_rows(all.coef.seven.features %>%
        filter(Approach == "Structure") %>%
        mutate(Features = 5)) %>%
    bind_rows(all.coef.five.features %>%
        filter(Approach == "Structure")%>%
        mutate(Features = 4))

all.coef.sequence.features <- all.coef.nine.features %>%
    filter(Approach == "Sequence") %>%
    mutate(Features = 2) %>%
    bind_rows(all.coef.five.features %>%
        filter(Approach == "Sequence") %>%
        mutate(Features = 1))

all.coef <- all.coef.combined.features %>%
    bind_rows(all.coef.structure.features) %>%
    bind_rows(all.coef.sequence.features)

# Save results
write.csv(all.metrics, "./analysis/classifier/with-cross-reactives/with-avg-eucdist/all_metrics.csv", row.names = FALSE)
write.csv(all.coef.combined.features, "./analysis/classifier/with-cross-reactives/with-avg-eucdist/all_coefs_combined_features.csv", row.names = FALSE)
write.csv(all.coef.structure.features, "./analysis/classifier/with-cross-reactives/with-avg-eucdist/all_coefs_structure_features.csv", row.names = FALSE)
write.csv(all.coef.sequence.features, "./analysis/classifier/with-cross-reactives/with-avg-eucdist/all_coefs_sequence_features.csv", row.names = FALSE)


# Plot results to identify best approaches, feature influence, etc
# Plotting accuracy (y axis), with bars colored by the values of approach + number of features, with a x axis being a combination of dataset + inclusion of NLV
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            Approach == "Combined" & Features == 5 ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset and approach", subtitle = "Model 1: Positive vs. Decoy pairs", y = "Accuracy") +
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/accuracy_by_dataset_and_approach_model1.png", width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            Approach == "Combined" & Features == 5 ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset and approach", subtitle = "Model 2: Positive vs. Unlikely and Decoy pairs", y = "Accuracy") +
    scale_fill_viridis_d(option = "viridis", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/accuracy_by_dataset_and_approach_model2.png", width = 20, height = 15, dpi = 300)

# Plotting ROC_AUC
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            Approach == "Combined" & Features == 5 ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC_AUC by dataset and approach", subtitle = "Model 1: Positive vs. Decoy pairs", y = "ROC_AUC") +
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/roc_auc_by_dataset_and_approach_model1.png", width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            Approach == "Combined" & Features == 5 ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC_AUC by dataset and approach", subtitle = "Model 2: Positive vs. Unlikely and Decoy pairs", y = "ROC_AUC") +
    scale_fill_viridis_d(option = "viridis", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/roc_auc_by_dataset_and_approach_model2.png", width = 20, height = 15, dpi = 300)

# Plotting PR_AUC
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            Approach == "Combined" & Features == 5 ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR_AUC by dataset and approach", subtitle = "Model 1: Positive vs. Decoy pairs", y = "PR_AUC") +
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.15,0.6),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/pr_auc_by_dataset_and_approach.png", width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            Approach == "Combined" & Features == 5 ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR_AUC by dataset and approach", subtitle = "Model 2: Positive vs. Unlikely and Decoy pairs", y = "PR_AUC") +
    scale_fill_viridis_d(option = "viridis", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.15,0.6),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/pr_auc_by_dataset_and_approach_model2.png", width = 20, height = 15, dpi = 300)

# Plotting recall
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            Approach == "Combined" & Features == 5 ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = recall, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Recall by dataset and approach", subtitle = "Model 1: Positive vs. Decoy pairs", y = "Recall") +
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.35,0.6),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/recall_by_dataset_and_approach_model1.png", width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            Approach == "Combined" & Features == 5 ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = recall, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Recall by dataset and approach", subtitle = "Model 2: Positive vs. Unlikely and Decoy pairs", y = "Recall") +
    scale_fill_viridis_d(option = "viridis", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.35,0.6),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/recall_by_dataset_and_approach_model2.png", width = 20, height = 15, dpi = 300)

# Precision
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            Approach == "Combined" & Features == 5 ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = precision, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Precision by dataset and approach", subtitle = "Model 1: Positive vs. Decoy pairs", y = "Precision") +
    scale_fill_viridis_d(option = "plasma", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.15,0.5),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/precision_by_dataset_and_approach_model1.png", width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            Approach == "Combined" & Features == 5 ~ paste(Approach, "approach (CDR3 only)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = precision, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Precision by dataset and approach", subtitle = "Model 2: Positive vs. Unlikely and Decoy pairs", y = "Precision") +
    scale_fill_viridis_d(option = "viridis", begin = 0.1, end = .9) +
    scale_y_continuous(limits=c(0.15,0.5),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/precision_by_dataset_and_approach_model2.png", width = 20, height = 15, dpi = 300)

# Plotting the same metrics based on the inclusion of full sequence similarity in the model
# Includes all structural methods (for reference), the relevant sequence method, and the relevant combined method(s)
# Set colors for consistency in plots excluding certain methods
model.approaches <- unique(all.metrics %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
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
            (Approach == "Sequence" & Features == 2) | 
            (Approach == "Combined" & (Features == 7 | Features == 9))) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset for approaches including full sequence similarity", subtitle = "Model 1: Positive vs Decoy pairs", y = "Accuracy") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/accuracy_by_dataset_full_sequence_model1.png", width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Features == 2) | 
            (Approach == "Combined" & (Features == 7 | Features == 9))) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset for approaches including full sequence similarity", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "Accuracy") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/accuracy_by_dataset_full_sequence_model2.png", width = 20, height = 15, dpi = 300)

# Without full sequence
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Features == 1) | 
            (Approach == "Combined" & Features == 5)) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset for approaches excluding full sequence similarity", subtitle = "Model 1: Positive vs Decoy pairs", y = "Accuracy") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/accuracy_by_dataset_cdr3_sequence_model1.png", width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Features == 1) | 
            (Approach == "Combined" & Features == 5)) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = accuracy, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Accuracy by dataset for approaches excluding full sequence similarity", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "Accuracy") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/accuracy_by_dataset_cdr3_sequence_model2.png", width = 20, height = 15, dpi = 300)

# ROC_AUC
# With full sequence
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Features == 2) | 
            (Approach == "Combined" & (Features == 7 | Features == 9))) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC_AUC by dataset for approaches including full sequence similarity", subtitle = "Model 1: Positive vs Decoy pairs", y = "ROC_AUC") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/roc_auc_by_dataset_full_sequence_model1.png", width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Features == 2) | 
            (Approach == "Combined" & (Features == 7 | Features == 9))) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC_AUC by dataset for approaches including full sequence similarity", subtitle = "Model 1: Positive vs Decoy pairs", y = "ROC_AUC") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/roc_auc_by_dataset_full_sequence_model2.png", width = 20, height = 15, dpi = 300)

# Without full sequence
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Features == 1) | 
            (Approach == "Combined" & Features == 5)) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC_AUC by dataset for approaches excluding full sequence similarity", subtitle = "Model 1: Positive vs Decoy pairs", y = "ROC_AUC") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/roc_auc_by_dataset_cdr3_sequence_model1.png", width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Features == 1) | 
            (Approach == "Combined" & Features == 5)) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = roc_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "ROC_AUC by dataset for approaches excluding full sequence similarity", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "ROC_AUC") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.45,0.65),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/roc_auc_by_dataset_cdr3_sequence_model2.png", width = 20, height = 15, dpi = 300)

# PR_AUC
# With full sequence
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Features == 2) | 
            (Approach == "Combined" & (Features == 7 | Features == 9))) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR_AUC by dataset for approaches including full sequence similarity", subtitle = "Model 1: Positive vs Decoy pairs", y = "PR_AUC") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.15,0.6),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/pr_auc_by_dataset_full_sequence_model1.png", width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Features == 2) | 
            (Approach == "Combined" & (Features == 7 | Features == 9))) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR_AUC by dataset for approaches including full sequence similarity", subtitle = "Model 1: Positive vs Decoy pairs", y = "PR_AUC") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.15,0.6),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/pr_auc_by_dataset_full_sequence_model2.png", width = 20, height = 15, dpi = 300)

# Without full sequence
# Model 1
ggplot(
    all.metrics %>% 
        filter(Model == "Model 1") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Features == 1) | 
            (Approach == "Combined" & Features == 5)) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR_AUC by dataset for approaches excluding full sequence similarity", subtitle = "Model 1: Positive vs Decoy pairs", y = "PR_AUC") +
    scale_fill_manual(values = model1.approach.colors) +
    scale_y_continuous(limits=c(0.15,0.6),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/pr_auc_by_dataset_cdr3_sequence_model1.png", width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.metrics %>% 
        filter(Model == "Model 2") %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        filter(Approach == "Structure" | 
            (Approach == "Sequence" & Features == 1) | 
            (Approach == "Combined" & Features == 5)) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Approach, y = pr_auc, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, strip.position = "bottom", scales = "free_x", nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "PR_AUC by dataset for approaches excluding full sequence similarity", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "PR_AUC") +
    scale_fill_manual(values = model2.approach.colors) +
    scale_y_continuous(limits=c(0.15,0.6),oob = rescale_none)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/pr_auc_by_dataset_cdr3_sequence_model2.png", width = 20, height = 15, dpi = 300)

# Feature coefficients
# Convert to long format for plotting, and plot features on X, Coefficients on Y, colored by Approach + Features
# Combined feature sets
ggplot(
    all.coef %>%
        filter(Approach == "Combined") %>%
        filter(Model == "Model 1") %>%
        pivot_longer(cols = starts_with("Value_"), names_to = "Feature", values_to = "Coefficient") %>%
        #rename_with(~ gsub("^Value_", "", .), starts_with("Value_")) %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Feature, y = Coefficient, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Feature coefficients by dataset and approach", subtitle = "Model 1: Positive vs Decoy pairs", y = "Coefficient") +
    scale_fill_manual(values = model1.approach.colors)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/combined_feature_coefficients_model1.png", width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.coef %>%
        filter(Approach == "Combined") %>%
        filter(Model == "Model 2") %>%
        pivot_longer(cols = starts_with("Value_"), names_to = "Feature", values_to = "Coefficient") %>%
        #rename_with(~ gsub("^Value_", "", .), starts_with("Value_")) %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Feature, y = Coefficient, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Feature coefficients by dataset and approach", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "Coefficient") +
    scale_fill_manual(values = model2.approach.colors)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/combined_feature_coefficients_model2.png", width = 20, height = 15, dpi = 300)

# Sequence feature sets
# Model 1
ggplot(
    all.coef %>%
        filter(Approach == "Sequence") %>%
        filter(Model == "Model 1") %>%
        pivot_longer(cols = c("Value_CDR3.similarity", "Value_full.similarity"), names_to = "Feature", values_to = "Coefficient") %>%
        #rename_with(~ gsub("^Value_", "", .), starts_with("Value_")) %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Feature, y = Coefficient, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Feature coefficients by dataset and approach", subtitle = "Model 1: Positive vs Decoy pairs", y = "Coefficient") +
    scale_fill_manual(values = model1.approach.colors)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/sequence_feature_coefficients_model1.png", width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.coef %>%
        filter(Approach == "Sequence") %>%
        filter(Model == "Model 2") %>%
        pivot_longer(cols = c("Value_CDR3.similarity", "Value_full.similarity"), names_to = "Feature", values_to = "Coefficient") %>%
        #rename_with(~ gsub("^Value_", "", .), starts_with("Value_")) %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Feature, y = Coefficient, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data, nrow = 1) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Feature coefficients by dataset and approach", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "Coefficient") +
    scale_fill_manual(values = model2.approach.colors)

# Structure feature sets
# Model 1
ggplot(
    all.coef %>%
        filter(Approach == "Structure") %>%
        filter(Model == "Model 1") %>%
        pivot_longer(cols = starts_with("Value_"), names_to = "Feature", values_to = "Coefficient") %>%
        #rename_with(~ gsub("^Value_", "", .), starts_with("Value_")) %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Feature, y = Coefficient, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Feature coefficients by dataset and approach", subtitle = "Model 1: Positive vs Decoy pairs", y = "Coefficient") +
    scale_fill_manual(values = model1.approach.colors)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/structure_feature_coefficients_model1.png", width = 20, height = 15, dpi = 300)

# Model 2
ggplot(
    all.coef %>%
        filter(Approach == "Structure") %>%
        filter(Model == "Model 2") %>%
        pivot_longer(cols = starts_with("Value_"), names_to = "Feature", values_to = "Coefficient") %>%
        #rename_with(~ gsub("^Value_", "", .), starts_with("Value_")) %>%
        mutate(Dataset = ifelse(Dataset == "CDR3", "<90% CDR3 Similarity\n", "<90% Full Seq. Similarity\n")) %>%
        mutate(NLV = ifelse(NLV == "Yes", "With NLV", "Without NLV")) %>%
        #mutate(`LogReg Approach` = paste(Approach, "approach, ", Features, " features")) %>%
        mutate(`LogReg Approach` = case_when(
            Approach == "Sequence" & Features == 1 ~ paste(Approach, "approach (CDR3 only)"),
            Approach == "Sequence" & Features == 2 ~ paste(Approach, "approach (CDR3 + full)"),
            .default = paste(Approach, "approach,", Features, "features"))) %>%
        mutate(Data = paste(Dataset, NLV)),
    aes(x = Feature, y = Coefficient, fill = `LogReg Approach`)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap( ~ Data) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Feature coefficients by dataset and approach", subtitle = "Model 2: Positive vs Unlikely and Decoy pairs", y = "Coefficient") +
    scale_fill_manual(values = model2.approach.colors)

ggsave("./analysis/classifier/with-cross-reactives/with-avg-eucdist/structure_feature_coefficients_model2.png", width = 20, height = 15, dpi = 300)

# Per-epitope performance metrics- for top 3 performing epitopes (GIL, GLC, YVL or GIL-YVL crossreactives)


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

