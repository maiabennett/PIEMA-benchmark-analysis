# PIEMA benchmark analysis

This repository contains the code and data for an analysis aimed at providing a benchmark of the PIEMA (Protein Interface Energy Map Alignment) tool for T cell receptor (TCR) structural similarity identification, found [here](https://github.com/maiabennett/PIEMA). The analysis also extends PIEMA through the implementation of multiple Logistic Regression classifiers incorporating sequence and structural features to make predictions on shared specificity between TCRs.

## Overview of the Analysis
### Data Selection
The analysis begins with the selection of high-confidence TCR binding data specific to 13 epitopes, each with more than 25 receptors available. This step is performed in the [process-inputs.R](process-inputs.R) script, which ensures that the dataset used for the analysis is robust and reliable.

![data-selection-overview](https://github.com/user-attachments/assets/8362bbb1-bd3a-4905-b3c4-ff516563c5b2)


### Data Analysis with PIEMA
The selected data is then analyzed using the PIEMA tool to generate structural similarity metrics. PIEMA's ability to identify structural similarities between TCRs is a key component of this analysis.

![PIEMA-overview](https://github.com/user-attachments/assets/5223bba0-8ee7-4ce9-a44e-dd0ff83d1272)


### Result Processing and Analysis
The results generated by PIEMA are processed and analyzed through various plots that assess the structural metrics and sequence similarity. This step is carried out in the process-results.R and plot-results.R scripts. The plots provide a visual representation of the data, making it easier to interpret the results and identify patterns.

### Logistic Classifier Implementation
To extend the capabilities of PIEMA, multiple Logistic Regression classifiers are implemented using the `tidymodels` R package suite. These classifiers incorporate both sequence and structural features to make predictions on shared specificity between TCRs. The implementation is detailed in the [logistic-classifier-structure.Rmd](./logistic-classifier-structure.Rmd), [logistic-classifier-combined.Rmd](./logistic-classifier-combined.Rmd), and [logistic-classifier-structure.Rmd](./logistic-classifier-sequence.Rmd) files. These documents describe the process of training the classifiers and generating preliminary analysis plots such as feature box plots, ROC curves, and PR curves.



### Logistic Classifier Efficacy Assessment
The efficacy of the Logistic Regression classifiers is assessed in the [logistic-classifier-analysis.R](./logistic-classifier-analysis.R) script. This script facilitates a more in-depth analysis of optimal model data curation and feature selection approaches. The assessment includes evaluating the performance metrics and coefficients of the classifiers to determine their effectiveness in predicting shared specificity.

### Case Study
A case study is conducted to examine receptor pairs where PIEMA works and doesn't work. This includes analyzing dissimilar sequences and similar sequences with high-scoring and low-scoring structural outputs from PIEMA. The case study provides insights into the strengths and limitations of PIEMA and the Logistic Regression classifiers and is found in the [case-study-analysis.R](./case-study-analysis.R) script.
