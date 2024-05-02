# Blood group prediction models

PLOS Computational Biology. 2024 Mar 21;20(3):e1011977. doi: 10.1371/journal.pcbi.1011977. PMID: 38512997

Published: March 21, 2024

[**A machine-learning method for biobank-scale genetic prediction of blood group antigens**](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011977)

K Hyvärinen, K Haimila, C Moslemi, Blood Service Biobank, ML Olsson, SR Ostrowski, OB Pedersen, C Erikstrup, J Partanen, J Ritari

## Prerequisites

1.  PLINK 1.90 <https://www.cog-genomics.org/plink/>

2.  bcftools <https://samtools.github.io/bcftools/bcftools.html>

3.  R packets: tidyverse, caret, ranger, data.table, e1071, ROCR, ggplotify, and ggbeeswarm

## 1. [Creating prediction models](https://github.com/FRCBS/Blood_group_prediction/tree/main/Creating_prediction_models)

Scripts for your own genotype and blood group data preprocessing and random forest model fitting.

Creates RF models, confusion matrix plots, boxplots with posterior probabilities, statistics tables, and tables of important variables.

## 2. [Applying the Finnish models](https://github.com/FRCBS/Blood_group_prediction/tree/main/Applying_Finnish_bloog_group_prediction_models)

The Finnish models and general blood group information are in `Finnish_Blood_group_prediction_models_genes.rds`

The script for blood group prediction: `001_Applying_prediction_models.R`
