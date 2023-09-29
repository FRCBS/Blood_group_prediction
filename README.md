# Blood_group_prediction

[**A machine-learning method for biobank-scale genetic prediction of blood group antigens**](https://www.medrxiv.org/content/10.1101/2023.09.18.23295700v2)**\
**K Hyvärinen, K Haimila, C Moslemi, Blood Service Biobank, ML Olsson, SR Ostrowski, OB Pedersen, C Erikstrup, J Partanen, J Ritari

## Prerequisites

1.  PLINK 1.90 <https://www.cog-genomics.org/plink/>

2.  bcftools <https://samtools.github.io/bcftools/bcftools.html>

3.  R packets:

    tidyverse, caret, ranger, data.table, e1071, and ggbeeswarm

## 1. Creating prediction models (./Creating_prediction_models/)

Scripts for your own genotype and blood group data preprocessing and random forest model fitting.

Creates RF models, confusion matrix plots, boxplots with posterior probabilities, statistics tables, and tables of important variables.

## 2. Applying the Finnish models (./Applying_Finnish_blood_group_prediction_models/)

The Finnish models and general blood group information are in `Finnish_Blood_group_prediction_models_genes.rds`

The script for blood group prediction: `001_Applying_prediction_models.R`
