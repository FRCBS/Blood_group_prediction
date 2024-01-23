###############################################################################
### Blood group imputation using random forest models
###############################################################################

library(tidyverse)
library(caret)
library(ranger)
library(data.table)
library(e1071)

############################################################################### 
### Import PLINK genotype dosage files

# Dosage files and BG_genes tibble are generated using 
# BG_genotype_preparations_replication.R script
  
dosage.files <-  list.files("data/BG_genes",      # List dosage files                     
                            pattern = c("*.raw"), 
                            full.names = TRUE) 
dosage.data <- lapply(dosage.files, read_table2)  # Read dosage files

# Name elements in the dosage data list
dosage.files_split <-  gsub("data/BG_genes/", "", # Create splitted list
                            dosage.files) %>% 
  str_split_fixed(., "_", 2) %>% .[,2] %>% 
  gsub("_2K.vcf.dosage.raw", '', ., fixed = T)
names(dosage.data) <- dosage.files_split     # Name elements with gene names
dosage.data %>% names()                           # Check the names

###############################################################################
### Merge genotype and phenotype data 

### Phenotype data from Phenotype_preparations_replication.R

str(BG_phenos)
# Structure specification example
  # tibble [1,192 × 54] (S3: tbl_df/tbl/data.frame), columns:
  # 1. $ releasedIdentifier: chr [1:1192] "Sample1" "Sample2" "Sample3" "Sample4" ...
  # 2. $ ABO               : chr [1:1192] "A" "B" "AB" "O" ...
  # 3. $ D                 : chr [1:1192] "D+" "D-" "D+" "D+" ...
  # 4. $ Kpa               : chr [1:1192] "Kpa-" "Kpa-" "Kpa-" "Kpa-" ...
  # 5. $ Kpb               : chr [1:1192] "Kpb+" "Kpb+" "Kpb+" "Kpb+" ...
  # 6. ...  n column, include all the blood group information into columns    

str(Pheno_info)
# Structure specification example
  # tibble [53 × 7] (S3: tbl_df/tbl/data.frame)
  # $ Phenotype    : chr [1:53] "ABO" "A1" "A2" "Yta" ...
  # $ BG_system    : chr [1:53] "ABO" "ABO" "ABO" "CartwrightYt" ...
  # $ Pheno_genes  : chr [1:53] "ABO" "ABO" "ABO" "ACHE" ...
  # $ Gene1        : chr [1:53] "ABO" "ABO" "ABO" "ACHE" ...
  # $ Gene2        : chr [1:53] NA NA NA NA ...
  # $ Gene3        : chr [1:53] NA NA NA NA ...
  # $ Plink_dataset: chr [1:53] "data/BG_genes/ABO_ABO_2K.vcf" "data/BG_genes/ABO_ABO_2K.vcf" ...

  # Column Pheno_genes contain all genes separated by _ e.g. GYPA_GYPB_GYPE, RHCE_RHD
  # Column Plink_dataset contain PLINK dataset prefixes generated using 001_BG_genotype_preparations_replication.R script

### Genotype data from 001_Genotype_preparation.R

# Merge data
Geno_pheno.data <- map(1:nrow(Pheno_info), function(i) {
  
  # Blood group phenotype 
  BG_PHENO <- Pheno_info[i, "Phenotype"] %>% unlist
  # Genes for the BG; char vector
  BG_GENES <- Pheno_info[i, "Pheno_genes"] %>% unlist
  # Extract dosages by gene names; list --> Reduce --> data frame
  BG_DOSAGE <- Reduce(function(x,y) left_join(x,y, by='IID'), 
                      dosage.data[BG_GENES]) %>%
    dplyr::select(., -matches('FID|PAT|MAT|SEX|PHENOTYPE'))
  # Extract phenotype for each sample ID and join to dosages
  OUT <- dplyr::select(BG_phenos, releasedIdentifier, all_of(BG_PHENO)) %>% 
    left_join(., BG_DOSAGE, by = c("releasedIdentifier" = "IID")) %>% data.frame
  return(OUT)
  
})

names(Geno_pheno.data) <-  unlist(Pheno_info[, 1])

# Remove rows with NAs
# removes all individual with at least one NA
Geno_pheno.data <- lapply(Geno_pheno.data, na.omit) 

# Note: If there are a lot of missing genotypes in the reference data set,
# skip the lines 82-84 and use missing value imputation presented 
# in lines 167-185. Training and test sets should be imputed separately.

###############################################################################
### Prepare phenotypes for data partitioning

### Remove phenotypes wo at least two level
# List of phenotype frequencies 
Pheno_freq <- lapply(Geno_pheno.data, function(x) { 
  table(x$Phenotype)        # Calculate table of frequencies for each phenotype
})

# Phenotype counts
Pheno_freq_omit <- lapply(Pheno_freq, function(x) {
  if(dim(x) < 2) {    # If dimension of freq table is < 2, there is only one phenotype level
    TRUE
  } else {
    FALSE
  }
}) %>% unlist         

Geno_pheno.data <- Geno_pheno.data[Pheno_freq_omit == FALSE] 
# FALSE means that dimension of freq table > 2

### Remove phenotypes with less than 4 cases in one level
# List of phenotype frequencies again 
Pheno_freq <- lapply(Geno_pheno.data, function(x) {
  table(x$Phenotype)        # Calculate table of frequencies for each phenotype                     
})

Pheno_freq_omit <- lapply(Pheno_freq, function(x) {
  if(any(x < 4)) {
    TRUE
  } else {
    FALSE
  }
}) %>% unlist 

Geno_pheno.data <- Geno_pheno.data[Pheno_freq_omit == FALSE]

###############################################################################
### Create data partitioning: split in training and test data 1:1, 
### 50% of each phenotype into both data sets

# Create data partitioning index
set.seed(3456)

# Change phenotype names - to _minus + to _plus
lapply(1:length(Geno_pheno.data), function(x) {
  Geno_pheno.data[[x]]$Phenotype <<- str_replace(Geno_pheno.data[[x]]$Phenotype,
                                                 "-", "_minus")
  Geno_pheno.data[[x]]$Phenotype <<- str_replace(Geno_pheno.data[[x]]$Phenotype,
                                                 "\\+", "_plus")
})

TrainIndex <- lapply(Geno_pheno.data, function(x) {
  createDataPartition(x$Phenotype, 
                      p = 0.5,              # 50% of data in the training set
                      list = FALSE,         # Avoids returning the data as a list
                      times = 1)            # Only one split 
})  
# The output is a list, elements are set of integers for 
# the rows of each phenotype that belong in the training set.                              

# Create train and test data
Train_data <- map(names(Geno_pheno.data), function(x) { 
  IDX <- TrainIndex[[x]]                         # x = Blood group name  
  DATA <- Geno_pheno.data[[x]][IDX, ]
  return(DATA)
})

Test_data <- map(names(Geno_pheno.data), function(x){
  IDX <- TrainIndex[[x]]                         # x = Blood group name
  DATA <- Geno_pheno.data[[x]][-IDX, ]           # Negative indexing TrainIndex
  return(DATA)
})

names(Train_data) <- names(TrainIndex)
names(Test_data) <- names(TrainIndex)

###############################################################################
### Impute missing values either by mean imputation
### Train and test sets separately

# # Mean imputation
# meanImpute <- function(x) {
#   # x cols 1-2 are character, the rest are numeric
#   x[, -c(1:2)] <- x[, -c(1:2)] %>% apply(., 2, function(y) {
#     y[is.na(y)] <- mean(y, na.rm=T)
#     return(y)
#   })
#   return(x)
# }
# Geno_pheno.data <- lapply(Geno_pheno.data, meanImpute) 
# 
# # Mean impute train and test sets
# Train_data <- lapply(Train_data, meanImpute) 
# Test_data <- lapply(Test_data, meanImpute) 

###############################################################################
### Create function including class weight generation and variant selection

Train_model <- function(x) {              # x = data frame element in Train data
  ID <- x$releasedIdentifier
  x <- select(x, -releasedIdentifier)
  x$Phenotype <- as.factor(x$Phenotype)          # Change Phenotype to factor 
  CLASS_WEIGHT <- 1 - ( table(x$Phenotype) / nrow(x) ) 
  MODEL <- ranger(Phenotype ~ ., 
                  data = x, 
                  num.trees = 2000,
                  mtry = 0.5 * ncol(x),                   # Half of the variants                
                  importance = "permutation",
                  write.forest = T, 
                  class.weights = CLASS_WEIGHT)
  VARIANTS <- as.data.frame(MODEL[["variable.importance"]]) %>%
    rownames_to_column() %>%
    rename(Variant = rowname, Importance = 'MODEL[["variable.importance"]]')
  BEST_VARIANTS <- filter(VARIANTS, Importance > 0)
  BEST_GENO_PHENO <- select(x, Phenotype, BEST_VARIANTS$Variant)
  MODEL <- ranger(Phenotype ~ ., 
                  data = BEST_GENO_PHENO, 
                  num.trees = 2000,
                  mtry = 0.5 * ncol(BEST_GENO_PHENO),                               
                  importance = "permutation",
                  probability = T, # Grow a probability forest as in Malley et al. (2012)
                  write.forest = T, 
                  class.weights = CLASS_WEIGHT)
  MODEL$predictions <- as.data.frame(MODEL$predictions)
  MODEL$predictions$releasedIdentifier <- ID  # Insert releasedIdentifier column 
  # to the OOB prediction results
  return(MODEL)
}

### Remove samples from Train_data with missing phenotype information
# expects the phenotype to be the first column
Train_data <- map(Train_data, function(x) filter(x, !is.na(x[, 1])))

### Fit model to Train data
Train_results <- lapply(Train_data, Train_model)

### Fit model to full data
Full_data_results <- lapply(Geno_pheno.data, Train_model)

### Write models to disk and include blood group gene information for applying the
### the prediction models

# Add blood group gene and phenotype information data to the large list Full_data_results
models_genes <- c(Full_data_results, BG_gene_rf = list(BG_gene_rf))
models_genes <- c(models_genes, Pheno_info = list(Pheno_info))
# Save in RDS format
saveRDS(models_genes, file = "results/Blood_group_prediction_models_genes.rds")

### Predict Test data
Test_results <- lapply(names(Test_data), function(x) {
  out <- predict(Train_results[[x]], data = Test_data[[x]])
  out$predictions <- 
    data.frame(releasedIdentifier = Test_data[[x]]$releasedIdentifier,
               out$predictions)
  return(out)
})

# Add names to the Test_results list
names(Test_results) <- names(Test_data)

###############################################################################
### Modify test data predictions and observed phenotypes into one data frame

# Create data frames with the test data prediction results and the correct value
Preds_obs <- lapply(1:length(Test_results), function(x) { 
  # Preds = predictions, obs = observations
  inner_join(Test_results[[x]]$predictions, 
             Test_data[[x]][c("releasedIdentifier", "Phenotype")], 
             by = "releasedIdentifier")
})

# Add names to the Preds_obs list
names(Preds_obs) <- names(Test_results)

# Categorical observed phenotype variables into numeric columns
lapply(1:length(Preds_obs), function(x) {
  Preds_obs[[x]] <<- data.frame(Preds_obs[[x]], 
                                model.matrix(~ 0 + Phenotype, data = Preds_obs[[x]]))
})                                              # ~ 0 means no intercept included
# Results columns e.g."PhenotypeCx_minus" and "PhenotypeCx_plus"

### Create class factors from test data predictions post probabilities
Preds_obs <- map(names(Preds_obs), function(x) {      
  DATA <- Preds_obs[[x]]
  PRED <- select(DATA, -releasedIdentifier, -contains("Phenotype")) 
  PRED_classes <- map_dfc(colnames(PRED), function(y) { 
    ifelse(PRED[, y] > 0.5, 1, 0)
  })                                                  
  colnames(PRED_classes) <- paste(colnames(PRED), "class", sep = "_")
  cbind(DATA, PRED_classes) %>% return()
})

# Add names to the Preds_obs list
names(Preds_obs) <- names(Test_results)

###############################################################################
### Modify train data OOB predictions and observed phenotypes into one data frame

# Create data frames with the train data OOB prediction results and the correct value
OOB_obs <- lapply(1:length(Train_results), function(x) {
  inner_join(Train_results[[x]]$predictions, 
             Train_data[[x]][c("releasedIdentifier", "Phenotype")], 
             by = "releasedIdentifier")
})

# Add names to the OOB_obs list
names(OOB_obs) <- names(Train_results)

# Categorical observed phenotype variables into numeric columns
lapply(1:length(OOB_obs), function(x) {
  OOB_obs[[x]] <<- data.frame(OOB_obs[[x]], 
                              model.matrix(~ 0 + Phenotype, data = OOB_obs[[x]]))
})                                              # ~ 0 means no intercept included
# Results columns e.g."PhenotypeCx_minus" and "PhenotypeCx_plus"

### Create class factors from train data OOB predictions
OOB_obs <- map(names(OOB_obs), function(x) {      
  DATA <- OOB_obs[[x]]
  PRED <- select(DATA, -releasedIdentifier, -contains("Phenotype")) 
  PRED_classes <- map_dfc(colnames(PRED), function(y) { 
    ifelse(PRED[, y] > 0.5, 1, 0)
  })                                                  
  colnames(PRED_classes) <- paste(colnames(PRED), "class", sep = "_")
  cbind(DATA, PRED_classes) %>% return()
})

# Add names to the OOB_obs list
names(OOB_obs) <- names(Train_results)

###############################################################################
### Modify full data OOB predictions and observed phenotypes into one data frame

# Create data frames with the full data OOB prediction results and the correct value
Full_OOB_obs <- lapply(1:length(Full_data_results), function(x) {
  inner_join(Full_data_results[[x]]$predictions, 
             Geno_pheno.data[[x]][c("releasedIdentifier", "Phenotype")], 
             by = "releasedIdentifier")
})

# Add names to the Full_OOB_obs list
names(Full_OOB_obs) <- names(Full_data_results)

# Categorical observed phenotype variables into numeric columns
lapply(1:length(Full_OOB_obs), function(x) {
  Full_OOB_obs[[x]] <<- data.frame(Full_OOB_obs[[x]], 
                                   model.matrix(~ 0 + Phenotype, 
                                                data = Full_OOB_obs[[x]]))
})                                              # ~ 0 means no intercept included
# Results columns e.g."PhenotypeCx_minus" and "PhenotypeCx_plus"

### Create class factors from full data OOB predictions
Full_OOB_obs <- map(names(Full_OOB_obs), function(x) {      
  DATA <- Full_OOB_obs[[x]]
  PRED <- select(DATA, -releasedIdentifier, -contains("Phenotype")) 
  PRED_classes <- map_dfc(colnames(PRED), function(y) { 
    ifelse(PRED[, y] > 0.5, 1, 0)
  })                                                  
  colnames(PRED_classes) <- paste(colnames(PRED), "class", sep = "_")
  cbind(DATA, PRED_classes) %>% return()
})

# Add names to the Full_OOB_obs list
names(Full_OOB_obs) <- names(Full_data_results)

###############################################################################
### Investigation of model accuracies and errors

### Test result confusion matrix and statistics

Conf_mat_test <- map(names(Preds_obs), function(x) {
  DATA <- Preds_obs[[x]]                         # Create Preds_obs data frame
  CONF <- select(DATA, -releasedIdentifier, -Phenotype) %>% 
    select(contains(c("Phenotype", "class"))) # Phenotype is reference, class is prediction
  PHENO_COLS <- which(grepl("Phenotype", colnames(CONF))) 
  # Grepl returns T/F for observed phenotype colums, 
  # which returns number of these columns
  Conf_mat_test <- map(PHENO_COLS, function(y) {
    CF <- confusionMatrix(data = CONF[, y + length(PHENO_COLS)] %>% factor, 
                          # Prediction class as factor
                          reference = CONF[, y] %>% factor)
    # Observed phenotype as factor
    return(CF)
  })
  names(Conf_mat_test) <- colnames(CONF)[PHENO_COLS] %>% gsub("Phenotype", "", .)           
  # Removes "Phenotype" from Conf_mat_test names
  return(Conf_mat_test)
})

# Add names to Conf_mat list and save as RDS
names(Conf_mat_test) <- names(Preds_obs)
saveRDS(Conf_mat_test, file = "results/Conf_mat_test.rds")

### Train result confusion matrix and statistics
Conf_mat_train <- map(names(OOB_obs), function(x) {
  print(x)
  DATA <- OOB_obs[[x]]                        
  CONF <- select(DATA, -releasedIdentifier, -Phenotype) %>% 
    select(contains(c("Phenotype", "class"))) 
  PHENO_COLS <- which(grepl("Phenotype", colnames(CONF))) 
  Conf_mat_train <- map(PHENO_COLS, function(y) {
    CF <- confusionMatrix(data = CONF[, y + length(PHENO_COLS)] %>% factor, 
                          reference = CONF[, y] %>% factor)                   
    return(CF)
  })
  names(Conf_mat_train) <- colnames(CONF)[PHENO_COLS] %>% gsub("Phenotype", "", .)          
  return(Conf_mat_train)
})

# Add names to Cof_mat_train list and save as RDS
names(Conf_mat_train) <- names(OOB_obs)
saveRDS(Conf_mat_train, file = "results/Conf_mat_train.rds")

### Full data result confusion matrix and statistics
Conf_mat_full <- map(names(Full_OOB_obs), function(x) {
  print(x)
  DATA <- Full_OOB_obs[[x]]                         
  CONF <- select(DATA, -releasedIdentifier, -Phenotype) %>%  
    select(contains(c("Phenotype", "class")))        
  PHENO_COLS <- which(grepl("Phenotype", colnames(CONF))) 
  Conf_mat_full <- map(PHENO_COLS, function(y) {
    CF <- confusionMatrix(data = CONF[, y + length(PHENO_COLS)] %>% factor, 
                          reference = CONF[, y] %>% factor)                 
    return(CF)
  })
  names(Conf_mat_full) <- colnames(CONF)[PHENO_COLS] %>% gsub("Phenotype", "", .)         
  return(Conf_mat_full)
})

# Add names to Cof_mat_full list and save as RDS
names(Conf_mat_full) <- names(Full_OOB_obs)
saveRDS(Conf_mat_full, file = "results/Conf_mat_full.rds")

###############################################################################

