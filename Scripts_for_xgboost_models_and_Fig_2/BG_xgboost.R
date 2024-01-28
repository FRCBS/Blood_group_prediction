
library(tidyverse)
library(data.table)
library(xgboost)
library(caret)
library(ggbeeswarm)


### Fit xgboost models for glood group antigen prediction
# Required objects are generated using "Creating_prediction_models" scripts

###############################################################################
### Functions
###############################################################################

Train_xgboost <- function(x) { # x = data frame element in Train data
  
  x <- na.omit(x)
  ID          <- x$releasedIdentifier
  x$Phenotype <- as.factor(x$Phenotype) %>% as.numeric # Change Phenotype to factor 
  x$Phenotype  <- x$Phenotype-1
  
  # how many boosting rounds before overfitting?
  mean.rounds <- map(1:100, function(i) {
    
    train.ind <- sample(1:nrow(x), nrow(x)*(2/3), replace = F)
    test.ind  <- (1:nrow(x))[-train.ind]
    
    Dtrain    <- xgb.DMatrix(data = dplyr::select(x, -c(releasedIdentifier, Phenotype)) %>% 
                               as.matrix %>% .[train.ind, ], 
                             label = x$Phenotype[train.ind])
    Dtest     <- xgb.DMatrix(data = dplyr::select(x, -c(releasedIdentifier, Phenotype)) %>% 
                               as.matrix %>% .[test.ind, ], 
                             label = x$Phenotype[test.ind])
    
    MODEL     <- xgb.train(Dtrain, watchlist = list(T1=Dtest), early_stopping_rounds = 4, 
                           maximize = F, nrounds = 200,
                           params = list(objective = "binary:logistic",
                                         subsample = 1,
                                         colsample_bytree = 1))
    MODEL$evaluation_log$T1_logloss %>% which.min() %>% return()
    
  }) %>% unlist %>% mean()     
  
  # model on full data
  xgboost(xgb.DMatrix(data = dplyr::select(x, -c(releasedIdentifier, Phenotype)) %>% as.matrix, 
                      label = x$Phenotype),
          nrounds = mean.rounds, objective = "binary:logistic") %>% return()
  
}



###############################################################################
### Phenotype data preparations for imputation of blood group antigens
###############################################################################

library(tidyverse)
library(data.table)

# BG_phenos and Pheno_info are generated using instructions 
# in 002_Phenotype_preparation.R

str(BG_phenos)
# Structure specification example
# tibble [1,192 × 54] (S3: tbl_df/tbl/data.frame), columns:
# 1. $ releasedIdentifier: chr [1:1192] Sample1" "Sample2" "Sample3" "Sample4" ...
# 2. $ ABO               : chr [1:1192] "A" "B" "AB" "O" ...
# 3. $ D                 : chr [1:1192] "D+" "D-" "D+" "D+" ...
# 4. $ Kpa               : chr [1:1192] "Kpa-" "Kpa-" "Kpa-" "Kpa-" ...
# 5. $ Kpb               : chr [1:1192] "Kpb+" "Kpb+" "Kpb+" "Kpb+" ...
# 6. ...  n. column, include all the blood group information into columns   

str(Pheno_info)
# Structure specification example
# tibble [53 × 7] (S3: tbl_df/tbl/data.frame)
# $ Phenotype    : chr [1:53] "ABO" "A1" "A2" "Yta" ...
# $ BG_system    : chr [1:53] "ABO" "ABO" "ABO" "Cartwright" ...
# $ Pheno_genes  : chr [1:53] "ABO" "ABO" "ABO" "ACHE" ...
# $ Gene1        : chr [1:53] "ABO" "ABO" "ABO" "ACHE" ...
# $ Gene2        : chr [1:53] NA NA NA NA ...
# $ Gene3        : chr [1:53] NA NA NA NA ...
# $ Plink_dataset: chr [1:53] "data/BG_genes/ABO_ABO_2K.vcf" "data/BG_genes/ABO_ABO_2K.vcf" ...

# Column Pheno_genes contains all genes separated by _ e.g. GYPA_GYPB_GYPE, RHCE_RHD
# Column Plink_dataset contains PLINK dataset prefixes generated using 001_Genotype_preparation.R script


###############################################################################
### Blood group imputation using gradient boosting models
###############################################################################

### Remove samples from Train_data with missing phenotype information

# expects the phenotype to be the first column
Train_data <- map(Train_data, function(x) filter(x, !is.na(x[, 1])))


### Fit model to Train data

# fit models for non-ABO
Train_results <- lapply(Train_data[-1], Train_xgboost)

# fit ABO separately
# creaste training data for ABO
Train_data_ABO <- map(c('A', 'O', 'B', 'AB'), function(z) {
  out <- Train_data[['ABO']]
  out[, 'Phenotype'] <- ifelse(out[, 'Phenotype']==z, 1, 0)
  return(out)
})
names(Train_data_ABO) <- c('A', 'O', 'B', 'AB')

# fit models for ABO
Train_results_ABO <- lapply(Train_data_ABO, Train_xgboost)
names(Train_results_ABO) <- c('A', 'O', 'B', 'AB')


### Predict Test data

Test_results <- lapply(names(Test_data)[-1], function(x) {
  # 
  out <- predict(Train_results[[x]], 
                 newdata = Test_data[[x]] %>% dplyr::select( -c(releasedIdentifier, Phenotype)) %>% as.matrix
                 ) # JRRI edit
  out$predictions <- data.frame(
    releasedIdentifier = Test_data[[x]]$releasedIdentifier,
    predictions=out)
  return(out)
})
# Add names to the Test_results list
names(Test_results) <- names(Test_data)[-1]

# predict ABO

# create test data for ABO
Test_data_ABO <- map(c('A', 'O', 'B', 'AB'), function(z) {
  out <- Test_data[['ABO']]
  out[, 'Phenotype'] <- ifelse(out[, 'Phenotype']==z, 1, 0)
  return(out)
})
names(Test_data_ABO) <- c('A', 'O', 'B', 'AB')

Test_results_ABO <- lapply(names(Test_data_ABO), function(x) {
  # 
  out <- predict(Train_results_ABO[[x]], 
                 newdata = Test_data_ABO[[x]] %>% dplyr::select( -c(releasedIdentifier, Phenotype)) %>% as.matrix
  ) # JRRI edit
  out$predictions <- data.frame(
    releasedIdentifier = Test_data_ABO[[x]]$releasedIdentifier,
    predictions=out)
  return(out)
})
# Add names to the Test_results list
names(Test_results_ABO) <- names(Test_data_ABO)


###############################################################################
### Modify test data predictions and observed phenotypes into one data frame

# Create data frames with the test data prediction results and the correct value
Preds_obs <- lapply(names(Test_results), function(x) { 
  # Preds = predictions, obs = observations
  inner_join(Test_results[[x]]$predictions, 
             Test_data[[x]][c("releasedIdentifier", "Phenotype")], 
             by = "releasedIdentifier") %>% na.omit()
})

# Add names to the Preds_obs list
names(Preds_obs) <- names(Test_results)

# Categorical observed phenotype variables into numeric columns
lapply(names(Preds_obs), function(x) {
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
### Modify test data predictions and observed phenotypes into one data frame
### ABO

# Create data frames with the test data prediction results and the correct value
Preds_obs_ABO <- lapply(names(Test_results_ABO), function(x) { 
  # Preds = predictions, obs = observations
  inner_join(Test_results_ABO[[x]]$predictions, 
             Test_data_ABO[[x]][c("releasedIdentifier", "Phenotype")], 
             by = "releasedIdentifier") %>% na.omit()
})

# Add names to the Preds_obs list
names(Preds_obs_ABO) <- names(Test_results_ABO)

# Categorical observed phenotype variables into numeric columns
lapply(names(Preds_obs_ABO), function(x) {
  Preds_obs_ABO[[x]] <<- data.frame(Preds_obs_ABO[[x]], 
                                model.matrix(~ 0 + Phenotype, data = Preds_obs_ABO[[x]]))
})                                              # ~ 0 means no intercept included
# Results columns e.g."PhenotypeCx_minus" and "PhenotypeCx_plus"

### Create class factors from test data predictions post probabilities
Preds_obs_ABO <- map(names(Preds_obs_ABO), function(x) {      
  DATA <- Preds_obs_ABO[[x]]
  PRED <- select(DATA, -releasedIdentifier, -contains("Phenotype")) 
  PRED_classes <- map_dfc(colnames(PRED), function(y) { 
    ifelse(PRED[, y] > 0.5, 1, 0)
  })                                                  
  colnames(PRED_classes) <- paste(colnames(PRED), "class", sep = "_")
  cbind(DATA, PRED_classes) %>% return()
})

# Add names to the Preds_obs list
names(Preds_obs_ABO) <- names(Test_results_ABO)


###############################################################################
### Investigation of model accuracies and errors

### Test result confusion matrix and statistics

Conf_mat_test <- map(names(Preds_obs), function(x) {
  # x <- 'Ytb'
  print(x)
  DATA <- Preds_obs[[x]]                         # Create Preds_obs data frame
  CONF <- select(DATA, -releasedIdentifier, -Phenotype) %>% 
    select(contains(c("Phenotype", "class"))) # Phenotype is reference, class is prediction
  PHENO_COLS <- which(grepl("Phenotype.*_plus", colnames(CONF))) 
  # Grepl returns T/F for observed phenotype colums, 
  # which returns number of these columns
  Conf_mat_test <- map(PHENO_COLS, function(y) {
    # y <- 2
    CF <- confusionMatrix(data = CONF[, y + 1] %>% factor, 
                          # Prediction class as factor
                          reference = CONF[, y] %>% factor)
    # Observed phenotype as factor
    return(CF)
  })
  names(Conf_mat_test) <- colnames(CONF)[PHENO_COLS] %>% gsub("Phenotype", "", .)           
  # Removes "Phenotype" from Conf_mat_test names
  return(Conf_mat_test)
})
# Add names to Conf_mat list
names(Conf_mat_test) <- names(Preds_obs)

# ABO conf mat
Conf_mat_test_ABO <- map(names(Preds_obs_ABO), function(x) {
  # x <- 'A'
  print(x)
  DATA <- Preds_obs_ABO[[x]]                         # Create Preds_obs data frame
  CONF <- select(DATA, -releasedIdentifier, -Phenotype) %>% 
    select(contains(c("Phenotype", "class"))) # Phenotype is reference, class is prediction
  PHENO_COLS <- which(grepl("Phenotype", colnames(CONF))) 
  # Grepl returns T/F for observed phenotype colums, 
  # which returns number of these columns
  Conf_mat_test <- map(PHENO_COLS, function(y) {
    # y <- 2
    CF <- confusionMatrix(data = CONF[, y + 1] %>% factor, 
                          # Prediction class as factor
                          reference = CONF[, y] %>% factor)
    # Observed phenotype as factor
    return(CF)
  })
  names(Conf_mat_test) <- colnames(CONF)[PHENO_COLS] %>% gsub("Phenotype", "", .)           
  # Removes "Phenotype" from Conf_mat_test names
  return(Conf_mat_test)
})
# Add names to Conf_mat list
names(Conf_mat_test_ABO) <- names(Preds_obs_ABO)


###############################################################################
### Plots of antigen-level accuracies


# collect CF from test results
test_CF_tables <- map2_dfr(c(Conf_mat_test_ABO, Conf_mat_test), 
                           c(Conf_mat_test_ABO, Conf_mat_test) %>% names, 
                           function(x, y) {
                             TABLE <- x[[1]][["table"]] %>% data.frame
                             TABLE <- mutate(TABLE, goodbad = ifelse(TABLE$Prediction == TABLE$Reference, 
                                                                     "good", "bad")) %>%
                               group_by(Reference) %>%
                               mutate(prop = Freq/sum(Freq))
                             TABLE$Antigen <- y
                             return(TABLE)
                           })
test_CF_tables$Antigen <- factor(test_CF_tables$Antigen, 
                                 levels = test_CF_tables$Antigen %>% unique)
fwrite(test_CF_tables, 'results/xgboost_test_CF_tables.tsv', sep='\t')
test_CF_tables <- fread('results/xgboost_test_CF_tables.tsv')

# Collect test metrics from results
test_metrics <- map2(c(Conf_mat_test_ABO, Conf_mat_test), 
                         c(Conf_mat_test_ABO, Conf_mat_test) %>% names, 
                         function(x, y) {
                           # x <- c(Conf_mat_test_ABO, Conf_mat_test)[[1]]
                           out <- cbind(x[[1]][3][[1]] %>% t %>% data.frame() %>% .[, c(1,3)],
                                        x[[1]][4][[1]] %>% t %>% data.frame()) %>% t
                           colnames(out) <- y
                           return(out)
                         }) %>% do.call(cbind, .)
test_metrics <- data.frame(Metric=rownames(test_metrics), test_metrics)
# write accuracy metrics table
test_metrics[, -1] %>% t %>% data.frame %>% rownames_to_column %>% rename(Antigen=rowname) %>% 
  fwrite(., 'results/xgboost_test_metrics_table.tsv', sep='\t')
# convert to long format
test_metrics <- test_metrics %>% pivot_longer(2:ncol(test_metrics))
test_metrics$name <- factor(test_metrics$name, test_metrics$name %>% unique)
fwrite(test_metrics, 'results/xgboost_test_metrics.tsv', sep='\t')
test_metrics <- fread('results/xgboost_test_metrics.tsv')


# Collect PPs from results
# post probs
PP_test <- map2(c(Preds_obs_ABO, Preds_obs),
                c(Preds_obs_ABO, Preds_obs) %>% names, 
                function(x, y) {
                  if(y %in% c('A', 'O', 'B', 'AB')) {
                    out <- data.frame(Antigen=y, x[, 2:3]) 
                  } else out <- data.frame(Antigen=y, x[, c(2, 5)]) 
                  colnames(out)[2:3] <- c('predictions', 'Phenotype')
                  out %>% return()
                }) %>% do.call(rbind, .)
PP_test$Antigen <- factor(PP_test$Antigen, 
                          levels = PP_test$Antigen %>% unique)
PP_test$Phenotype <- factor(PP_test$Phenotype, levels = c('0', '1'))
fwrite(PP_test, 'results/xgboost_test_PP_test.tsv', sep='\t')
PP_test <- fread('results/xgboost_test_PP_test.tsv')



# Plot CF tables 

xgboost_test_data_CF_table <- test_CF_tables %>% 
  mutate(BG_system = case_when(grepl("^A|^B|AB|O|^A1|^A2", Antigen) ~ "ABO",
                               grepl("Coa|Cob", Antigen) ~ "Colton",
                               grepl("Doa|Dob", Antigen) ~ "Dombrock",                     
                               grepl("^D|^C|^c|^E|^e|Cw|Cx|hrS|hrB", Antigen) ~ "Rh",
                               grepl("K|Kpa|Ula", Antigen) ~"Kell",
                               grepl("Jka|Jkb", Antigen) ~ "Kidd",
                               grepl("Fya|Fyb", Antigen) ~ "Duffy",
                               grepl("^M|^N|^S|^s", Antigen) ~ "MNS",
                               grepl("Ytb", Antigen) ~ "Cartwright",
                               grepl("Lua", Antigen) ~ "Lutheran",
                               grepl("LWb", Antigen) ~ "Landsteiner-Wiener",
                               grepl("^HPA1a|^HPA1b", Antigen) ~ "HPA",
                               grepl("P1", Antigen) ~ "P1PK",
                               grepl("Lea|Leb", Antigen) ~ "Lewis",
                               grepl("Lsa", Antigen) ~ "Gerbich" ))
xgboost_test_data_CF_table$Reference <- xgboost_test_data_CF_table$Reference %>% factor
xgboost_test_data_CF_table$Prediction <- xgboost_test_data_CF_table$Prediction %>% factor

p.full_data_CF_plots <- ggplot(data = xgboost_test_data_CF_table %>% tibble, 
                               aes(x = Reference, y = Prediction, 
                                   fill = goodbad, alpha = 0.5)) +
  geom_tile() + 
  geom_text(aes(label = Freq), vjust = .5, alpha = 1, size = 4) +
  scale_fill_manual(values = c(good = "green3", bad = "firebrick2")) +
  theme_bw() +
  xlim(rev(levels(xgboost_test_data_CF_table$Reference))) + 
  facet_wrap(BG_system~Antigen) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12)) +
  theme(panel.grid=element_blank()) + 
  theme(legend.position = "none")

jpeg('results/Confusion_matrices_xgboost_test_data_by_BG.jpeg', 
     width=15, height=15, res=600, units='in')
print(p.full_data_CF_plots)
dev.off()



# Plot Post probs

xgboost_test_BG_PP_table <- PP_test %>% 
  mutate(BG_system = case_when(grepl("^A|^B|AB|O|^A1|^A2", Antigen) ~ "ABO",
                               grepl("Coa|Cob", Antigen) ~ "Colton",
                               grepl("Doa|Dob", Antigen) ~ "Dombrock",                     
                               grepl("^D|^C|^c|^E|^e|Cw|Cx|hrS|hrB", Antigen) ~ "Rh",
                               grepl("K|Kpa|Ula", Antigen) ~"Kell",
                               grepl("Jka|Jkb", Antigen) ~ "Kidd",
                               grepl("Fya|Fyb", Antigen) ~ "Duffy",
                               grepl("^M|^N|^S|^s", Antigen) ~ "MNS",
                               grepl("Ytb", Antigen) ~ "Cartwright",
                               grepl("Lua", Antigen) ~ "Lutheran",
                               grepl("LWb", Antigen) ~ "Landsteiner-Wiener",
                               grepl("^HPA1a|^HPA1b", Antigen) ~ "HPA",
                               grepl("P1", Antigen) ~ "P1PK",
                               grepl("Lea|Leb", Antigen) ~ "Lewis",
                               grepl("Lsa", Antigen) ~ "Gerbich" ))
xgboost_test_BG_PP_table <- xgboost_test_BG_PP_table %>% rename(Posterior_probability=predictions, Reference=Phenotype)
xgboost_test_BG_PP_table$Reference <- factor(xgboost_test_BG_PP_table$Reference)

p.test_BG_pp_boxplot <- ggplot(xgboost_test_BG_PP_table, 
                               aes(Reference, Posterior_probability)) +
  geom_quasirandom(dodge.width=0.8, shape=21) +
  geom_hline(yintercept=0.5, color="grey", linetype="dashed", size=0.4) +
  ylab('Posterior probability') + xlab('Reference') +
  facet_wrap(BG_system ~ Antigen) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12)) +
  theme(panel.grid=element_blank()) + 
  theme(legend.position = "none")+
  theme(panel.border=element_rect(fill=NA, size=0.3), 
        panel.grid=element_blank())

jpeg('./results/PP_boxplots_xgboost_test_data_by_BG.jpeg', 
     width=15, height=15, res=600, units='in')
print(p.test_BG_pp_boxplot)
dev.off()



# Plot Metrics

colnames(test_metrics)[2] <- 'Antigen'
test_metrics <- filter(test_metrics, 
                       Metric %in% c("Accuracy","AccuracyLower","Sensitivity","Specificity",
                                     "Pos.Pred.Value","Neg.Pred.Value","F1","Prevalence","Balanced.Accuracy"))
test_metrics <- test_metrics %>% dplyr::select(-Antigen) %>% group_by(., Metric) %>% 
  summarise(Mmean=mean(value)) %>% left_join(., test_metrics) %>% arrange(Mmean)
test_metrics$Metric <- factor(test_metrics$Metric, levels = test_metrics$Metric %>% unique)
test_metrics <- rename(test_metrics, Antigen=name)

xgboost_test_BG_Metric_table <- test_metrics %>% 
  mutate(BG_system = case_when(grepl("^A|^B|AB|O|^A1|^A2", Antigen) ~ "ABO",
                               grepl("Coa|Cob", Antigen) ~ "Colton",
                               grepl("Doa|Dob", Antigen) ~ "Dombrock",                     
                               grepl("^D|^C|^c|^E|^e|Cw|Cx|hrS|hrB", Antigen) ~ "Rh",
                               grepl("K|Kpa|Ula", Antigen) ~"Kell",
                               grepl("Jka|Jkb", Antigen) ~ "Kidd",
                               grepl("Fya|Fyb", Antigen) ~ "Duffy",
                               grepl("^M|^N|^S|^s", Antigen) ~ "MNS",
                               grepl("Ytb", Antigen) ~ "Cartwright",
                               grepl("Lua", Antigen) ~ "Lutheran",
                               grepl("LWb", Antigen) ~ "Landsteiner-Wiener",
                               grepl("^HPA1a|^HPA1b", Antigen) ~ "HPA",
                               grepl("P1", Antigen) ~ "P1PK",
                               grepl("Lea|Leb", Antigen) ~ "Lewis",
                               grepl("Lsa", Antigen) ~ "Gerbich" ))


p.test_BG_metrics <- xgboost_test_BG_Metric_table %>% ggplot(aes(value, Metric, color=value)) +
  geom_point(alpha=0.8, size = 2) +
  xlab('Value') +
  scale_color_gradient(high="green4", low="firebrick2") +
  facet_wrap(BG_system ~ Antigen) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12)) +
  theme(panel.grid=element_blank()) + 
  theme(legend.position = "none")+
  theme(panel.border=element_rect(fill=NA, size=0.3), 
        panel.grid=element_blank())

jpeg('./results/Metric_plots_xgboost_test_data_by_BG.jpeg', 
     width=15, height=15, res=600, units='in')
print(p.test_BG_metrics)
dev.off()


