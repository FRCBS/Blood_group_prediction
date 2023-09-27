###############################################################################
### Blood group prediction figures and table
###############################################################################

library(tidyverse)
library(caret)
library(ranger)
library(data.table)
library(e1071)
library(ggbeeswarm)

###############################################################################
## Variables for this script are generated 
## in the script 003_Blood_group_prediction.R and they are assumed to be 
## available in the global environment.
###############################################################################
### Confusion matrix plots

### Test data: All blood groups confusion matrix plots

# Create test data ABO confusion matrix table
ABO_CF_tables <- map_dfr(1:length(Conf_mat_test[["ABO"]]), function(x) {
  TABLE <- Conf_mat_test[["ABO"]][[x]][["table"]] %>% data.frame
  TABLE <- mutate(TABLE, goodbad = ifelse(TABLE$Prediction == TABLE$Reference, 
                                          "good", "bad")) %>%
    group_by(Reference) %>%
    mutate(prop = Freq/sum(Freq))
  TABLE$Antigen <- names(Conf_mat_test[["ABO"]][x])
  return(TABLE)
})

# Create test data other blood group confusion matrix table
BG_CF_tables <- map_dfr(2:length(Conf_mat_test), function(x) {
  TABLE <- Conf_mat_test[[x]][[2]][["table"]] %>% data.frame
  TABLE <- mutate(TABLE, goodbad = ifelse(TABLE$Prediction == TABLE$Reference, 
                                          "good", "bad")) %>%
    group_by(Reference) %>%
    mutate(prop = Freq/sum(Freq))
  TABLE$Antigen <- names(Conf_mat_test[[x]][2]) %>% 
    str_replace("_plus", "+")
  return(TABLE)
})

# Join CF tables and create blood antigen system column
# Please update according your own antigen data
Test_data_CF_table <-  bind_rows(ABO_CF_tables, BG_CF_tables) %>% 
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
                               grepl("HPA1a|HPA1b", Antigen) ~ "HPA1",
                               grepl("P1", Antigen) ~ "P1PK",
                               grepl("Lea|Leb", Antigen) ~ "Lewis",
                               grepl("Lsa", Antigen) ~ "Gerbich" ))

p.test_data_CF_plots <- ggplot(data = Test_data_CF_table, 
                               mapping = aes(x = Reference, y = Prediction, 
                                             fill = goodbad, alpha = 0.5)) +
  geom_tile() + 
  geom_text(aes(label = Freq), vjust = .5, alpha = 1, size = 4) +
  scale_fill_manual(values = c(good = "green3", bad = "firebrick2")) +
  theme_bw() +
  xlim(rev(levels(BG_CF_tables$Reference))) + 
  facet_wrap(BG_system~Antigen) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12)) +
  theme(panel.grid=element_blank()) + 
  theme(legend.position = "none")

jpeg('./results/Confusion_matrices_test_by_BG.jpeg', 
     width=15, height=15, res=600, units='in')
print(p.test_data_CF_plots)
dev.off()

### Train data: All blood groups confusion matrix plots

# Create train data ABO confusion matrix table
ABO_CF_tables <- map_dfr(1:length(Conf_mat_train[["ABO"]]), function(x) {
  TABLE <- Conf_mat_train[["ABO"]][[x]][["table"]] %>% data.frame
  TABLE <- mutate(TABLE, goodbad = ifelse(TABLE$Prediction == TABLE$Reference, 
                                          "good", "bad")) %>%
    group_by(Reference) %>%
    mutate(prop = Freq/sum(Freq))
  TABLE$Antigen <- names(Conf_mat_train[["ABO"]][x])
  return(TABLE)
})

# Create train data other blood group confusion matrix plots
BG_CF_tables <- map_dfr(2:length(Conf_mat_train), function(x) {
  TABLE <- Conf_mat_train[[x]][[2]][["table"]] %>% data.frame
  TABLE <- mutate(TABLE, goodbad = ifelse(TABLE$Prediction == TABLE$Reference, 
                                          "good", "bad")) %>%
    group_by(Reference) %>%
    mutate(prop = Freq/sum(Freq))
  TABLE$Antigen <- names(Conf_mat_train[[x]][2]) %>% 
    str_replace("_plus", "+")
  return(TABLE)
})

# Join CF tables and create blood antigen system column
# Please update according your own antigen data
Train_data_CF_table <-  bind_rows(ABO_CF_tables, BG_CF_tables) %>% 
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

p.train_data_CF_plots <- ggplot(data = Train_data_CF_table, 
                                mapping = aes(x = Reference, y = Prediction, 
                                              fill = goodbad, alpha = 0.5)) +
  geom_tile() + 
  geom_text(aes(label = Freq), vjust = .5, alpha = 1, size = 4) +
  scale_fill_manual(values = c(good = "green3", bad = "firebrick2")) +
  theme_bw() +
  xlim(rev(levels(BG_CF_tables$Reference))) + 
  facet_wrap(BG_system~Antigen) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12)) +
  theme(panel.grid=element_blank()) + 
  theme(legend.position = "none")

jpeg('./results/Confusion_matrices_train_by_BG.jpeg', 
     width=15, height=15, res=600, units='in')
print(p.train_data_CF_plots)
dev.off()

### Full data: All blood groups confusion matrix plots

# Create full data ABO confusion matrix table
ABO_CF_tables <- map_dfr(1:length(Conf_mat_full[["ABO"]]), function(x) {
  TABLE <- Conf_mat_full[["ABO"]][[x]][["table"]] %>% data.frame
  TABLE <- mutate(TABLE, goodbad = ifelse(TABLE$Prediction == TABLE$Reference, 
                                          "good", "bad")) %>%
    group_by(Reference) %>%
    mutate(prop = Freq/sum(Freq))
  TABLE$Antigen <- names(Conf_mat_full[["ABO"]][x])
  return(TABLE)
})

# Create full data other blood group confusion matrix plots
BG_CF_tables <- map_dfr(2:length(Conf_mat_full), function(x) {
  TABLE <- Conf_mat_full[[x]][[2]][["table"]] %>% data.frame
  TABLE <- mutate(TABLE, goodbad = ifelse(TABLE$Prediction == TABLE$Reference, 
                                          "good", "bad")) %>%
    group_by(Reference) %>%
    mutate(prop = Freq/sum(Freq))
  TABLE$Antigen <- names(Conf_mat_full[[x]][2]) %>% 
    str_replace("_plus", "+")
  return(TABLE)
})

# Join CF tables and create blood antigen system column
# Please update according your own antigen data
Full_data_CF_table <-  bind_rows(ABO_CF_tables, BG_CF_tables) %>% 
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

p.full_data_CF_plots <- ggplot(data = Full_data_CF_table, 
                               mapping = aes(x = Reference, y = Prediction, 
                                             fill = goodbad, alpha = 0.5)) +
  geom_tile() + 
  geom_text(aes(label = Freq), vjust = .5, alpha = 1, size = 4) +
  scale_fill_manual(values = c(good = "green3", bad = "firebrick2")) +
  theme_bw() +
  xlim(rev(levels(BG_CF_tables$Reference))) + 
  facet_wrap(BG_system~Antigen) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12)) +
  theme(panel.grid=element_blank()) + 
  theme(legend.position = "none")

jpeg('./results/Confusion_matrices_full_data_by_BG.jpeg', 
     width=15, height=15, res=600, units='in')
print(p.full_data_CF_plots)
dev.off()

###############################################################################
### Boxplots with posterior probabilities (PP)

### Test data: All blood groups boxplots

# Create test data ABO PP table
ABO_PP_table <- Preds_obs[["ABO"]] %>% 
  select(-Phenotype, -A_class, -AB_class, -B_class, -O_class)
TEMP_1 <- pivot_longer(ABO_PP_table, cols = 2:5, 
                       names_to = "Antigen", 
                       values_to = "Posterior_probability") %>% 
  select(-PhenotypeA, -PhenotypeAB, -PhenotypeB, -PhenotypeO)
TEMP_2 <- pivot_longer(ABO_PP_table, cols = 6:9, 
                       names_to = "Phenotype", 
                       values_to = "Reference") %>% 
  select(-A, -AB, -B, -O)
ABO_PP_table <- cbind(TEMP_1, TEMP_2[2:3])

# Create test data other blood group PP table
BG_PP_table <- map_dfr(2:length(Preds_obs), function(x) {
  TABLE <- Preds_obs[[x]]
  TABLE <- select(TABLE, -contains("class"), -contains("minus"))
  TABLE$Antigen <- colnames(TABLE)[2] %>% 
    str_replace("_plus", "+")
  colnames(TABLE)[2] <- "Posterior_probability"
  colnames(TABLE)[4] <- "Reference"
  return(TABLE)
})

# Combine test data PP tables
# Please update according your own antigen data
All_BG_PP_table <- bind_rows(ABO_PP_table, BG_PP_table) %>% 
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

p.test_BG_pp_boxplot <- ggplot(All_BG_PP_table, 
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

jpeg('./results/PP_boxplots_test_data_by_BG.jpeg', 
     width=15, height=15, res=600, units='in')
print(p.test_BG_pp_boxplot)
dev.off()

### Train data: All blood groups boxplots

# Create train data ABO PP table

ABO_PP_table <- OOB_obs[["ABO"]] %>% 
  select(-Phenotype, -A_class, -AB_class, -B_class, -O_class)
TEMP_1 <- pivot_longer(ABO_PP_table, cols = 1:4, 
                       names_to = "Antigen", 
                       values_to = "Posterior_probability") %>% 
  select(-PhenotypeA, -PhenotypeAB, -PhenotypeB, -PhenotypeO)
TEMP_2 <- pivot_longer(ABO_PP_table, cols = 6:9, 
                       names_to = "Phenotype", 
                       values_to = "Reference") %>% 
  select(-A, -AB, -B, -O)
ABO_PP_table <- cbind(TEMP_1, TEMP_2[2:3])

# Create train data other blood group PP table
BG_PP_table <- map_dfr(2:length(OOB_obs), function(x) {
  TABLE <- OOB_obs[[x]]
  TABLE <- select(TABLE, -contains("class"), -contains("minus"))
  TABLE$Antigen <- colnames(TABLE)[1] %>% 
    str_replace("_plus", "+")
  colnames(TABLE)[1] <- "Posterior_probability"
  colnames(TABLE)[4] <- "Reference"
  return(TABLE)
})

# Combine train data PP tables
# Please update according your own antigen data
All_BG_PP_table <- bind_rows(ABO_PP_table, BG_PP_table) %>% 
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

p.train_BG_pp_boxplot <- ggplot(All_BG_PP_table, 
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

jpeg('./results/PP_boxplots_train_data_by_BG.jpeg', 
     width=15, height=15, res=600, units='in')
print(p.train_BG_pp_boxplot)
dev.off()

### Full data: All blood groups boxplots

# Create all data ABO PP table

ABO_PP_table <- Full_OOB_obs[["ABO"]] %>% 
  select(-Phenotype, -A_class, -AB_class, -B_class, -O_class)
TEMP_1 <- pivot_longer(ABO_PP_table, cols = 1:4, 
                       names_to = "Antigen", 
                       values_to = "Posterior_probability") %>% 
  select(-PhenotypeA, -PhenotypeAB, -PhenotypeB, -PhenotypeO)
TEMP_2 <- pivot_longer(ABO_PP_table, cols = 6:9, 
                       names_to = "Phenotype", 
                       values_to = "Reference") %>% 
  select(-A, -AB, -B, -O)
ABO_PP_table <- cbind(TEMP_1, TEMP_2[2:3])

# Create all data other blood group PP table
BG_PP_table <- map_dfr(2:length(Full_OOB_obs), function(x) {
  TABLE <- Full_OOB_obs[[x]]
  TABLE <- select(TABLE, -contains("class"), -contains("minus"))
  TABLE$Antigen <- colnames(TABLE)[1] %>% 
    str_replace("_plus", "+")
  colnames(TABLE)[1] <- "Posterior_probability"
  colnames(TABLE)[4] <- "Reference"
  return(TABLE)
})

# Combine all data PP tables
# Please update according your own antigen data
All_BG_PP_table <- bind_rows(ABO_PP_table, BG_PP_table) %>% 
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

p.all_data_BG_pp_boxplot <- ggplot(All_BG_PP_table, 
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

jpeg('./results/PP_boxplots_full_data_by_BG.jpeg', 
     width=15, height=15, res=600, units='in')
print(p.all_data_BG_pp_boxplot)
dev.off()

############################################################
### Statistics tables

### Test data: All blood groups statistics

# Create test data ABO stats table
ABO_stat_table <- map_dfr(1:length(Conf_mat_test[["ABO"]]), function(x) {
  TABLE <- Conf_mat_test[["ABO"]][[x]][["byClass"]][["Sensitivity"]] %>%
    data.frame()
  colnames(TABLE) <- "Sensitivity"
  TABLE$Specificity <- Conf_mat_test[["ABO"]][[x]][["byClass"]][["Specificity"]]
  TABLE$PosPredValue <- Conf_mat_test[["ABO"]][[x]][["byClass"]][["Pos Pred Value"]]
  TABLE$NegPredValue <- Conf_mat_test[["ABO"]][[x]][["byClass"]][["Neg Pred Value"]]
  TABLE$Balanced_accuracy <- Conf_mat_test[["ABO"]][[x]][["byClass"]][["Balanced Accuracy"]]
  TABLE$Antigen <- names(Conf_mat_test[["ABO"]][x])
  return(TABLE)
})

# Create test data other blood groups stats table
BG_stat_table <- map_dfr(2:length(Conf_mat_test), function(x) {
  TABLE <- Conf_mat_test[[x]][[2]][["byClass"]][["Sensitivity"]] %>%
    data.frame()
  colnames(TABLE) <- "Sensitivity"
  TABLE$Specificity <- Conf_mat_test[[x]][[2]][["byClass"]][["Specificity"]]
  TABLE$PosPredValue <- Conf_mat_test[[x]][[2]][["byClass"]][["Pos Pred Value"]]
  TABLE$NegPredValue <- Conf_mat_test[[x]][[2]][["byClass"]][["Neg Pred Value"]]
  TABLE$Balanced_accuracy <- Conf_mat_test[[x]][[2]][["byClass"]][["Balanced Accuracy"]]
  TABLE$Antigen <- names(Conf_mat_test[[x]][2]) %>% 
    str_replace("_plus", "+")
  return(TABLE)
})

# Combine stat tables
# Please update according your own antigen data
Test_stat_table <- bind_rows(ABO_stat_table, BG_stat_table)%>% 
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
                               grepl("Lsa", Antigen) ~ "Gerbich" )) %>% 
  relocate(BG_system) %>% relocate(Antigen, .after = BG_system)

# Write table
write.table(Test_stat_table, "results/Test_stat_table", 
            quote = F, row.names = F, col.names = T)

### Train data: All blood groups statistics

# Create train data ABO stats table
ABO_stat_table <- map_dfr(1:length(Conf_mat_train[["ABO"]]), function(x) {
  TABLE <- Conf_mat_train[["ABO"]][[x]][["byClass"]][["Sensitivity"]] %>%
    data.frame()
  colnames(TABLE) <- "Sensitivity"
  TABLE$Specificity <- Conf_mat_train[["ABO"]][[x]][["byClass"]][["Specificity"]]
  TABLE$PosPredValue <- Conf_mat_train[["ABO"]][[x]][["byClass"]][["Pos Pred Value"]]
  TABLE$NegPredValue <- Conf_mat_train[["ABO"]][[x]][["byClass"]][["Neg Pred Value"]]
  TABLE$Balanced_accuracy <- Conf_mat_train[["ABO"]][[x]][["byClass"]][["Balanced Accuracy"]]
  TABLE$Antigen <- names(Conf_mat_train[["ABO"]][x])
  return(TABLE)
})

# Create test data other blood groups stats table
BG_stat_table <- map_dfr(2:length(Conf_mat_train), function(x) {
  TABLE <- Conf_mat_train[[x]][[2]][["byClass"]][["Sensitivity"]] %>%
    data.frame()
  colnames(TABLE) <- "Sensitivity"
  TABLE$Specificity <- Conf_mat_train[[x]][[2]][["byClass"]][["Specificity"]]
  TABLE$PosPredValue <- Conf_mat_train[[x]][[2]][["byClass"]][["Pos Pred Value"]]
  TABLE$NegPredValue <- Conf_mat_train[[x]][[2]][["byClass"]][["Neg Pred Value"]]
  TABLE$Balanced_accuracy <- Conf_mat_train[[x]][[2]][["byClass"]][["Balanced Accuracy"]]
  TABLE$Antigen <- names(Conf_mat_train[[x]][2]) %>% 
    str_replace("_plus", "+")
  return(TABLE)
})

# Combine stat tables
# Please update according your own antigen data
Train_stat_table <- bind_rows(ABO_stat_table, BG_stat_table)%>% 
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
                               grepl("Lsa", Antigen) ~ "Gerbich" )) %>% 
  relocate(BG_system) %>% relocate(Antigen, .after = BG_system)

# Write table
write.table(Train_stat_table, "results/Train_stat_table", 
            quote = F, row.names = F, col.names = T)

### Full data: All blood groups statistics

# Create full data ABO stats table
ABO_stat_table <- map_dfr(1:length(Conf_mat_full[["ABO"]]), function(x) {
  TABLE <- Conf_mat_full[["ABO"]][[x]][["byClass"]][["Sensitivity"]] %>%
    data.frame()
  colnames(TABLE) <- "Sensitivity"
  TABLE$Specificity <- Conf_mat_full[["ABO"]][[x]][["byClass"]][["Specificity"]]
  TABLE$PosPredValue <- Conf_mat_full[["ABO"]][[x]][["byClass"]][["Pos Pred Value"]]
  TABLE$NegPredValue <- Conf_mat_full[["ABO"]][[x]][["byClass"]][["Neg Pred Value"]]
  TABLE$Balanced_accuracy <- Conf_mat_full[["ABO"]][[x]][["byClass"]][["Balanced Accuracy"]]
  TABLE$Antigen <- names(Conf_mat_full[["ABO"]][x])
  return(TABLE)
})

# Create full data other blood groups stats table
BG_stat_table <- map_dfr(2:length(Conf_mat_full), function(x) {
  TABLE <- Conf_mat_full[[x]][[2]][["byClass"]][["Sensitivity"]] %>%
    data.frame()
  colnames(TABLE) <- "Sensitivity"
  TABLE$Specificity <- Conf_mat_full[[x]][[2]][["byClass"]][["Specificity"]]
  TABLE$PosPredValue <- Conf_mat_full[[x]][[2]][["byClass"]][["Pos Pred Value"]]
  TABLE$NegPredValue <- Conf_mat_full[[x]][[2]][["byClass"]][["Neg Pred Value"]]
  TABLE$Balanced_accuracy <- Conf_mat_full[[x]][[2]][["byClass"]][["Balanced Accuracy"]]
  TABLE$Antigen <- names(Conf_mat_full[[x]][2]) %>% 
    str_replace("_plus", "+")
  return(TABLE)
})

# Combine stat tables
# Please update according your own antigen data
Full_stat_table <- bind_rows(ABO_stat_table, BG_stat_table)%>% 
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
                               grepl("Lsa", Antigen) ~ "Gerbich" )) %>% 
  relocate(BG_system) %>% relocate(Antigen, .after = BG_system)

# Write table
write.table(Full_stat_table, "results/Full_stat_table", 
            quote = F, row.names = F, col.names = T)

###############################################################################
### Important variables in blood group RF models may be found in 
### Full_data_results[[x]][["variable.importance"]], x = list number or name
###############################################################################