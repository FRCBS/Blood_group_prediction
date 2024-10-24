###############################################################################
### Applying DF10/R12 blood group imputation models for
### FinnGen DF10/R12 genotype data
### Kati HyvÃärinen & Jarmo Ritari
###############################################################################

library(tidyverse)
library(caret)
library(ranger)
library(e1071)
library(data.table)

###############################################################################
### Prerequisites
# 1. PLINK 1.90 https://www.cog-genomics.org/plink/
# 2. bcftools https://samtools.github.io/bcftools/bcftools.html 
# 3. Folders data and results (and src)
# 4. data/Blood_group_prediction_models_genes.rds
# 5. Chromosomal genotype vcfs (gzipped) in folder data_location/Genotype_input_data
#    example: data/Genotype_input/your_data_chr1.vcf.gz <- EDIT!
# 6. Folder data/BG_genes/

###############################################################################
### Import models and blood group information

models_info <- readRDS("data/DF10_Blood_group_prediction_models_genes.rds")

###############################################################################
# ### Format genotype data

### Object BG_genes
BG_genes <- as_tibble(read_table("data/Example_BG_genes"))
str(BG_genes)
# Structure specification
# tibble [23 Ã— 3] (S3: tbl_df/tbl/data.frame)
#$ BG_system: chr [1:23] "ABO" "CartwrightYt" "Colton" "Cromer" ...
#$ Gene1_3  : chr [1:23] "Gene1" "Gene1" "Gene1" "Gene1" ...
#$ Genes    : chr [1:23] "ABO" "ACHE" "AQP1" "CD55" ...

### Object BG_gene_rf
BG_gene_rf <- as.data.frame(read_table("data/Example_BG_gene_rf"))
str(BG_gene_rf)
# Structure specification
#'data.frame':	23 obs. of  9 variables:
#$ ensembl_gene_id  : chr  "ENSG00000128274" "ENSG00000175164" ...
#$ hgnc_symbol      : chr  "A4GALT" "ABO" "ACHE" "ACKR1" ...
#$ chromosome_name  : chr  "22" "9" "7" "1" ...
#$ start_position   : int  42692121 133233278 100889994 159203307 ...
#$ end_position     : int  42721298 133276024 100896974 159206500  ...
#$ BG_system        : chr  "P1PK" "ABO" "CartwrightYt" "Duffy" ...
#$ Gene1_3          : chr  "Gene1" "Gene1" "Gene1" "Gene1" ...
#$ rf_start_position: num  4.27e+07 1.33e+08 1.01e+08 1.59e+08 ...
#$ rf_end_position  : num  4.27e+07 1.33e+08 1.01e+08 1.59e+08 ...

###############################################################################
### Extract blood group genes (+-2K) from chromosomal vcf files

# 1. Chromosomal genotype vcfs
# 2. Output folder data/BG_genes
# 3. Object BG_gene_rf

for (i in 1:nrow(BG_gene_rf)) {
  CHR <- BG_gene_rf[i, "chromosome_name"]
  BEG <- BG_gene_rf[i, "rf_start_position"]
  END <- BG_gene_rf[i, "rf_end_position"]
  BG <- BG_gene_rf[i, "BG_system"]
  GENE <- BG_gene_rf[i, "hgnc_symbol"]
  system(command =
           paste0("bcftools view -r ", "chr", CHR, ":", BEG, "-", END,
                  " data_location/Genotype_input_data/your_data", CHR,
                  ".vcf.gz -Ov -o data/BG_genes/",
                  BG, "_", GENE, "_2K.vcf"))
}

####################################################
### Blood group gene vfc to PLINK files

# List vcf files
vcf_files <- list.files(path = "data/BG_genes", pattern = "*.vcf", 
                        full.names = TRUE)

# Convert vcf to plink binary data set
for (file in vcf_files) {
  system(command = paste0("plink --vcf ", file, 
                          " --double-id --make-bed --out ", file))
}

####################################################
### Merged binary files for multi-gene phenotypes
# Merge MNS binary files
system(command = paste0("plink --bfile data/BG_genes/MNS_GYPA_2K.vcf --bmerge data/BG_genes/MNS_GYPB_2K.vcf --make-bed --out data/BG_genes/MNS_GYPA_GYPB_2K.vcf"))
system(command = paste0("plink --bfile data/BG_genes/MNS_GYPA_GYPB_2K.vcf --bmerge data/BG_genes/MNS_GYPE_2K.vcf --make-bed --out data/BG_genes/MNS_GYPA_GYPB_GYPE_2K.vcf"))

# Merge Rh binary files
system(command = paste0("plink --bfile data/BG_genes/Rh_RHCE_2K.vcf --bmerge data/BG_genes/Rh_RHD_2K.vcf --make-bed --out data/BG_genes/Rh_RHCE_RHD_2K.vcf"))

# Merge Lewis binary files
system(command = paste0("plink --bfile data/BG_genes/Lewis_FUT2_2K.vcf --bmerge data/BG_genes/Lewis_FUT3_2K.vcf --make-bed --out data/BG_genes/Lewis_FUT2_FUT3_2K.vcf"))

# Merge P1PK binary files
system(command = paste0("plink --bfile data/BG_genes/P1PK_A4GALT_2K.vcf --bmerge data/BG_genes/P1PK_B3GALNT1_2K.vcf --make-bed --out data/BG_genes/P1PK_A4GALT_B3GALNT1_2K.vcf"))

####################################################
### Genotype dosage files

# List of plink binary data set names
plink_files <- list.files("data/BG_genes", 
                          pattern = c("*.bed"), full.names = TRUE)
plink_files <- str_replace_all(plink_files, ".bed", "") 
# Remove ".bed" in order to create binary data set name

# Remove multigenic single gene files from plink_files
plink_files <- plink_files[!grepl('MNS_GYPA_2K|MNS_GYPB_2K|MNS_GYPE_2K|MNS_GYPA_GYPB_2K|Rh_RHCE_2K|Rh_RHD_2K|Lewis_FUT2_2K|Lewis_FUT3_2K|P1PK_B3GALNT1_2K|P1PK_A4GALT_2K', 
                                  plink_files)]
# Create genotype dosage files
for (file in plink_files) {
  system(command = paste0("plink --bfile ", file, 
                          " --recodeA --out ", file, ".dosage"))
}

###############################################################################
### Harmonize dosage files

### Import PLINK genotype dosage files
# List dosage files
input.dosage.files <-  list.files("data/BG_genes",
                            pattern = c("*.raw$"),
                            full.names = TRUE)

# Read dosage files
# The allele calculated in dosage is the final letter in the variant name
input_dosage <- lapply(input.dosage.files, read_table) 

# Remove unnecessary PLINK columns and rename IID to releasedIdentifier
input_dosage <- lapply(input_dosage, select, -FID, -PAT, -MAT, -SEX, -PHENOTYPE)
input_dosage <- lapply(input_dosage, rename, releasedIdentifier = IID)

# Name elements in the dosage data list
dosage.files_split <-  gsub("data/BG_genes/", "", # Create splitted list
                            input.dosage.files) %>% 
  str_split_fixed(., "_", 2) %>% .[,2] %>% 
  gsub("_2K.vcf.dosage.raw", '', ., fixed = T)
names(input_dosage) <- dosage.files_split     

### Model variants
# Load important variables for each model
important_variables <- map(models_info[1:36], function(x) {
  return(names(x[["variable.importance"]]))
})

# Selecting phenotypes in "important_variables" matching existing dosage data
input_phenotypes <- filter(models_info[["Pheno_info"]], Pheno_genes %in%
                             names(input_dosage))$Phenotype
important_variables <- important_variables[input_phenotypes]
important_variables <- important_variables[!is.na(names(important_variables))]

### Input .bim files for harmonization
# List input .bim files
bim.files <- plink_files %>% paste0(., ".bim")

# Read bim.files
input_bim <- lapply(bim.files, read_table2, col_names = F)

# Name elements in the input_bim 
names(input_bim) <- names(input_dosage)   

##############################################
### Modify dosage file

modified_dosage <- map(names(important_variables), function(x) {
  print(x) # x = Blood group phenotype

  # Extract blood group gene names listing elements in input_dosage and input_bim
  GENE_NAME <- filter(models_info[["Pheno_info"]], Phenotype == x)[1, "Pheno_genes"] %>% 
    unlist # To vector
  # Select dosage element by gene name
  DOSAGE_DATA <- input_dosage[[GENE_NAME]]
  
  ## Important variants
  # Important variants as data.frame
  IMPORTANT_VARS <- str_split_fixed(important_variables[[x]], "_", 5) %>% 
    data.frame %>% 
    unite(., "Variant_position", X1, X2)
  IMPORTANT_VARS <- rename(IMPORTANT_VARS, 
                           Model_ref = X3, 
                           Model_alt = X4, 
                           Model_dosage_allele = X5)
  IMPORTANT_VARS$Model_variant <- important_variables[[x]]
  IMPORTANT_VARS$Variant_ID <- unite(IMPORTANT_VARS, "Variant_ID", 
                                        Variant_position, 
                                        Model_ref, Model_alt)$Variant_ID
  
  ## Input data
  # Input data dosage allele
  DOSAGE_ALLELE <- str_split(colnames(DOSAGE_DATA)[-1], "_") %>% 
    map(., function(y) {
      y[length(y)]
    }) %>% unlist
  
  # Extended input bim
  EXT_BIM <- input_bim[[GENE_NAME]] %>% rename(CHR = X1, 
                                             Variant_ID = X2,
                                             Morgans = X3,
                                             BP = X4,
                                             Input_A1 = X5,
                                             Input_A2 = X6)
  # Add input dosage allele information
  EXT_BIM$Input_dosage_allele <- DOSAGE_ALLELE
  # Modify CHR column to include chr and number, e.g. chr9
  EXT_BIM$CHR <- tolower(EXT_BIM$CHR)
  if (grepl("chr", EXT_BIM$CHR[1]) == F) {
    EXT_BIM$CHR <- paste0("chr", EXT_BIM$CHR)
  }
 
  # Create "Input_dosage_variant" ID column based on input data
  EXT_BIM$Input_dosage_variant <- unite(EXT_BIM, "Input_dosage_variant", 
                                        Variant_ID, 
                                        Input_dosage_allele)$Input_dosage_variant
  

  # Filter out non-standard variant coding
  keep.variant <- (cbind(grepl("A|T|C|G", EXT_BIM$Input_A1), 
                         grepl("A|T|C|G", EXT_BIM$Input_A2)) %>% 
                     rowSums(na.rm = T)) == 2
  EXT_BIM <- EXT_BIM[keep.variant, ]
 
  ## Missing variants 
  # Identify important variants from input data
  IMPORTANT_BIM <- left_join(IMPORTANT_VARS, EXT_BIM, 
                             by = "Variant_ID")
  
  ### Identify dosage allele discrepancy
  INPUT_DOSAGE_DISCR <- IMPORTANT_BIM$Input_dosage_variant[IMPORTANT_BIM$Model_dosage_allele != 
                                                             IMPORTANT_BIM$Input_dosage_allele]
  # Identify discrepancy alleles from dosage data and change orientation
  map(INPUT_DOSAGE_DISCR, function(y) {
    DOSAGE_DATA[,y] <<- recode(DOSAGE_DATA[,y] %>% unlist, "0" = 2, "1" = 1, "2" = 0) 
  })
  
  # Replace input dosage variant names with model variants
  map(2:ncol(DOSAGE_DATA), function(y) {
    COLUMN_NAME <- colnames(DOSAGE_DATA)[y]
    MODEL_VAR <- filter(IMPORTANT_BIM, 
                        Input_dosage_variant == COLUMN_NAME)$Model_variant
    if(length(MODEL_VAR)>0) colnames(DOSAGE_DATA)[y] <<- MODEL_VAR   
  })

  # Identify missing variants from input data
  MISSING_VAR <- setdiff(IMPORTANT_VARS$Model_variant, IMPORTANT_BIM$Model_variant)
  # Remove missing variants from input dosage data
  DOSAGE_DATA <- select(DOSAGE_DATA, releasedIdentifier,
                        IMPORTANT_BIM$Model_variant)
  # Add value 1 to missing variant dosage and change format to data.frame
  MISSING_VAR_MAT <- matrix(1, nrow = nrow(DOSAGE_DATA), 
                            ncol = length(MISSING_VAR)) %>% data.frame()
  colnames(MISSING_VAR_MAT) <- MISSING_VAR
  # Add missing variants to input dosage data
  DOSAGE_DATA <- cbind(DOSAGE_DATA, MISSING_VAR_MAT)
  
  ## Missing data
  # Calculate median dosage values based on input data and replace missing 
  # values by median value
  DOSAGE_DATA[DOSAGE_DATA == -9] <- NA 
  map(2:ncol(DOSAGE_DATA), function(y) {
    DOSAGE_DATA[,y][is.na(DOSAGE_DATA[,y])] <<- median(DOSAGE_DATA[, y], 
                                                       na.rm = T) 
  })  
  
  return(DOSAGE_DATA)
})

names(modified_dosage) <- names(important_variables)

###############################################################################
### Apply models

# Selecting Finnish models and Pheno_info in "models_info" matching existing input data
models_info <- models_info[c(names(important_variables), 'Pheno_info')]

# Filter Pheno_info to match existing input data
models_info[["Pheno_info"]] <- filter(models_info[["Pheno_info"]], Phenotype %in% 
                                        names(important_variables))

# Apply models
Prediction_results <- lapply(names(models_info)[1:(length(models_info) -1)], function(x) {
  print(x)
  out <- predict(models_info[[x]], data = modified_dosage[[x]])
  out$predictions <- 
    data.frame(releasedIdentifier = modified_dosage[[x]]$releasedIdentifier,
               out$predictions)
  return(out)
})

names(Prediction_results) <- names(models_info)[1:(length(models_info) -1)]

### Write out blood group prediction results
# Save in RDS format
saveRDS(Prediction_results, file = "results/Blood_group_prediction_results.rds")

# Combined prediction results
Prediction_tables <- lapply(1:length(Prediction_results), function(x) {
  TABLE <- Prediction_results[[x]][["predictions"]]
  return(TABLE)
})
# Add names to the Prediction_tables
names(Prediction_tables) <- names(Prediction_results)

# Combined result table with posterior probabilities
Prediction_table <- Prediction_tables %>% reduce(left_join, 
                                                 by = "releasedIdentifier")
fwrite(Prediction_table, "results/Blood_group_predictions_PP",
       sep = " ", col.names = T)

## Antigen result table
# Select antigen/phenotype positive results
Antigen_table <- Prediction_table %>% select(-matches('minus'))
# Posterior probabilities to 1 and 0
Antigen_table <- Antigen_table %>% mutate_if(is.numeric, ~if_else(. > 0.5, 1, 0))
fwrite(Antigen_table, "results/BG_prediction_antigen_pos", 
       sep = " ", col.names = T)

###############################################################################
