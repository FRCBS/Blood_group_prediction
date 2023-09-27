###############################################################################
### Genotype data preparations for prediction of blood group antigens
###############################################################################

library(tidyverse)

###############################################################################
### General information

# Files provided: "Example_BG_genes" and "Example_BG_gene_rf".

# Prerequisites
# 1. PLINK 1.90 https://www.cog-genomics.org/plink/
# 2. bcftools https://samtools.github.io/bcftools/bcftools.html 
# 3. Folders data and results (and src) 
# 4. Chromosomal genotype vcfs (gzipped) in folder data/Reference_input
#    example: data/Reference_input/your_data_chr1.vcf.gz
# 5. Folder data/BG_genes/
# 6. File "Example_BG_genes" in folder data/BG_genes/
# 7. File "Example/BG_gene_rf" in folder data/BG_genes/

### Object BG_genes
# The following script requires a tibble "BG_genes" describing blood group 
# gene information data. Please find an example "Example_BG_genes" included. 
# Import the example file and remove or add information according your own 
# blood group data.
BG_genes <- as_tibble(read_table("data/BG_genes/Example_BG_genes"))
str(BG_genes)
# Structure specification
  # tibble [23 × 3] (S3: tbl_df/tbl/data.frame)
    #$ BG_system: chr [1:23] "ABO" "CartwrightYt" "Colton" "Cromer" ...
    #$ Gene1_3  : chr [1:23] "Gene1" "Gene1" "Gene1" "Gene1" ...
    #$ Genes    : chr [1:23] "ABO" "ACHE" "AQP1" "CD55" ...
  
    # There may be more than one gene for BG_system named in column Gene1_3:
    # "Gene1", "Gene2", "Gene3"
  # Order alphabetically by BG_system

### Object BG_gene_rf
# The following script requires a data frame "BG_gene_rf" describing blood group 
# gene information data.
# Import the example file and remove or add information according your own 
# blood group data.
BG_gene_rf <- as.data.frame(read_table("data/BG_genes/Example_BG_gene_rf"))
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

# Prerequisites: 
# 1. Chromosomal genotype vcfs in folder data/Reference_input
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
                  " data/Reference_input/your_data_chr", CHR, 
                  "_.vcf.gz -Ov -o data/BG_genes/", 
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