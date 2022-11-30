region_args = commandArgs(trailingOnly=TRUE)

## v4 updates build_CS function 
## 1. use SNP-specific yty for Pr(Model)
## 2. provides two options for Pr(Model): one based on BF+gprior, another based on Wald-type+gprior

## --- Load Packages 
.libPaths("/project/dconti_624/Users/shenjiay/RPackages/")
library(data.table)
library(susieR)
library(tidyverse)
library("BayesFactor")
library("Rmpfr") 


source('/project/dconti_624/Users/shenjiay/Anqi/scripts/mJAM_fn_in_progress/mJAM_SuSiE_Forward.R')
source('/project/dconti_624/Users/shenjiay/Anqi/scripts/mJAM_fn_in_progress/mJAM_LDpruning.R')
source('/project/dconti_624/Users/shenjiay/Anqi/scripts/mJAM_fn_in_progress/mJAM_condp_v3.R')
source('/project/dconti_624/Users/shenjiay/Anqi/scripts/mJAM_fn_in_progress/mJAM_Forward.R')
muffled_cov2cor = function (x) {
  withCallingHandlers(cov2cor(x),
                      warning = function(w) {
                        if (grepl("had 0 or NA entries; non-finite result is doubtful",
                                  w$message))
                          invokeRestart("muffleWarning")
                      }) }
expand.grid.unique <- function(x, y, include.equals=TRUE){
  x <- unique(x)
  y <- unique(y)
  g <- function(i){
    z <- setdiff(y, x[seq_len(i-include.equals)])
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  do.call(rbind, lapply(seq_along(x), g))
}

#######################################################################
########       Section 1: Load data and check      
#######################################################################

# N_GWAS_EUR <- 726828
# N_GWAS_AA <- 80999
# N_GWAS_Asian <- 106599
# N_GWAS_LA <- 30336

N_cases_EUR <- 122188; N_ctrl_EUR <- 604640
N_cases_AA <- 19391; N_ctrl_AA <- 61608
N_cases_LA <- 3931; N_ctrl_LA <- 26405
N_cases_Asian <- 10809; N_ctrl_Asian <- 95790

N_GWAS_EUR <- 4/(1/N_cases_EUR + 1/N_ctrl_EUR)
N_GWAS_AA <- 4/(1/N_cases_AA + 1/N_ctrl_AA)
N_GWAS_LA <- 4/(1/N_cases_LA + 1/N_ctrl_LA)
N_GWAS_Asian <- 4/(1/N_cases_Asian + 1/N_ctrl_Asian)


## Find the MarkerNames within this region, according to CHR and POS
## MarkerName = minimac3_Name

##################
CHR_center = as.integer(region_args[1])
POS_lower = as.integer(region_args[2])
POS_upper = as.integer(region_args[3])

condp_cut <- 1e-7
within_pop_threshold <- as.numeric(region_args[4])
across_pop_threshold <- as.numeric(region_args[5])
o_pre <- as.character(region_args[6])
selected_old_names <- as.character(region_args[7:length(region_args)])
##################

CHR_center = 4
POS_lower = 110010499
POS_upper = 111610499
condp_cut <- 1e-7
within_pop_threshold <- 0.1
across_pop_threshold <- 0.1
o_pre <- "v4.1b_ME_"
selected_old_names <- c("4:110810499:G:T", "4:110848652:A:G")


if(o_pre == "v4.1a_ME_"){
  source('/project/dconti_624/Users/shenjiay/Anqi/scripts/mJAM_fn_in_progress/mJAM_build_CS_v2_med_yty_test.R')
}else if(o_pre == "v4.1b_ME_"){
  source('/project/dconti_624/Users/shenjiay/Anqi/scripts/mJAM_fn_in_progress/mJAM_build_CS_v2_med_yty_test.R')
}else{
  source('/project/dconti_624/Users/shenjiay/Anqi/scripts/mJAM_fn_in_progress/mJAM_build_CS_v2_med_yty_test.R')
}

Forward_Dir <- paste0("/project/haiman_625/Users/wanganqi/Biobank/META2021/JAM/results_final451/" , CHR_center, "_", POS_lower, "_", POS_upper, "/")
if (dir.exists(Forward_Dir)==FALSE){
  system(paste('mkdir -p ', Forward_Dir, sep=''))
}
COJO_Dir = paste0(Forward_Dir, "COJO/")
if (dir.exists(COJO_Dir)==FALSE){
  system(paste('mkdir -p ', COJO_Dir, sep=''))
}
PlinkDosage_Dir = paste0(Forward_Dir, "Plink/")
if (dir.exists(PlinkDosage_Dir)==FALSE){
  system(paste('mkdir -p ', PlinkDosage_Dir, sep=''))
}
susie_Dir = paste0(Forward_Dir, "susie/")
if (dir.exists(susie_Dir)==FALSE){
  system(paste('mkdir -p ', susie_Dir, sep=''))
}


sumstats <- read.table("/project/haiman_625/Users/wanganqi/Biobank/META2021/JAM/ME_4_study_Nov2021_full.harmonzied.txt.P_lt_1e-3", 
                       header = T) %>% 
  separate(MarkerName, into = c("CHR", "POS", "A1", "A2"), sep = ":", remove = F) %>% 
  mutate(POS = as.integer(POS)) %>% 
  dplyr::select(-c(A1, A2))


sumstats1 <- sumstats %>% 
  mutate(Effect = as.numeric(Effect), P.value = as.numeric(P.value), 
         Eur_17_study_Nov2021_1.txt.Freq1 = as.numeric(Eur_17_study_Nov2021_1.txt.Freq1), 
         AAPC_10_study_Oct2021_1.txt.Freq1 = as.numeric(AAPC_10_study_Oct2021_1.txt.Freq1), 
         Asian_5_study_Oct2021_1.txt.Freq1 = as.numeric(Asian_5_study_Oct2021_1.txt.Freq1), 
         Hispanic_5_study_Oct2021_1.txt.Freq1 = as.numeric(Hispanic_5_study_Oct2021_1.txt.Freq1),
         num_missing = is.na(Eur_17_study_Nov2021_1.txt.Freq1)+is.na(AAPC_10_study_Oct2021_1.txt.Freq1)
         +is.na(Asian_5_study_Oct2021_1.txt.Freq1)+is.na(Hispanic_5_study_Oct2021_1.txt.Freq1), 
         EUR_rare = (is.na(Eur_17_study_Nov2021_1.txt.Freq1)|Eur_17_study_Nov2021_1.txt.Freq1<0.05|Eur_17_study_Nov2021_1.txt.Freq1>0.95),
         AA_rare=  (is.na(AAPC_10_study_Oct2021_1.txt.Freq1)|AAPC_10_study_Oct2021_1.txt.Freq1<0.05|AAPC_10_study_Oct2021_1.txt.Freq1>0.95),
         LA_rare= (is.na(Hispanic_5_study_Oct2021_1.txt.Freq1)|Hispanic_5_study_Oct2021_1.txt.Freq1<0.05|Hispanic_5_study_Oct2021_1.txt.Freq1>0.95),
         ASN_rare= (is.na(Asian_5_study_Oct2021_1.txt.Freq1)|Asian_5_study_Oct2021_1.txt.Freq1<0.05|Asian_5_study_Oct2021_1.txt.Freq1>0.95), 
         rare = ((EUR_rare+AA_rare+LA_rare+ASN_rare)-num_missing)/(4-num_missing)) %>% 
  filter(CHR == CHR_center,POS>=POS_lower,POS<=POS_upper) %>% 
  rename(MarkerName_old = MarkerName) %>% 
  mutate(Allele1 = toupper(Allele1), Allele2 = toupper(Allele2)) %>%  
  unite("MarkerName", c(MarkerName_old,Allele1, Allele2), sep= ":", remove = FALSE)

sumstats2 <- sumstats1 %>% 
  filter(!is.na(P.value), MarkerName_old != "1:16936052:C:T", MarkerName != "1:152535918:C:G",
         (is.na(Eur_17_study_Nov2021_1.txt.Freq1)|(Eur_17_study_Nov2021_1.txt.Freq1>0.02&Eur_17_study_Nov2021_1.txt.Freq1<0.98)), 
         (is.na(AAPC_10_study_Oct2021_1.txt.Freq1)|(AAPC_10_study_Oct2021_1.txt.Freq1>0.02&AAPC_10_study_Oct2021_1.txt.Freq1<0.98)), 
         (is.na(Asian_5_study_Oct2021_1.txt.Freq1)|(Asian_5_study_Oct2021_1.txt.Freq1>0.02&Asian_5_study_Oct2021_1.txt.Freq1<0.98)), 
         (is.na(Hispanic_5_study_Oct2021_1.txt.Freq1)|(Hispanic_5_study_Oct2021_1.txt.Freq1>0.02&Hispanic_5_study_Oct2021_1.txt.Freq1<0.98)))


#########################
## output filtered SNPs and simplified sum stats
filtered_SNPs <- tibble(SNP = character(), reason = character())
if(nrow(sumstats2) < nrow(sumstats1)){
  filtered_SNPs <- filtered_SNPs %>% add_row(SNP = setdiff(sumstats1$MarkerName, sumstats2$MarkerName), 
                                             reason = "Freq1 exactly 0 or 1")
}
#########################

## Subset dosage files
Asian_in_filename <- "/project/haiman_625/Users/wanganqi/Biobank/META2021/JAM/data/ONCO-Asian/ONCO-Asian_ME_4_study_Nov2021_P_lt_1e-3_extract_0KB_Rsq_gteq0.0_dosage.A1dosage.metaName"
Asian_out_filename <- paste0("/project/haiman_625/Users/wanganqi/Biobank/META2021/JAM/data/ONCO-Asian/", CHR_center ,".", POS_lower, ".", POS_upper,".txt")
if(!file.exists(Asian_out_filename)){
  system(paste0("awk -F ' ' -v OFS='\t' ",
                "'{if ($1==",CHR_center, " && $3>=",POS_lower," && $3<=",POS_upper, ") print}' ",
                Asian_in_filename, "  > ", Asian_out_filename))
}

AAPC_in_filename <- "/project/haiman_625/Users/wanganqi/Biobank/META2021/JAM/data/ONCO-AAPC/ONCO-AAPC_ME_4_study_Nov2021_P_lt_1e-3_extract_0KB_Rsq_gteq0.0_dosage.valid8298.metaName"
AAPC_out_filename <- paste0("/project/haiman_625/Users/wanganqi/Biobank/META2021/JAM/data/ONCO-AAPC/", CHR_center ,".", POS_lower, ".", POS_upper, ".txt")
if(!file.exists(AAPC_out_filename)){
  system(paste0("awk -F ' ' -v OFS='\t' ",
                "'{if ($2==",CHR_center, " && $3>=",POS_lower," && $3<=",POS_upper, ") print}' ",
                AAPC_in_filename, "  > ", AAPC_out_filename))
}
# AAPC_dos <- fread(AAPC_in_filename, nrow = 10)

LAPC_in_filename <- "/project/haiman_625/Users/wanganqi/Biobank/META2021/JAM/data/ONCO-LAPC/ONCO-LAPC_ME_4_study_Nov2021_P_lt_1e-3_extract_0KB_Rsq_gteq0.0_dosage.valid2244.metaName"
LAPC_out_filename <- paste0("/project/haiman_625/Users/wanganqi/Biobank/META2021/JAM/data/ONCO-LAPC/", CHR_center ,".", POS_lower, ".", POS_upper,".txt")
if(!file.exists(LAPC_out_filename)){
  system(paste0("awk -F ' ' -v OFS='\t' ",
                "'{if ($2==",CHR_center, " && $3>=",POS_lower," && $3<=",POS_upper, ") print}' ",
                LAPC_in_filename, "  > ", LAPC_out_filename))
}
# LAPC_dos <- fread(LAPC_in_filename, nrow = 10)

EUR_in_filename <- "/project/haiman_625/Users/wanganqi/Biobank/META2021/JAM/data/ONCO-EUR/ONCO-EUR_ME_4_study_Nov2021_P_lt_1e-3_extract_0KB_Rsq_gteq0.0_dosage.onco89622.metaName"
EUR_out_filename <- paste0("/project/haiman_625/Users/wanganqi/Biobank/META2021/JAM/data/ONCO-EUR/", CHR_center ,".", POS_lower, ".", POS_upper, ".txt")
if(!file.exists(EUR_out_filename)){
  system(paste0("awk -F ' ' -v OFS='\t' ",
                "'{if ($1==",CHR_center, " && $3>=",POS_lower," && $3<=",POS_upper, ") print}' ",
                EUR_in_filename, "  > ", EUR_out_filename))
}


## Pull out summary statistics & dosage using identified MarkerNames 
Asian_dosage <- fread(Asian_out_filename)
EUR_dosage <- fread(EUR_out_filename)
AAPC_dosage <- fread(AAPC_out_filename)
LAPC_dosage <- fread(LAPC_out_filename)


#########################
## sumstat allele1 = dosage 1st column 
Asian_match <- sumstats2 %>% select(c(CHR, POS, Allele1, Allele2, MarkerName, MarkerName_old)) %>% 
  left_join(select(Asian_dosage, c(V1, V2, V3, V4, V5)), by = c("POS" = "V3")) %>% 
  mutate(noflip = (Allele1 == V4 & Allele2 == V5 & MarkerName_old == V2), 
         flip = (Allele1 == V5 & Allele2 == V4 & MarkerName_old == V2), 
         nomatch = (noflip+flip!=1)) %>% 
  filter(!nomatch)
Asian_flip_SNPs <- filter(Asian_match, flip) %>% pull(V2)
Asian_dosage1 <- Asian_dosage
Asian_dosage1[Asian_dosage$V2 %in% Asian_flip_SNPs, 6:ncol(Asian_dosage)] = 2-Asian_dosage[Asian_dosage$V2 %in% Asian_flip_SNPs, 6:ncol(Asian_dosage)]
Asian_dosage10 <- Asian_dosage1 %>% 
  mutate(MAF = rowMeans(.[,6:ncol(Asian_dosage1)])/2,
         V4_new = ifelse(V2 %in% Asian_flip_SNPs, V5, V4), 
         V5_new = ifelse(V2 %in% Asian_flip_SNPs, V4, V5)) %>% 
  unite("MarkerName", c(V2, V4_new, V5_new), sep= ":", remove = FALSE) %>% 
  select(-c(V2, V4, V5)) %>% 
  select(V1, MarkerName, V3, V4_new, V5_new, MAF,everything())
Asian_dosage11 <- Asian_dosage10 %>% 
  filter(MAF>0 & MAF<1) %>% 
  select(-MAF) 
### 
EUR_match <- sumstats2 %>% select(c(CHR, POS, Allele1, Allele2, MarkerName,MarkerName_old)) %>% 
  left_join(select(EUR_dosage, c(V1, V2, V3, V4, V5)), by = c("POS" = "V3")) %>% 
  mutate(noflip = (Allele1 == V4 & Allele2 == V5 & MarkerName_old == V2), 
         flip = (Allele1 == V5 & Allele2 == V4 & MarkerName_old == V2), 
         nomatch = (noflip+flip!=1)) %>% 
  filter(!nomatch)
EUR_flip_SNPs <- filter(EUR_match, flip) %>% pull(V2)
EUR_dosage1 <- EUR_dosage
EUR_dosage1[EUR_dosage$V2 %in% EUR_flip_SNPs, 6:ncol(EUR_dosage)] = 2-EUR_dosage[EUR_dosage$V2 %in% EUR_flip_SNPs, 6:ncol(EUR_dosage)]
EUR_dosage10 <- EUR_dosage1 %>% 
  mutate(MAF = rowMeans(.[,6:ncol(EUR_dosage1)])/2, 
         V4_new = ifelse(V2 %in% EUR_flip_SNPs, V5, V4), 
         V5_new = ifelse(V2 %in% EUR_flip_SNPs, V4, V5)) %>% 
  unite("MarkerName", c(V2, V4_new, V5_new), sep= ":", remove = FALSE) %>% 
  select(-c(V2, V4, V5)) %>% 
  select(V1, MarkerName, V3, V4_new, V5_new, everything())
EUR_dosage11 <- EUR_dosage10 %>% 
  filter(MAF>0 & MAF<1) %>% 
  select(-MAF) 
### 


AAPC_match <- sumstats2 %>% select(c(CHR, POS, Allele1, Allele2, MarkerName, MarkerName_old)) %>% 
  left_join(select(AAPC_dosage, c(V1, V2, V3, V4, V5)), by = c("POS" = "V3")) %>% 
  mutate(noflip = (Allele1 == V4 & Allele2 == V5 & MarkerName_old == V1), 
         flip = (Allele1 == V5 & Allele2 == V4 & MarkerName_old == V1), 
         nomatch = (noflip+flip!=1)) %>% 
  filter(!nomatch)
AAPC_flip_SNPs <- filter(AAPC_match, flip) %>% pull(V1)
AAPC_dosage1 <- AAPC_dosage
AAPC_dosage1[AAPC_dosage$V1 %in% AAPC_flip_SNPs, 6:ncol(AAPC_dosage)] = 2-AAPC_dosage[AAPC_dosage$V1 %in% AAPC_flip_SNPs, 6:ncol(AAPC_dosage)]
AAPC_dosage10 <- AAPC_dosage1 %>% 
  mutate(MAF = rowMeans(.[,6:ncol(AAPC_dosage1)])/2,
         V4_new = ifelse(V1 %in% AAPC_flip_SNPs, V5, V4), 
         V5_new = ifelse(V1 %in% AAPC_flip_SNPs, V4, V5)
  ) %>% 
  unite("MarkerName", c(V1, V4_new, V5_new), sep= ":", remove = FALSE) %>% 
  select(-c(V1, V4, V5)) %>% 
  select(MarkerName, V2, V3, V4_new, V5_new, everything())
AAPC_dosage11 <- AAPC_dosage10 %>% 
  filter(MAF>0 & MAF<1) %>% 
  select(-MAF) 
### 
LAPC_match <- sumstats2 %>% select(c(CHR, POS, Allele1, Allele2, MarkerName, MarkerName_old)) %>% 
  left_join(select(LAPC_dosage, c(V1, V2, V3, V4, V5)), by = c("POS" = "V3")) %>% 
  mutate(noflip = (Allele1 == V4 & Allele2 == V5& MarkerName_old == V1), 
         flip = (Allele1 == V5 & Allele2 == V4& MarkerName_old == V1), 
         nomatch = (noflip+flip!=1)) %>% 
  filter(!nomatch)
LAPC_flip_SNPs <- filter(LAPC_match, flip) %>% pull(V1)
LAPC_dosage1 <- LAPC_dosage
LAPC_dosage1[LAPC_dosage$V1 %in% LAPC_flip_SNPs, 6:ncol(LAPC_dosage)] = 2-LAPC_dosage[LAPC_dosage$V1 %in% LAPC_flip_SNPs, 6:ncol(LAPC_dosage)]
LAPC_dosage10 <- LAPC_dosage1 %>% 
  mutate(MAF = rowMeans(.[,6:ncol(LAPC_dosage1)])/2,
         V4_new = ifelse(V1 %in% LAPC_flip_SNPs, V5, V4), 
         V5_new = ifelse(V1 %in% LAPC_flip_SNPs, V4, V5)) %>% 
  unite("MarkerName", c(V1, V4_new, V5_new), sep= ":", remove = FALSE) %>% 
  select(-c(V1, V4, V5)) %>% 
  select(MarkerName, V2, V3, V4_new, V5_new, everything())
LAPC_dosage11 <- LAPC_dosage10 %>% 
  filter(MAF>0 & MAF<1) %>% 
  select(-MAF) 
#########################



## missing = either dosage missing or sumstats missing 
EUR_missing <- union(setdiff(sumstats2$MarkerName, EUR_dosage11$MarkerName), 
                     sumstats2 %>% filter(is.na(Eur_17_study_Nov2021_1.txt.Freq1)) %>% pull(MarkerName))
AA_missing <- union(setdiff(sumstats2$MarkerName, AAPC_dosage11$MarkerName), 
                    sumstats2 %>% filter(is.na(AAPC_10_study_Oct2021_1.txt.Freq1)) %>% pull(MarkerName))
LA_missing <- union(setdiff(sumstats2$MarkerName, LAPC_dosage11$MarkerName), 
                    sumstats2 %>% filter(is.na(Hispanic_5_study_Oct2021_1.txt.Freq1)) %>% pull(MarkerName))
ASN_missing <- union(setdiff(sumstats2$MarkerName, Asian_dosage11$MarkerName), 
                     sumstats2 %>% filter(is.na(Asian_5_study_Oct2021_1.txt.Freq1)) %>% pull(MarkerName))
missing_names <- Reduce(intersect, list(EUR_missing,AA_missing,LA_missing,ASN_missing))

## remove SNPs that are missing in all ethnic groups
CommonMarkerNames <- setdiff(sumstats2$MarkerName,missing_names)

forced_snps <- sumstats2 %>% filter(MarkerName_old %in% selected_old_names, 
                                    MarkerName %in% CommonMarkerNames) %>% 
  arrange(as.numeric(P.value)) %>% pull(MarkerName)

if(length(forced_snps)==0){
  write.table("index SNP(s) all missing", paste0(Forward_Dir, o_pre,"no_index.txt"))
  q()
}

#########################
## output filtered SNPs and simplified sum stats & which ethnic is missing. 
EUR_dosage_miss <- tibble(SNP = setdiff(sumstats2$MarkerName, EUR_dosage10$MarkerName), reason = "EUR dosage missing") %>% 
  add_row(SNP = setdiff(EUR_dosage10$MarkerName, EUR_dosage11$MarkerName), reason = "EUR dosage freq 0 or 1")
EUR_sumstats_miss <- tibble(SNP = sumstats2 %>% filter(is.na(Eur_17_study_Nov2021_1.txt.Freq1)) %>% pull(MarkerName), reason = "EUR sumstat missing")
AA_dosage_miss <- tibble(SNP = setdiff(sumstats2$MarkerName, AAPC_dosage10$MarkerName), reason = "AA dosage missing")%>% 
  add_row(SNP = setdiff(AAPC_dosage10$MarkerName, AAPC_dosage11$MarkerName), reason = "AA dosage freq 0 or 1")
AA_sumstats_miss <- tibble(SNP = sumstats2 %>% filter(is.na(AAPC_10_study_Oct2021_1.txt.Freq1)) %>% pull(MarkerName), reason = "AA sumstat missing")
LA_dosage_miss <- tibble(SNP = setdiff(sumstats2$MarkerName, LAPC_dosage10$MarkerName), reason = "LA dosage missing")%>% 
  add_row(SNP = setdiff(LAPC_dosage10$MarkerName, LAPC_dosage11$MarkerName), reason = "LA dosage freq 0 or 1")
LA_sumstats_miss <- tibble(SNP = sumstats2 %>% filter(is.na(Hispanic_5_study_Oct2021_1.txt.Freq1)) %>% pull(MarkerName), reason = "LA sumstat missing")
ASN_dosage_miss <- tibble(SNP = setdiff(sumstats2$MarkerName, Asian_dosage10$MarkerName), reason = "ASN dosage missing")%>% 
  add_row(SNP = setdiff(Asian_dosage10$MarkerName, Asian_dosage11$MarkerName), reason = "ASN dosage freq 0 or 1")
ASN_sumstats_miss <- tibble(SNP = sumstats2 %>% filter(is.na(Asian_5_study_Oct2021_1.txt.Freq1)) %>% pull(MarkerName), reason = "ASN sumstat missing")
ethnic_missing <- bind_rows(EUR_dosage_miss, EUR_sumstats_miss, AA_dosage_miss, AA_sumstats_miss, 
                            LA_dosage_miss, LA_sumstats_miss, ASN_dosage_miss, ASN_sumstats_miss) %>% 
  filter(SNP %in% missing_names) %>% 
  arrange(SNP)
if(length(missing_names)>0){
  filtered_SNPs <- filtered_SNPs %>% bind_rows(ethnic_missing)
}
if(nrow(filtered_SNPs)>0){
  filtered_SNPs <- filtered_SNPs %>% 
    left_join(select(sumstats1, c(MarkerName, P.value)), by = c("SNP" = "MarkerName")) %>% 
    arrange(as.numeric(P.value))
}


######### Run this for COJO analysis ######### 
## filter out any dosage missing SNPs
any_dosage_miss = ethnic_missing %>% 
  filter(reason %in% c("EUR dosage missing", "LA dosage missing", "AA dosage missing","ASN dosage missing")) %>% 
  pull(SNP) %>% unique()
CommonMarkerNames <- setdiff(CommonMarkerNames,any_dosage_miss)

# write.table(filtered_SNPs, paste0(Forward_Dir,o_pre,"filtered_snps.txt"), quote = F, row.names = F, col.names = T, sep = '\t')
#########################



## Filter sumstats and dosage based on CommonMarkerNames
sumstats3 <- sumstats2[match(CommonMarkerNames, sumstats2$MarkerName),]
Asian_dosage2 <- Asian_dosage11[match(CommonMarkerNames, Asian_dosage11$MarkerName),]
EUR_dosage2 <- EUR_dosage11[match(CommonMarkerNames, EUR_dosage11$MarkerName),]
AAPC_dosage2 <- AAPC_dosage11[match(CommonMarkerNames, AAPC_dosage11$MarkerName),]
LAPC_dosage2 <- LAPC_dosage11[match(CommonMarkerNames, LAPC_dosage11$MarkerName),]

dim(sumstats3);dim(Asian_dosage2);dim(AAPC_dosage2); dim(EUR_dosage2); dim(LAPC_dosage2)

## Transpose dosage 
EUR_dosage3 <- t(EUR_dosage2[, 6:ncol(EUR_dosage2)]); colnames(EUR_dosage3) <- CommonMarkerNames
Asian_dosage3 <- t(Asian_dosage2[, 6:ncol(Asian_dosage2)]); colnames(Asian_dosage3) <- CommonMarkerNames
AAPC_dosage3 <- t(AAPC_dosage2[, 6:ncol(AAPC_dosage2)]); colnames(AAPC_dosage3) <- CommonMarkerNames
LAPC_dosage3 <- t(LAPC_dosage2[, 6:ncol(LAPC_dosage2)]); colnames(LAPC_dosage3) <- CommonMarkerNames


## Check MAF of summary statistics vs MAF of dosage 
Meta_EUR_MAF <- as.numeric(sumstats3$Eur_17_study_Nov2021_1.txt.Freq1)
Dosage_EUR_MAF <- colMeans(EUR_dosage3)/2
flip_snps <- which( (Dosage_EUR_MAF>0.5 & Meta_EUR_MAF<0.5) | (Dosage_EUR_MAF<0.5 & Meta_EUR_MAF>0.5) )
EUR_dosage3_flipped <- EUR_dosage3 

Meta_Asian_MAF <- as.numeric(sumstats3$Asian_5_study_Oct2021_1.txt.Freq1)
Dosage_Asian_MAF <- colMeans(Asian_dosage3)/2
flip_snps <- which( (Dosage_Asian_MAF>0.5 & Meta_Asian_MAF<0.5) | (Dosage_Asian_MAF<0.5 & Meta_Asian_MAF>0.5) )
Asian_dosage3_flipped <- Asian_dosage3 

Meta_LAPC_MAF <- as.numeric(sumstats3$Hispanic_5_study_Oct2021_1.txt.Freq1)
Dosage_LAPC_MAF <- colMeans(LAPC_dosage3)/2
flip_snps <- which( (Dosage_LAPC_MAF>0.5 & Meta_LAPC_MAF<0.5) | (Dosage_LAPC_MAF<0.5 & Meta_LAPC_MAF>0.5) )
LAPC_dosage3_flipped <- LAPC_dosage3 

Meta_AAPC_MAF <- as.numeric(sumstats3$AAPC_10_study_Oct2021_1.txt.Freq1)
Dosage_AAPC_MAF <- colMeans(AAPC_dosage3)/2
flip_snps <- which( (Dosage_AAPC_MAF>0.5 & Meta_AAPC_MAF<0.5) | (Dosage_AAPC_MAF<0.5 & Meta_AAPC_MAF>0.5) )
AAPC_dosage3_flipped <- AAPC_dosage3 



#######################################################################
########       Section 2: Run mJAM - SandP
#######################################################################

## --- Prepare input for mJAM index SNP selection 
Input_Dosage <- vector(mode = "list", length = 4)
N_GWAS <- c(N_GWAS_EUR, N_GWAS_AA, N_GWAS_LA, N_GWAS_Asian)

## EUR
subset_EUR <- CommonMarkerNames
chr1_EUR_sub <- filter(sumstats3, MarkerName %in% subset_EUR)
beta_EUR <- as.numeric(chr1_EUR_sub$Eur_17_study_Nov2021_1.txt.Effect)
names(beta_EUR) <- subset_EUR
sebeta_EUR <- as.numeric(chr1_EUR_sub$Eur_17_study_Nov2021_1.txt.StdErr)
names(sebeta_EUR) <- subset_EUR
maf_EUR <- as.numeric(chr1_EUR_sub$Eur_17_study_Nov2021_1.txt.Freq1)
names(maf_EUR) <- subset_EUR
Input_Dosage[[1]] <- EUR_dosage3_flipped[,subset_EUR]


## AA
beta_AA <- as.numeric(chr1_EUR_sub$AAPC_10_study_Oct2021_1.txt.Effect)
names(beta_AA) <- subset_EUR
sebeta_AA <- as.numeric(chr1_EUR_sub$AAPC_10_study_Oct2021_1.txt.StdErr)
names(sebeta_AA) <- subset_EUR
maf_AA <- as.numeric(chr1_EUR_sub$AAPC_10_study_Oct2021_1.txt.Freq1)
names(maf_AA) <- subset_EUR
Input_Dosage[[2]] <- AAPC_dosage3_flipped[,subset_EUR]



## LA
beta_LA <- as.numeric(chr1_EUR_sub$Hispanic_5_study_Oct2021_1.txt.Effect)
names(beta_LA) <- subset_EUR
sebeta_LA <- as.numeric(chr1_EUR_sub$Hispanic_5_study_Oct2021_1.txt.StdErr)
names(sebeta_LA) <- subset_EUR
maf_LA <- as.numeric(chr1_EUR_sub$Hispanic_5_study_Oct2021_1.txt.Freq1)
names(maf_LA) <- subset_EUR
Input_Dosage[[3]] <- LAPC_dosage3_flipped[,subset_EUR]


## Asian
beta_Asian <- as.numeric(chr1_EUR_sub$Asian_5_study_Oct2021_1.txt.Effect)
names(beta_Asian) <- subset_EUR
sebeta_Asian <- as.numeric(chr1_EUR_sub$Asian_5_study_Oct2021_1.txt.StdErr)
names(sebeta_Asian) <- subset_EUR
maf_Asian <- as.numeric(chr1_EUR_sub$Asian_5_study_Oct2021_1.txt.Freq1)
names(maf_Asian) <- subset_EUR
Input_Dosage[[4]] <- Asian_dosage3_flipped[,subset_EUR]


################################################################################
Marg_Result <- data.frame(SNP = subset_EUR, 
                          beta_EUR = beta_EUR,sebeta_EUR = sebeta_EUR,
                          beta_AA = beta_AA, sebeta_AA = sebeta_AA, 
                          beta_LA = beta_LA, sebeta_LA = sebeta_LA,
                          beta_Asian = beta_Asian, sebeta_Asian = sebeta_Asian)
MAF_Result <- data.frame(SNP = subset_EUR,
                         MAF_EUR = maf_EUR,
                         MAF_AA = maf_AA,
                         MAF_LA = maf_LA,
                         MAF_Asian = maf_Asian)

N_SNP <- numSNPs <- nrow(Marg_Result)


Original_Input_Dosage <- Input_Dosage
Input_MarglogOR <- Input_MargSElogOR <- Input_MAF <- vector("list", length(Input_Dosage))
for(pop in 1:length(Input_Dosage)){
  Input_MarglogOR[[pop]] <- Marg_Result[,2*pop]
  Input_MargSElogOR[[pop]] <- Marg_Result[,2*pop+1]
  Input_MAF[[pop]] <- MAF_Result[,1+pop]
  names(Input_MarglogOR[[pop]]) <- Marg_Result[,1]
  names(Input_MargSElogOR[[pop]]) <- Marg_Result[,1]
  names(Input_MAF[[pop]]) <- MAF_Result[,1]
}


####################################
# COJO
####################################


PopName = paste0(CHR_center, "_", POS_lower, "_", POS_upper)
# Output files preparation
MAP_FileName = paste0(PopName,".map")
PED_FileName = paste0(PopName,".ped")
MA_FileName = paste0(PopName,".ma")
N_SNP = length(CommonMarkerNames)
# Construct design matrix for MAP files, all .map files are the same: dimension N_SNP
MAP_Chr = rep(CHR_center,N_SNP)
MAP_RS = colnames(Input_Dosage[[1]])
MAP_Dist = rep(0,N_SNP)
MAP_Pos = sumstats3$POS
MAP_Design = cbind(MAP_Chr,MAP_RS,MAP_Dist,MAP_Pos)
# Construct design matrix for MA files: dimension N_SNP
MA_SNP = CommonMarkerNames
# A1 and A2
Default_A1 = "A"
Default_A2 = "C"
MA_A1 = rep(Default_A1, N_SNP)
MA_A2 = rep(Default_A2, N_SNP)
MA_Design = cbind(MA_SNP,MA_A1,MA_A2)

# Initialize the output for the PED file for all populations
PED_AllPop = NULL
N_Sample = sapply(Input_Dosage, nrow)
for (e in 1:length(Input_Dosage)){
  # Construct design matrix for PED files: dimension N_Sample
  PED_FamilyID = c(1:N_Sample[e])
  PED_SampleID = rep(1,N_Sample[e])
  PED_PaternalID = rep(0,N_Sample[e])
  PED_MaternalID = rep(0,N_Sample[e])
  # Sex: 1==male 2==female other==unknown
  PED_Sex = rep(1,N_Sample[e])
  PED_Design = cbind(PED_FamilyID,PED_SampleID,PED_PaternalID,PED_MaternalID,PED_Sex)
  # Construct PED file
  # Dosage matrix to A and C
  temp_Dosage = Input_Dosage[[e]]
  temp_1A = ifelse(round(temp_Dosage) >= 1,Default_A1,Default_A2)
  temp_2A = ifelse(round(temp_Dosage) >= 2,Default_A1,Default_A2)
  # Combine the above two matrix in a way that fits the format requirement
  temp_PED = matrix(rbind(temp_1A,temp_2A),nrow=N_Sample[e],ncol=2*N_SNP)
  # Column bind the PED_Design, Y and temp_PED to output
  PED_Affection = rep(0, N_Sample[e])
  temp_PED = cbind(PED_Design,PED_Affection,temp_PED)
  # Append temp_PED to PED_AllPop
  PED_AllPop = rbind(PED_AllPop,temp_PED)
}
# Output MAP file for all population
write.table(MAP_Design,paste(PlinkDosage_Dir,MAP_FileName,sep=""),quote=F,sep="\t",row.name=F,col.name=F,append=F)
# For PED_AllPop, adjust the Family_ID, the first column
PED_AllPop[,1] = c(1:sum(N_Sample))
# Output PED file for all populations
write.table(PED_AllPop,paste(PlinkDosage_Dir,PED_FileName,sep=""),quote=F,sep="\t",row.name=F,col.name=F,append=F)


# For MA_AllPop, obtain the required info from metasoft results
# Get the mean Freq from all populations
# MA_AllPop_Freq = rep(0,N_SNP)
# for (e in 1:length(N_Sample)){
#   MA_AllPop_Freq = MA_AllPop_Freq + Input_MAF[[e]]*N_Sample[e]
# }
# MA_AllPop_Freq = MA_AllPop_Freq/sum(N_Sample)
MA_AllPop_Freq = sumstats3$Freq1
  
# BETA_FE is the 4th column in MetaSoft output
MA_AllPop_Beta =  sumstats3$Effect
# STD_FE is the 5th column in MetaSoft output
MA_AllPop_SE = sumstats3$StdErr
# PVALUE_FE is the 3rd column in MetaSoft output
MA_AllPop_P = sumstats3$P.value
# Total number of individuals in all populations
MA_AllPop_N = rep(round(sum(N_GWAS)),N_SNP)
# Combine the results for MA_AllPop
MA_AllPop = cbind(MA_Design,MA_AllPop_Freq,MA_AllPop_Beta,MA_AllPop_SE,MA_AllPop_P,MA_AllPop_N)
colnames(MA_AllPop) = c("SNP","A1","A2","freq","b","se","p","N")
# Output MA file for all population
write.table(MA_AllPop,paste0(COJO_Dir,MA_FileName),
            quote=F,sep="\t",row.name=F,col.name=T,append=F)


# Write local plink command to one file for all replications in current folder
PlinkFolder = "/project/dconti_624/Users/shenjiay/Softwares/plink_linux/"
Plink_Command = paste("./plink --file ",PlinkDosage_Dir,PopName," --make-bed --out ",COJO_Dir,PopName,sep="")

# Run plink commands
system(paste("cd", PlinkFolder, "&& chmod +x plink &&",Plink_Command))

# Write HPC COJO command to one file for all replications in current folder
GCTAFolder = "/project/dconti_624/Users/shenjiay/Softwares/gcta_1.93.2beta/"
BF_Alpha = condp_cut
CojoInput_Files = PopName
CojoOutput_Files = "AllPopRes"

# Output Cojo commands
# Bonferroni alpha
Cojo_Command = paste("./gcta64  --bfile ",COJO_Dir,CojoInput_Files,
                     " --cojo-file ",COJO_Dir,CojoInput_Files,
                     ".ma --cojo-slct --cojo-p ",BF_Alpha," --out ",COJO_Dir,CojoOutput_Files,sep="")
system(paste("cd", GCTAFolder, "&& chmod +x gcta64 &&",Cojo_Command))



# ===========================================================================
# Analyse COJO results 
COJO_res <- read.table(paste0(COJO_Dir, CojoOutput_Files, ".jma.cojo"), header = T, sep = '\t')


sumstats3 %>% 
  ggplot(aes(x = POS, y = -log10(P.value)))+
  geom_point(color = "grey")+
  geom_point(data = sumstats3[sumstats3$MarkerName %in% COJO_res$SNP,], 
             color = 'black')+
  geom_point(data = sumstats3[sumstats3$MarkerName %in% COJO_res$SNP,], 
             colour="red", shape=1, size=2.5, stroke=1.5)+
  geom_text(
    data = sumstats3[match(COJO_res$SNP, CommonMarkerNames),],
    label= COJO_res$SNP,
    nudge_x = 0.25, nudge_y = 0.3, check_overlap = T, size = 2.5
  )

ggsave(paste0(COJO_Dir,o_pre,"MHT.png"), width = 9, height = 5)




###########################################################################
##############      susie       ################
###########################################################################
source('/project/dconti_624/Users/shenjiay/Anqi/scripts/mJAM_fn_in_progress/mJAM_SuSiE_Forward.R')
Input_p_cases = c(N_cases_EUR/(N_cases_EUR+N_ctrl_EUR), 
                  N_cases_AA/(N_cases_AA+N_ctrl_AA), 
                  N_cases_LA/(N_cases_LA+N_ctrl_LA), 
                  N_cases_Asian/(N_cases_Asian+N_ctrl_Asian) )

mJAM.susie.results <- mJAM_SuSiE(marginal.betas = Input_MarglogOR,
                                 marginal.se = Input_MargSElogOR,
                                 EAFs = Input_MAF,
                                 N_GWAS = N_GWAS,
                                 X_ref = Input_Dosage,
                                 SNP_names = CommonMarkerNames,
                                 logORToBeta = FALSE,
                                 p_cases = Input_p_cases,
                                 SuSiE_num_comp = 10,
                                 SuSiE_coverage = 0.95,
                                 SuSiE_min_abs_corr = sqrt(across_pop_threshold),
                                 max_iter = 500)

# Output the posterior probabilities
mJAM_SuSiE_Output <- mJAM.susie.results$summary
mJAM_SuSiE_Post_FileName = paste("mJAM_SuSiE_Post_",".txt",sep="")
write.table(mJAM_SuSiE_Output,paste(susie_Dir,mJAM_SuSiE_Post_FileName,sep=""),
            quote=F,sep="\t",row.name=T,col.name=T,append=F)

# Output credible sets 0.95 and 0.8
if(is.null(mJAM.susie.results$fit$sets$cs)){
  mJAM_SuSiE_CS_0.95_Output <- mJAM_SuSiE_get_cs(mJAM.susie.results$fit,0.95)
}else{
  mJAM_SuSiE_CS_0.95_Output <- cbind(mJAM_SuSiE_get_cs(mJAM.susie.results$fit,0.95), 
                                     mJAM.susie.results$fit$sets$purity)
}
mJAM_SuSiE_CS_0.95_FileName = paste("mJAM_SuSiE_CS_0.95_",".txt",sep="")
write.table(mJAM_SuSiE_CS_0.95_Output,
            paste(susie_Dir,mJAM_SuSiE_CS_0.95_FileName,sep=""),
            quote=F,sep="\t",row.name=F,col.name=T,append=F)


## plot susie results
temp_max_cs_size <- max(mJAM_SuSiE_CS_0.95_Output$cs.size)
temp_susie_pip <- mJAM_SuSiE_Output %>% mutate(cs.element = 1:length(CommonMarkerNames))

mJAM_SuSiE_CS_0.95_Output2 = mJAM_SuSiE_CS_0.95_Output %>% 
  separate(credible.set, 
           into = paste0("cs.element", 1:temp_max_cs_size), 
           sep = ",") %>% 
  pivot_longer(cols = 3:(temp_max_cs_size+2),
               values_to = "cs.element") %>% 
  filter(!is.na(cs.element)) %>% select(-name) %>% mutate(cs.element = as.integer(cs.element)) %>% 
  left_join(temp_susie_pip, by = "cs.element") %>%
  group_by(index) %>%
  mutate(max_pip = (pip == max(pip)),
         index_SNP = ifelse(max_pip,cs.element, NA)) 

plot_df = sumstats3 %>% 
  mutate(cs.element = 1:length(CommonMarkerNames)) %>% 
  left_join(mJAM_SuSiE_CS_0.95_Output2, by = "cs.element") 

plot_df %>% 
  ggplot(aes(x = POS, y = -log10(as.numeric(P.value))))+
  geom_point(color ="grey")+
  geom_point(aes(color = index, group = index), data = filter(plot_df, !is.na(index)))+
  geom_point(data = filter(plot_df, !is.na(index_SNP)),
             colour="red", shape=1, size=2.5, stroke=1.5)+
  geom_text(
    data = filter(plot_df, !is.na(index_SNP)),
    label= filter(plot_df, !is.na(index_SNP)) %>% pull(MarkerName_old),
    nudge_x = 0.25, nudge_y = 0.3, check_overlap = T, size = 2.5
  )+
  geom_hline(yintercept = -log10(5e-8), linetype="dashed", color = "red")+
  labs(x = "Position", y = "-log10(meta.p.value)", title = paste0(CHR_center, ":",POS_lower, "-", POS_upper))

ggsave(paste0(susie_Dir,o_pre,"plot_cs.png"), width = 9, height = 5)
  




