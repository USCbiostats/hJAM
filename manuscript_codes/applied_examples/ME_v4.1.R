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

mJAM_LDpruning = function(target, testing, R, within_thre = 0.95, across_thre = 0.80){
  
  numEthnic <- length(R)
  high_corr_within <- high_corr_across <- vector(mode = "list", length = numEthnic)
  
  for(i in 1:numEthnic){
    high_corr_within[[i]] <- which(abs(R[[i]][target,]) > within_thre)
    high_corr_within[[i]] <- intersect(testing, high_corr_within[[i]])
    high_corr_across[[i]] <- which(abs(R[[i]][target,]) > across_thre)
    high_corr_across[[i]] <- intersect(testing, high_corr_across[[i]])
  }
  
  remove_within <- Reduce(union, high_corr_within)
  remove_across <- Reduce(intersect, high_corr_across)
  
  return(list(remove_within = remove_within, remove_across = remove_across))
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
# CHR_center = as.integer(region_args[1])
# POS_lower = as.integer(region_args[2])
# POS_upper = as.integer(region_args[3])
# 
# condp_cut <- 1e-7
# within_pop_threshold <- as.numeric(region_args[4])
# across_pop_threshold <- as.numeric(region_args[5])
# o_pre <- as.character(region_args[6])
# selected_old_names <- as.character(region_args[7:length(region_args)])
##################

CHR_center = 6
POS_lower = 116400434
POS_upper = 118000434
condp_cut <- 1e-7
within_pop_threshold <- 0.1
across_pop_threshold <- 0.1
o_pre <- "v4.1b_ME_"
selected_old_names <- c("6:117208509:G:T")

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
write.table(filtered_SNPs, paste0(Forward_Dir,o_pre,"filtered_snps.txt"), quote = F, row.names = F, col.names = T, sep = '\t')
#########################



######### Run this for COJO analysis ######### 
## filter out any dosage missing SNPs
any_dosage_miss = bind_rows(EUR_dosage_miss, EUR_sumstats_miss, AA_dosage_miss, AA_sumstats_miss, 
                            LA_dosage_miss, LA_sumstats_miss, ASN_dosage_miss, ASN_sumstats_miss) %>% 
  filter(reason %in% c("EUR dosage missing", "LA dosage missing", "AA dosage missing","ASN dosage missing")) %>% 
  pull(SNP) %>% unique()
any_dosage_miss
CommonMarkerNames <- setdiff(CommonMarkerNames,any_dosage_miss)
forced_snps

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

## --- Check missing data
numEthnic <- length(Input_MarglogOR)
missing_gwas_tbl <- tibble(missing_ethnic_idx = integer(), missing_snp_idx = integer())
missing_dosage_tbl <- tibble(missing_ethnic_idx = integer(), missing_snp_idx = integer())

for(i in 1:numEthnic){
  temp_missing_gwas <- union(which(is.na(Input_MarglogOR[[i]])), which(is.na(Input_MargSElogOR[[i]])))
  temp_missing_dosage <- which(colSums(is.na(Input_Dosage[[i]]))>0)
  if(length(temp_missing_gwas)>0){
    missing_gwas_tbl <- missing_gwas_tbl %>% add_row(missing_ethnic_idx =i, missing_snp_idx = temp_missing_gwas)
  }
  if(length(temp_missing_dosage)>0){
    missing_dosage_tbl <- missing_dosage_tbl %>% add_row(missing_ethnic_idx =i, missing_snp_idx = temp_missing_dosage)
  }
}
missing_tbl <- bind_rows(missing_dosage_tbl, missing_gwas_tbl) %>% distinct()

dim(missing_dosage_tbl)
dim(missing_gwas_tbl)
dim(missing_tbl)
missing_tbl %>% group_by(missing_snp_idx) %>% summarize(n_miss=n_distinct(missing_ethnic_idx)) %>% arrange(desc(n_miss))


## --- Calculate the LD matrix of Input_Dosage
Dosage_cor <- vector("list",length(Input_Dosage))
for(i in 1:length(Input_Dosage)){
  if (nrow(missing_dosage_tbl) > 0 && i%in%missing_dosage_tbl$missing_ethnic_idx){
    ## --- Get missing SNP index
    temp_missing_snp_idx <- filter(missing_dosage_tbl, missing_ethnic_idx == i) %>% pull(missing_snp_idx)
    ## --- Get Dosage_cor with complete SNP
    temp_Dosage_cor <- cor(Input_Dosage[[i]][,-temp_missing_snp_idx])^2
    ## --- Fill in missing SNPs with zeros
    Dosage_cor[[i]] <- diag(1, nrow = numSNPs, ncol = numSNPs)
    Dosage_cor[[i]][-temp_missing_snp_idx, -temp_missing_snp_idx] <- temp_Dosage_cor
  }else{
    Dosage_cor[[i]] <- cor(Input_Dosage[[i]])^2
  }
  colnames(Dosage_cor[[i]]) <- rownames(Dosage_cor[[i]]) <- Marg_Result[,1]
}
cor(Input_Dosage[[1]][,forced_snps])^2
cor(Input_Dosage[[2]][,forced_snps])^2
cor(Input_Dosage[[3]][,forced_snps])^2
cor(Input_Dosage[[4]][,forced_snps])^2 

## -- Calculate sufficient statistics 
GItGI <- GIty <- yty <- yty_med <- vector("list", numEthnic)
for (i in 1:numEthnic){
  if (nrow(missing_tbl) > 0 && i%in%missing_tbl$missing_ethnic_idx){
    ## --- Get missing SNP index
    temp_missing_snp_idx <- filter(missing_tbl, missing_ethnic_idx == i) %>% pull(missing_snp_idx)
    ## --- Get GItGI, GIty, yty with complete SNP
    temp_X_ref <- Input_Dosage[[i]][,-temp_missing_snp_idx]
    temp_MAFs <- Input_MAF[[i]][-temp_missing_snp_idx]
    temp_GItGI <- get_XtX(N_outcome = N_GWAS[i], Gl = temp_X_ref, maf = temp_MAFs)
    ##
    temp.marginal.betas <- Input_MarglogOR[[i]][-temp_missing_snp_idx]
    temp.marginal.se <- Input_MargSElogOR[[i]][-temp_missing_snp_idx]
    temp.GIty <- get_z(maf = temp_MAFs, betas = temp.marginal.betas, N_outcome = N_GWAS[i])
    ##
    # yty[[i]] <- get_yty(maf = temp_MAFs, N_outcome = N_GWAS[i], betas = temp.marginal.betas, betas.se = temp.marginal.se)
    ## --- Fill in missing SNPs with zeros
    GItGI[[i]] <- matrix(0, nrow = numSNPs, ncol = numSNPs)
    GItGI[[i]][-temp_missing_snp_idx, -temp_missing_snp_idx] <- temp_GItGI
    GIty[[i]] <- rep(0, numSNPs)
    GIty[[i]][-temp_missing_snp_idx] <- temp.GIty
  }else{
    GItGI[[i]] <- get_XtX(N_outcome = N_GWAS[i], Gl = Input_Dosage[[i]], 
                          maf = Input_MAF[[i]])
    GIty[[i]] <- get_z(maf = Input_MAF[[i]], 
                       betas = Input_MarglogOR[[i]], N_outcome = N_GWAS[i])
  }
  temp_yty <- get_yty(maf = Input_MAF[[i]], N_outcome = N_GWAS[i], 
                      betas = Input_MarglogOR[[i]], 
                      betas.se = Input_MargSElogOR[[i]])
  yty[[i]] <- temp_yty$yTy.all
  yty_med[[i]] <- temp_yty$yTy.est
  colnames(GItGI[[i]]) <- rownames(GItGI[[i]]) <- Marg_Result[,1] 
  names(GIty[[i]]) <- names(yty[[i]]) <- Marg_Result[,1] 
  print(i)
}
rare_SNPs <- chr1_EUR_sub %>% filter(rare >= 0.5) %>% pull(MarkerName)

## Run Forward selection
iter_count <- 0
selected_ids <- NULL
prev_index_snp <- NULL
pruned_snps <- NULL
condp_list <- NULL
curr_min_condp <- 0
prev_min_condp <- 0
all_CS <- all_CS2 <- all_99CS <- NULL
subset_EUR <- has_dosage_SNP <- Marg_Result$SNP
pruned_snp_df <- tibble(pruned_SNP = character(), index_SNP = character(), reason = character())
chr1_EUR_sub <- filter(sumstats3, MarkerName %in% subset_EUR)
GItGI_curr <- GItGI
GIty_curr <- GIty
yty_curr <- yty

while(iter_count >= 0){
  ## step 1: select top SNPs in the remaining list
  ## selected_id should be the ID in Input_MarglogOR
  if(length(unique(pruned_snps))==N_SNP){break}
  if(iter_count == 0){Input_id = NULL}else{Input_id = match(selected_ids, subset_EUR)}
  
  if(sum(rare_SNPs %in% subset_EUR)>0){
    rare_id = match(rare_SNPs, subset_EUR)
    rare_id = rare_id[!is.na(rare_id)]
  }else{
    rare_id = NULL
  }
  ## for round 1, use the top SNP of meta-analysis
  newFS_RES <- mJAM_get_condp(GItGI = GItGI_curr, GIty = GIty_curr, yty = yty_curr, 
                              yty_med = yty_med,
                              N_GWAS = N_GWAS, selected_id = Input_id,
                              rare_id = rare_id)
  
  if(iter_count ==0){
    curr_index_snp <- forced_snps[1]
    selected_ids <- c(selected_ids, curr_index_snp)
    condp_list <- c(condp_list, sumstats3 %>% filter(MarkerName == curr_index_snp) %>% pull(P.value))
    ## output mJAM marginal p and meta marginal p 
    marginal_est <- tibble(SNP = subset_EUR, mJAM_effect = signif(newFS_RES$effect_est, digits = 3), 
                           meta_effect = chr1_EUR_sub$Effect,
                           mJAM_se = signif(newFS_RES$se_est, digits = 3),meta_se = chr1_EUR_sub$StdErr, 
                           mJAM_p = signif(newFS_RES$condp, digits = 3), meta_p = chr1_EUR_sub$P.value) %>% 
      arrange(meta_p, mJAM_p)
    write.table(marginal_est, paste0(Forward_Dir,o_pre, "mJAM_vs_meta_est.txt"), quote = F, row.names = F, col.names = T, sep = '\t')
    troubled_est <- filter(marginal_est, abs(-log10(mJAM_p)+log10(meta_p))>2)
    med_probs_old <- select(marginal_est, c(SNP,mJAM_p,mJAM_effect)) %>% 
      rename(prev_cond_p = mJAM_p, prev_cond_effect = mJAM_effect)
    # write.table(troubled_est, paste0(Forward_Dir,o_pre, "troubled_est.txt"), quote = F, row.names = F, col.names = T, sep = '\t')
    ##
    p_pval <- tibble(mJAM_pval = marginal_est$mJAM_p, meta_pval = marginal_est$meta_p) %>%
      ggplot(aes(x = -log10(mJAM_pval), y = -log10(meta_pval))) +
      geom_point(color="grey")+
      scale_x_continuous(limits = c(0, max(c(-log10(marginal_est$mJAM_p), -log10(marginal_est$meta_p)))*1.01))+
      scale_y_continuous(limits = c(0, max(c(-log10(marginal_est$mJAM_p), -log10(marginal_est$meta_p)))*1.01))+
      geom_abline()+labs(title = "pval")
    p_b <- tibble(mJAM_b = marginal_est$mJAM_effect, meta_b = marginal_est$meta_effect) %>%
      ggplot(aes(x = mJAM_b, y = meta_b)) +
      geom_point(color="grey")+
      scale_x_continuous(limits = c(0, max(c(marginal_est$mJAM_effect, marginal_est$meta_effect))*1.01))+
      scale_y_continuous(limits = c(0, max(c(marginal_est$mJAM_effect, marginal_est$meta_effect))*1.01))+
      geom_abline()+labs(title = "beta")
    p_se <- tibble(mJAM_se = marginal_est$mJAM_se, meta_se = marginal_est$meta_se) %>%
      ggplot(aes(x = mJAM_se, y = meta_se)) +
      geom_point(color="grey")+
      scale_x_continuous(limits = c(0, max(c(marginal_est$mJAM_se, marginal_est$meta_se))*1.01))+
      scale_y_continuous(limits = c(0, max(c(marginal_est$mJAM_se, marginal_est$meta_se))*1.01))+
      geom_abline()+labs(title = "se")
    p_diag <- gridExtra::grid.arrange(p_pval, p_b, p_se, nrow = 1)
    ggsave(paste0(Forward_Dir,o_pre,"mJAM_vs_meta_diag.png"), plot = p_diag, width = 12, height =5)
  }else{
    if(iter_count>=length(forced_snps)) {break}
    ## output mJAM marginal p and meta marginal p 
    chr1_EUR_sub <- chr1_EUR_sub[match(subset_EUR,chr1_EUR_sub$MarkerName),]
    condp_est <- tibble(SNP = subset_EUR, 
                        mJAM_effect = signif(newFS_RES$effect_est, digits = 3), 
                        meta_effect = chr1_EUR_sub$Effect,
                        mJAM_se = signif(newFS_RES$se_est, digits = 3),meta_se = chr1_EUR_sub$StdErr, 
                        mJAM_p = signif(newFS_RES$condp, digits = 3), meta_p = chr1_EUR_sub$P.value,
                        selected_SNPs = paste(selected_ids, collapse = ",")) %>% 
      arrange(mJAM_p, meta_p)
    curr_index_snp <- forced_snps[iter_count+1]
    selected_ids <- c(selected_ids, curr_index_snp)
    condp_list <- c(condp_list, newFS_RES$condp[match(curr_index_snp, subset_EUR)])
    # write.table(condp_est, paste0(Forward_Dir,o_pre,"index", iter_count+1,"_mJAM_vs_meta_est.txt"), quote = F, row.names = F, col.names = T, sep = '\t')
  }
  message(paste("No.", iter_count+1,"selected index SNP:", curr_index_snp))
  
  ## step 2: construct credible sets based on the selected SNP
  ###############################################################################################
  All_id = subset_EUR
  X_id = curr_index_snp # X_id = C_id = "2:242480862:C:T:T:C"
  prev_X_list = prev_index_snp # prev_X_list = "2:242157241:A:G:A:G"
  PrCS_weights = "Pr(M_C)"
  coverage = 0.95
  Pr_Med_cut = 0.1
  ## --- Get "sum" of posterior probability of mediation test
  Testing_id <- All_id[!All_id %in% c(X_id, prev_X_list)]
  Post_Med_Prob <-  rep(NA, length = length(All_id)); names(Post_Med_Prob) <- All_id
  Med_Effect_Size <- rep(NA, length(All_id)); names(Med_Effect_Size) <- All_id
  Post_Model_Prob <- rep(NA, length = length(All_id)); names(Post_Model_Prob) <- All_id
  Post_Model_Prob_Ratio <- rep(NA, length = length(All_id)); names(Post_Model_Prob_Ratio) <- All_id
  Post_CS_Prob <- rep(NA, length = length(All_id)); names(Post_CS_Prob) <- All_id
  R2_est <- rep(NA, length = length(All_id)); names(R2_est) <- All_id
  n_miss <- rep(NA, length = length(All_id)); names(n_miss) <- All_id
  
  ## --- Fill post probs for index SNP
  Post_Med_Prob[X_id] <- 1
  Med_Effect_Size[X_id] <- 0
  if(PrCS_weights == "Pr(M_C)"){
    temp_PrM_res <- mJAM_get_PrM(GItGI = GItGI_curr, GIty = GIty_curr,
                 yty = yty_curr, yty_med = yty_med, 
                 N_GWAS = N_GWAS, C_id = X_id, 
                 prev_X_list = prev_X_list, g = sum(N_GWAS),
                 rare_SNPs = rare_SNPs)
    Post_Model_Prob[X_id] <- Post_CS_Prob[X_id] <- sum_Post_CS_Porb <- Index_Model_Prob <- temp_PrM_res$post_prob
    Post_Model_Prob_Ratio[X_id] <- 1
    R2_est[X_id] <- temp_PrM_res$R2_est
    n_miss[X_id] <- temp_PrM_res$n_miss
  }
  if(PrCS_weights == "Pr(Wald)"){
    Post_Model_Prob[X_id] <- Post_CS_Prob[X_id] <- sum_Post_CS_Porb <- Index_Model_Prob <- 
      mJAM_get_PrM_Wald(GItGI = GItGI_curr, GIty = GIty_curr, yty = yty_curr, yty_med = yty_med, 
                        N_GWAS = N_GWAS, C_id = X_id, rare_SNPs = rare_SNPs, 
                        prev_X_list = NULL, g = sum(N_GWAS))
    Post_Model_Prob_Ratio[X_id] <- 1
  }
  
  for(C_id in Testing_id){
    temp_med <- mJAM_get_PrMed(GItGI = GItGI_curr, GIty = GIty_curr,yty = yty_curr,
                               yty_med = yty_med, 
                               N_GWAS  = N_GWAS, g = sum(N_GWAS),
                               C_id = C_id, 
                               X_id = X_id, 
                               prev_X_list = NULL)
    Post_Med_Prob[C_id] <- temp_med$Post_Med_Prob
    Med_Effect_Size[C_id] <- temp_med$Med_Effect_Size
    
    if(PrCS_weights == "Pr(M_C)"){
      temp_PrM_res <- mJAM_get_PrM(GItGI = GItGI_curr, GIty = GIty_curr,
                              yty = yty_curr, yty_med = yty_med, 
                              N_GWAS = N_GWAS, C_id = C_id, 
                              prev_X_list = prev_X_list, g = sum(N_GWAS),
                              rare_SNPs = rare_SNPs)
      temp_r2 <- temp_PrM_res$post_prob
      R2_est[C_id] <- temp_PrM_res$R2_est
      n_miss[C_id] <- temp_PrM_res$n_miss
    }
    if(PrCS_weights == "Pr(Wald)"){
      temp_r2 <- mJAM_get_PrM_Wald(GItGI = GItGI_curr, GIty = GIty_curr, N_GWAS = N_GWAS, 
                                   yty = yty_curr, yty_med = yty_med, rare_SNPs = rare_SNPs, 
                                   C_id = C_id, prev_X_list = NULL, g = sum(N_GWAS))
    }
    Post_Model_Prob[C_id] <- temp_r2
    Post_Model_Prob_Ratio[C_id] <- as.double(temp_r2 / Index_Model_Prob)
    Post_CS_Prob[C_id] <- temp_r2*Post_Med_Prob[C_id]
    sum_Post_CS_Porb <- sum_Post_CS_Porb + temp_r2*Post_Med_Prob[C_id]
  }
  
  ### --------------------------------------------------------- ##############
  
  ## --- Get "standardized" posterior probability of mediation test
  Post_Model_Prob <- rep(NA, length = length(All_id)); names(Post_Model_Prob) <- All_id
  # Post_CS_Prob <- rep(NA, length = length(All_id)); names(Post_CS_Prob) <- All_id
  Std_CS_Prob <- rep(NA, length(All_id)); names(Std_CS_Prob) <- All_id
  
  ## --- Fill post probs for index SNP
  if(PrCS_weights == "Pr(Wald)"){
    Post_Med_Prob[X_id] <- 1
    Post_Model_Prob[X_id] <- Post_CS_Prob[X_id] <- 
      mJAM_get_PrM_Wald(GItGI = GItGI_curr, GIty = GIty_curr, 
                        yty = yty_curr, yty_med = yty_med, rare_SNPs = rare_SNPs, 
                        N_GWAS = N_GWAS, C_id = X_id, prev_X_list = NULL, 
                        g = sum(N_GWAS))
    Std_CS_Prob[X_id] <- Post_Model_Prob[X_id]/sum_Post_CS_Porb
  }
  if(PrCS_weights == "Pr(M_C)"){
    Post_Med_Prob[X_id] <- 1
    Post_Model_Prob[X_id] <- Post_CS_Prob[X_id] <- 
      mJAM_get_PrM(GItGI = GItGI_curr, GIty = GIty_curr, 
                   yty = yty_curr, yty_med = yty_med, rare_SNPs = rare_SNPs, 
                   N_GWAS = N_GWAS, C_id = X_id, prev_X_list = prev_X_list, 
                   g = sum(N_GWAS))$post_prob
    Std_CS_Prob[X_id] <-
      as.numeric(mJAM_get_PrM(GItGI = GItGI_curr, GIty = GIty_curr, 
                              yty = yty_curr, yty_med = yty_med,rare_SNPs = rare_SNPs, 
                              N_GWAS = N_GWAS, C_id = X_id, prev_X_list = prev_X_list, 
                              g = sum(N_GWAS))$post_prob/sum_Post_CS_Porb)
  }
  
  for(C_id in Testing_id){
    
    if(PrCS_weights == "Pr(Wald)"){
      temp_r2 <-  mJAM_get_PrM_Wald(GItGI = GItGI_curr, GIty = GIty_curr, 
                                    yty = yty_curr, yty_med = yty_med, rare_SNPs = rare_SNPs, 
                                    N_GWAS = N_GWAS, C_id = C_id, prev_X_list = NULL, 
                                    g = sum(N_GWAS))
    }
    if(PrCS_weights == "Pr(M_C)"){
      temp_r2 <- mJAM_get_PrM(GItGI = GItGI_curr, GIty = GIty_curr, 
                              yty = yty_curr, yty_med = yty_med, rare_SNPs = rare_SNPs, 
                              N_GWAS = N_GWAS, C_id = C_id, prev_X_list = prev_X_list, 
                              g = sum(N_GWAS))$post_prob
    }
    
    Post_Model_Prob[C_id] <- temp_r2
    Std_CS_Prob[C_id] <- as.numeric(temp_r2*Post_Med_Prob[C_id]/sum_Post_CS_Porb)
  }

  if(PrCS_weights == "Pr(M_C)"){
    Post_Model_Prob_str <- rep(NA, length(Post_Model_Prob))
    R2_str <- rep(NA, length(R2_est))
    for(j in 1:length(Post_Model_Prob)){
      Post_Model_Prob_str[j] <- sub("\'mpfr1\' ", "", capture.output(Post_Model_Prob[[names(Post_Model_Prob[j])]]))
      R2_str[j] <- as.numeric(sub("\'mpfr1\' ", "", capture.output(R2_est[[names(R2_est[j])]])))
    }
  }
  if(PrCS_weights == "Pr(Wald)"){
    Post_Model_Prob_str <- Post_Model_Prob
  }
  
  Post_Model_Prob_Ratio2 <- Post_Model_Prob_Ratio
  Post_Model_Prob_Ratio2[which(Post_Model_Prob_Ratio>1)] <- 1
  Post_Med_Prob2 <- Post_Med_Prob
  Post_Med_Prob2[which(Post_Med_Prob<Pr_Med_cut)] <- 0
  SD_Post_CS_Prob2 <- (Post_Model_Prob_Ratio2*Post_Med_Prob2)/sum(Post_Model_Prob_Ratio2*Post_Med_Prob2, na.rm = T)
  
  Post_Med_Prob_df <-
    tibble(SNP_names = All_id,
           R2 = R2_str, n_miss = n_miss,
           Post_Model_Prob = Post_Model_Prob_str,
           Post_Model_Prob_Ratio = Post_Model_Prob_Ratio,
           Post_Model_Prob_Ratio2 = Post_Model_Prob_Ratio2,
           Med_Effect_Size = Med_Effect_Size,
           Post_Med_Prob = Post_Med_Prob,
           Post_Med_Prob2 = Post_Med_Prob2,
           SD_Post_CS_Prob = SD_Post_CS_Prob2) %>%
    arrange(desc(SD_Post_CS_Prob)) %>%
    mutate(CumSum_Porb = cumsum(SD_Post_CS_Prob),
           EmpiricalCut = min(CumSum_Porb[CumSum_Porb >= coverage], na.rm = T),
           duplicated_one = duplicated(CumSum_Porb),
           CS_in = ((CumSum_Porb <= EmpiricalCut) & !duplicated_one))%>%
    select(-duplicated_one)

  curr_CS <- Post_Med_Prob_df %>% 
    rename(CS_SNP = SNP_names) %>% 
    mutate(index_SNP = X_id)
  ###############################################################################################
  write.table(curr_CS, paste0(Forward_Dir,o_pre,"index",iter_count+1,"_CS.txt"),quote = F, 
              row.names = F, col.names = T, sep = '\t')
  
  ## step 3.2: Create CS using BF-based model prob  
  curr_CS_snp <- filter(curr_CS, CS_in == T) %>% pull(CS_SNP)
  all_CS <- rbind(all_CS, filter(curr_CS, CS_in == T))
  message(paste(length(curr_CS_snp), "CS SNPs selected for", curr_index_snp))
  
  ## step 3.3: prune out CS SNPs and SNPs in LD with index SNP; subset input statistics
  # curr_CS_snp <- curr_index_snp
  
  # pruned_snps <- c(pruned_snps, curr_CS_snp)
  # pruned_snp_df <- pruned_snp_df %>% add_row(pruned_SNP = curr_CS_snp, 
  #                                            index_SNP = curr_index_snp, 
  #                                            reason = "inCS")
  pruned_snps <- c(pruned_snps, curr_index_snp)
  curr_LD_snp <- mJAM_LDpruning(target = match(curr_index_snp,has_dosage_SNP),
                                testing = match(has_dosage_SNP[-match(pruned_snps,has_dosage_SNP)],has_dosage_SNP),
                                R = Dosage_cor,
                                within_thre = within_pop_threshold, across_thre = across_pop_threshold)
  pruned_snps <- c(pruned_snps,has_dosage_SNP[c(curr_LD_snp$remove_within, curr_LD_snp$remove_across)])
  pruned_snp_df <- pruned_snp_df %>% add_row(pruned_SNP = has_dosage_SNP[c(curr_LD_snp$remove_within, curr_LD_snp$remove_across)], 
                                             index_SNP = curr_index_snp, 
                                             reason = "highLD")
  
  ## step 4: update input statistics
  if(length(pruned_snps)>0){
    subset_EUR <- has_dosage_SNP[-match(pruned_snps, has_dosage_SNP)]
  }
  subset_EUR <- union(subset_EUR, selected_ids)
  subset_EUR <- has_dosage_SNP[has_dosage_SNP%in%subset_EUR]
  message(paste(length(unique(pruned_snps)), "SNPs got pruned.", length(subset_EUR), "SNPs left."))
  
  if(length(subset_EUR)>1){
    for(e in 1:numEthnic){
      GItGI_curr[[e]] <- GItGI[[e]][subset_EUR, subset_EUR] 
      GIty_curr[[e]] <- GIty[[e]][subset_EUR] 
      yty_curr[[e]] <- yty[[e]][subset_EUR] 
      colnames(GItGI_curr[[e]]) <- rownames(GItGI_curr[[e]]) <- subset_EUR
      names(GIty_curr[[e]]) <- names(yty_curr[[e]]) <- subset_EUR
    }
    ## step 5: continue to the next iteration
    iter_count <- iter_count+1
    prev_index_snp <- c(prev_index_snp, curr_index_snp)
    message("Continue to next round of index SNP selection.")
  }else{
    ## step 5: continue to the next iteration
    iter_count <- iter_count+1
    prev_index_snp <- c(prev_index_snp, curr_index_snp)
    message("Analysis ended.")
    break
  }
  message("Analysis ended.")
}

#########################

if(length(selected_ids)==1){
  MULTI_index <- tibble(SNP = selected_ids,MarkerName_old = sumstats3$MarkerName_old[match(selected_ids,sumstats3$MarkerName)],
                        condp = condp_list, finalp = min(chr1_EUR_sub$P.value), pcut = as.character(condp_cut))
}else{
  final_condp_selected <- mJAM_get_condp_selected(GItGI = GItGI, GIty = GIty, yty = yty,
                                                  N_GWAS = N_GWAS,yty_med = yty_med,
                                                  selected_id = selected_ids,
                                                  rare_SNPs = rare_SNPs)

  MULTI_index <- tibble(SNP = selected_ids, MarkerName_old = sumstats3$MarkerName_old[match(selected_ids,sumstats3$MarkerName)],
                        condp = condp_list, finalp = final_condp_selected$condp[,1], pcut = as.character(condp_cut)) %>%
    left_join(select(marginal_est, c(SNP, meta_p, mJAM_p)), by = "SNP") %>%
    select(c(SNP, MarkerName_old,meta_p, mJAM_p, everything()))
}
MULTI_index <- MULTI_index %>% mutate(condp = format(condp, scientific = T, digits = 3),
                                      finalp = format(finalp, scientific = T, digits = 3))
# write.table(MULTI_index, paste0(Forward_Dir, o_pre,"index_snps.txt"),quote = F, row.names = F, col.names = T, sep = '\t')

output_CS <- all_CS %>%
  select(c(index_SNP, CS_SNP,Post_Model_Prob_Ratio,Post_Model_Prob_Ratio2,
           Med_Effect_Size,Post_Med_Prob, Post_Med_Prob2, SD_Post_CS_Prob,CumSum_Porb)) %>%
  left_join(select(sumstats3, c(MarkerName, MarkerName_old)), by = c("CS_SNP" = "MarkerName")) %>%
  mutate(Post_Model_Prob_Ratio = format(Post_Model_Prob_Ratio, digits = 3),
         Post_Model_Prob_Ratio2 = format(Post_Model_Prob_Ratio2, digits = 3),
         Med_Effect_Size = format(Med_Effect_Size, digits = 3),
         Post_Med_Prob = format(Post_Med_Prob, digits = 3),
         Post_Med_Prob2 = format(Post_Med_Prob2, digits = 3),
         SD_Post_CS_Prob = format(SD_Post_CS_Prob, digits = 3),
         CumSum_Porb = format(CumSum_Porb, digits = 3))
write.table(output_CS, paste0(Forward_Dir,o_pre, "95CS.txt"),quote = F, row.names = F, col.names = T, sep = '\t')

###
# write.table(pruned_snp_df, paste0(Forward_Dir, o_pre,"pruned_SNPs.txt"),quote = F, row.names = F, col.names = T, sep = '\t')
###


ColNames <- colnames(sumstats3)
ColNames <- sub("Eur_17_study_Nov2021_1.txt","EUR", ColNames)
ColNames <- sub("AAPC_10_study_Oct2021_1.txt","AA", ColNames)
ColNames <- sub("Hispanic_5_study_Oct2021_1.txt","LA", ColNames)
ColNames <- sub("Asian_5_study_Oct2021_1.txt","ASN", ColNames)
sumstats4 <- filter(sumstats3, MarkerName %in% all_CS$CS_SNP)
troubled_sumstats <- filter(sumstats3, MarkerName %in% troubled_est$SNP)
colnames(sumstats4) <- ColNames
colnames(troubled_sumstats) <- ColNames
sumstats4_simp <- sumstats4 %>%
  left_join(all_CS, by = c("MarkerName" = "CS_SNP")) %>% 
  select(c(index_SNP, MarkerName, MarkerName_old, Effect, P.value, EUR.Freq1, EUR.Effect, EUR.P.value, AA.Freq1, AA.Effect, AA.P.value,
           ASN.Freq1, ASN.Effect, ASN.P.value, LA.Freq1, LA.Effect, LA.P.value)) %>%
  rowwise() %>%
  mutate(Freq.common = sum(EUR.Freq1>0.1&EUR.Freq1<0.9,AA.Freq1>0.1&AA.Freq1<0.9,
                           ASN.Freq1>0.1&ASN.Freq1<0.9,LA.Freq1>0.1&LA.Freq1<0.9 ,na.rm = T),
         Pos.effect = sum(EUR.Effect>0&EUR.P.value<0.001,AA.Effect>0&AA.P.value<0.001,
                          ASN.Effect>0&ASN.P.value<0.001,LA.Effect>0&LA.P.value<0.001 ,na.rm = T),
         Neg.effect = sum(EUR.Effect<0&EUR.P.value<0.001,AA.Effect<0&AA.P.value<0.001,
                          ASN.Effect<0&ASN.P.value<0.001,LA.Effect<0&LA.P.value<0.001 ,na.rm = T),
         Same.dir = ((Pos.effect == 0)|(Neg.effect == 0))) %>%
  select(-c(Pos.effect, Neg.effect))
troubled_sumstats <- troubled_sumstats %>%
  select(c(MarkerName,Effect, P.value, EUR.Freq1, EUR.Effect, EUR.P.value, AA.Freq1, AA.Effect, AA.P.value,
           ASN.Freq1, ASN.Effect, ASN.P.value, LA.Freq1, LA.Effect, LA.P.value)) %>%
  rowwise() %>%
  mutate(Freq.common = sum(EUR.Freq1>0.1&EUR.Freq1<0.9,AA.Freq1>0.1&AA.Freq1<0.9,
                           ASN.Freq1>0.1&ASN.Freq1<0.9,LA.Freq1>0.1&LA.Freq1<0.9 ,na.rm = T),
         Pos.effect = sum(EUR.Effect>0&EUR.P.value<0.001,AA.Effect>0&AA.P.value<0.001,
                          ASN.Effect>0&ASN.P.value<0.001,LA.Effect>0&LA.P.value<0.001 ,na.rm = T),
         Neg.effect = sum(EUR.Effect<0&EUR.P.value<0.001,AA.Effect<0&AA.P.value<0.001,
                          ASN.Effect<0&ASN.P.value<0.001,LA.Effect<0&LA.P.value<0.001 ,na.rm = T),
         Same.dir = ((Pos.effect == 0)|(Neg.effect == 0))) %>%
  select(-c(Pos.effect, Neg.effect))
# write.table(sumstats4, paste0(Forward_Dir, o_pre,"CS_sumstats.txt"),quote = F, row.names = F, col.names = T, sep = "\t")
write.table(sumstats4_simp, paste0(Forward_Dir,o_pre, "CS_sumstats_simplified.txt"),quote = F, row.names = F, col.names = T, sep = "\t")
# write.table(troubled_sumstats, paste0(Forward_Dir,o_pre, "troubled_sumstats_simplified.txt"),quote = F, row.names = F, col.names = T, sep = "\t")

## plot multi-ethnic CS
cs_snp <- all_CS %>%
  left_join(sumstats3, by = c("CS_SNP"="MarkerName")) %>%
  mutate(is_index = (index_SNP == CS_SNP)) %>%
  left_join(tibble(index_SNP = selected_ids, oldname = sumstats3$MarkerName_old[match(selected_ids,sumstats3$MarkerName)],
                   order = 1:length(selected_ids)), by = "index_SNP") %>%
  mutate(index_SNP_legend = paste0("(#",order,")",oldname))

EUR_cor_cs <- AA_cor_cs <- LA_cor_cs <- Asian_cor_cs <- rep(NA, nrow(cs_snp))
for(i in 1:nrow(cs_snp)){
  temp_cs_cnp <- as.character(cs_snp[i,1])
  temp_index_snp <- as.character(cs_snp[i,14])
  EUR_cor_cs[i] <- Dosage_cor[[1]][temp_cs_cnp, temp_index_snp]
  AA_cor_cs[i] <- Dosage_cor[[2]][temp_cs_cnp, temp_index_snp]
  LA_cor_cs[i] <- Dosage_cor[[3]][temp_cs_cnp, temp_index_snp]
  Asian_cor_cs[i] <- Dosage_cor[[4]][temp_cs_cnp, temp_index_snp]
}
cs_snp$EUR_cor_cs = EUR_cor_cs; cs_snp$AA_cor_cs = AA_cor_cs
cs_snp$LA_cor_cs = LA_cor_cs; cs_snp$Asian_cor_cs = Asian_cor_cs
CS_corr = select(cs_snp, c(index_SNP, CS_SNP,EUR_cor_cs, AA_cor_cs,LA_cor_cs,Asian_cor_cs)) %>%
  mutate(EUR_cor_cs = format(EUR_cor_cs, digits = 3),
         AA_cor_cs = format(AA_cor_cs, digits = 3),
         LA_cor_cs = format(LA_cor_cs, digits = 3),
         Asian_cor_cs = format(Asian_cor_cs, digits = 3)) %>% 
  group_by(index_SNP) %>% 
  arrange(desc(EUR_cor_cs), desc(AA_cor_cs),desc(Asian_cor_cs), desc(LA_cor_cs))
write.table(CS_corr,
            paste0(Forward_Dir, o_pre,"CS_corr.txt"),quote = F, row.names = F, col.names = T, sep = '\t')

sumstats3 %>%
  ggplot(aes(x = POS, y = -log10(as.numeric(P.value))))+
  geom_point(color="grey")+
  geom_point(aes(color = index_SNP_legend, group = index_SNP_legend), data = cs_snp)+
  geom_point(data = filter(cs_snp, is_index == T),
             colour="red", shape=1, size=2.5, stroke=1.5)+
  geom_text(
    data = sumstats3[match(selected_ids, CommonMarkerNames),],
    label= 1:length(selected_ids),
    nudge_x = 0.25, nudge_y = 0.3, check_overlap = T, size = 2.5
  )+
  geom_hline(yintercept = -log10(5e-8), linetype="dashed", color = "red")+
  geom_hline(yintercept = -log10(condp_cut), linetype="dotted", color = "blue")+
  labs(x = "Position", y = "-log10(meta.p.value)", title = paste0(CHR_center, ":",POS_lower, "-", POS_upper))
ggsave(paste0(Forward_Dir,o_pre,"plot_cs.png"), width = 9, height = 5)

## plot ethnic-specific CS
p_vals <- as.numeric(c(sumstats3$Eur_17_study_Nov2021_1.txt.P.value, sumstats3$Hispanic_5_study_Oct2021_1.txt.P.value,
                       sumstats3$AAPC_10_study_Oct2021_1.txt.P.value, sumstats3$Asian_5_study_Oct2021_1.txt.P.value))
if(min(p_vals, na.rm = T)<1e-20){
  nonEUR_p_vals <- as.numeric(c(sumstats3$Hispanic_5_study_Oct2021_1.txt.P.value,
                                sumstats3$AAPC_10_study_Oct2021_1.txt.P.value, sumstats3$Asian_5_study_Oct2021_1.txt.P.value))
}else{
  nonEUR_p_vals <- p_vals
}

p_EUR <- sumstats3 %>%
  ggplot(aes(x = POS, y = -log10(as.numeric(Eur_17_study_Nov2021_1.txt.P.value))))+
  geom_point(color="grey", size  = 1)+
  geom_point(aes(color = index_SNP, group = index_SNP), data = cs_snp, size  = 1)+
  geom_point(aes(group=index_SNP, alpha=abs(EUR_cor_cs)), data = cs_snp, size  = 1)+
  geom_point(data = filter(cs_snp, is_index == T),
             colour="red", shape=1, size=1.5, stroke=1)+
  geom_hline(yintercept = -log10(5e-8), linetype="dashed", color = "red")+
  geom_hline(yintercept = -log10(condp_cut), linetype="dotted", color = "blue")+
  geom_text(
    data = sumstats3[match(selected_ids, CommonMarkerNames),],
    label= 1:length(selected_ids),
    nudge_x = 0.25, nudge_y = 0.9, check_overlap = T, size = 2.5
  )+
  scale_y_continuous(limits = c(0, 1.05*(-log10(min(p_vals, na.rm = T)))))+
  labs(x = "Position", y = "-log10(EUR.p.value)",
       title = paste0(CHR_center, ":",POS_lower, "-", POS_upper, " -EUR"))+
  theme(legend.position="none")
p_AA <- sumstats3 %>%
  ggplot(aes(x = POS, y = -log10(as.numeric(AAPC_10_study_Oct2021_1.txt.P.value))))+
  geom_point(color="grey", size  = 1)+
  geom_point(aes(color = index_SNP, group = index_SNP), data = cs_snp, size  = 1)+
  geom_point(aes(group=index_SNP, alpha=abs(AA_cor_cs)), data = cs_snp, size  = 1)+
  geom_point(data = filter(cs_snp, is_index == T),
             colour="red", shape=1, size=1.5, stroke=1)+
  geom_hline(yintercept = -log10(5e-8), linetype="dashed", color = "red")+
  geom_hline(yintercept = -log10(condp_cut), linetype="dotted", color = "blue")+
  geom_text(
    data = sumstats3[match(selected_ids, CommonMarkerNames),],
    label= 1:length(selected_ids),
    nudge_x = 0.25, nudge_y = 0.9, check_overlap = T, size = 2.5
  )+
  scale_y_continuous(limits = c(0, 1.05*(-log10(min(nonEUR_p_vals, na.rm = T)))))+
  labs(x = "Position", y = "-log10(AA.p.value)",
       title = paste0(CHR_center, ":",POS_lower, "-", POS_upper, " -AA"))+
  theme(legend.position="none")
p_LA <- sumstats3 %>%
  ggplot(aes(x = POS, y = -log10(as.numeric(Hispanic_5_study_Oct2021_1.txt.P.value))))+
  geom_point(color="grey", size  = 1)+
  geom_point(aes(color = index_SNP, group = index_SNP), data = cs_snp, size  = 1)+
  geom_point(aes(group=index_SNP, alpha=abs(LA_cor_cs)), data = cs_snp, size  = 1)+
  geom_point(data = filter(cs_snp, is_index == T),
             colour="red", shape=1, size=1.5, stroke=2)+
  geom_hline(yintercept = -log10(5e-8), linetype="dashed", color = "red")+
  geom_hline(yintercept = -log10(condp_cut), linetype="dotted", color = "blue")+
  geom_text(
    data = sumstats3[match(selected_ids, CommonMarkerNames),],
    label= 1:length(selected_ids),
    nudge_x = 0.25, nudge_y = 0.9, check_overlap = T, size = 2.5
  )+
  scale_y_continuous(limits = c(0, 1.05*(-log10(min(nonEUR_p_vals, na.rm = T)))))+
  labs(x = "Position", y = "-log10(LA.p.value)",
       title = paste0(CHR_center, ":",POS_lower, "-", POS_upper, " -LA"))+
  theme(legend.position="none")
p_Asian <- sumstats3 %>%
  ggplot(aes(x = POS, y = -log10(as.numeric(Asian_5_study_Oct2021_1.txt.P.value))))+
  geom_point(color="grey", size  = 1)+
  geom_point(aes(color = index_SNP, group = index_SNP), data = cs_snp, size  = 1)+
  geom_point(aes(group=index_SNP, alpha=abs(Asian_cor_cs)), data = cs_snp, size  = 1)+
  geom_point(data = filter(cs_snp, is_index == T),
             colour="red", shape=1, size=1.5, stroke=1)+
  geom_hline(yintercept = -log10(5e-8), linetype="dashed", color = "red")+
  geom_hline(yintercept = -log10(condp_cut), linetype="dotted", color = "blue")+
  geom_text(
    data = sumstats3[match(selected_ids, CommonMarkerNames),],
    label= 1:length(selected_ids),
    nudge_x = 0.25, nudge_y = 0.9, check_overlap = T, size = 2.5
  )+
  scale_y_continuous(limits = c(0, 1.05*(-log10(min(nonEUR_p_vals, na.rm = T)))))+
  labs(x = "Position", y = "-log10(Asian.p.value)",
       title = paste0(CHR_center, ":",POS_lower, "-", POS_upper, " -Asian"))+
  theme(legend.position="none")
p_4 <- gridExtra::grid.arrange(p_EUR, p_AA, p_LA, p_Asian, nrow = 2)
ggsave(paste0(Forward_Dir,o_pre,"plot_cs_ethnic_spec.png"), plot = p_4, width = 8, height =6)


