# Run analysis for NUCKS1, PM20D1, and RP11
rm(list=ls())

# Import summary statistics of GTEx -----------------------------------------------
load("CHR1_eQTL.Rdata")
library(hJAM); library(MendelianRandomization); library(ggplot2); library(ggpubr); library(dplyr)
library(reshape2); library(tidyr)
source("MetaXcan_fun.R")

## Obtain significant signals for each exposure of interest
load("NUCKS1_weights.Rdata")
load("PM20D1_weights.Rdata")

eQTL_weights = rbind(NUCKS1_weights, PM20D1_weights)
eQTL_weights = eQTL_weights[!duplicated(eQTL_weights$variant_id), ]
eQTL_weights = merge(eQTL_weights, eQTL_rsID[, c('variant_id', 'rsID')], by = 'variant_id')
eQTL_weights$rsID = as.character(eQTL_weights$rsID)
eQTL_weights = eQTL_weights[!duplicated(eQTL_weights$rsID), ]

# Change the number of X
numX = 2

# Import reference data ---------------------------------------------------
ref_dta = read.table("1000G_EUR_chr1.raw", header = T)
ref_id = ref_dta[, 1]
ref_dta = ref_dta[, -c(1:6)]

## Check the heatmap of the correlation matrix
rho = cor(ref_dta)
melted_cormat = melt(rho)
head(melted_cormat)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

# Set the eQTL order based on the reference data
eQTL_order = as.data.frame(matrix(unlist(strsplit(colnames(ref_dta), "_")), ncol = 2, byrow = T))
colnames(eQTL_order) = c("rsID", "A1")
eQTL_order$indx = 1:nrow(eQTL_order)
eQTL_order$rsID = as.character(eQTL_order$rsID)

r2 = c(3, 4, 5, 6)
r.jam.linear = mr.mvivw.coef = 
  mr.mvivw.noCor.coef = mr.ivw.coef =  mr.ivw.noCor.coef = 
  s.predixcan.coefR = list()

# Generate output files
outfile = list()
for(i_X in 1:2){
  outfile[[i_X]] = as.data.frame(matrix(0, nrow = 6*length(r2), ncol = 9))
  outfile[[i_X]][, 1] = rep(r2, each = 6)
  colnames(outfile[[i_X]]) = c("r2", "Methods", "num_eQTLs", "est", "se", "p-value", "lower", "upper", "95% CI")
  
  outfile[[i_X]][, 2] = rep(c("hJAM", "MVIVW MR", "MVIVW MR (no Corr)","IVW MR", "IVW MR (no Corr)",
                              "S-PrediXcan"), times = length(r2))
  outfile[[i_X]]$Methods = factor(outfile[[i_X]]$Methods, c("S-PrediXcan", "IVW MR (no Corr)", "IVW MR", 
                                                            "MVIVW MR (no Corr)", "MVIVW MR", 
                                                            "hJAM"))
}


for(i_r2 in 1:length(r2)){
  
  # Import the priority prunning results
  pp_eQTL = read.table(paste0("out_r", r2[i_r2], "_chr1.results"),
                       header = T)
  pp_eQTL_rsID = as.character(pp_eQTL[pp_eQTL$selected==1, 'name'])
  
  # Get the selected eQTLs
  eQTL_order_pp = eQTL_order[eQTL_order$rsID %in% pp_eQTL_rsID, 'indx']
  
  ## Order the order of eQTLs in the weight matrix
  eQTL_weights_order = merge(eQTL_weights, eQTL_order, by = 'rsID')
  eQTL_weights_order = eQTL_weights_order[order(eQTL_weights_order$indx), ]
  
  eQTL_weights_order[, c('CHR', 'POS', 'var_A1', 'var_A2', 'build')] = strsplit(as.character(eQTL_weights$variant_id), "_") %>%
    unlist() %>% matrix(ncol = 5, byrow = T) %>% as.data.frame() 
  
  all.equal(eQTL_weights_order$A1, eQTL_weights_order$var_A1)
  
  # Import prca data ---------------------------------------------------
  load("eQTL_CHR1_prca.Rdata")
  eQTL_weights_prca = merge(eQTL_weights_order, 
                            eQTL_prca[, c('variant_id', 'Freq1', 'FreqSE', 'Effect', 'StdErr', 'prca_A1', 'prca_A2')], 
                            by = 'variant_id')
  eQTL_weights_prca = eQTL_weights_prca[!duplicated(eQTL_weights_prca$rsID), ]
  eQTL_weights_prca = eQTL_weights_prca[order(eQTL_weights_prca$indx), ]
  
  all.equal(as.character(eQTL_weights_prca$A1), as.character(eQTL_weights_prca$prca_A1))
  
  # Construct Z matrix ------------------------------------------------------
  Z = seZ = matrix(0, nrow = length(pp_eQTL_rsID), ncol = 2)
  colnames(Z) = c("NUCKS1", "PM20D1")
  Z_indx = eQTL_weights_prca$indx[eQTL_order_pp]
  Z[, 1] = eQTL_weights_prca$alpha.nucks1[eQTL_order_pp]
  Z[, 2] = eQTL_weights_prca$alpha.pm20d1[eQTL_order_pp]
  
  seZ[, 1] = eQTL_weights_prca$se.nucks1[eQTL_order_pp]
  seZ[, 2] = eQTL_weights_prca$se.pm20d1[eQTL_order_pp]
  
  betas.gwas = eQTL_weights_prca$Effect[eQTL_order_pp]
  se.betaG = eQTL_weights_prca$StdErr[eQTL_order_pp]
  
  # Check the correlation
  rho = cor(ref_dta[, eQTL_order_pp])
  cond.Z = get_cond_A(marginal_A = Z, Gl = ref_dta[, eQTL_order_pp], N.Gx = 503, ridgeTerm = T)
  for(i in 1:nrow(cond.Z)){
    if(cond.Z[i, 1]<0){
      i_eQTL_order_pp = eQTL_order_pp[i]
      ref_dta[, i_eQTL_order_pp] = ifelse(ref_dta[, i_eQTL_order_pp]==2, 0, 
                                          ifelse(ref_dta[, i_eQTL_order_pp]==0, 2, ref_dta[, i_eQTL_order_pp]))
      cond.Z[i, ] = (-1)*cond.Z[i, ]
      Z[i, ] = (-1)*Z[i, ]
      betas.gwas[i] = (-1)*betas.gwas[i]
    }
  }
  
  # Run hJAM
  r.jam.linear[[i_r2]] = hJAM_lnreg(betas.Gy = betas.gwas, 
                                    N.Gy = nrow(ref_dta), 
                                    Gl = ref_dta[, eQTL_order_pp], 
                                    A = cond.Z, ridgeTerm = T)
  
  # Run MVIVW MR
  mr.mvivw.coef[[i_r2]] = mr_mvivw(mr_mvinput(bx = Z, bxse = seZ, by = betas.gwas,
                                              byse = se.betaG, correlation = cor(ref_dta[, eQTL_order_pp])))
  
  mr.mvivw.noCor.coef[[i_r2]] = mr_mvivw(mr_mvinput(bx = Z, bxse = seZ, by = betas.gwas,
                                                    byse = se.betaG))
  
  # Run MR Egger and IVW MR
  # Run marginal analysis
  l_snp = list()
  NUCKS1_weights_rsID = merge(NUCKS1_weights, eQTL_rsID, by = 'variant_id')
  PM20D1_weights_rsID = merge(PM20D1_weights, eQTL_rsID, by = 'variant_id')
  
  # Import the priority prunning results
  l_snp[[1]] = as.character(NUCKS1_weights_rsID[as.character(NUCKS1_weights_rsID$rsID) %in% pp_eQTL_rsID, 'variant_id'])
  l_snp[[2]] = as.character(PM20D1_weights_rsID[as.character(PM20D1_weights_rsID$rsID) %in% pp_eQTL_rsID, 'variant_id'])
  
  # Calculate varX
  ## Calculate rho
  ref_dta_rho = cor(ref_dta)
  varG = apply(ref_dta[, eQTL_order_pp], 2, var)
  sigma_l_all = sqrt(varG)
  varp = function(p) sqrt(2*p*(1-p))
  
  for(i_X in 1:numX){
    
    # hJAM
    outfile[[i_X]][(i_r2-1)*6+1, 3:8] = 
      c(r.jam.linear[[i_r2]]$numSNP, r.jam.linear[[i_r2]]$Estimate[i_X], r.jam.linear[[i_r2]]$StdErr[i_X],
        r.jam.linear[[i_r2]]$Pvalue[i_X], r.jam.linear[[i_r2]]$Lower.CI[i_X], r.jam.linear[[i_r2]]$Upper.CI[i_X])
    
    # mvivw
    outfile[[i_X]][(i_r2-1)*6+2, 3:8] = 
      c(r.jam.linear[[i_r2]]$numSNP, mr.mvivw.coef[[i_r2]]$Estimate[i_X], mr.mvivw.coef[[i_r2]]$StdError[i_X],
        mr.mvivw.coef[[i_r2]]$Pvalue[i_X], mr.mvivw.coef[[i_r2]]$CILower[i_X], mr.mvivw.coef[[i_r2]]$CIUpper[i_X])
    
    # mvivw no correlation
    outfile[[i_X]][(i_r2-1)*6+3, 3:8] = 
      c(r.jam.linear[[i_r2]]$numSNP, mr.mvivw.noCor.coef[[i_r2]]$Estimate[i_X], mr.mvivw.noCor.coef[[i_r2]]$StdError[i_X], 
        mr.mvivw.noCor.coef[[i_r2]]$Pvalue[i_X], mr.mvivw.noCor.coef[[i_r2]]$CILower[i_X], mr.mvivw.noCor.coef[[i_r2]]$CIUpper[i_X])
    
    # Use top SNPs only
    # Construct the data: Z and betas
    topSNPs = eQTL_weights_order[eQTL_weights_order$variant_id %in% l_snp[[i_X]], ]
    cat(paste0("numX = ", i_X, ": observed ", nrow(topSNPs)), "\n")
    cat(paste0("numX = ", i_X, ": expected ", length(l_snp[[i_X]])), "\n\n")
    
    ref_dta_topSNP = as.matrix(ref_dta[, topSNPs$indx], nrow = nrow(ref_dta))
    Z_topSNP = Z[which(Z_indx %in% topSNPs$indx), i_X]
    condZ_topSNP = cond.Z[which(Z_indx %in% topSNPs$indx), i_X]
    seZ_topSNP = seZ[which(Z_indx %in% topSNPs$indx), i_X]
    betas.gwas_topSNP = betas.gwas[which(Z_indx %in% topSNPs$indx)]
    se.betaG_topSNP = se.betaG[which(Z_indx %in% topSNPs$indx)]
    
    mr.ivw.coef[[i_X]] = mr_ivw(mr_input(bx = Z_topSNP, bxse = seZ_topSNP, by = betas.gwas_topSNP,
                                         byse = se.betaG_topSNP, correlation = cor(ref_dta_topSNP)))
    mr.ivw.noCor.coef[[i_X]] = mr_ivw(mr_input(bx = Z_topSNP, bxse = seZ_topSNP, by = betas.gwas_topSNP,
                                               byse = se.betaG_topSNP))
    
    # Calculate for metaxcan
    ## gamma_g
    rhox_topSNP = ref_dta_rho[topSNPs$indx, topSNPs$indx]
    px_topSNP = eQTL_weights_prca[topSNPs$indx, "Freq1"]
    varx_topSNP = varp(px_topSNP)
    Gamma_g = (varx_topSNP%*%t(varx_topSNP))*rhox_topSNP
    
    ## sigma_l
    sigma_l = sigma_l_all[which(Z_indx %in% topSNPs$indx)]
    
    ## sigma_g
    sigma_g = sqrt(t(Z_topSNP)%*%Gamma_g%*%Z_topSNP)
    
    s.predixcan.coefR[[i_X]] = metaXcan(betas = betas.gwas_topSNP, se.betaG = se.betaG_topSNP, Gamma_G = Gamma_g,
                                        Z = Z_topSNP, sigma_g, sigma_l)
    
    # ivw with top SNPs
    outfile[[i_X]][(i_r2-1)*6+4, 3:8] = 
      c(length(rhox_topSNP), mr.ivw.coef[[i_X]]$Estimate, mr.ivw.coef[[i_X]]$StdError,
        mr.ivw.coef[[i_X]]$Pvalue, mr.ivw.coef[[i_X]]$CILower, mr.ivw.coef[[i_X]]$CIUpper)
    
    # ivw with top SNPs, no correlation
    outfile[[i_X]][(i_r2-1)*6+5, 3:8] = 
      c(length(rhox_topSNP), mr.ivw.noCor.coef[[i_X]]$Estimate, mr.ivw.noCor.coef[[i_X]]$StdError, 
        mr.ivw.noCor.coef[[i_X]]$Pvalue,mr.ivw.noCor.coef[[i_X]]$CILower, mr.ivw.noCor.coef[[i_X]]$CIUpper)
    
    # S-PrediXcan
    outfile[[i_X]][(i_r2-1)*6+6, 3:8] = 
      c(length(rhox_topSNP), s.predixcan.coefR[[i_X]]$betas, s.predixcan.coefR[[i_X]]$ses, 
        s.predixcan.coefR[[i_X]]$p, s.predixcan.coefR[[i_X]]$beta-1.96*s.predixcan.coefR[[i_X]]$ses, s.predixcan.coefR[[i_X]]$beta+1.96*s.predixcan.coefR[[i_X]]$ses)
  }
}

# Import S-PrediXcan
s.predixcan.coef = list(c("betas" = 0.138788061, "ses" = 0.031331262, "p" = 9.44E-06),
                        c("betas" = 0.315447412, "ses" = 0.067053816, "p" = 2.55E-06))

# PLOT 
traits = c("NUCKS1", "PM20D1")
yaddon = c(1, 0.8)
yaddon_low = c(0.7, 0.8)
plist = estlist = list()
for(i_X in 1:numX){
  
  plot.dta = outfile[[i_X]] %>% mutate(gene = traits[i_X])
  
  plot.dta$r2 = plot.dta$r2/10
  plot.dta$est = exp(plot.dta$est)
  plot.dta$lower = exp(plot.dta$lower)
  plot.dta$upper = exp(plot.dta$upper)
  plot.dta$num_eQTLs = factor(plot.dta$num_eQTLs)
  
  plot.dta$`95% CI` = paste0(round(plot.dta$est, 3), " (", round(plot.dta$lower, 3), ", ", round(plot.dta$upper, 3), ")")
  
  maxy = max(plot.dta$est)
  
  annotation.dta = plot.dta %>% select(r2, Methods, '95% CI') %>% mutate(annotate_y = maxy+0.6)
  
  estlist[[i_X]] = 
    ggplot(plot.dta, aes(x=r2, y=est)) + #theme_bw() +
    geom_hline(yintercept=1, color='red2', size = 0.4)+
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper), width= 0.01, position=position_dodge(0.05)) +
    coord_flip() + ylim(yaddon_low[i_X], maxy+yaddon[i_X]) +
    xlab("R squared") + ylab('Odds ratio and 95% CI') +
    geom_label(data = annotation.dta, aes(x = r2, y = annotate_y, label = `95% CI`), 
               size = 3) +
    facet_wrap(Methods~.) + 
    ggtitle(traits[i_X]) +
    theme(plot.title = element_text(size = 14, face = "bold"),
    ) 
  
  annotation.dta = plot.dta %>% 
    mutate(p.show = ifelse(`p-value`<0.001, 'p < 0.001', paste0("p = ", round(`p-value`, 3)))) %>%
    select(r2, Methods, p.show) %>% mutate(annotate_y = maxy+0.6)
  
  plist[[i_X]] = 
    ggplot(plot.dta, aes(x=r2, y=est)) + #theme_bw() +
    geom_hline(yintercept=1, color='red2', size = 0.4)+
    geom_point()+
    geom_errorbar(aes(ymin=lower, ymax=upper), width= 0.01, position=position_dodge(0.05)) +
    coord_flip() + ylim(yaddon_low[i_X], maxy+yaddon[i_X]) +
    xlab("R squared") + ylab('Odds ratio and 95% CI') +
    geom_label(data = annotation.dta, aes(x = r2, y = annotate_y, label = p.show), 
               size = 3) +
    facet_wrap(Methods~.) + 
    ggtitle(traits[i_X]) +
    theme(plot.title = element_text(size = 14, face = "bold"),
    ) 
  
  cat("\n============", traits[i_X],"===========\n\n")
  print(plot.dta)
  
  write.table(plot.dta, paste0("univariate_result_", traits[i], ".txt"), sep = '\t', quote = F, row.names = F, col.names = T)
}

library(gridExtra)
png("GTEx.plot.png", width = 8, height = 9, res = 200, units = "in")
grid.arrange(plist[[1]], plist[[2]], ncol = 1)
dev.off()

png("GTEx.est.plot.png", width = 8, height = 9, res = 200, units = "in")
grid.arrange(estlist[[1]], estlist[[2]], ncol = 1)
dev.off()
