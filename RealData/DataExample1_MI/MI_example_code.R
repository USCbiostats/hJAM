rm(list=ls())

# Import summary statistics -----------------------------------------------
source("MetaXcan_fun.R")
library(hJAM)
load("20190219_asthma_smoke_bmi_t2d.RData")
load("20190420_i_jam_snp_eaf.RData")
colnames(bmi)
colnames(t2d)

# Filter the SNPs
nrow(bmi)
nrow(t2d)

combined = rbind(bmi, t2d)
combined = combined[!duplicated(combined$SNP), ]

rsid = combined$SNP
write.table(rsid, "snp_asked_all_X.txt", quote = F, row.names = F, col.names = F)
snp_ask = as.character(combined[, 1])

# Change the number of X
numX = 2

# Import reference data ---------------------------------------------------
load("frq_out.rdata")
load("ref_dta.rdata")

frq_out$ind = 1:nrow(frq_out)
combined = merge(combined, frq_out, by = 'SNP')
combined = combined[order(combined$ind), ] ## Order the order of SNPs in combined data

# Harmonize the data
all.equal(combined$REF, as.character(combined$Major_A))
diff_a1 = which(combined$REF != frq_out$Major_A)
for(i_snp in diff_a1){
  combined[i_snp, 'BETA.bmi'] = (-1)*combined[i_snp, 'BETA.bmi']
  combined[i_snp, 'BETA.t2d'] = (-1)*combined[i_snp, 'BETA.t2d']
  combined[i_snp, 'BETA.cvd'] = (-1)*combined[i_snp, 'BETA.cvd']
  tmp.alt = combined[i_snp, 'REF']
  tmp.ref = combined[i_snp, 'ALT']
  
  combined[i_snp, 'ALT'] = tmp.alt
  combined[i_snp, 'REF'] = tmp.ref
}
all.equal(combined$REF, as.character(combined$Major_A)) # check if the ref_allele are the same

## Obtain significant signals for each exposure of interest
traits = c("BMI", "Type 2 Diabetes")
traits_dataset = c(nrow(bmi), nrow(t2d))

l_snp = list()
l_snp[[1]] = bmi$SNP
l_snp[[2]] = t2d$SNP

# Construct Z matrix ------------------------------------------------------
Z = seZ = Z.jam = matrix(0, nrow = nrow(combined), ncol = 2)
colnames(Z) = c("bmi", "t2d")
Z[, 1] = combined$BETA.bmi
Z[, 2] = combined$BETA.t2d

seZ[, 1] = combined$SE.bmi
seZ[, 2] = combined$SE.t2d

Z.jam = get_cond_A(marginal_A = Z, Gl=ref_dta, N.Gx = 339224, ridgeTerm = T)

# Joint X's
betas.gwas = combined$BETA.cvd
se.betaG = combined$SE.cvd

r.jam.linear = hJAM_lnreg(betas.Gy = betas.gwas, A = Z.jam, N.Gy = 459324, Gl = ref_dta, ridgeTerm = T)
r.jam.linear

r.jam.bias = hJAM_egger(betas.Gy = betas.gwas, A = Z.jam, N.Gy = 459324, Gl = ref_dta, ridgeTerm = T)
r.jam.bias

library(MendelianRandomization)
rho = cor(ref_dta)
mr.mvivw.coef = mr_mvivw(mr_mvinput(bx = Z, bxse = seZ, by = betas.gwas,
                                    byse = se.betaG, correlation = rho))
mr.mvivw.coef
mr.mvivw.noCor.coef = mr_mvivw(mr_mvinput(bx = Z, bxse = seZ, by = betas.gwas,
                                          byse = se.betaG))
mr.mvivw.noCor.coef

# MVMR Egger
mr_mvegger(mr_mvinput(bx = Z, bxse = seZ, by = betas.gwas,
                      byse = se.betaG, correlation = rho))

# Marginal X's
mr.ivw.coef = mr.ivw.noCor.coef = list()
mr.egger.coef = mr.egger.top.coef = list()
r.jam.marginal = r.jam.top.marginal = list()
spred.coef = list()

# Calculate varX
## Calculate rho
ref_dta_rho = cor(ref_dta)
varG = apply(ref_dta, 2, var)
sigma_l_all = sqrt(varG)

varp = function(p) sqrt(2*p*(1-p))

for(i in 1:ncol(Z)){
  
  # Use top SNPs only
  # Construct the data: Z and betas
  topSNPs = frq_out[frq_out$SNP %in% l_snp[[i]], ]
  ref_dta_topSNP = as.matrix(ref_dta[, topSNPs$ind], nrow = nrow(ref_dta))
  Z_topSNP = Z[topSNPs$ind, i]
  seZ_topSNP = seZ[topSNPs$ind, i]
  betas.gwas_topSNP = betas.gwas[topSNPs$ind]
  se.betaG_topSNP = se.betaG[topSNPs$ind]
  
  # Change the direction
  negative_ind = which(Z_topSNP<0)
  ref_dta_topSNP[, negative_ind] = ifelse(ref_dta_topSNP[, negative_ind] == 0, 2,
                                          ifelse(ref_dta_topSNP[, negative_ind] == 2, 0, 1))
  Z_topSNP[negative_ind] = (-1)*Z_topSNP[negative_ind]
  betas.gwas_topSNP[negative_ind] = (-1)*betas.gwas_topSNP[negative_ind]
  
  # Calculate for metaxcan
  ## gamma_g
  rhox_topSNP = cor(ref_dta_topSNP)
  px_topSNP = frq_out[topSNPs$ind, "ref_frq"]
  px_topSNP[negative_ind] = 1 - px_topSNP[negative_ind]
  varx_topSNP = varp(px_topSNP)
  Gamma_g = (varx_topSNP%*%t(varx_topSNP))*rhox_topSNP
  
  ## sigma_l
  sigma_l = sigma_l_all[topSNPs$ind]
  
  ## sigma_g
  sigma_g = sqrt(t(Z_topSNP)%*%Gamma_g%*%Z_topSNP)
  
  spred.coef[[i]] = metaXcan(betas = betas.gwas_topSNP, se.betaG = se.betaG_topSNP, Gamma_G = Gamma_g,
                             Z = Z_topSNP, sigma_g, sigma_l)
  
  mr.ivw.noCor.coef[[i]] = mr_ivw(mr_input(bx = Z_topSNP, bxse = seZ_topSNP, by = betas.gwas_topSNP,
                                           byse = se.betaG_topSNP))
  
  mr.ivw.coef[[i]] = mr_ivw(mr_input(bx = Z_topSNP, bxse = seZ_topSNP, by = betas.gwas_topSNP,
                                     byse = se.betaG_topSNP, correlation = rhox_topSNP))
  
  mr.egger.top.coef[[i]] = mr_egger(mr_input(bx = Z_topSNP, bxse = seZ_topSNP, by = betas.gwas_topSNP,
                                             byse = se.betaG_topSNP))
}
