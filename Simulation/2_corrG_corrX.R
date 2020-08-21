rm(list=ls())
# import functions
source("sourceFunc.R")

# require library
library(mvtnorm); library(plyr); library(ggplot2); library(ggpubr); library(reshape2)
library(hJAM); library(MendelianRandomization); library(AER)

# simulate data : parameter set up
set.seed(123)
N_Gl = 500; N_Gx = 1000; N_Gy = 5000 # population set up
numSNP = 30; numGenes = 2  # genotype set up
pi1 = pi2 = c(0,0.1)
pigroup = length(pi1)*length(pi2)
iterationTimes = 1000
SIGMA_sim = getSIGMA(numSNP/numGenes, numGenes, NULL, presetcorrG = 0.6)
SIGMA_sim[1:10, 11:30] = 0
SIGMA_sim[11:20, c(1:10, 21:30)] = 0
SIGMA_sim[21:30, 1:20] = 0

MAF = runif(numSNP, 0.05,0.3)

# Create matrix to store output summaries
out.cond.1 = matrix(NA, nrow=iterationTimes*numGenes*pigroup, ncol = 24)
colnames(out.cond.1) = 
  c('iterTimes', 'pi.indicator','truePi', 'tsls.pi', 'tsls.se',
    'hjAM.pi','hjAM.se', 'metaXcan.pi', 'metaXcan.se',
    'mr.ivw.pi','mr.ivw.se', 'mr.ivw.noCor.pi', 'mr.ivw.noCor.se',
    'mr.mvivw.pi','mr.mvivw.se', 'mr.mvivw.noCor.pi','mr.mvivw.noCor.se', 
    'tsls.pvalue', 'hJAM.pvalue', 'metaXcan.pvalue','mr.ivw.pvalue',
    'mr.ivw.noCor.pvalue', 'mr.mvivw.pvalue','mr.mvivw.noCor.pvalue')
out.cond.1[,1] = rep(1:iterationTimes, each=2)
out.cond.1[,2] = rep(1:2,times=iterationTimes)

currentIter = 0
for (g1 in 1:length(pi1)){
  for (g2 in 1:length(pi2)){
    currentIter = currentIter + 1
    pi = c(pi1[g1], pi2[g2])
    rowStart = (currentIter-1)*iterationTimes*numGenes+1
    rowEnd = currentIter*iterationTimes*numGenes
    out.cond.1[rowStart:rowEnd,3] = rep(pi,times=iterationTimes)
    
    cat(c(g1,g2),"\n")
    for (k in 1:iterationTimes){
      
      # ------------------------  Simulate data : G, E, Y -------------------------------# 
      
      Gl = createGeno(N_Gl, numSNP, MAF, SIGMA_sim)
      Gx = createGeno(N_Gx, numSNP, MAF, SIGMA_sim)
      Gy = createGeno(N_Gy, numSNP, MAF, SIGMA_sim)
      
      A = matrix(0, nrow=numSNP, ncol=numGenes)
      A[4:6,1] = simBetabyR2(r2 = 0.1, Nrow=1, Ncol=3, varY=1, X=Gx[,4:6], sd=1) #alpha 1
      A[14:16,2] = simBetabyR2(r2 = 0.1, Nrow=1, Ncol=3, varY=1, X=Gx[,14:16], sd=1) #alpha 2
      A[24:26,] = simBetabyR2(r2 = 0.5, Nrow=2, Ncol=3, varY=1, X=Gx[,24:26], sd=1)
      
      Ex = matrix(NA, nrow=N_Gx, ncol=numGenes)
      Ey = matrix(NA, nrow=N_Gy, ncol=numGenes)
      
      for(i in 1:numGenes){
        Ex[,i] = rnorm(N_Gx, mean=Gx%*%A[,i], sd=1)
        Ey[,i] = rnorm(N_Gy, mean=Gy%*%A[,i], sd=1)
      }
      
      # simulate Y data
      Yy = rnorm(N_Gy, mean=Ey%*%pi, sd=1)
      
      # row indicator
      s1 = rowStart+(k-1)*numGenes
      s2 = rowStart+k*numGenes-1
      
      if(k %% 100==0)
        cat(k,"\n")
      
      # ------------------------  Individual data -------------------------------# 
      Gy.reg = ivreg(Yy ~ Ey | Gy) 
      out.cond.1[s1:s2,c(4:5,18)] = summary(Gy.reg)$coef[2:(numGenes+1),c(1:2,4)]
      
      # ------------------------  Obtain marginal Z matrix -------------------------------# 
      # Obtain Z matrix: Gx-->X : dataset 1
      Z = seZ = pZ = matrix(0,nrow=numSNP, ncol=numGenes)   # dim = SNP by numGenes
      for(i in 1:numGenes){
        var.Ex = var(Ex[, i])
        for(margin.Z in 1:numSNP){
          Ex.reg = lm(Ex[,i] ~ Gx[,margin.Z])
          Z[margin.Z,i] = summary(Ex.reg)$coef[2,1]
          seZ[margin.Z,i] = summary(Ex.reg)$coef[2,2]
          pZ[margin.Z,i] = summary(Ex.reg)$coef[2,4]
        }
      }
      
      # get summary statistics from Gy -> Y
      betas.gwas = se.betaG = p.betaG = rep(0,ncol(Gy)) #beta_x
      for (i in 1:ncol(Gy)){
        Yy.reg = lm(Yy ~ Gy[,i])
        betas.gwas[i] = summary(Yy.reg)$coef[2,1]
        se.betaG[i] = summary(Yy.reg)$coef[2,2]
        p.betaG[i] = summary(Yy.reg)$coef[2,4]
      }
      
      #-------------------- screening step for alpha
      
      lower.5.p.alpha = ifelse(pZ < 0.2, T, F)
      lower.5.pZ = ifelse(lower.5.p.alpha[,1] + lower.5.p.alpha[,2] > 0, T, F)
      #
      # # remove rows with p_beta > 0.1 in conditional Z matrix
      Z1 = Z[lower.5.pZ, ]
      seZ1 = seZ[lower.5.pZ, ]
      #
      # # remove rows with p_beta > 0.1 in beta list
      betas.gwas = betas.gwas[lower.5.pZ]
      se.betaG = se.betaG[lower.5.pZ]
      #
      # # remove rows with p_beta > 0.1 in Gl
      Gl = Gl[, lower.5.pZ]
      Gx = Gx[, lower.5.pZ]
      #
      
      #-------------------- Obtain conditional Z matrix
      Z.jam = get_cond_Z(marginal_Z = Z1, Gl = Gl, N.Gx = N_Gx, ridgeTerm = T)
      colnames(Z.jam) = c('pi1', 'pi2')
      # ------------------- Run hJAM 
      
      r.jam.linear = hJAM_lnreg(betas.Gy = betas.gwas, N.Gy = N_Gy, Gl = Gl, Z = Z.jam)
      out.cond.1[s1:s2, 6] = r.jam.linear$Estimate
      out.cond.1[s1:s2, 7] = r.jam.linear$StdErr
      out.cond.1[s1:s2, 19] = r.jam.linear$Pvalue
      
      rho_Gl = cor(Gl)
      mr.mvivw.coef = mr_mvivw(mr_mvinput(bx = Z1, bxse = seZ1, by = betas.gwas,
                                          byse = se.betaG, correlation = rho_Gl))
      mr.mvivw.noCor.coef = mr_mvivw(mr_mvinput(bx = Z1, bxse = seZ1, by = betas.gwas,
                                                byse = se.betaG))
      out.cond.1[s1:s2,14] = mr.mvivw.coef$Estimate
      out.cond.1[s1:s2,15] = mr.mvivw.coef$StdError
      out.cond.1[s1:s2,23] = mr.mvivw.coef$Pvalue
      
      out.cond.1[s1:s2,16] = mr.mvivw.noCor.coef$Estimate
      out.cond.1[s1:s2,17] = mr.mvivw.noCor.coef$StdError
      out.cond.1[s1:s2,24] = mr.mvivw.noCor.coef$Pvalue
      
      #-------------------- Marginal X analysis
      
      for(i in 1:numGenes){
        s = ifelse(i==1,rowStart+(k-1)*numGenes,rowStart+k*numGenes-1)
        
        # analysis strategy 1
        # Reference for MetaXcan: sigma^2_g, Gamma_g, w_lg
        ref.MX = metaXcanRef(Gl=Gl, Gx=Gx, Z.jam[,i])
        Gamma_G = ref.MX$Gamma_G
        sigma_g = ref.MX$sigma_g
        sigma_l = ref.MX$sigma_l
        
        metaXcan.coef = metaXcan(betas = betas.gwas, se.betaG = se.betaG, 
                                 Z=Z.jam[,i], Gamma_G, sigma_g = sigma_g, sigma_l = sigma_l)
        out.cond.1[s, 8] = metaXcan.coef$betas
        out.cond.1[s, 9] = metaXcan.coef$ses
        out.cond.1[s, 20] = metaXcan.coef$p
        
        mr.ivw.coef = mr_ivw(mr_input(bx = Z1[,i], bxse = seZ1[,i], by = betas.gwas, 
                                      byse = se.betaG, correlation = rho_Gl))
        
        mr.ivw.noCor.coef = mr_ivw(mr_input(bx = Z1[,i], bxse = seZ1[,i], by = betas.gwas, 
                                            byse = se.betaG))
        
        out.cond.1[s,10] = mr.ivw.coef$Estimate
        out.cond.1[s,11] = mr.ivw.coef$StdError
        out.cond.1[s,21] = mr.ivw.coef$Pvalue
        
        out.cond.1[s,12] = mr.ivw.noCor.coef$Estimate
        out.cond.1[s,13] = mr.ivw.noCor.coef$StdError
        out.cond.1[s,22] = mr.ivw.noCor.coef$Pvalue
      }
      
      rm(Z1, seZ1)
    }
  }       
}

# PLOT for conditional A matrix -----------------------------------------------------------------------
setwd("/home/rcf-40/jian848/lai/MRJAM/1SCRIPTS/simulationNov2019_noCor_paper/PLOTS")
#write.table(corrE, "corrG-indepX-corrE.txt",quote=F, row.names = F, col.names = T)

out.cond.1 = as.data.frame(out.cond.1)
out.cond.1$group = rep(1:pigroup, each=iterationTimes*numGenes)

write.table(out.cond.1, paste0("numSNP_",numSNP,"_corrG-corrX-simData1.txt"), quote=F, row.names = F, col.names = T)
