#####
sim_args = commandArgs(trailingOnly=TRUE)

sim_scenario = sim_args[1]
n_study_per_pop = as.integer(sim_args[2])
sample_size = sim_args[3]
n_snps = as.integer(sim_args[4])
n_causal = as.integer(sim_args[5])
r2 = as.numeric(sim_args[6])
effect_size = as.numeric(sim_args[7])
#####
# sim_scenario = "B"
# n_study_per_pop = 3
# sample_size = "Balanced"
# n_snps = 50
# n_causal = 1
# r2 = 0.6
# effect_size = 0.03
#####

.libPaths("/project/dconti_624/Users/shenjiay/RPackages/")
library(tidyverse)
source("/project/dconti_624/Users/shenjiay/mJAM/HPCSimu/CompareJAMEstimates/JAM.functions.v2.R")

########################################################
# Section 0: Specify parameters for the simulation 
########################################################
# Define causal IDs and correlation structure IDs to iterate
# Total number of SNPs
N_SNP = n_snps
block_size = ifelse(N_SNP == 50, 10, 50)
if(!(N_SNP %in% c(50,1000))){stop("Please specify N_SNP and its blocksize.")}

# Initiate the Causal index
N_Causal = n_causal
if(N_Causal != 0){
  Causal_ID = seq(1,N_SNP,by = block_size)[1:N_Causal]
}else{
  Causal_ID = NULL
}

# Initiate the correlation IDs = N of populations
N_Study_PerPop = n_study_per_pop
Cor_ID = c(rep(1,N_Study_PerPop),rep(2,N_Study_PerPop),rep(3,N_Study_PerPop))

# Specify number of LD and the correlation 
# N_LD = 9
LD_corr = r2

# Number of subjects of each population
Sample_Size_Scenario = sample_size

if(Sample_Size_Scenario == "Balanced"){
  N_Sample_Per_Study = as.integer(15000/N_Study_PerPop)
  N_Sample = c(rep(N_Sample_Per_Study,N_Study_PerPop),rep(N_Sample_Per_Study,N_Study_PerPop),rep(N_Sample_Per_Study,N_Study_PerPop))
}else{
  if(N_Study_PerPop == 3){
    N_Sample = c(rep(11000,3), rep(2000,3),  rep(2000,3))
  }else{
    stop("N_Study_PerPop != 3. Please change N_Sample so that they sum up to 45,000.")
  }
}


# Number of replications 
N_Rep = 500
# Starting point of replication
Start_Rep = 1

# Set the effect size 
Effect_size = effect_size

# Set the p-value cutoff for index selection 
## final cutoff would be Bonp/N_SNP
Bonp <-  0.05

# The name of populations 
N_Study = length(Cor_ID)
StudyName = paste("S",c(1:N_Study),sep="")

########################################################
# Section 1: Create folders
########################################################

# Automatically generate all desired directories
GlobalPath = "/project/dconti_624/Users/shenjiay/mJAM/HPCSimu/"
Causal_Folder = paste0(N_SNP,"SNPs_",length(Causal_ID),"Causal_",paste(Causal_ID,collapse="_"),"_Cor", LD_corr)
Dosage_Folder = paste("Beta_",Effect_size,"_N_", sum(N_Sample),"_",Sample_Size_Scenario,"_", N_Study_PerPop,"S",sep="")
Dosage_Dir = paste(GlobalPath,"Dosage/",Dosage_Folder,"/",Causal_Folder,"/",sep="")
Marg_Dir = paste(GlobalPath,"MarginalEffects/",Dosage_Folder,"/",Causal_Folder,"/",sep="")
Meta_Dir = paste(GlobalPath,"MetaSoftResults/",Dosage_Folder,"/",Causal_Folder,"/",sep="")
SuSiE_Dir = paste(GlobalPath,"SuSiEv2Results/",Dosage_Folder,"/",Causal_Folder,"/",sep="")
SandP_Dir = paste(GlobalPath,"SandPResults/",Dosage_Folder,"/",Causal_Folder,"/p",Bonp, "/",sep="")
MsCAVIAR_Dir = paste(GlobalPath,"MsCAVIARResults/",Dosage_Folder,"/",Causal_Folder,"/",sep="")
CojoInput_Dir = paste(GlobalPath,"COJO/CojoInput/",Dosage_Folder,"/",Causal_Folder,"/",sep="")
CojoResults_Dir = paste(GlobalPath,"COJO/CojoResults/",Dosage_Folder,"/",Causal_Folder,"/",sep="")
PlinkDosage_Dir = paste(GlobalPath,"COJO/PlinkDosage/",Dosage_Folder,"/",Causal_Folder,"/",sep="")
Analysis_Output_Dir = paste0(GlobalPath, "Analysis/", Dosage_Folder,"/",Causal_Folder,"/")

# To created nested folders, use mkdir from CMD
if (dir.exists(Dosage_Dir)==FALSE){
  system(paste('mkdir -p ', Dosage_Dir, sep=''))
}
if (dir.exists(Marg_Dir)==FALSE){
  system(paste('mkdir -p ','"',Marg_Dir,'"',sep=''))
}
if (dir.exists(Meta_Dir)==FALSE){
  system(paste('mkdir -p ','"',Meta_Dir,'"',sep=''))
}
if (dir.exists(CojoInput_Dir)==FALSE){
  system(paste('mkdir -p ','"',CojoInput_Dir,'"',sep=''))
}
if (dir.exists(CojoResults_Dir)==FALSE){
  system(paste('mkdir -p ','"',CojoResults_Dir,'"',sep=''))
}
if (dir.exists(PlinkDosage_Dir)==FALSE){
  system(paste('mkdir -p ','"',PlinkDosage_Dir,'"',sep=''))
}
if (dir.exists(Analysis_Output_Dir)==FALSE){
  system(paste('mkdir -p ','"',Analysis_Output_Dir,'"',sep=''))
}
if (dir.exists(SandP_Dir)==FALSE){
  system(paste('mkdir -p ','"',SandP_Dir,'"',sep=''))
}
if (dir.exists(SuSiE_Dir)==FALSE){
  system(paste('mkdir -p ','"',SuSiE_Dir,'"',sep=''))
}


##############################################################################
# Section 1: This script is to simulate dosage and generate marginal results
##############################################################################
# Set desired MAF and compute the cut-off from a standard normal distribution
MAF = c(rep(0.2,N_Study_PerPop), rep(0.4,N_Study_PerPop), rep(0.6, N_Study_PerPop))
Cut_L = qnorm((1-MAF)^2)
Cut_U = qnorm(1-MAF^2)
# Define the three correlation structures
# CAUTION: the length of each correlation vector should be N_SNP and the first element should be 1.
Cor1 = rep(c(0,rep(LD_corr,block_size-1)),N_SNP/block_size)
Cor2 = rep(c(0,rep(LD_corr,block_size-1)),N_SNP/block_size)
Cor3 = rep(c(0,rep(LD_corr,block_size-1)),N_SNP/block_size)
# Cor1 = Cor2 = Cor3 =rep(0,N_SNP)
# if(N_LD>0){
#   for(i in 1:length(Causal_ID)){
#     Cor1[(Causal_ID[i]+1):(Causal_ID[i]+N_LD)] = LD_corr
#     Cor2[(Causal_ID[i]+1):(Causal_ID[i]+N_LD)] = LD_corr
#     Cor3[(Causal_ID[i]+1):(Causal_ID[i]+N_LD)] = LD_corr
#   }
# }


# CAUTION!
# Construct a conversion vector that stores the beta required in the linear regression to generate SNPs with desired correlation
# Desired correlationis a vector from 0 to 0.9 with 0.1 space between each pair
CortoBeta_Cor = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
CortoBeta_Beta = c(0, 0.25 ,0.4 ,0.6 ,0.84 ,1.1 ,1.4 ,1.8 , 2.5, 4, 8)
set.seed(2022)
##################################
# Iterate through all replications
##################################
# r2_cont_matrix <- matrix(NA, nrow = N_Rep, ncol = N_Study)
for (k in Start_Rep:(Start_Rep+N_Rep-1)){
  print(paste("Now running: marginal results replication ",k,sep=""))
  Dosage_FileName = paste("Dosage_",StudyName,"_",k,".txt",sep="")
  Y_FileName = paste("Y_",StudyName,"_",k,".txt",sep="")
  Marg_FileName = paste0("Marg_",k,".txt")
  MAF_FileName = paste0("MAF_",k,".txt")
  # Initialize the marginal output with SNP names
  SNP_Name = paste("rs",c(1:N_SNP),sep="")
  Marg_Output <-  MAF_Output <-  SNP_Name
  #################################
  # Iterate through all populations
  #################################
  for (e in 1:N_Study){
    ##### METHOD 1 of simulating data ######
    ## Repeatedly generate dosages until the cholesky decomposition is possible
    itr = 0
    repeat{
      itr = itr + 1
      # print(itr)
      # Generate SNP_Ref from normal distribution for this specific populations
      SNP_Ref = rnorm(N_Sample[e],mean=0,sd=1)
      SNP_Ref = (SNP_Ref >= Cut_U[e]) + (SNP_Ref >= Cut_L[e])
      # Generate all other SNPs with desired correlation structure
      Dosage = SNP_Ref
      for (i in 2:N_SNP){
        if (eval(parse(text = paste("Cor",Cor_ID[e],sep="")))[i] == 0){
          # Re-generate reference SNP
          SNP_Ref = rnorm(N_Sample[e],mean=0,sd=1)
          SNP_Ref = (SNP_Ref >= Cut_U[e]) + (SNP_Ref >= Cut_L[e])
          Dosage = cbind(Dosage,SNP_Ref)
        }
        if (eval(parse(text = paste("Cor",Cor_ID[e],sep="")))[i]> 0){
          temp_SNP = rnorm(N_Sample[e],mean=CortoBeta_Beta[which(CortoBeta_Cor==eval(parse(text = paste("Cor",Cor_ID[e],sep="")))[i])]*(SNP_Ref-mean(SNP_Ref)),1)
          temp_SNP = (temp_SNP >= Cut_U[e]) + (temp_SNP >= Cut_L[e])
          Dosage = cbind(Dosage,temp_SNP)
        }
      }
      colnames(Dosage) = SNP_Name
      temp_Chol = NULL
      temp_Chol = try(chol(t(Dosage)%*%Dosage))
      # temp_corr = stats::cor(Dosage)
      # if (!is.null(dim(temp_Chol)) & !is.singular.matrix(temp_corr)){
      if (!is.null(dim(temp_Chol))){
        #print("Great! Found Dosage with valid cholesky decomposition!")
        break
      }
    }
    ##### METHOD 2 of simulating data ######
    # sigma <- matrix(0,ncol = N_SNP,nrow = N_SNP)
    # sigma[Causal_ID:(Causal_ID+N_LD),Causal_ID:(Causal_ID+N_LD)] <- matrix(0.75, nrow = N_LD+1, ncol = N_LD+1)
    # diag(sigma) <- rep(1, N_SNP)
    # AllData <- simulateData(M = N_SNP, N.ref = N_Sample[e], N.ipd = N_Sample[e],
    #                         CausalEffects = Effect_size,SIGMA = sigma)
    # Dosage <- AllData$W
    #######
    # Output the dosage to desired location
    write.table(Dosage,paste(Dosage_Dir,Dosage_FileName[e],sep=""),quote=F,sep="\t",row.name=F,col.name=T,append=F)

    ######################
    ## --- Simulate continuous outcome based on selected causal SNP
    Beta = rep(0, N_SNP)
    Beta[Causal_ID] = rep(Effect_size, length(Causal_ID))
    X_cent <- scale(Dosage, center = TRUE, scale = FALSE)
    Y <- rnorm(N_Sample[e], mean=X_cent%*%Beta, sd = 1)
    Y_cent <- scale(Y, center = TRUE, scale = FALSE)

    # Output the dosage to desired location
    write.table(Y_cent,paste(Dosage_Dir,Y_FileName[e],sep=""),quote=F,sep="\t",row.name=F,col.name=F,append=F)

    ## --- Get summary statistics
    temp_Marg_Output <- susieR::univariate_regression(X_cent, Y)
    temp_MAF_Output <- colMeans(Dosage)/2
    Marg_Output <- cbind(Marg_Output, temp_Marg_Output$betahat, temp_Marg_Output$sebetahat)
    MAF_Output <- cbind(MAF_Output, temp_MAF_Output)

    # Marg_p = 2*pnorm(abs( temp_Marg_Output$betahat/temp_Marg_Output$sebetahat), lower.tail = F)
  }
  # Output the Marginal results for all populations
  write.table(Marg_Output,paste0(Marg_Dir,Marg_FileName),quote=F,sep="\t",row.name=F,col.name=F,append=F)
  write.table(MAF_Output,paste0(Marg_Dir,MAF_FileName),quote=F,sep="\t",row.name=F,col.name=F,append=F)
}
# write.table(r2_cont_matrix,paste0(Analysis_Output_Dir,"r2_whole_matrix.txt"),quote=F,sep="\t",row.name=F,col.name=F,append=F)


####################################
# Section 2: Run metasoft from CMD
####################################
Metasoft_Folder = "/project/dconti_624/Users/shenjiay/Softwares/Metasoft/"
for (k in Start_Rep:(Start_Rep+N_Rep-1)){
  Input = paste(Marg_Dir,"Marg_",k,".txt",sep="")
  Output = paste(Meta_Dir,"Meta_",k,".txt",sep="")
  MetaCMD_temp = paste("/project/dconti_624/Users/shenjiay/Softwares/jre1.8.0_271/bin/java",
                       "-jar", paste0(Metasoft_Folder, "Metasoft.jar"),
                       "-input",Input,
                       "-output",Output,
                       "-pvalue_table", paste0(Metasoft_Folder, "HanEskinPvalueTable.txt"),
                       sep=" ")
  system(MetaCMD_temp)
}


##############################################################################
# Section 3: Run mJAM-SnP
##############################################################################
source('/project/dconti_624/Users/shenjiay/Anqi/scripts/mJAM_fn_in_progress/mJAM_SuSiE_Forward.R')
source('/project/dconti_624/Users/shenjiay/Anqi/scripts/mJAM_fn_in_progress/mJAM_build_CS.R')
source('/project/dconti_624/Users/shenjiay/Anqi/scripts/mJAM_fn_in_progress/mJAM_LDpruning.R')
source('/project/dconti_624/Users/shenjiay/Anqi/scripts/mJAM_fn_in_progress/mJAM_condp.R')
source('/project/dconti_624/Users/shenjiay/Anqi/scripts/mJAM_fn_in_progress/mJAM_Forward.R')


for (k in Start_Rep:(Start_Rep+N_Rep-1)){
  print(paste("Now running: mJAM-SnP Replication ",k,sep=""))
  Dosage_FileName = paste("Dosage_",StudyName,"_",k,".txt",sep="")
  p_cases_FileName = paste("Y_",StudyName,"_",k,".txt",sep="")
  Marg_FileName = paste("Marg_",k,".txt",sep="")
  MAF_FileName = paste("MAF_",k,".txt",sep="")

  # Load the Marginal results: the 1st column is rs id, the rest are beta and SE for each of the ethnic group
  Marg_Result = read.table(paste(Marg_Dir,Marg_FileName,sep=""),header=F)
  MAF_Result = read.table(paste(Marg_Dir,MAF_FileName,sep=""),header=F)
  names(Marg_Result)[1] <- "SNP"
  names(MAF_Result)[1] <- "SNP"

  # Separate the marginal results by ethnic group
  Input_Dosage <- vector("list", N_Study)
  # Input_p_cases <-  rep(0, N_Study)
  has_dosage_SNP <- paste0("rs",1:N_SNP)
  for (e in 1:N_Study){
    Input_Dosage[[e]] = data.matrix(read.table(paste(Dosage_Dir,Dosage_FileName[e],sep=""),header=T))
    colnames(Input_Dosage[[e]]) <- has_dosage_SNP
    # Input_p_cases[e] = mean(read.table(paste(Dosage_Dir,p_cases_FileName[e],sep=""),header=T)[,1])
  }

  mJAM_Forward_res <- mJAM_Forward(N_GWAS = N_Sample,
                                   X_ref = Input_Dosage,
                                   Marg_Result = Marg_Result,
                                   EAF_Result = MAF_Result,
                                   condp_cut = 0.05/N_SNP,
                                   within_pop_threshold = 0.50,
                                   across_pop_threshold = 0.20,
                                   Pr_Med_cut = 0)

  write.table(mJAM_Forward_res$index, paste0(SandP_Dir, k, "_index.txt"),quote = F, row.names = F, col.names = T)
  write.table(mJAM_Forward_res$cs, paste0(SandP_Dir, k, "_CS.txt"),quote = F, row.names = F, col.names = T)

}


#########################################################################
# Section 3: This script is to run mJAM from our simulated datasets
# Caution: The code in this section assumes there are three populations
#########################################################################
# Load R package

# devtools::install_github("pjnewcombe/R2BGLiMS")
# devtools::install_github("stephenslab/susieR")
# install.packages("matrixcalc")

source('/project/dconti_624/Users/shenjiay/Anqi/scripts/mJAM_fn_in_progress/mJAM_SuSiE_Forward.R')
##################################
# Iterate through all replications
##################################
for (k in Start_Rep:(Start_Rep+N_Rep-1)){
  print(paste("Now running: mJAM-SuSiE Replication ",k,sep=""))
  Dosage_FileName = paste("Dosage_",StudyName,"_",k,".txt",sep="")
  p_cases_FileName = paste("Y_",StudyName,"_",k,".txt",sep="")
  Marg_FileName = paste("Marg_",k,".txt",sep="")
  MAF_FileName = paste("MAF_",k,".txt",sep="")
  
  
  # Load the Marginal results: the 1st column is rs id, the rest are beta and SE for each of the ethnic group
  Marg_Result = read.table(paste(Marg_Dir,Marg_FileName,sep=""),header=F)
  MAF_Result = read.table(paste(Marg_Dir,MAF_FileName,sep=""),header=F)
  # names(Marg_Result)[1] <- "SNP"
  # names(MAF_Result)[1] <- "SNP"
  
  # Separate the marginal results by ethnic group
  Input_MarglogOR <- Input_MargSElogOR <- Input_MAF <-  Input_Dosage <- vector("list", N_Study)
  Input_p_cases <-  rep(0, N_Study)
  for (e in 1:N_Study){
    Input_Dosage[[e]] = data.matrix(read.table(paste(Dosage_Dir,Dosage_FileName[e],sep=""),header=T))
    Input_p_cases[e] = mean(read.table(paste(Dosage_Dir,p_cases_FileName[e],sep=""),header=T)[,1])
    Input_MarglogOR[[e]] <- Marg_Result[,2*e]
    Input_MargSElogOR[[e]] <- Marg_Result[,2*e+1]
    Input_MAF[[e]] <- MAF_Result[,e+1]
  }
  
  ######### SuSiE #################
  mJAM.susie.results <- mJAM_SuSiE(marginal.betas = Input_MarglogOR,
                                   marginal.se = Input_MargSElogOR,
                                   EAFs = Input_MAF,
                                   N_GWAS = N_Sample,
                                   X_ref = Input_Dosage,
                                   SNP_names = as.character(Marg_Result$V1),
                                   logORToBeta = FALSE,
                                   p_cases = Input_p_cases,
                                   SuSiE_num_comp = 10,
                                   SuSiE_coverage = 0.95,
                                   SuSiE_min_abs_corr = 0.5,
                                   max_iter = 500)
  
  # Output the posterior probabilities
  mJAM_SuSiE_Output <- mJAM.susie.results$summary
  mJAM_SuSiE_Post_FileName = paste("mJAM_SuSiE_Post_",k,".txt",sep="")
  write.table(mJAM_SuSiE_Output,paste(SuSiE_Dir,mJAM_SuSiE_Post_FileName,sep=""),
              quote=F,sep="\t",row.name=T,col.name=T,append=F)
  # Output credible sets 0.95 and 0.8
  if(is.null(mJAM.susie.results$fit$sets$cs)){
    mJAM_SuSiE_CS_0.95_Output <- mJAM_SuSiE_get_cs(mJAM.susie.results$fit,0.95)
  }else{
    mJAM_SuSiE_CS_0.95_Output <- cbind(mJAM_SuSiE_get_cs(mJAM.susie.results$fit,0.95),
                                       mJAM.susie.results$fit$sets$purity)
  }
  mJAM_SuSiE_CS_0.95_FileName = paste("mJAM_SuSiE_CS_0.95_",k,".txt",sep="")
  write.table(mJAM_SuSiE_CS_0.95_Output,
              paste(SuSiE_Dir,mJAM_SuSiE_CS_0.95_FileName,sep=""),
              quote=F,sep="\t",row.name=F,col.name=T,append=F)
  
  # Output the alpha matrix at convergence
  alpha_output <- rbind(mJAM.susie.results$fit$alpha,
                        susie_get_pip(mJAM.susie.results$fit)) %>%
    as.data.frame() %>%
    mutate(L = c(1:10, 0),  ## L = 0 means the overall pip
           index_SNP = c(apply(mJAM.susie.results$fit$alpha,1,which.max), NA),
           resid_var = 1,
           niter = mJAM.susie.results$fit$niter,
           converged = mJAM.susie.results$fit$converged)
  write.table(alpha_output,
              paste0(SuSiE_Dir,paste0("mJAM_SuSiE_alpha_",k,".txt")),
              quote=F,sep="\t",row.name=F,col.name=T,append=F)
}


#################################################################################
# Section 5: Prepare the data to be ready for plink and cojo stepwise selection
#################################################################################
##################################
# Iterate through all replications
##################################
for (k in Start_Rep:(Start_Rep+N_Rep-1)){
  print(paste("Now running: plink and cojo replication ",k,sep=""))
  Dosage_FileName = paste("Dosage_",StudyName,"_",k,".txt",sep="")
  Y_FileName = paste("Y_",StudyName,"_",k,".txt",sep="")
  Marg_FileName = paste("Marg_",k,".txt",sep="")
  Meta_FileName = paste("Meta_",k,".txt",sep="")
  # Load the dosage of all ethnic groups
  Dosage_VarName = paste("Dosage_",StudyName,sep="")
  Y_VarName = paste("Y_",StudyName,sep="")
  Freq_VarName = paste("Freq_",StudyName,sep="")
  for (e in 1:N_Study){
    temp_Dosage = data.matrix(read.table(paste(Dosage_Dir,Dosage_FileName[e],sep=""),header=T))
    temp_Y = data.matrix(read.table(paste(Dosage_Dir,Y_FileName[e],sep=""),header=F))
    assign(Dosage_VarName[e], temp_Dosage)
    assign(Y_VarName[e], temp_Y)
    assign(Freq_VarName[e], colMeans(temp_Dosage)/2)
  }
  # Load the Marginal results: the 1st column is rs id, the rest are beta and SE for each of the ethnic group
  NameMarg_Result = read.table(paste(Marg_Dir,Marg_FileName,sep=""),header=F)
  # Separate the marginal results by ethnic group and compute the p values
  MarglogOR_VarName = paste("MarglogOR_",StudyName,sep="")
  MargSElogOR_VarName = paste("MargSElogOR_",StudyName,sep="")
  P_VarName = paste("P_",StudyName,sep="")
  for (e in 1:N_Study){
    assign(MarglogOR_VarName[e], NameMarg_Result[,2*e])
    assign(MargSElogOR_VarName[e], NameMarg_Result[,2*e+1])
    # Compute P values
    assign(P_VarName[e], 2*pnorm(abs(eval(parse(text = MarglogOR_VarName[e]))/(eval(parse(text = MargSElogOR_VarName[e])))),lower.tail=FALSE) )
  }
  # Output files preparation
  MAP_FileName = paste(k,"_",StudyName,".map",sep="")
  PED_FileName = paste(k,"_",StudyName,".ped",sep="")
  MA_FileName = paste(k,"_",StudyName,".ma",sep="")
  # Construct design matrix for MAP files, all .map files are the same: dimension N_SNP
  MAP_Chr = rep(7,N_SNP)
  MAP_RS = paste("rs",c(1:N_SNP),sep="")
  MAP_Dist = rep(0,N_SNP)
  MAP_Pos = c(1:N_SNP)
  MAP_Design = cbind(MAP_Chr,MAP_RS,MAP_Dist,MAP_Pos)
  # Construct design matrix for MA files: dimension N_SNP
  MA_SNP = paste("rs",c(1:N_SNP),sep="")
  # Use A as the default risk allele, C as the default alternative allele
  Default_A1 = "A"
  Default_A2 = "C"
  MA_A1 = rep(Default_A1,N_SNP)
  MA_A2 = rep(Default_A2,N_SNP)
  MA_Design = cbind(MA_SNP,MA_A1,MA_A2)
  # Initialize the output for the PED file for all populations
  PED_AllPop = NULL
  for (e in 1:N_Study){
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
    temp_Dosage = eval(parse(text = Dosage_VarName[e]))
    temp_Y = eval(parse(text = Y_VarName[e]))
    temp_1A = ifelse(temp_Dosage >= 1,Default_A1,Default_A2)
    temp_2A = ifelse(temp_Dosage >= 2,Default_A1,Default_A2)
    # Combine the above two matrix in a way that fits the format requirement
    temp_PED = matrix(rbind(temp_1A,temp_2A),nrow=N_Sample[e],ncol=2*N_SNP)
    # Column bind the PED_Design, Y and temp_PED to output
    PED_Affection = temp_Y + 1
    temp_PED = cbind(PED_Design,PED_Affection,temp_PED)
    # Append temp_PED to PED_AllPop
    PED_AllPop = rbind(PED_AllPop,temp_PED)
  }
  # Output MAP file for all population
  write.table(MAP_Design,paste(PlinkDosage_Dir,paste(k,"_AllPop.map",sep=""),sep=""),quote=F,sep="\t",row.name=F,col.name=F,append=F)
  # For PED_AllPop, adjust the Family_ID, the first column
  PED_AllPop[,1] = c(1:sum(N_Sample))
  # Output PED file for all populations
  write.table(PED_AllPop,paste(PlinkDosage_Dir,paste(k,"_AllPop.ped",sep=""),sep=""),quote=F,sep="\t",row.name=F,col.name=F,append=F)
  
  # For MA_AllPop, obtain the required info from metasoft results
  Meta_Result = read.table(paste(Meta_Dir,Meta_FileName,sep=""),header=F, skip=1)
  # Ge the mean Freq from all populations
  MA_AllPop_Freq = rep(0,N_SNP)
  for (e in 1:N_Study){
    MA_AllPop_Freq = MA_AllPop_Freq + eval(parse(text = Freq_VarName[e]))
  }
  MA_AllPop_Freq = MA_AllPop_Freq/N_Study
  # BETA_FE is the 4th column in MetaSoft output
  MA_AllPop_Beta =  Meta_Result$V4
  # STD_FE is the 5th column in MetaSoft output
  MA_AllPop_SE = Meta_Result$V5
  # PVALUE_FE is the 3rd column in MetaSoft output
  MA_AllPop_P = Meta_Result$V3
  # Total number of individuals in all populations
  MA_AllPop_N = rep(sum(N_Sample),N_SNP)
  # Combine the results for MA_AllPop
  MA_AllPop = cbind(MA_Design,MA_AllPop_Freq,MA_AllPop_Beta,MA_AllPop_SE,MA_AllPop_P,MA_AllPop_N)
  colnames(MA_AllPop) = c("SNP","A1","A2","freq","b","se","p","N")
  # Output MA file for all population
  write.table(MA_AllPop,paste(CojoInput_Dir,paste(k,"_AllPop.ma",sep=""),sep=""),quote=F,sep="\t",row.name=F,col.name=T,append=F)
  
  # Write local plink command to one file for all replications in current folder
  PlinkFolder = "/project/dconti_624/Users/shenjiay/Softwares/plink_linux/"
  PlinkDosage_Files = paste(k,"_AllPop",sep="")
  Plink_Command = paste("./plink --file ",PlinkDosage_Dir,PlinkDosage_Files," --make-bed --out ",CojoInput_Dir,PlinkDosage_Files,sep="")
  # Run plink commands
  system(paste("cd", PlinkFolder, "&& chmod +x plink &&",Plink_Command))
  
  # Write HPC COJO command to one file for all replications in current folder
  GCTAFolder = "/project/dconti_624/Users/shenjiay/Softwares/gcta_1.93.2beta/"
  BF_Alpha = 0.05/N_SNP
  CojoInput_Files = paste(k,"_AllPop",sep="")
  
  # Output Cojo commands
  # Nominal alpha
  Cojo_Command = paste("./gcta64  --bfile ",CojoInput_Dir,CojoInput_Files,
                       " --cojo-file ",CojoInput_Dir,CojoInput_Files,
                       ".ma --cojo-slct --cojo-p 0.05 --out ",CojoResults_Dir,CojoInput_Files,"NoBF",sep="")
  system(paste("cd", GCTAFolder, "&& chmod +x gcta64 &&",Cojo_Command))
  
  # Bonferroni alpha
  Cojo_Command = paste("./gcta64  --bfile ",CojoInput_Dir,CojoInput_Files,
                       " --cojo-file ",CojoInput_Dir,CojoInput_Files,
                       ".ma --cojo-slct --cojo-p ",BF_Alpha," --out ",CojoResults_Dir,CojoInput_Files,sep="")
  system(paste("cd", GCTAFolder, "&& chmod +x gcta64 &&",Cojo_Command))
}


