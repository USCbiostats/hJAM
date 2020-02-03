## ----setup---------------------------------------------------------------
library(hJAM)

## ----data_check----------------------------------------------------------
data("conditional_A")
data("betas.Gy")
data("SNPs_info")
conditional_A[1:10, ]
betas.Gy[1:10]
SNPs_info[1:10, ]

## ----graphic_check-------------------------------------------------------
scatter_plot_p = SNPs_scatter_plot(A = conditional_A, betas.Gy = betas.Gy, num_X = 2)
heatmap_p = SNPs_heatmap(Gl)

## ----Amatrix-------------------------------------------------------------
data("marginal_A")
cond_A = get_cond_A(marginal_A = marginal_A, Gl = Gl, N.Gx = 339224, ridgeTerm = T)
cond_A[1:10, ]

## ----hjam_lnreg----------------------------------------------------------
hJAM::hJAM_lnreg(betas.Gy = betas.Gy, N.Gy = 459324, A = conditional_A, Gl = Gl, ridgeTerm = TRUE) # 459324 is the sample size of the UK Biobank GWAS of MI

## ----hjam_egger----------------------------------------------------------
hJAM::hJAM_egger(betas.Gy = betas.Gy, N.Gy = 459324, A = conditional_A, Gl = Gl, ridgeTerm = TRUE) # 459324 is the sample size of the UK Biobank GWAS of MI

