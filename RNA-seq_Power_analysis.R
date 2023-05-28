##################################################################################
# TITLE - RNA-seq. Power analysis                                                #
# AUTHOR - Mukhtar Sadykov                                                       #
# DATE - December 1st, 2021                                                      #
# DESCRIPTION - This R script performs power analysis for RNA-seq experiments.   #
# It estimates the sample size needed to achieve a desired level of statistical. #
# power. By simulating data based on user-defined parameters, the script assesses#
# the ability of the experiment to detect differential gene expression.          #
##################################################################################

############################## IMPORTANT NOTES BEFORE STARTING ###############################
# 1. A "Data" folder has been included into the repository, this will include the processed  #
#    and normalized gene expression data, such as count matrices and metadata, once the work #
#    is ready for publication.                                                               #
# 2. A "Results" folder has been included into the repository. We have used this path to have#
#    all "write.csv and "save" commands for our different pipelines                          #
##############################################################################################

##############################
## POWER ANALYSIS ##
##############################

## based on ssizeRNA package
library(ssizeRNA)
library(edgeR)

# check.power(m = 14, mu = 10, disp = 0.1, fc = 2, sims = 10)
cts  <- read.csv("./Data/SKOs.csv", row.names = 1) # 60663
study <- read.delim("./Data/meta_SKOs.csv", row.names=1)

# Calculate dispersion for each gene
d <- DGEList(counts=cts, samples=study)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
disp <- d$tagwise.dispersion

##calculate the power to resolve the difference between H1_0 and H1_1 for the DKO H1_0_1
mu <- apply(cts[, c("H1_0_1", "H1_0_2", "H1_0_3", "H1_0_4",
"H1_1_1", "H1_1_2", "H1_1_3", "H1_1_4")], 1, mean)

size1 <- ssizeRNA_single(nGenes = 30000, pi0 = 0.80, m = 200, mu = mu,
disp = disp, fc = 2, fdr = 0.05,
power = 0.8, maxN = 15)

