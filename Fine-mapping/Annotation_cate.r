#=============================================================#
# # # generate SNP annotation of  Number of 1-SNP credible sets
#=============================================================#
# read local credible set
setwd("/public/home/shilulu/project_hearing-loss/new_run/all_meta/gctb/fine_mapping/")
lcs = "ARHL_MVP_AJHG_BBJ_reformatMETAL_gctb_SbayesRC_finemapping.lcs"
library(data.table)
library(dplyr)

lcs = fread('ARHL_MVP_AJHG_BBJ_reformatMETAL_gctb_SbayesRC_finemapping.lcs')[, .(size, SNP)][, SNP := gsub("[\t\r\n]", "", SNP)]; 
lcs_snp = lcs[size==1, SNP]
baseline2 = fread('/public/share/wchirdzhq2022/Wulab_share/GCTB_ref/annot_baseline2.2.txt')
lcs_snp_annot = baseline2[SNP %in% lcs_snp]
fwrite(lcs_snp_annot, 'finemapping_snp.txt', col.names=TRUE, sep='\t')