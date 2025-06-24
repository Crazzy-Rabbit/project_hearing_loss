# # # 2) generate homologous gene in R
# if (!require("BiocManager", quietly = TRUE)) 
#     install.packages("BiocManager")
# BiocManager::install("biomaRt")
library(dplyr)
library(data.table)
library(biomaRt)
listEnsembl()
# connect Ensembl
mouse <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

homologs <- getBM(attributes = c("chromosome_name","ensembl_gene_id", "external_gene_name",
                                "mmusculus_homolog_ensembl_gene", "mmusculus_homolog_orthology_type"),
                  filters = "with_mmusculus_homolog", 
                  values = TRUE, 
                  mart = human)
# head(homologs)
homologs_1to1 <- subset(homologs, mmusculus_homolog_orthology_type == "ortholog_one2one")
# rm MT X Y chromosome genes
homGene_1to1_chr1_22 <- homologs_1to1 %>% 
  filter(!chromosome_name %in% c("MT", "X", "Y")) %>%
  rename(human_ensembl_gene = ensembl_gene_id)

fwrite(homGene_1to1_chr1_22, file="HomGene1to1_mouse2human.txt", sep="\t", row.names=FALSE)