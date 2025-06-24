# generate ENSG.coord.txt
library(biomaRt)
library(data.table)
human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
human_gene <- getBM(attributes = c("chromosome_name","ensembl_gene_id", "external_gene_name"),
                    filters = "", values = TRUE, mart = human)
colnames(human_gene) = c("CHR", "ENSG_ID", "GENE")

glist_hg19 <- fread("/public/home/shilulu/script/plot_smr/glist-hg19")
colnames(glist_hg19) = c("CHR", "START", "END", "GENE")
merge_df <- merge(glist_hg19, human_gene, by="GENE")
glist_hg19_ENSG <- merge_df[, c("ENSG_ID", "CHR.y", "START", "END")]
colnames(glist_hg19_ENSG) = c("GENE", "CHR", "START", "END")

fwrite(glist_hg19_ENSG, file="ENSG_coord.txt", sep="\t")