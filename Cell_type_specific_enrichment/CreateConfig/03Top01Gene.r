# # # 3) extract 1:1 homologs scRNA genes from mouse cochlea
library(dplyr)
library(data.table)
setwd("/public/share/wchirdzhq2022/Wulab_share/LDSC/Mouse_cochleae/")
scGene = fread("Sc_P8_12_20.txt")
homologs_1to1 = fread("HomGene1to1_mouse2human.txt")[,2:3,with=FALSE]
colnames(homologs_1to1) = c("human_gene", "gene")

merge_df = merge(scGene, homologs_1to1, by="gene")
merge_df = subset(merge_df, select = -c(gene))
colnames(merge_df) = gsub("[[:space:]/-]", "_", colnames(merge_df))
colnames(merge_df) = gsub("'", "_", colnames(merge_df))
merge_df = merge_df %>% select(human_gene, everything())

control = merge_df[, 1]
fwrite(control, file="Sc_P8_12_20_control.GeneSet", sep="\t", col.names=FALSE)

lapply(2:ncol(merge_df), function(i){
  temp_dt = data.frame(merge_df[[1]], merge_df[[i]])
  # colnames(df1) = names(merge_df)[c(37, i)]
  sorted_dt = temp_dt[order(temp_dt[[2]], decreasing = TRUE), ]
  top_10 = sorted_dt[1:floor(0.1*nrow(sorted_dt)), ]
  top_10_gene = top_10[, 1, drop=FALSE]

  cell_type = colnames(merge_df)[i]
  fwrite(top_10_gene, file=paste0("Sc_P8_12_20_", cell_type, ".GeneSet"), sep="\t", col.names=FALSE)
})