#===============================================================#
####     VLOOKUP in R and calculate FDR and plot LDSC cts    ####
#===============================================================#
#! /public/home/shilulu/anaconda3/envs/R4.2.0/bin/Rscript

data <- read.table('ARHL_MVP_AJHG_BBJ.cell_type_results.txt', header=TRUE, sep='\t')

data$FDR <- p.adjust(data$Coefficient_P_value, method='BH')
# save the FDR file
write.table(data, file='ARHL_MVP_AJHG_BBJ.cell_type_results.FDR.txt', sep='\t', row.name=FALSE, quote=FALSE)

library(ggplot2)
library(data.table)

data <- fread("ARHL_MVP_AJHG_BBJ.cell_type_results.FDR_plot.txt")
if (any(data$FDR <= 0.05)) {
  # 计算FDR <= 0.05的值中，对应的max p值，作为 fig 中的阈值线
  p_value_threshold <- max(data$Coefficient_P_value[data$FDR <= 0.05], na.rm=TRUE)
  
  log10_p_value_threshold <- -log10(p_value_threshold)
} else {
  log10_p_value_threshold <- NaN
}
# color for dif tissue
custom_colors <- c("Hair cells" = "#DD7694", 
                "Supporting cells" = "#BCD4E7", 
                "Surrounding structures" = "#056E83", 
                "Lateral wall" = "#E9D9BF", 
                "Circulating cells" = "#D4920A", 
                "Glial cells" = "#5AA4AE", 
                "Neurons" = "#65472F")   
# set Name as factor
data$Name <- factor(data$Name, levels=unique(data$Name))
data$Region <- factor(data$Region, levels=unique(data$Region))

p=ggplot(data, aes(x=Name, y=-log10(Coefficient_P_value), fill=Region)) +
  geom_bar(stat="identity", position="dodge", width=0.8)+
  geom_hline(yintercept=log10_p_value_threshold, linetype="dashed", color="red")+
  scale_fill_manual(values=custom_colors) +
  scale_y_continuous(expand=c(0,0), limits=c(0,3))+
  labs(x=NULL, y="-log10(P)", title="Cells from Mouse cochlear")+
  theme_bw()+
  theme(axis.text.x=element_text(size=9, angle=60, hjust=1, vjust=1),
        axis.text.y=element_text(size=9),
        legend.position="left",
        legend.title=element_blank(),
        # panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        strip.background=element_rect(colour=NA))

ggsave("Mouse_cochlear_ldsc-seg.png", p, width=8, height=4, dpi=500)
