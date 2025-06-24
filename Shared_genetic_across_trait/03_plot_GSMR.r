#========================================================#
### delete row contained "nan"
#========================================================#
#! /public/home/shilulu/anaconda3/envs/R4.2.0/bin/Rscript
library(data.table)
library(dplyr)

gsmr <- fread("ARHL_gsmr_new.gsmr")
gsmr[gsmr == "nan"] <- NA

# p < 0.05 / row number / 2 for 2 dirction gsmr
# remove row contained "nan" (all)
final <- gsmr %>%
  filter(p < 0.05 / (nrow(gsmr))) %>%
  filter(complete.cases(.))

fwrite(final, file="ARHL_gsmr_new_noNA_anno.gsmr", sep="\t", na="nan", quote=FALSE)

#### 2. forest plot
setwd("F:/Github/PHD_job/2_project_hearing loss/process/badgers/gsmr")

library(ggplot2)
library(data.table)
library(patchwork)
data = fread("ARHL_gsmr_new.gsmr")

data$lower=data$bxy-1.96*data$se
data$upper=data$bxy+1.96*data$se


data$FDR <- p.adjust(data$p, method='BH')
final = data[FDR < 0.05, ]

final$fill = ifelse(final$p < 0.05 / 24, "type1", "type2")

final$Exposure = factor(final$Exposure, levels = rev(unique(final$Exposure)))
final$Outcome = factor(final$Outcome, levels = rev(unique(final$Outcome)))

data1 = final[Exposure != "ARHL",]
data2 = final[Exposure == "ARHL",]


p1 = ggplot(data=data1, aes(y=Exposure, x=bxy)) +
  geom_point(aes(fill = fill), shape=21, size=3, color="#1B9E77") +
  scale_fill_manual(values = c(type1 = "#1B9E77", type2 = "white")) + 
  geom_errorbar(aes(xmin=lower, xmax=upper), width=0.2, linewidth=0.8, color="#1B9E77")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"), 
        panel.background = element_blank(),
        axis.title = element_blank(),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.position = "none")
p1

p2 = ggplot(data=data2, aes(y=Outcome, x=bxy)) +
  geom_point(aes(fill = fill), shape=21, size=3, color="#D95F02") +
  scale_fill_manual(values = c(type1 = "#D95F02", type2 = "white")) + 
  geom_errorbar(aes(xmin=lower, xmax=upper), width=0.2, linewidth=0.8, color="#D95F02")+
  geom_vline(xintercept = 0, linetype = "dashed")+
  theme_bw() + 
  theme(axis.text.x = element_text( size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"), 
        panel.background = element_blank(),
        axis.title = element_blank(),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),
        legend.position = "none")
p2


p = p1 + p2 + plot_layout(widths = c(2, 2))
p

ggsave("gsmr_plot.png", p, dpi=600, width=9, height=3)
ggsave("gsmr_plot.pdf", p, width=9, height=3)