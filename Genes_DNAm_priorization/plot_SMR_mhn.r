#------------- plotMhnSMR ----//
library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)

dt = fread("ARHL_MVP_AJHG_BBJ_reformatMETAL_eQTLGen.smr")[, .(Gene, ProbeChr, Probe_bp, p_SMR, b_SMR, p_HEIDI)]
plt = dt[,.(Outco_Gene, Outco_Chr, Expo_bp, p_SMR)]
colnames(plt) = c("Gene", "CHR", "BP", "P")

Gene = dt[p_SMR < threshold & p_HEIDI > 0.01, ]
Gene$color = Gene[, ifelse(b_SMR > 0, brewer.pal(8,"Dark2")[1:2][2], brewer.pal(8,"Dark2")[1:2][1])]
Gene$pch = Gene[, ifelse(b_SMR > 0, 24, 25)]


plotMhnSMR <- function(plt, highlight=NULL, threshold=0.05/nrow(plt), ) {
  
  chr_lengths <- plt %>%
    group_by(CHR) %>%
    summarise(chr_len = max(BP))
  data <- plt %>%
    arrange(CHR, BP) %>%
    mutate(chr_cumsum = cumsum(c(0, head(chr_lengths$chr_len, -1)))[CHR]) %>%
    mutate(BPcum = BP + chr_cumsum)
  
  axis_set <- data %>%
    group_by(CHR) %>%
    summarise(center = (max(BPcum) + min(BPcum)) / 2)
  
  if (!is.null(highlight)){
    hlt = merge(highlight, data, by=("Gene"))
  }
  
}


hlt_dt = merge(Gene, data, by.x=c("Outco_Gene", "Expo_bp"), by.y=c("Gene", "BP"))
hlt_dt$pch = as.factor(hlt_dt$pch)

threshold = -log10(0.05 / nrow(dt))
p = ggplot() +
  geom_point(data = data, aes(x = BPcum, y = -log10(P), color = factor(CHR)), alpha = 0.9) +
  geom_point(data=hlt_dt, aes(x=BPcum, y=-log10(P), shape=pch),fill=hlt_dt$color, color=hlt_dt$color, alpha=0.9, size=2.5) +
  geom_text_repel(
    data = hlt_dt,
    aes(x=BPcum, y=-log10(P), label = Outco_Gene),
    color=hlt_dt$color,
    size = 3,
    fontface="italic",
    max.overlaps = Inf,
    max.iter = 10000, 
    force = 2,
    force_pull = 0.01,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50") + 
  geom_hline(yintercept=-log10(threshold), color="grey68", linetype="dashed") +
  scale_color_manual(values=rep(c("#bfbfbf", "#b3dfe7"), 22)) + 
  scale_shape_manual(values = c(24, 25)) + 
  scale_x_continuous(label=axis_set$CHR, breaks=axis_set$center, expand=expansion(mult = c(0, 0.01))) +
  scale_y_continuous(expand=expansion(mult = c(0, 0.01))) + 
  labs(x="Chromosome", y=expression(-log[10] (italic(P_SMR)))) +
  theme_classic() +
  theme(legend.position="none",
        axis.text.x = element_text(size=10, color="black"),
        axis.text.y = element_text(size=10, color="black"),
        axis.title.x = element_text(size=15, margin=margin(t=10), color="black"),
        axis.title.y = element_text(size=15, margin=margin(t=10), color="black"),
        axis.line = element_line(color="black"), 
        axis.ticks = element_line(color="black"))
p
ggsave("mQTL_mnh.png", p, height = 4, width = 10, dpi=600)
ggsave("mQTL_mnh.pdf", p, height = 4, width = 10)