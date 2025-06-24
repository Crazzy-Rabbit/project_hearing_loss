#===========================================#
#                                           #
# use opera to generate the omics plot file #
#                                           #
#===========================================#
opera="/public/home/wuyang/program/OPERA-main/opera_Linux"
bfile="/public/home/lijunpeng/smr-1.3.1-linux-x86_64/g1000_eur/g1000_eur"
mQTL="/public/home/gaikai/Age/data/mQTL/LBC_BSGS_meta/bl_mqtl_chr"
eQTLGen="/public/share/wchirdzhq2022/Wulab_share/SMR_xQTL_besd/xQTL/eQTL/besd_eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense"
Genelist="/public/home/shilulu/script/plot_smr/glist_hg19_refseq.txt"
gwas="/public/home/shilulu/project_hearing-loss/new_run/all_meta/ARHL_MVP_AJHG_BBJ_reformatMETAL"


# 17 ENSG00000004975 DVL2
# 17 ENSG00000170291 ELP5
# 19 ENSG00000178093 TSSK6
# 19 ENSG00000186010 NDUFA13
# 19 ENSG00000250067 YJEFN3
# 1 ENSG00000230896 RP11-767N6.7
# 6 ENSG00000220614 RP11-480N24.4

cat plot_probe | while read -r id; do
  arr=($id)
  chrs=${arr[0]}
  probe=${arr[1]}
  gene=${arr[2]}
  echo "${mQTL}${chrs}" > ${gene}_mlist
  cmd="$opera --bfile $bfile \
    --gwas-summary $gwas \
    --beqtl-summary $eQTLGen \
    --besd-flist ${gene}_mlist \
    --plot \
    --probe ${probe} \
    --probe-wind 2000 \
    --gene-list $Genelist \
    --out ${gene}_plot > ${gene}.opera.log 2>&1"
  qsubshcom "$cmd" 10 20G opera 90:00:00 ""
done

source("E:/Shi/ALL_of_my_Job/WCH_script/SMRplot/plot_OmicsSMR_xQTL.r")
setwd("E:/Shi/ALL_of_my_Job/24-28川大华西/2_project_hearing loss/process/SMR")
SMRData = ReadomicSMRData("ACADVL_plot.ENSG00000072778.txt")
pdf("ACADVL_opera.pdf", height=12, width=10)
omicSMRLocusPlot(data=SMRData, esmr_thresh=3.19e-06, msmr_thresh=5.37e-07, m2esmr_thresh=5.75e-04, m2esmr_heidi=0.01,
                window=200, anno_methyl=TRUE, highlight="rs507506", annoSig_only=TRUE, max_anno_probe=6,
                eprobeNEARBY="ENSG00000072778",mprobeNEARBY=c("cg03613822", "cg12805420"), 
                epi_plot=TRUE, funcAnnoFile="E:/Shi/ALL_of_my_Job/WCH_script/SMRplot/funcAnno.RData")
dev.off()
# epi_plot=TRUE if you want plot epigenome
# anno_methyl=TRUE if you want the mQTL can annote in plot
pdf("ACADVL_DVL2_ELP5.pdf", height=12, width=10)
omicSMRLocusPlot(data=SMRData, esmr_thresh=3.19e-06, msmr_thresh=5.37e-07, m2esmr_thresh=5.75e-04, m2esmr_heidi=0.01,
                window=200, anno_methyl=TRUE, annoSig_only=TRUE, max_anno_probe=6,
                eprobeNEARBY=c("ENSG00000072778","ENSG00000170291","ENSG00000004975"),mprobeNEARBY=c("cg12805420"), 
                )
dev.off()