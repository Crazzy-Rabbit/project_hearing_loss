#=================================#
# # # start run smr and filter
# whole blood eQTL #
cd /public/home/shilulu/project_hearing-loss/new_run/all_meta/SMR/whole_blood

#! /bin/bash
eQTLGen="/public/share/wchirdzhq2022/Wulab_share/SMR_xQTL_besd/xQTL/eQTL/besd_eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense"
gwas="/public/home/shilulu/project_hearing-loss/new_run/all_meta/ARHL_MVP_AJHG_BBJ_reformatMETAL"


qsubshcom "Rscript run_SMRinR_filter.r --gwas ${gwas} --qtl ${eQTLGen} --out eQTLGen" 10 100G SMR 60:00:00 ""

#=====================================#
# whole blood methQTL #
#! /bin/bash 
mQTL="/public/share/wchirdzhq2022/Wulab_share/SMR_xQTL_besd/xQTL/mQTL/LBC_BSGS_meta/bl_mqtl_chr"
gwas="/public/home/shilulu/project_hearing-loss/new_run/all_meta/ARHL_MVP_AJHG_BBJ_reformatMETAL"
outdir="/public/home/shilulu/project_hearing-loss/new_run/all_meta/SMR/whole_blood/mQTL"
outprx=$(basename -- "$gwas")"_mQTL_chr"
smr=${outprx}_all.smr
# smr each chr
cmd1="Rscript run_SMRinR_filter.r --func run_smr_qtl --gwas ${gwas} --qtl ${mQTL}{TASK_ID} --out ${outdir}/${outprx}{TASK_ID}"
sub1=`qsubshcom "$cmd1" 10 100G mqtl_wb 60:00:00 "-array=1-22"`

# combined all chr result
sub2=`qsubshcom "awk 'NR==1 || FNR >1' ${outdir}/${outprx}{1..22}.smr > ${outprx}_all.smr" 1 20G awk 1:00:00 "-wait=${sub1}"`

# filter p_smr & p_heidi
cmd2="Rscript run_SMRinR_filter.r --func run_filter --smr ${smr} --out ${outprx}"
qsubshcom "$cmd2" 10 100G mqtl_flt 60:00:00 "-wait=${sub2}"
#=====================================#


#=====================================================================================#
#                                                                                     #
# two molecular traits test of whole blood methQTL and eQTL, only for significant QTL #
#                                                                                     #
#=====================================================================================#
# --beqtl-summary the first one reads mQTL summary data as the exposure. The second one reads eQTL summary data from as the outcome.

# step1: extract significant probeID from smr result
for file in ARHL_MVP_AJHG_BBJ_reformatMETAL_*smr_0.01heidi.txt; do
    qtl=$(echo $file | sed -n 's/.*reformatMETAL_\([^_]*\)_.*/\1/p')
    Rscript -e "suppressWarnings(suppressMessages(library(data.table))); 
                file=fread('${file}')[,.(probeID)]; 
                setwd('/public/home/shilulu/project_hearing-loss/new_run/all_meta/SMR/whole_blood/');
                fwrite(file, 'eQTLGen_probeID.txt', col.names=FALSE, sep='\t')"
done

# step2: run smr of mqtl to eqtl 
#! /bin/bash
SMR="/public/home/lijunpeng/smr-1.3.1-linux-x86_64/smr-1.3.1"
bfile="/public/home/lijunpeng/smr-1.3.1-linux-x86_64/g1000_eur/g1000_eur"
mQTL="/public/share/wchirdzhq2022/Wulab_share/SMR_xQTL_besd/xQTL/mQTL/LBC_BSGS_meta/bl_mqtl_chr"
eQTLGen="/public/share/wchirdzhq2022/Wulab_share/SMR_xQTL_besd/xQTL/eQTL/besd_eQTLGen/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense"
outdir="/public/home/shilulu/project_hearing-loss/new_run/all_meta/SMR/whole_blood/m2eqtl"

cmd="${SMR} --bfile ${bfile} \
    --beqtl-summary ${mQTL}{TASK_ID} \
    --beqtl-summary ${eQTLGen} \
    --extract-exposure-probe mqtl_probeID.txt \
    --extract-outcome-probe ${probe} \
    --diff-freq-prop 0.5 \
    --thread-num 10 \
    --out ${outdir}/m2e_chr{TASK_ID} > ${outdir}/m2e_chr{TASK_ID}.log 2>&1"
qsubshcom "$cmd" 1 100G smr 60:00:00 "-array=1-22"
  # combined all chr 
awk 'NR==1 || FNR >1' ${outdir}/m2e_chr*.smr > m2e_all_chr.smr