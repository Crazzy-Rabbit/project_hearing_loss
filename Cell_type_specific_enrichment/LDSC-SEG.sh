conda activate ld
#***************************************************************#
## LDSC-SEG for Mouse cochleae in Philippe Jean et al. PNAS
#***************************************************************#
#! /bin/bash
SUMS="/public/home/shilulu/software/ldsc/munge_sumstats.py"
SNPlist="/public/home/lijunpeng/LDSC_LD_data/eur_w_ld_chr/w_hm3.snplist"
meta_file="/public/home/shilulu/project_hearing-loss/new_run/all_meta/ARHL_MVP_AJHG_BBJ_reformatMETAL.gz"
outdir="/public/home/shilulu/project_hearing-loss/new_run/all_meta/ldsc"

# step 1: covert to .sumstats.gz
cmd1="${SUMS} --sumstats ${meta_file} \
--merge-alleles ${SNPlist} \
--chunksize 500000 \
--a1 A1 \
--a2 A2 \
--out ${outdir}/ARHL_MVP_AJHG_BBJ_ldsc"
sub_id1=`qsubshcom "$cmd1" 1 100G sums_ldsc 90:00:00 ""`


#! /bin/bash                     
LDSC="/public/home/shilulu/software/ldsc/ldsc.py"
REF_LD_CHR="/public/share/wchirdzhq2022/Wulab_share/LDSC/1000G_EUR_Phase3_baseline/baseline."
REF_LD_CTS="/public/share/wchirdzhq2022/Wulab_share/LDSC/Mouse_cochleae/Mouse_cochleae_gene_expr.ldcts"
W_LD_CHR="/public/share/wchirdzhq2022/Wulab_share/LDSC/weights_hm3_no_hla/weights."
out_sums="/public/home/shilulu/project_hearing-loss/new_run/all_meta/ldsc"
out_ldscseg="/public/home/shilulu/project_hearing-loss/new_run/all_meta/ldsc/Mouse_cochleae"

cmd="${LDSC} --h2-cts ${out_sums}/ARHL_MVP_AJHG_BBJ_ldsc.sumstats.gz \
             --ref-ld-chr ${REF_LD_CHR} \
             --out ${out_ldscseg}/ARHL_MVP_AJHG_BBJ \
             --ref-ld-chr-cts ${REF_LD_CTS} \
             --w-ld-chr ${W_LD_CHR}"
qsubshcom "$cmd" 1 50G ldsc_cts 10:00:00 ""