#--------------------------- GSMR ----//
#
# GSMR for some traits select from BADGERS
# 
#--------------------------- GSMR ----//
dir="/public/home/shilulu/Wulab_project/ARHL/GSMR/cleanGWAS"
echo "Loneliness $dir/clean_Loneliness_GCST006924.txt.gz
Insomnia $dir/clean_Sleeplessness-insomnia.fastGWA.gz
LongStandIllness $dir/clean_Long-standing_illness.fastGWA.gz
Overall_health_rating $dir/clean_Overall_health_rating.fastGWA.gz
Snoring $dir/clean_Snoring.fastGWA.gz
Neuroticism $dir/clean_Neuroticism_score.fastGWA.gz
Tiredness $dir/clean_tiredness_lethargy2080.fastGWA.gz
Leg_fat_percentage $dir/clean_leg_fat_percentage_right.fastGWA.gz
Noisy_workplace $dir/clean_Noisy_workplace.fastGWA.gz
Taking_other_prescription_medications $dir/clean_Taking_other_prescription_medications.fastGWA.gz
derpess_mood $dir/clean_depressed_mood.fastGWA.gz
Ever_smoke $dir/clean_ever_smoke.fastGWA.gz" > exposure_new.txt

echo "ARHL /public/home/shilulu/Wulab_project/ARHL/14May2025/Demeta_raw_hm3_Meta_rescale_keeponly_SNP_14May2025.txt" > outcome.txt

#=============================================#
#! /bin/bash 
gcta="/public/home/wuyang/bin/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1"
ref="/public/share/wchirdzhq2022/Wulab_share/1000GenomePhase3_Ref_hg37/g1000_eur/g1000_eur"
exposure="/public/home/shilulu/Wulab_project/ARHL/GSMR/cleanGWAS/exposure_new.txt"
outcome="/public/home/shilulu/Wulab_project/ARHL/GSMR/cleanGWAS/outcome.txt"

cmd="${gcta} --bfile ${ref} \
--gsmr-file ${exposure} ${outcome} \
--gsmr-direction 2 \
--gwas-thresh 5e-08 \
--diff-freq 0.5 \
--clump-r2 0.05 \
--gsmr2-beta \
--gsmr-snp-min 10 \
--heidi-thresh 0.01 \
--effect-plot \
--out result/ARHL_gsmr_new \
--thread-num 20"
qsubshcom "$cmd" 20 100G gsmr 90:00:00 ""