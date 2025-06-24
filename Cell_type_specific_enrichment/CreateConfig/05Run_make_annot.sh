# use make_annot.py to generate annot.gz file 
# using 100K widow as Hilary K. Finucane et al. Nat Genet
conda activate ldsc
bfile="/public/share/wchirdzhq2022/Wulab_share/LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC"
ldsc="/public/home/shilulu/software/ldsc"

ls *.GeneSet | while read id; do
  annot=$(basename -- "${id}" ".GeneSet")
  # for chr in {1..22}; do 
  cmd="python ${ldsc}/make_annot.py \
      --gene-set-file ${annot}.GeneSet \
      --gene-coord-file ${ldsc}/ENSG_coord.txt \
      --bimfile ${bfile}.{TASK_ID}.bim \
      --windowsize 100000 \
      --annot-file ${annot}.{TASK_ID}.annot.gz"
  # done 
  qsubshcom "$cmd" 1 10G ldsc_anot 1:00:00 "-array=1-22"
done

## Step 2: Computing LD scores with an annot file
bfile="/public/share/wchirdzhq2022/Wulab_share/LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC"
ldsc="/public/home/shilulu/software/ldsc"
ldsc_dir="/public/share/wchirdzhq2022/Wulab_share/LDSC"

awk '{if ($1!="SNP") {print $1} }' ${ldsc_dir}/w_hm3.snplist > ${ldsc_dir}/listHM3.txt

ls *.GeneSet | while read id; do
  annot=$(basename -- "${id}" ".GeneSet")
  cmd="python ${ldsc}/ldsc.py --l2 \
    --bfile ${bfile}.{TASK_ID} \
    --print-snps ${ldsc_dir}/listHM3.txt \
    --ld-wind-cm 1 \
    --annot ${annot}.{TASK_ID}.annot.gz \
    --thin-annot \
    --out ${annot}.{TASK_ID}"
  qsubshcom "$cmd" 1 10G ldsc_l2 1:00:00 "-array=1-22"
done

# generate .ldcts file for my sc data
dir=`pwd`
ls *GeneSet | while read id; do
  tissue=$(basename -- ${id} | sed 's/^Sc_P20_//;s/\.GeneSet$//')
  if [[ "$tissue" == "control" ]]; then
    echo "control not to be the tissue or cell type"
  else
    cell=${dir}/$(basename -- ${id} "GeneSet")
    control=${dir}/"Sc_P20_control."
    echo "${tissue} ${cell},${control}" >> Mouse_cochleae_P20_gene_expr.ldcts
  fi
done
# rm control file in cell 
