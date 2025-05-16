#!/bin/bash
#
#SBATCH --job-name=PySCENIC_all_NETs
#SBATCH --output=slurm_output_all_NETs.txt
#
#SBATCH --ntasks=12
#SBATCH --nodes=1
#SBATCH --time=80:00:00
#SBATCH --mem=120G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hilmar.berger@charite.de

OUTPUT_FOLDER=../Data/PySCENIC
INPUT_DATA=../Data/Processed/NET_all_NET_with_Tosti_complete.loom
RESOURCES=../Data/PySCENIC_refs

mkdir $OUTPUT_FOLDER

pyscenic grn --num_workers 8 \
         -o $OUTPUT_FOLDER/expr_mat.adjacencies.tsv \
          $INPUT_DATA \
          $RESOURCES/lambert2018.txt --sparse -m grnboost2

pyscenic ctx $OUTPUT_FOLDER/expr_mat.adjacencies.tsv \
         $RESOURCES//hg19-500bp-upstream-7species.mc9nr.feather \
         $RESOURCES/hg19-tss-centered-5kb-7species.mc9nr.feather \
         $RESOURCES/hg19-tss-centered-10kb-7species.mc9nr.feather \
        --annotations_fname $RESOURCES/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
        --expression_mtx_fname $INPUT_DATA \
        --mode "dask_multiprocessing" \
        --output $OUTPUT_FOLDER/regulons.csv \
        --num_workers 8

pyscenic aucell $INPUT_DATA \
                $OUTPUT_FOLDER/regulons.csv \
                -o $OUTPUT_FOLDER/auc_mtx.csv \
                --num_workers 8

