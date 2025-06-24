#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH --job-name=DS_leafcutter_Amygdala_DPseq
#SBATCH -c 1
#SBATCH -o ../logs/o_leafcutter_step_5_DPseq_Amygdala.txt
#SBATCH -e ../logs/e_leafcutter_step_5_DPseq_Amygdala.txt
#SBATCH --mail-type=ALL
#SBATCH --time=3-00:00:00


ml leafcutter

./prepare_results.R -o ../leafcutter/data/covariates_results/step_5_DPseq_Amygdala.Rdata \
             -m ../data/modSep_Amygdala_for_leafcutter_DPseq.txt -f 0.05 -c Amygdala_DPseq  \
             ../leafcutter/data/covariates/leafcutter_perind_numers.counts.gz \
             ../leafcutter/data/covariates_results/Amygdala_covariates_DPseq_cluster_significance.txt \
             ../leafcutter/data/covariates_results/Amygdala_covariates_DPseq_effect_sizes.txt \
             ../leafcutter/data/covariates_results/gencode_hg19
