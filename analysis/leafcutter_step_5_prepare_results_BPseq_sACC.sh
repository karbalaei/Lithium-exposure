#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=20G
#SBATCH --job-name=DS_leafcutter_sACC_BPseq
#SBATCH -c 1
#SBATCH -o ../logs/o_leafcutter_step_5_Bpseq_sACC.txt
#SBATCH -e ../logs/e_leafcutter_step_5_Bpseq_sACC.txt
#SBATCH --mail-type=ALL
#SBATCH --time=3-00:00:00


ml leafcutter

./prepare_results.R -o ../leafcutter/data/covariates_results/step_5_BPseq_sACC.Rdata \
             -m ../data/modSep_sACC_for_leafcutter_BPseq.txt -f 0.05 -c sACC_BPseq  \
             ../leafcutter/data/covariates/leafcutter_perind_numers.counts.gz \
             ../leafcutter/data/covariates_results/sACC_covariates_BPseq_cluster_significance.txt \
             ../leafcutter/data/covariates_results/sACC_covariates_BPseq_effect_sizes.txt \
             ../leafcutter/data/covariates_results/gencode_hg19
