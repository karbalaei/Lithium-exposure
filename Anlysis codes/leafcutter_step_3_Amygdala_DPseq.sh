#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=50G
#SBATCH --job-name=DS_leafcutter_Amygdala_DPseq
#SBATCH -c 1
#SBATCH -o ../logs/o_leafcutter_step_3_DS_Amygdala_DPseq.txt
#SBATCH -e ../logs/e_leafcutter_step_3_DS_Amygdala_DPseq.txt
#SBATCH --mail-type=ALL
#SBATCH --time=3-00:00:00


set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"



ml leafcutter

echo "DPseq_model Amygdala"


leafcutter_ds.R --num_threads 64  ../leafcutter/data/covariates/leafcutter_perind_numers.counts.gz ../data/modSep_Amygdala_for_leafcutter_DPseq.txt -o ../leafcutter/data/covariates_results/Amygdala_covariates_DPseq 



