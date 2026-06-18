#!/bin/bash
#SBATCH -p shared
#SBATCH --mem=10G
#SBATCH --job-name=DS_leafcutter_step_4
#SBATCH -c 1
#SBATCH -o ../logs/o_leafcutter_step_4.txt
#SBATCH -e ../logs/e_leafcutter_step_4.txt
#SBATCH --mail-type=ALL
#SBATCH --time=3-00:00:00

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz  

mv gencode.v19.annotation.gtf.gz ../leafcutter/data/covariates_results/

chmod u+x gtf2leafcutter.pl
chmod u+x prepare_results.R


./gtf2leafcutter.pl -o ../leafcutter/data/covariates_results/gencode_hg19 ../leafcutter/data/covariates_results/gencode.v19.annotation.gtf.gz





