#!/bin/sh

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.annotation.gtf.gz  

mv gencode.v48.annotation.gtf.gz ../leafcutter/data/covariates_results/

chmod u+x gtf2leafcutter.pl
chmod u+x prepare_results.R


./gtf2leafcutter.pl -o ../leafcutter/data/covariates_results/gencode_hg19 ../leafcutter/data/covariates_results/gencode.v48.annotation.gtf.gz





