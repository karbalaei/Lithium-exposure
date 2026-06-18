# Differentially expressed features

To elucidate region-specific molecular profiles within this cohort, we conducted transcriptomic differential expression analyses on postmortem tissue specimens from individuals with Bipolar Disorder, stratified by documented lithium exposure. Analyses were performed across two anatomically distinct brain regions crucial to emotional regulation: the amygdala and the subgenual anterior cingulate cortex (sACC). .

After loading preprocessed data,  low expressed featured were removed from the analysis. Then using model.matrix function and some qSVs data from previous study , model was created to compare two groups of data.

The comparison was done by applying few functions which more important ones and their application are as follow:

[calcNormFactors](https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/calcNormFactors): Calculates scaling factors to normalize for differences in library sizes between samples.

[voom](https://rdrr.io/bioc/limma/man/voom.html): Converts raw RNA-seq counts to log-CPM values and estimates precision weights to handle the mean-variance relationship.

[lmFit](https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/lmFit): Fits a linear model to the expression data for each gene/exon to prepare for differential expression testing.