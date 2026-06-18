# Differentially Expressed Features

To elucidate region-specific molecular profiles within this cohort, we conducted transcriptomic differential expression analyses on postmortem tissue specimens from individuals with Bipolar Disorder, stratified by documented lithium exposure. Analyses were performed across two anatomically distinct brain regions crucial to emotional regulation: the amygdala and the subgenual anterior cingulate cortex (sACC).

## Methodological Workflow

Following the retrieval of preprocessed count data, lowly expressed features were filtered out to improve statistical power. A design matrix was then constructed using the `model.matrix()` function, incorporating quality surrogate variables (qSVs) derived from previous benchmarks to rigorously control for unobserved technical covariates and batch effects.

The differential expression pipeline was executed using the standard `edgeR`/`limma` framework, leveraging the following core functions:

* [**`calcNormFactors`**](https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/calcNormFactors): Calculates scaling factors to normalize for differences in library sizes between samples.
* [**`voom`**](https://rdrr.io/bioc/limma/man/voom.html): Converts raw RNA-seq counts to $\text{log-CPM}$ values and estimates observational-level precision weights to model the mean-variance relationship.
* [**`lmFit`**](https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/lmFit): Fits a feature-wise linear model to the expression data, preparing the dataset for empirical Bayes moderation and differential expression testing.