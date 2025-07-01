# Lithium-exposure

## Workflow:

![Flowchart](https://github.com/karbalaei/Lithium-exposure/blob/main/graphs/Flowchart.jpg)

--- 
## Analysis Replication Steps:

### 1. Preprocessing

**File**: `./analysis/preprocessing/preprocessing_bpseq_model.R`

**Inputs**:

-   RSE gene, exon, junction, and transcript data

    ``` r
    load("../../../zandiHyde_bipolar_rnaseq/data/zandiHypde_bipolar_rseGene_n511.rda") load("../../../zandiHyde_bipolar_rnaseq/data/zandiHypde_bipolar_rseExon_n511.rda") load("../../../zandiHyde_bipolar_rnaseq/data/zandiHypde_bipolar_rseJxn_n511.rda") load("../../../zandiHyde_bipolar_rnaseq/data/zandiHypde_bipolar_rseTx_n511.rda")
    ```

-   qSV matrix from bpseq paper

    ``` r
    load("../../../zandiHyde_bipolar_rnaseq/case_control/qSV_mat.rda")
    ```

**Outputs:**

-   Preprocessed RSE gene, exon, junction, and transcript data

````         
``` R
save(rse_gene_lithium, file="../preprocessed_data/rse_gene_lithium_bpseq_model_processed_data.Rdata")
save(rse_exon_lithium, file="../preprocessed_data/rse_exon_lithium_bpseq_model_processed_data.Rdata")
save(rse_jxn_lithium, file="../preprocessed_data/rse_jxn_lithium_bpseq_model_processed_data.Rdata")
save(rse_tx_lithium, file="../preprocessed_data/rse_tx_lithium_bpseq_model_processed_data.Rdata")
```
````

-   qSV matrix specific to this lithium study

    ``` r
    save(qSV_mat_lithium, file="../preprocessed_data/qSV_mat_lithium_bpseq_model_processed_data.Rdata")
    ```

### 2. Differential expression analysis

**File**: `./analysis/de_analysis/de_qsv_BPseq_model.R`

**Inputs**:

-   Preprocessed and filtered RSE gene, exon, junction, and transcript data

    ``` R
    load("../preprocessed_data/rse_gene_lithium_bpseq_model_processed_data.Rdata")
    load("../preprocessed_data/rse_exon_lithium_bpseq_model_processed_data.Rdata")
    load("../preprocessed_data/rse_jxn_lithium_bpseq_model_processed_data.Rdata")
    load("../preprocessed_data/rse_tx_lithium_bpseq_model_processed_data.Rdata")
    ```

-   qSV matrix created from preprocessing

    ``` R
    load("../preprocessed_data/qSV_mat_lithium_bpseq_model_processed_data.Rdata")
    ```

**Outputs**:

-   P values of features, separated by gene/exon/junction/transcript level and brain region. Dated by run time.

    ``` R
    save(
      Lithium_BPseqmodel, 
      file = paste0("../results/results_bbpseqmodel_de_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
      )
    ```
---
## Plotting Replication Steps:

### 1. Differential expression plots

**File**: `./analysis/de_analysis/de_qsv_BPseq_model.R`

**Inputs**:

-   Preprocessed and filtered RSE gene, exon, junction, and transcript data

    ``` R
    load("../preprocessed_data/rse_gene_lithium_bpseq_model_processed_data.Rdata")
    load("../preprocessed_data/rse_exon_lithium_bpseq_model_processed_data.Rdata")
    load("../preprocessed_data/rse_jxn_lithium_bpseq_model_processed_data.Rdata")
    load("../preprocessed_data/rse_tx_lithium_bpseq_model_processed_data.Rdata")
    ```

-   qSV matrix created from preprocessing

    ``` R
    load("../preprocessed_data/qSV_mat_lithium_bpseq_model_processed_data.Rdata")
    ```

**Outputs**:

-   P values of features, separated by gene/exon/junction/transcript level and brain region. Dated by run time.

    ``` R
    save(
      Lithium_BPseqmodel, 
      file = paste0("../results/results_bbpseqmodel_de_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
      )
    ```