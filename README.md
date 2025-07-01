# Lithium-exposure

## Workflow:

![Flowchart](https://github.com/karbalaei/Lithium-exposure/blob/main/graphs/Flowchart.jpg)

## Replication Steps:

### 1. Preprocessing

**File**: `./analysis/preprocessing/preprocessing_bpseq_model.R`

**Inputs**:

-   RSE gene, exon, junction, and transcript data

    ``` R
    load("../../../zandiHyde_bipolar_rnaseq/data/zandiHypde_bipolar_rseGene_n511.rda") load("../../../zandiHyde_bipolar_rnaseq/data/zandiHypde_bipolar_rseExon_n511.rda") load("../../../zandiHyde_bipolar_rnaseq/data/zandiHypde_bipolar_rseJxn_n511.rda") load("../../../zandiHyde_bipolar_rnaseq/data/zandiHypde_bipolar_rseTx_n511.rda")
    ```

-   qSV matrix from bpseq paper

    ``` R
    load("../../../zandiHyde_bipolar_rnaseq/case_control/qSV_mat.rda")
    ```

**Outputs:**

-    Preprocessed RSE gene, exon, junction, and transcript data

    ``` R
    save(rse_gene_lithium, file="../preprocessed_data/rse_gene_lithium_bpseq_model_processed_data.Rdata")
    save(rse_exon_lithium, file="../preprocessed_data/rse_exon_lithium_bpseq_model_processed_data.Rdata")
    save(rse_jxn_lithium, file="../preprocessed_data/rse_jxn_lithium_bpseq_model_processed_data.Rdata")
    save(rse_tx_lithium, file="../preprocessed_data/rse_tx_lithium_bpseq_model_processed_data.Rdata")
    ```

-   qSV matrix specific to this lithium study

    ``` R
    save(qSV_mat_lithium, file="../preprocessed_data/qSV_mat_lithium_bpseq_model_processed_data.Rdata")
    ```

### 2. Differential expression analysis

**File**: `./analysis/de_analysis/`

**Inputs**: 

**Outputs**:


