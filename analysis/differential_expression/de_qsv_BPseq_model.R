### !!! README: RUN THIS FILE FROM THE Lithium-exposure/ DIRECTORY !!! ###

### Library ####

library(SummarizedExperiment)
library(jaffelab)
library(recount)
library(sva)
library(edgeR)
library(dplyr)
library(readr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(purrr)
library(here)
library(enrichplot)
library(ggnewscale)
library(sessioninfo)
library("KEGGREST")
library(vctrs)
library(tibble)
library(stringr)

### Load preprocessed data ###
### (output of preprocessing/preprocessing_bpseq_model.R) ###

load("../preprocessed_data/rse_gene_lithium_bpseq_model_processed_data.Rdata")
load("../preprocessed_data/rse_exon_lithium_bpseq_model_processed_data.Rdata")
load("../preprocessed_data/rse_jxn_lithium_bpseq_model_processed_data.Rdata")
load("../preprocessed_data/rse_tx_lithium_bpseq_model_processed_data.Rdata")
load("../preprocessed_data/qSV_mat_lithium_bpseq_model_processed_data.Rdata")


### Differentially expressed feature analysis by region #####

modSep_lithium  = model.matrix(
  ~lithium_group + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
    mitoRate + rRNA_rate + totalAssignedGene + RIN + abs(ERCCsumLogErr), 
  data=colData(rse_gene_lithium)
) 


## split by region 

sACC_Index_lithium = which(colData(rse_gene_lithium)$BrainRegion == "sACC")
mod_sACC_lithium = cbind(modSep_lithium[sACC_Index_lithium,], qSV_mat_lithium[sACC_Index_lithium, ])


Amyg_Index_lithium = which(colData(rse_gene_lithium)$BrainRegion == "Amygdala")
mod_Amyg_lithium = cbind(modSep_lithium[Amyg_Index_lithium,], qSV_mat_lithium[Amyg_Index_lithium, ])

####  Gene ####

###### sACC #######
dge_sACC = DGEList(counts = assays(rse_gene_lithium[ ,sACC_Index_lithium])$counts, 
                   genes = rowData(rse_gene_lithium))

dge_sACC = calcNormFactors(dge_sACC)

vGene_sACC = voom(dge_sACC,mod_sACC_lithium, plot=F)


fitGene_sACC = lmFit(vGene_sACC)
eBGene_sACC = eBayes(fitGene_sACC)
outGene_sACC = topTable(eBGene_sACC,coef=2,
                        p.value = 1,number=nrow(rse_gene_lithium) )
outGene_sACC = outGene_sACC[rownames(rse_gene_lithium),]

sum(outGene_sACC$adj.P.Val < 0.05) # 0
sum(outGene_sACC$P.Value < 0.005) # 50

outGene_sACC <-  outGene_sACC %>% dplyr::rename("common_gene_id" = gencodeID)


##### Amygdala #####
dge_Amyg = DGEList(counts = assays(rse_gene_lithium[,Amyg_Index_lithium])$counts, 
                   genes = rowData(rse_gene_lithium))


dge_Amyg = calcNormFactors(dge_Amyg)

vGene_Amyg = voom(dge_Amyg,mod_Amyg_lithium, plot=F)

fitGene_Amyg = lmFit(vGene_Amyg)
eBGene_Amyg = eBayes(fitGene_Amyg)
outGene_Amyg = topTable(eBGene_Amyg,coef=2,
                        p.value = 1,number=nrow(rse_gene_lithium))
outGene_Amyg = outGene_Amyg[rownames(rse_gene_lithium),]

sum(outGene_Amyg$adj.P.Val < 0.05) #0
sum(outGene_Amyg$P.Value < 0.005) #237

outGene_Amyg <-  outGene_Amyg %>% dplyr::rename("common_gene_id" = gencodeID)


##### core output ####

nam = c("common_gene_id", "logFC", "AveExpr","t", "P.Value", "adj.P.Val", "B")

geneOut = merge(outGene_Amyg[,nam], outGene_sACC[,nam], by ='row.names', all = TRUE) %>% column_to_rownames("Row.names")

colnames(geneOut) <-  colnames(geneOut) %>% str_replace_all("\\.x" , "\\_Amygdala") %>% str_replace_all("\\.y" , "\\_sACC")

write.csv(geneOut ,  here("results" , "Lithium_BPseq_model_gene.csv"))

geneOut_list <-  list(outGene_Amyg , outGene_sACC)

names(geneOut_list) <-  c("Amygdala" , "sACC")


#### Exon ####

##### sACC ######
dee_sACC = DGEList(counts = assays(rse_exon_lithium[,sACC_Index_lithium])$counts, 
                   genes = rowData(rse_exon_lithium))
dee_sACC = calcNormFactors(dee_sACC)

vExon_sACC = voom(dee_sACC,mod_sACC_lithium, plot=F)

fitExon_sACC = lmFit(vExon_sACC)
eBExon_sACC = eBayes(fitExon_sACC)
outExon_sACC = topTable(eBExon_sACC,coef=2,
                        p.value = 1,number=nrow(rse_exon_lithium))
outExon_sACC = outExon_sACC[rownames(rse_exon_lithium),]
sum(outExon_sACC$adj.P.Val < 0.05) #1
sum(outExon_sACC$P.Value < 0.005) #995

outExon_sACC <-  outExon_sACC %>% dplyr::rename("common_gene_id" = gencodeID)

##### Amygdala ######

dee_Amyg = DGEList(counts = assays(rse_exon_lithium[,Amyg_Index_lithium])$counts, 
                   genes = rowData(rse_exon_lithium))
dee_Amyg = calcNormFactors(dee_Amyg)

vExon_Amyg = voom(dee_Amyg,mod_Amyg_lithium, plot=F)

fitExon_Amyg = lmFit(vExon_Amyg)
eBExon_Amyg = eBayes(fitExon_Amyg)
outExon_Amyg = topTable(eBExon_Amyg,coef=2,
                        p.value = 1,number=nrow(rse_exon_lithium))
outExon_Amyg = outExon_Amyg[rownames(rse_exon_lithium),]
sum(outExon_Amyg$adj.P.Val < 0.05) #743
sum(outExon_Amyg$P.Value < 0.005) #4062

outExon_Amyg <-  outExon_Amyg %>% dplyr::rename("common_gene_id" = gencodeID)


##### core output #####

nam = c("common_gene_id" , "logFC", "AveExpr","t", "P.Value", "adj.P.Val", "B")

exonOut = merge(outExon_Amyg[,nam], outExon_sACC[,nam], by ='row.names', all = TRUE) %>% column_to_rownames("Row.names")

colnames(exonOut) <-  colnames(exonOut) %>% str_replace_all("\\.x" , "\\_Amygdala") %>% str_replace_all("\\.y" , "\\_sACC")

write.csv(exonOut , here("results", "Lithium_BPseq_model_exon.csv"))


exonOut_list <-  list(outExon_Amyg , outExon_sACC)

names(exonOut_list) <-  c("Amygdala" , "sACC")


####  Junction ####


##### sACC ######
dje_sACC = DGEList(counts = assays(rse_jxn_lithium[,sACC_Index_lithium])$counts, 
                   genes = rowData(rse_jxn_lithium))
dje_sACC = calcNormFactors(dje_sACC)

vJxn_sACC = voom(dje_sACC,mod_sACC_lithium, plot=F)

fitJxn_sACC = lmFit(vJxn_sACC)
eBJxn_sACC = eBayes(fitJxn_sACC)
outJxn_sACC = topTable(eBJxn_sACC,coef=2,
                       p.value = 1,number=nrow(rse_jxn_lithium))
outJxn_sACC = outJxn_sACC[rownames(rse_jxn_lithium),]
sum(outJxn_sACC$adj.P.Val < 0.05) #  3 
sum(outJxn_sACC$P.Value < 0.005) #  1238


outJxn_sACC <-  outJxn_sACC %>% dplyr::rename("common_gene_id" = newGeneID , Symbol_General=newGeneSymbol)

##### Amygdala ######
dje_Amyg = DGEList(counts = assays(rse_jxn_lithium[,Amyg_Index_lithium])$counts, 
                   genes = rowData(rse_jxn_lithium))
dje_Amyg = calcNormFactors(dje_Amyg)

vJxn_Amyg = voom(dje_Amyg,mod_Amyg_lithium, plot=F)

fitJxn_Amyg = lmFit(vJxn_Amyg)
eBJxn_Amyg = eBayes(fitJxn_Amyg)
outJxn_Amyg = topTable(eBJxn_Amyg,coef=2,
                       p.value = 1,number=nrow(rse_jxn_lithium))
outJxn_Amyg = outJxn_Amyg[rownames(rse_jxn_lithium),]
sum(outJxn_Amyg$adj.P.Val < 0.05) # 479
sum(outJxn_Amyg$P.Value < 0.005) # 2894


outJxn_Amyg <-  outJxn_Amyg %>% dplyr::rename("common_gene_id" = newGeneID , Symbol_General=newGeneSymbol)


##### core output #####
nam = c("common_gene_id", "Symbol_General",  "logFC", "AveExpr","t", "P.Value", "adj.P.Val", "B")

jxnOut = merge(outJxn_Amyg[,nam], outJxn_sACC[,nam], by ='row.names', all = TRUE) %>% column_to_rownames("Row.names")

colnames(jxnOut) <-  colnames(jxnOut) %>% str_replace_all("\\.x" , "\\_Amygdala") %>% str_replace_all("\\.y" , "\\_sACC")

write.csv(jxnOut , here("results", "Lithium_BPseq_model_Jxn.csv"))


jxnOut_list <-  list(outJxn_Amyg , outJxn_sACC)

names(jxnOut_list) <-  c("Amygdala" , "sACC")


####  Transcript ####

txExprs = log2(assays(rse_tx_lithium)$tpm+ 1)

##### sACC ######

fitTx_sACC = lmFit(txExprs[,sACC_Index_lithium], mod_sACC_lithium)

eBTx_sACC = eBayes(fitTx_sACC)
outTx_sACC = topTable(eBTx_sACC,coef=2,
                      p.value = 1,number=nrow(rse_tx_lithium), 
                      genelist = rowRanges(rse_tx_lithium))
outTx_sACC = outTx_sACC[rownames(rse_tx_lithium),]
sum(outTx_sACC$adj.P.Val < 0.05)  ## 0
sum(outTx_sACC$P.Value < 0.005)  ## 316


outTx_sACC <-  outTx_sACC %>% dplyr::rename("common_gene_id" = ID.gene_id , Symbol_General= "ID.gene_name")


##### Amygdala ######

fitTx_Amyg = lmFit(txExprs[,Amyg_Index_lithium], mod_Amyg_lithium)
eBTx_Amyg = eBayes(fitTx_Amyg)
outTx_Amyg = topTable(eBTx_Amyg,coef=2,
                      p.value = 1,number=nrow(rse_tx_lithium),
                      genelist = rowRanges(rse_tx_lithium))
outTx_Amyg = outTx_Amyg[rownames(rse_tx_lithium),]
sum(outTx_Amyg$adj.P.Val < 0.05)  ## 0
sum(outTx_Amyg$P.Value < 0.005)  ## 291


outTx_Amyg <-  outTx_Amyg %>% dplyr::rename("common_gene_id" = ID.gene_id , Symbol_General= "ID.gene_name")

##### core output #####
nam = c("common_gene_id" , "Symbol_General" ,   "logFC", "AveExpr","t", "P.Value", "adj.P.Val", "B")


txOut = merge(outTx_Amyg[,nam], outTx_sACC[,nam], by ='row.names', all = TRUE) %>% column_to_rownames("Row.names")

colnames(txOut) <-  colnames(txOut) %>% str_replace_all("\\.x" , "\\_Amygdala") %>% str_replace_all("\\.y" , "\\_sACC")


write.csv(txOut , here("results", "Lithium_BPseq_model_Transcript.csv"))

txOut_list <-  list(outTx_Amyg , outTx_sACC)

names(txOut_list) <-  c("Amygdala" , "sACC")


Lithium_BPseqmodel <-  list("gene" = geneOut_list, "exon" = exonOut_list, "jxn" = jxnOut_list ,"tx" = txOut_list)


save(
  Lithium_BPseqmodel, 
  file = paste0("../results/results_bbpseqmodel_de_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
  )

## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()


