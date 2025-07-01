library(jaffelab)
library(SummarizedExperiment)
library(sva)
library(readxl)
library(edgeR)
library(recount)

#load rse data
load("../../../zandiHyde_bipolar_rnaseq/data/zandiHypde_bipolar_rseGene_n511.rda")
load("../../../zandiHyde_bipolar_rnaseq/data/zandiHypde_bipolar_rseExon_n511.rda")
load("../../../zandiHyde_bipolar_rnaseq/data/zandiHypde_bipolar_rseJxn_n511.rda")
load("../../../zandiHyde_bipolar_rnaseq/data/zandiHypde_bipolar_rseTx_n511.rda")
load("../../../zandiHyde_bipolar_rnaseq/case_control/qSV_mat.rda")


#load("../zandiHyde_bipolar_rnaseq/data/degradation_rse_BipSeq_BothRegions.rda") #load cov_rse

## add ancestry 
load("../../../zandiHyde_bipolar_rnaseq/genotype_data/zandiHyde_bipolar_MDS_n511.rda")
mds = mds[rse_gene$BrNum,1:5]
colnames(mds) = paste0("snpPC", 1:5)
colData(rse_gene) = cbind(colData(rse_gene), mds)

# load lithium exposure data
lithiumdf = read.csv("../../lithium_data/lithium_confirmed_exposure_statuses.csv")

###########################
### BASIC PREPROCESSING ###
###########################

# factorize diagnoses
identical(colnames(rse_gene), colnames(cov_rse)) # TRUE
rse_gene$Dx = factor(ifelse(rse_gene$PrimaryDx == "Control", "Control","Bipolar"), 
                     levels = c("Control", "Bipolar"))

# get snp PCs
mds = mds[rse_gene$BrNum,1:5]
colnames(mds) = paste0("snpPC", 1:5)
colData(rse_gene) = cbind(colData(rse_gene), mds)

# get lithium exposure statuses 
lithiumdf$lithium = factor(lithiumdf$confirmed_lithium_exposed, levels = c(0, 1))
lithiumdf = lithiumdf[, c("BrNum", "lithium")]
rse_gene$lithium_group = NA
rse_gene$lithium_group = lithiumdf$lithium[match(rse_gene$BrNum, lithiumdf$BrNum)]

table(rse_gene$lithium_group , rse_gene$BrainRegion)
# Amygdala sACC
# 0       25   27
# 1        9   10

amyggroup = rse_gene[, rse_gene$BrainRegion == "Amygdala"]

table(amyggroup$Sex, amyggroup$lithium_group)
# 0  1
# F  9  4
# M 16  5
# amyggroup[,amyggroup$lithium_group == 0]

saccgroup = rse_gene[, rse_gene$BrainRegion == "sACC"]

table(saccgroup$Sex, saccgroup$lithium_group)
# 0  1
# F  9  5
# M 18  5


###########
# filter ##
## gene
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')
geneIndex = rowMeans(assays(rse_gene)$rpkm) > 0.25  ## both regions
rse_gene = rse_gene[geneIndex,]

## exon
assays(rse_exon)$rpkm = recount::getRPKM(rse_exon, 'Length')
exonIndex = rowMeans(assays(rse_exon)$rpkm) > 0.3
rse_exon = rse_exon[exonIndex,]

## junction
rowRanges(rse_jxn)$Length <- 100
assays(rse_jxn)$rp10m = recount::getRPKM(rse_jxn, 'Length')
jxnIndex = rowMeans(assays(rse_jxn)$rp10m) > 0.35 & rowData(rse_jxn)$Class != "Novel"
rse_jxn = rse_jxn[jxnIndex,]

## transcript
txIndex = rowMeans(assays(rse_tx)$tpm) > 0.4 
rse_tx = rse_tx[txIndex,]


###########################
### Generate qSV matrix ###
###########################
modJoint = model.matrix(~Dx*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 + 
                          mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr, 
                        data=colData(rse_gene))

degExprs = log2(assays(cov_rse)$count+1)
k = num.sv(degExprs, modJoint) # 19
qSV_mat = prcomp(t(degExprs))$x[,1:k]
varExplQsva = getPcaVars(prcomp(t(degExprs)))
varExplQsva[1:k]
sum(varExplQsva[1:k]) # 87.731%

# model w/o interaction to subset by region
modSep_lithium = model.matrix(~lithium_group + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 +
                                          mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr, 
                                        data=colData(rse_gene)) 


###############################################################
### Get rows from qSV matrix relevant to our patient sample ###
###############################################################


##### extract qSVs related to Lithium exposures #####

qSV_mat_lithium = qSV_mat[rownames(modSep_lithium),]

##### extract samples related to Lithium exposures #####

rse_gene_lithium = rse_gene[,rownames(modSep_lithium)]
rse_exon_lithium = rse_exon[,rownames(modSep_lithium)]
rse_jxn_lithium = rse_jxn[,rownames(modSep_lithium)]
rse_tx_lithium = rse_tx[,rownames(modSep_lithium)]


### Filtering data ####
##### Gene ######
assays(rse_gene_lithium)$rpkm = recount::getRPKM(rse_gene_lithium, 'Length')
geneIndex = rowMeans(assays(rse_gene_lithium)$rpkm) > 0.25  ## both regions
rse_gene_lithium = rse_gene_lithium[geneIndex,]

#####  Exon #####
assays(rse_exon_lithium)$rpkm = recount::getRPKM(rse_exon_lithium, 'Length')
exonIndex = rowMeans(assays(rse_exon_lithium)$rpkm) > 0.3
rse_exon_lithium = rse_exon_lithium[exonIndex,]


#### Junction #####
rowRanges(rse_jxn_lithium)$Length <- 100
assays(rse_jxn_lithium)$rp10m = recount::getRPKM(rse_jxn_lithium, 'Length')
jxnIndex = rowMeans(assays(rse_jxn_lithium)$rp10m) > 0.35 & rowData(rse_jxn_lithium)$Class != "Novel"
rse_jxn_lithium = rse_jxn_lithium[jxnIndex,]


#### Transcript ####
txIndex = rowMeans(assays(rse_tx_lithium)$tpm) > 0.4
rse_tx_lithium = rse_tx_lithium[txIndex,]


#####################
### SAVE OBJECTS ###
####################
#save(rse_gene_lithium, file="./processed_data/rse_gene_lithium_bpseq_model.Rdata")
#save(rse_exon_lithium, file="./processed_data/rse_exon_lithium_bpseq_model.Rdata")
#save(rse_jxn_lithium, file="./processed_data/rse_jxn_lithium_bpseq_model.Rdata")
#save(rse_tx_lithium, file="./processed_data/rse_tx_lithium_bpseq_model.Rdata")
#save(qSV_mat_lithium, file="./processed_data/qSV_mat_lithium_bpseq_model.Rdata")
