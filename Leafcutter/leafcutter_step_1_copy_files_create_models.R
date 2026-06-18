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
library(tidyverse)

options(scipen=999)


load(here::here("data", "zandiHypde_bipolar_rseGene_n511.rda"))
load(here("data","degradation_rse_BipSeq_BothRegions.rda")) #load cov_rse

identical(colnames(rse_gene), colnames(cov_rse)) # TRUE

rse_gene$Dx = factor(ifelse(rse_gene$PrimaryDx == "Control", "Control","Bipolar"), 
                     levels = c("Control", "Bipolar"))



##### add ancestry  ####
load(here("data","zandiHyde_bipolar_MDS_n511.rda")) # load mds


mds = mds[rse_gene$BrNum,1:5]
colnames(mds) = paste0("snpPC", 1:5)
colData(rse_gene) = cbind(colData(rse_gene), mds)

#### load lithium exposure ####

drug_info = read.delim(here::here("lithium","drug_info.csv"),as.is=TRUE , sep = ",")


rse_gene$lithium = drug_info$lithium[match(rse_gene$BrNum, drug_info$BrNum)]
rse_gene$lifetime_lithium = drug_info$life_lithium[match(rse_gene$BrNum, drug_info$BrNum)]
rse_gene$lithium_group = NA
rse_gene$lithium_group[which(rse_gene$lithium==1)] ="Lithium"
rse_gene$lithium_group[which(rse_gene$lithium==0 & rse_gene$lifetime_lithium==0)] ="Control"
rse_gene$lithium_group <-  factor(rse_gene$lithium_group)
table(rse_gene$lithium_group , rse_gene$BrainRegion)

# Amygdala sACC
# 0       25   27
# 1        9   10


#### qSVs #####

regions <- list(sacc = "sACC", amyg = "Amygdala")

rse_gene_lithium <-  rse_gene[ , rse_gene$lithium_group %in% c("Control" , "Lithium")  ]

rse_gene_lithium$lithium_group <-  factor(rse_gene_lithium$lithium_group)


bam_files <- gsub("dcl01/lieber/ajaffe/lab","dcs04/lieber/lcolladotor/hydeZandi_LIBD2005", rse_gene_lithium$bamFile)

leafcutter_ids = gsub("_accepted_hits.sorted.bam" , "_junctions_primaryOnly_regtools.bed" , 
		gsub("/dcs04/lieber/lcolladotor/hydeZandi_LIBD2005/zandiHyde_bipolar_rnaseq/preprocessed_data/HISAT2_out/", "" , bam_files))



rse_gene_lithium$leafcutter_ids = leafcutter_ids

pd = colData(rse_gene_lithium) %>% as.data.frame() %>%  rownames_to_column("id")

###### DPseq model #####


modJoint_DPseq = model.matrix(~Dx*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
                          mitoRate + rRNA_rate + totalAssignedGene + RIN +  abs(ERCCsumLogErr), 
                        data=colData(rse_gene))


degExprs = log2(assays(cov_rse)$count+1)
k_DPseq = num.sv(degExprs, modJoint_DPseq) # 19
qSV_mat_DPseq = prcomp(t(degExprs))$x[,1:k_DPseq]
varExplQsva_DPseq = getPcaVars(prcomp(t(degExprs)))
varExplQsva_DPseq[1:k_DPseq]
sum(varExplQsva_DPseq[1:k_DPseq]) # 87.976%


qSV_mat_DPseq = qSV_mat_DPseq[pd$id ,] 

####  model w/o interaction to subset by region and extracting qSVs and samples related to Lithium exposures #### 

modSep_DPseq <- map(regions, ~cbind(model.matrix(~lithium_group + AgeDeath + Sex  + mitoRate + 
                                           rRNA_rate +totalAssignedGene + RIN + abs(ERCCsumLogErr) +
                                           snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5,
                                           data= pd[pd$BrainRegion == .x,]),
                              qSV_mat_DPseq[pd$BrainRegion == .x,] , id =  pd[pd$BrainRegion == .x , "leafcutter_ids"]))

map(modSep_DPseq, colnames)

write.table(modSep_DPseq$sacc, file= here("data" , "modSep_sACC_Lithium_DPseq.txt") , row.names = F , sep = "\t" , quote=F)
write.table(modSep_DPseq$amyg, file= here("data" , "modSep_Amygdala_Lithium_DPseq.txt") , row.names = T , sep = "\t" ,quote=F)

Amygdala_model_DPseq = read.delim(file= here("data" , "modSep_Amygdala_Lithium_DPseq.txt") , sep = "\t")
sACC_model_DPseq = read.delim(file= here("data" , "modSep_sACC_Lithium_DPseq.txt") , sep = "\t")


sACC_model_DPseq$lithium_groupLithium  = if_else(sACC_model_DPseq$lithium_groupLithium  =="1" , "Lithium " , "Control")
sACC_model_DPseq$SexM = if_else(sACC_model_DPseq$SexM =="0" , "Female" , "Male")


Amygdala_model_DPseq$lithium_groupLithium  = if_else(Amygdala_model_DPseq$lithium_groupLithium  =="1" , "Lithium " , "Control")
Amygdala_model_DPseq$SexM = if_else(Amygdala_model_DPseq$SexM =="0" , "Female" , "Male")



Amygdala_model_DPseq[ , -1] %>% relocate(id) %>% 
write.table( file= here("data" , "modSep_Amygdala_for_leafcutter_DPseq.txt") , row.names = F , sep = "\t" , col.names = F , quote=F)

sACC_model_DPseq[ , -1] %>% relocate(id) %>% 
write.table(file= here("data" , "modSep_sACC_for_leafcutter_DPseq.txt") , row.names = F , sep = "\t" , col.names = F , quote=F)


###### BPseq model #####


modJoint_BPseq = model.matrix(~Dx*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5 +
                          mitoRate + rRNA_rate + totalAssignedGene + RIN +  abs(ERCCsumLogErr), 
                        data=colData(rse_gene))


degExprs = log2(assays(cov_rse)$count+1)
k_BPseq = num.sv(degExprs, modJoint_BPseq) # 19
qSV_mat_BPseq = prcomp(t(degExprs))$x[,1:k_BPseq]
varExplQsva_BPseq = getPcaVars(prcomp(t(degExprs)))
varExplQsva_BPseq[1:k_BPseq]
sum(varExplQsva_BPseq[1:k_BPseq]) # 87.976%


qSV_mat_BPseq = qSV_mat_BPseq[pd$id ,] 

####  model w/o interaction to subset by region and extracting qSVs and samples related to Lithium exposures #### 

modSep_BPseq <- map(regions, ~cbind(model.matrix(~lithium_group + AgeDeath + Sex  + mitoRate + 
                                           rRNA_rate +totalAssignedGene + RIN + abs(ERCCsumLogErr) +
                                           snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5,
                                           data= pd[pd$BrainRegion == .x,]),
                              qSV_mat_BPseq[pd$BrainRegion == .x,] , id =  pd[pd$BrainRegion == .x , "leafcutter_ids"]))

map(modSep_BPseq, colnames)

write.table(modSep_BPseq$sacc, file= here("data" , "modSep_sACC_Lithium_BPseq.txt") , row.names = F , sep = "\t" , quote=F)
write.table(modSep_BPseq$amyg, file= here("data" , "modSep_Amygdala_Lithium_BPseq.txt") , row.names = T , sep = "\t" ,quote=F)

Amygdala_model_BPseq = read.delim(file= here("data" , "modSep_Amygdala_Lithium_BPseq.txt") , sep = "\t")
sACC_model_BPseq = read.delim(file= here("data" , "modSep_sACC_Lithium_BPseq.txt") , sep = "\t")


sACC_model_BPseq$lithium_groupLithium  = if_else(sACC_model_BPseq$lithium_groupLithium  =="1" , "Lithium " , "Control")
sACC_model_BPseq$SexM = if_else(sACC_model_BPseq$SexM =="0" , "Female" , "Male")


Amygdala_model_BPseq$lithium_groupLithium  = if_else(Amygdala_model_BPseq$lithium_groupLithium  =="1" , "Lithium " , "Control")
Amygdala_model_BPseq$SexM = if_else(Amygdala_model_BPseq$SexM =="0" , "Female" , "Male")



Amygdala_model_BPseq[ , -1] %>% relocate(id) %>% 
write.table( file= here("data" , "modSep_Amygdala_for_leafcutter_BPseq.txt") , row.names = F , sep = "\t" , col.names = F , quote=F)

sACC_model_BPseq[ , -1] %>% relocate(id) %>% 
write.table(file= here("data" , "modSep_sACC_for_leafcutter_BPseq.txt") , row.names = F , sep = "\t" , col.names = F , quote=F)


### make a copy of bed files which are actually in BED12 format and are junc files to junc_LIBD  

### Method D which used in this analysis
junc_files <- gsub("_accepted_hits.sorted.bam","_junctions_primaryOnly_regtools.bed", gsub("HISAT2_out","Counts/junction",
                                                                                           gsub("dcl01/lieber/ajaffe/lab","dcs04/lieber/lcolladotor/hydeZandi_LIBD2005", rse_gene_lithium$bamFile)))


# identify the folders
new.folder <- "/dcs04/lieber/lcolladotor/hydeZandi_LIBD2005/Lithium_project/leafcutter/data/junc_LIBD"

# find the files that you want

# copy the files to the new folder
file.copy(junc_files, new.folder)


### Method A in codes; do not need to run

# bam_files <- gsub("dcl01/lieber/ajaffe/lab","dcs04/lieber/lcolladotor/hydeZandi_LIBD2005", rse_gene_lithium$bamFile)
# 
# bam_files_ids = gsub("_accepted_hits.sorted.bam" , "_junctions_primaryOnly_regtools.bed" , 
#                      gsub("/dcs04/lieber/lcolladotor/hydeZandi_LIBD2005/zandiHyde_bipolar_rnaseq/preprocessed_data/HISAT2_out/", "" , bam_files))
# 
# # identify the folders
# new.folder <- "/dcs04/lieber/lcolladotor/hydeZandi_LIBD2005/Lithium_project/leafcutter/data/bam_files"
# 
# # find the files that you want
# 
# # copy the files to the new folder
# file.copy(bam_files, new.folder)

#R13892_H7JLCBBXX_junctions_primaryOnly_regtools.bed
#R13892_H7JLCBBXX_accepted_hits.sorted.bam
#junc_files <- gsub("dcl01/lieber/ajaffe/lab","dcs04/lieber/lcolladotor/hydeZandi_LIBD2005", rse_gene_lithium$bamFile)





### Method C in codes; do not need to run

# count_files <- gsub("_accepted_hits.sorted.bam","_junctions_primaryOnly_regtools.count", gsub("HISAT2_out","Counts/junction",
#                                                                                               gsub("dcl01/lieber/ajaffe/lab","dcs04/lieber/lcolladotor/hydeZandi_LIBD2005", rse_gene_lithium$bamFile)))
# 
# 
# # identify the folders
# new.folder <- "/dcs04/lieber/lcolladotor/hydeZandi_LIBD2005/Lithium_project/leafcutter/data/count_LIBD"
# 
# # find the files that you want
# 
# # copy the files to the new folder
# file.copy(count_files, new.folder)







