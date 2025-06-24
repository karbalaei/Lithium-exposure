##
library(sva)
library(lmerTest)
library(SummarizedExperiment)
library(jaffelab)
library(WGCNA)
library(broom)
library(clusterProfiler)
library(readxl)
library(RColorBrewer)
library(here)
library(gridExtra)
library(grid)
library(purrr)
library(dplyr)
library(tidyverse)

## load data


load(here::here("data", "zandiHypde_bipolar_rseGene_n511.rda"))
load(here::here("data","degradation_rse_BipSeq_BothRegions.rda")) #load cov_rse

## checks
identical(colnames(rse_gene), colnames(cov_rse)) # TRUE
rse_gene$Dx = factor(ifelse(rse_gene$PrimaryDx == "Control", "Control","Bipolar"), 
				levels = c("Control", "Bipolar"))
			
## add ancestry 
load(here::here("data","zandiHyde_bipolar_MDS_n511.rda")) # load mds


mds = mds[rse_gene$BrNum,1:5]
colnames(mds) = paste0("snpPC", 1:5)
colData(rse_gene) = cbind(colData(rse_gene), mds)

#### Lithium exposure info ####


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

#########
## model #
##########


modJoint = model.matrix(~Dx*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 +
	mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr, 
	data=colData(rse_gene))


## qsva
degExprs = log2(assays(cov_rse)$count+1)
k = num.sv(degExprs, modJoint) # 19
qSV_mat = prcomp(t(degExprs))$x[,1:k]


## join and move around region, dx and interaction for cleaning
#modQsva = cbind(modJoint[,c(1:4,16,5:15)], qSV_mat)

rse_gene_lithium <-  rse_gene[ , rse_gene$lithium_group %in% c("Control" , "Lithium")  ]

rse_gene_lithium$lithium_group <-  factor(rse_gene_lithium$lithium_group)

modJoint_lithium  = model.matrix(~lithium_group*BrainRegion + AgeDeath + Sex + snpPC1 + snpPC2 + snpPC3 +
                                 mitoRate + rRNA_rate + totalAssignedGene + RIN + ERCCsumLogErr, 
                               data=colData(rse_gene))


qSV_mat_Lithium <-  qSV_mat[colnames(rse_gene_lithium) , ]

modQsvaLithium = cbind(modJoint_lithium[,c(1:4,14,5:13)], qSV_mat_Lithium)

  
identical(colnames(rse_gene_lithium) , rownames(modJoint_lithium))

#########################
## load wgcna output ####
#########################

## load
load(here::here("data", "constructed_network_signed_bicor.rda"))


# get colors
net$colorsLab = labels2colors(net$colors)
colorDat = data.frame(num = net$colors, col = net$colorsLab, 
	stringsAsFactors=FALSE)
colorDat$Label = paste0("ME", colorDat$num)
colorDat = colorDat[order(colorDat$num),]
colorDat = colorDat[!duplicated(colorDat$num),]
colorDat$numGenes = table(net$colors)[as.character(colorDat$num)]


### Lithium 

#m = modQsva[,1:5]  # this is what was protected

m_Lithium <-  modQsvaLithium [, 1:5]


MEs = net$MEs

colnames(MEs) = colorDat$col[match(colnames(MEs), colorDat$Label)]

MEs = MEs[,colorDat$col]

rownames(MEs) <-   colnames(rse_gene) 

MEs_Lithium <-  MEs[colnames(rse_gene_lithium) , ]


## check
statList_Lithium = lapply(MEs_Lithium, function(x) summary(lmer(x ~ m_Lithium  + (1|rse_gene_lithium$BrNum) - 1))$coef)

# modified
bpdEffect_Lithium = as.data.frame(t(sapply(statList_Lithium, function(x) x[2,])))
regionEffect_Lithium = as.data.frame(t(sapply(statList_Lithium, function(x) x[3,])))
ageEffect_Lithium = as.data.frame(t(sapply(statList_Lithium, function(x) x[4,])))
intEffect_Lithium = as.data.frame(t(sapply(statList_Lithium, function(x) x[5,])))
colnames(bpdEffect_Lithium)= colnames(regionEffect_Lithium) = colnames(ageEffect_Lithium) = colnames(intEffect_Lithium) = c(
  "slope", "se", "df", "t", "pvalue")

signif(bpdEffect_Lithium, 3)
signif(regionEffect_Lithium, 3)
signif(ageEffect_Lithium, 3)
signif(intEffect_Lithium, 3)


##################
# make boxplots ##
lab = paste0(substr(rse_gene_lithium$BrainRegion,1,4), ":", rse_gene_lithium$lithium_group) 
lab = factor(lab, levels = c("sACC:Control", "sACC:Lithium", "Amyg:Control", "Amyg:Lithium"))
 

pdf(here::here("graphs" , "BPseq_model" ,  "wgcna" ,  "MEs_vs_Lithium.pdf"),w=8,h=6)
palette(brewer.pal(4,"Paired"))
par(mar=c(3,6,4,2), cex.axis=0.8 ,cex.lab=1 ,cex.main = 1.4 )
for(i in 1:ncol(MEs)) {
	boxplot(MEs_Lithium[,i] ~ lab, outline = FALSE, xlab="",
		ylim = quantile(unlist(MEs),c(0.002,0.998)),main = colnames(MEs)[i],
		names = gsub(":", "\n", levels(lab)), ylab = "Module Eigengene")
	points(MEs_Lithium[,i] ~ jitter(as.numeric(lab),amount=0.1), pch=21, bg=lab)
	legend("top", c(paste0("Region p=", signif(regionEffect_Lithium[i,5],3)), 
		paste0("Lithium p=", signif(bpdEffect_Lithium[i,5],3))),cex=1)
}
dev.off()
 
#rse_gene_lithium$lithium_group_named <-  str_replace_all(rse_gene_lithium$lithium_group , '0' , "Control") %>%   str_replace_all( '1' , "Lithium")
  
clean_ME = t(cleaningY(t(MEs_Lithium), m_Lithium, P=2))
pdf(here::here("graphs" , "BPseq_model" , "wgcna" ,"clean_MEs_vs_Lithium.pdf"),w=5,h=6)
palette(brewer.pal(4,"Paired"))
par(mar=c(3,6,4,2), cex.axis=0.8,cex.lab=1,cex.main = 1.4)
for(i in 1:ncol(clean_ME)) {
	boxplot(clean_ME[,i] ~ rse_gene_lithium$lithium_group, outline = FALSE, xlab="",
		ylim = quantile(unlist(clean_ME),c(0.002,0.998)),main = colnames(MEs_Lithium)[i],
		ylab = "Module Eigengene (Adj)")
	points(clean_ME[,i] ~ jitter(as.numeric(rse_gene_lithium$lithium_group),amount=0.1), pch=21, bg=lab)
	legend("top", paste0("Lithium p=", signif(bpdEffect_Lithium[i,5],3)),cex=1)
}
dev.off()
	


##############################
## enrichment of DEGs #####
#############################

#### Lithium ####

load(here::here("results" ,  "Lithium results BPseqmodel.RDS"))

#rm(outExon_bothRegion, outJxn_bothRegion, outTx_bothRegion)

identical(rownames(Lithium_BPseqmodel$gene$Amygdala), fNames)

length(rownames(Lithium_BPseqmodel$gene$Amygdala))  #24784
length(fNames) #25136

sum(fNames %in% rownames(Lithium_BPseqmodel$gene$Amygdala)  ) #24784


sum(rownames(Lithium_BPseqmodel$gene$Amygdala) %in% fNames) #24784

net$fNames <- fNames

#### Genes #####

#rm(list = c(tt , tt_df , deModEnrich , numGenes ,  prop.table_tt , prop.table_tt_df))

Enrichment_fun_gene <- function( i , j , name) {

Avi_id <-   which(net$fNames %in% rownames(Lithium_BPseqmodel[[i]][[j]]))

net_Lithium <-  list()

net_Lithium$colorsLab <-  net$colorsLab[Avi_id]

tt = table(net_Lithium$colorsLab, Lithium_BPseqmodel[[i]][[j]]$P.Value < 0.005)
tt = tt[rownames(bpdEffect_Lithium),]

prop.table_tt  = prop.table(tt,1)


tt_df =  as.data.frame(tt)%>% setNames(c("Module" , "Avail" , "Num")) %>% 
  pivot_wider( names_from = Avail , values_from = Num) 

prop.table_tt_df =  as.data.frame(prop.table_tt)%>% setNames(c("Module" , "Avail" , "Num")) %>% 
  pivot_wider( names_from = Avail , values_from = Num) %>% mutate(across(where(is.numeric), round, 3)) 
  

# manual chi-sq?
deModEnrich = as.data.frame(t(sapply(colorDat$col, function(cc) {
  tab = table(net_Lithium$colorsLab == cc,Lithium_BPseqmodel[[i]][[j]]$P.Value < 0.005)
  c(chisq.test(tab)$p.value, getOR(tab))
})))
colnames(deModEnrich) = c("Pvalue", "OR")

numGenes <-  table(net_Lithium$colorsLab)  %>%  as.data.frame() %>% setNames(c("Module" , "Size" )) 

deModEnrich <-  deModEnrich %>% rownames_to_column("Module") %>%
  mutate(across(where(is.numeric), signif, 3)) %>%
  merge(numGenes , by = "Module") %>% 
  merge(tt_df , by = "Module") %>% 
  merge(prop.table_tt_df , by = "Module") %>% 
  setNames(c("Module", "Pvalue",   "OR" ,      "numGenes" , "#Not DEfs" ,  "#DEfs" , "%Not DEfs" ,  "%DEfs" ))


write.csv(prop.table_tt , here::here("results" , "wgcna_Enrichment" , "BPseq_model" , paste0("prop.table of  availability of DEf in WGCNA modules" , name , ".csv" ) ), row.names = T)

write.csv(deModEnrich , here::here("results" , "wgcna_Enrichment" , "BPseq_model" , paste0("Enrichment on WGCNA modules stat graphs " , name , ".csv" ) ), row.names = T)


jpeg(here::here("graphs" , "BPseq_model" , "wgcna" ,  paste0("Enrichment on WGCNA modules stat graphs " , name , ".jpg" ) ) , res = 300 , width = 5000 , height = 4000)

par( mfrow= c(2,1)  , mar=c(4,4,4,2))

plot(prop.table(tt,1)[,2], -log10(bpdEffect_Lithium$pvalue) ,
     xlab="prop of availability of DEf in modules", ylab="-log10 pvalue of bpdEffect")

plot(-log10(deModEnrich$Pvalue), -log10(bpdEffect_Lithium$pvalue ) ,
     xlab="-log10 pvalue from Chi-square test", ylab="-log10 pvalue of bpdEffect")


dev.off()


jpeg(here::here("graphs" , "BPseq_model" , "wgcna" ,  paste0("Enrichment DEfs on WGCNA modules stat graphs " , name , ".jpg" ) ) , res = 200 , width = 3000 , height = 4000)

p<-tableGrob(deModEnrich)
grid.arrange(p)

dev.off()

jpeg(here::here("graphs" , "BPseq_model" , "wgcna" ,  paste0("prop.table of  availability of DEf in WGCNA modules " , name , ".jpg" ) ) , res = 200 , width = 3000 , height = 4000)

p<-tableGrob(prop.table(tt,1))
grid.arrange(p)

dev.off()


}

gene_plan <-  data.frame( i = 1 , j = c(1 , 2) , name = c("Gene_Amygdala" , "Gene_sACC"))

pmap(gene_plan , Enrichment_fun_gene)

#### Exon #####

#rm(list = tt , tt_df , deModEnrich , numGenes ,  prop.table_tt , prop.table_tt_df)


name = "Exon_Amygdala" # i =2 j = 1
name = "Exon_sACC" # i =1 j = 2

i =2

j=1

Enrichment_fun_exon <- function( i , j , name) {

common_id <-   intersect(net$fNames , Lithium_BPseqmodel[[i]][[j]]$common_gene_id)

Avi_id <-   which(net$fNames %in%  common_id)

Lithium_BPseqmodel_common <-  Lithium_BPseqmodel[[i]][[j]][Lithium_BPseqmodel[[i]][[j]]$common_gene_id %in% common_id ,]

Lithium_BPseqmodel_DEf <-  Lithium_BPseqmodel_common[Lithium_BPseqmodel_common$P.Value < 0.005 , 'common_gene_id' ] %>% unique()

univers <-  Lithium_BPseqmodel_common$common_gene_id %>% unique() %>% as.data.frame() %>%
  setNames("Ids") %>% mutate("DEf" = if_else(Ids %in% Lithium_BPseqmodel_DEf , TRUE , FALSE))

net_Lithium <-  list()

net_Lithium$colorsLab <-  net$colorsLab[Avi_id]

tt = table(net_Lithium$colorsLab, univers$DEf )
tt = tt[rownames(bpdEffect_Lithium),]

prop.table_tt  = prop.table(tt,1)


tt_df =  as.data.frame(tt)%>% setNames(c("Module" , "Avail" , "Num")) %>% 
  pivot_wider( names_from = Avail , values_from = Num) 

prop.table_tt_df =  as.data.frame(prop.table_tt)%>% setNames(c("Module" , "Avail" , "Num")) %>% 
  pivot_wider( names_from = Avail , values_from = Num) %>% mutate(across(where(is.numeric), round, 3)) 


# manual chi-sq?

deModEnrich = as.data.frame(t(sapply(colorDat$col, function(cc) {
  tab = table(net_Lithium$colorsLab == cc,univers$DEf)
  c(chisq.test(tab)$p.value, getOR(tab))
})))
colnames(deModEnrich) = c("Pvalue", "OR")

numGenes <-  table(net_Lithium$colorsLab)  %>%  as.data.frame() %>% setNames(c("Module" , "Size" )) 

deModEnrich <-  deModEnrich %>% rownames_to_column("Module") %>%
  mutate(across(where(is.numeric), signif, 3)) %>%
  merge(numGenes , by = "Module") %>% 
  merge(tt_df , by = "Module") %>% 
  merge(prop.table_tt_df , by = "Module") %>% 
  setNames(c("Module", "Pvalue",   "OR" ,      "numGenes" , "#Not DEfs" ,  "#DEfs" , "%Not DEfs" ,  "%DEfs" ))



write.csv(prop.table_tt , here::here("results" , "wgcna_Enrichment" , "BPseq_model" , paste0("prop.table of  availability of DEf in WGCNA modules" , name , ".csv" ) ), row.names = T)

write.csv(deModEnrich , here::here("results" , "wgcna_Enrichment" , "BPseq_model" , paste0("Enrichment on WGCNA modules stat graphs " , name , ".csv" ) ), row.names = T)


jpeg(here::here("graphs" , "BPseq_model" , "wgcna" ,  paste0("Enrichment on WGCNA modules stat graphs " , name , ".jpg" ) ) , res = 300 , width = 5000 , height = 4000)

par( mfrow= c(2,1)  , mar=c(4,4,4,2))

plot(prop.table(tt,1)[,2], -log10(bpdEffect_Lithium$pvalue) ,
     xlab="prop of availability of DEf in modules", ylab="-log10 pvalue of bpdEffect")

plot(-log10(deModEnrich$Pvalue), -log10(bpdEffect_Lithium$pvalue ) ,
     xlab="-log10 pvalue from Chi-square test", ylab="-log10 pvalue of bpdEffect")


dev.off()


jpeg(here::here("graphs" , "BPseq_model" , "wgcna" ,  paste0("Enrichment DEfs on WGCNA modules stat graphs " , name , ".jpg" ) ) , res = 200 , width = 3000 , height = 4000)

p<-tableGrob(deModEnrich)
grid.arrange(p)

dev.off()

jpeg(here::here("graphs" , "BPseq_model" , "wgcna" ,  paste0("prop.table of  availability of DEf in WGCNA modules " , name , ".jpg" ) ) , res = 200 , width = 3000 , height = 4000)

p<-tableGrob(prop.table(tt,1))
grid.arrange(p)

dev.off()

}

exon_plan <-  data.frame( i = 2 , j = c(1 , 2) , name = c("Exon_Amygdala" , "Exon_sACC"))

pmap(exon_plan , Enrichment_fun_exon)

#### Junction #####

#rm(list = tt , tt_df , deModEnrich , numGenes ,  prop.table_tt , prop.table_tt_df)

name = "Junction_Amygdala" # i =2 j = 1
name = "Junction_sACC" # i =1 j = 2

i =3

j=1

Enrichment_fun_junction <- function( i , j , name) {
  
  common_id <-   intersect(net$fNames , unique(Lithium_BPseqmodel[[i]][[j]]$common_gene_id))
  
  Avi_id <-   which(net$fNames %in%  common_id)
  
  Lithium_BPseqmodel_common <-  Lithium_BPseqmodel[[i]][[j]][Lithium_BPseqmodel[[i]][[j]]$common_gene_id %in% common_id ,]
  
  Lithium_BPseqmodel_DEf <-  Lithium_BPseqmodel_common[Lithium_BPseqmodel_common$P.Value < 0.005 , 'common_gene_id' ] %>% unique()
  
  
  univers <-  Lithium_BPseqmodel_common$common_gene_id %>% unique() %>% as.data.frame() %>%
    setNames("Ids") %>% mutate("DEf" = if_else(Ids %in% Lithium_BPseqmodel_DEf , TRUE , FALSE)) 
  net_Lithium <-  list()
  
  net_Lithium$colorsLab <-  net$colorsLab[Avi_id]
  
  
  
  #Avi_id2 <-   which(net$fNames %in% univers)
  
  tt = table(net_Lithium$colorsLab, univers$DEf )
  tt = tt[rownames(bpdEffect_Lithium),]
  
  tt_df =  as.data.frame(tt)%>% setNames(c("Module" , "Avail" , "Num")) %>% 
    pivot_wider( names_from = Avail , values_from = Num) 
  
  
  prop.table_tt  = prop.table(tt,1)
  
  prop.table_tt_df =  as.data.frame(prop.table_tt)%>% setNames(c("Module" , "Avail" , "Num")) %>% 
    pivot_wider( names_from = Avail , values_from = Num) %>% mutate(across(where(is.numeric), round, 3)) 
  

  # manual chi-sq?
  
  deModEnrich = as.data.frame(t(sapply(colorDat$col, function(cc) {
    tab = table(net_Lithium$colorsLab == cc,univers$DEf)
    c(chisq.test(tab)$p.value, getOR(tab))
  })))
  colnames(deModEnrich) = c("Pvalue", "OR")
  
  numGenes <-  table(net_Lithium$colorsLab)  %>%  as.data.frame() %>% setNames(c("Module" , "Size" )) 
  
  deModEnrich <-  deModEnrich %>% rownames_to_column("Module") %>%
    mutate(across(where(is.numeric), signif, 3)) %>%
    merge(numGenes , by = "Module") %>% 
    merge(tt_df , by = "Module") %>% 
    merge(prop.table_tt_df , by = "Module") %>% 
    setNames(c("Module", "Pvalue",   "OR" ,      "numGenes" , "#Not DEfs" ,  "#DEfs" , "%Not DEfs" ,  "%DEfs" ))
  
  write.csv(prop.table_tt , here::here("results" , "wgcna_Enrichment" , "BPseq_model" , paste0("prop.table of  availability of DEf in WGCNA modules" , name , ".csv" ) ), row.names = T)
  
  write.csv(deModEnrich , here::here("results" , "wgcna_Enrichment" , "BPseq_model" , paste0("Enrichment on WGCNA modules stat graphs " , name , ".csv" ) ), row.names = T)
  
  
  jpeg(here::here("graphs" , "BPseq_model" , "wgcna" ,  paste0("Enrichment on WGCNA modules stat graphs " , name , ".jpg" ) ) , res = 300 , width = 5000 , height = 4000)
  
  par( mfrow= c(2,1)  , mar=c(4,4,4,2))
  
  plot(prop.table(tt,1)[,2], -log10(bpdEffect_Lithium$pvalue) ,
       xlab="prop of availability of DEf in modules", ylab="-log10 pvalue of bpdEffect")
  
  plot(-log10(deModEnrich$Pvalue), -log10(bpdEffect_Lithium$pvalue ) ,
       xlab="-log10 pvalue from Chi-square test", ylab="-log10 pvalue of bpdEffect")
  
  
  dev.off()
  
  
  jpeg(here::here("graphs" , "BPseq_model" , "wgcna" ,  paste0("Enrichment DEfs on WGCNA modules stat graphs " , name , ".jpg" ) ) , res = 400 , width = 3000 , height = 4000)
  
  p<-tableGrob(deModEnrich)
  grid.arrange(p)
  
  dev.off()
  
  jpeg(here::here("graphs" , "BPseq_model" , "wgcna" ,  paste0("prop.table of  availability of DEf in WGCNA modules " , name , ".jpg" ) ) , res = 400 , width = 3000 , height = 4000)
  
  p<-tableGrob(prop.table(tt,1))
  grid.arrange(p)
  
  dev.off()
  
}

junction_plan <-  data.frame( i = 3 , j = c(1 , 2) , name = c("Junction_Amygdala" , "Junction_sACC"))

pmap(junction_plan , Enrichment_fun_junction)

#### Transcript #####
#rm(list = tt , tt_df , deModEnrich , numGenes ,  prop.table_tt , prop.table_tt_df)

# name = "Junction_Amygdala" # i =2 j = 1
# name = "Junction_sACC" # i =1 j = 2
# 
# i =4
# 
# j=1

Enrichment_fun_transcript <- function( i , j , name) {
  
  common_id <-   intersect(net$fNames , unique(Lithium_BPseqmodel[[i]][[j]]$common_gene_id))
  
  Avi_id <-   which(net$fNames %in%  common_id)
  
  Lithium_BPseqmodel_common <-  Lithium_BPseqmodel[[i]][[j]][Lithium_BPseqmodel[[i]][[j]]$common_gene_id %in% common_id ,]
  
  Lithium_BPseqmodel_DEf <-  Lithium_BPseqmodel_common[Lithium_BPseqmodel_common$P.Value < 0.005 , 'common_gene_id' ] %>% unique()
  
  univers <-  Lithium_BPseqmodel_common$common_gene_id %>% unique() %>% as.data.frame() %>%
    setNames("Ids") %>% mutate("DEf" = if_else(Ids %in% Lithium_BPseqmodel_DEf , TRUE , FALSE))
  net_Lithium <-  list()
  
  net_Lithium$colorsLab <-  net$colorsLab[Avi_id]
  
  tt = table(net_Lithium$colorsLab, univers$DEf )
  tt = tt[rownames(bpdEffect_Lithium),]
  
  prop.table_tt  = prop.table(tt,1)
    
  tt_df =  as.data.frame(tt)%>% setNames(c("Module" , "Avail" , "Num")) %>% 
    pivot_wider( names_from = Avail , values_from = Num) 
  
  prop.table_tt_df =  as.data.frame(prop.table_tt)%>% setNames(c("Module" , "Avail" , "Num")) %>% 
    pivot_wider( names_from = Avail , values_from = Num) %>% mutate(across(where(is.numeric), round, 3)) 
    
  # manual chi-sq?
    deModEnrich = as.data.frame(t(sapply(colorDat$col, function(cc) {
    tab = table(net_Lithium$colorsLab == cc,univers$DEf)
    c(chisq.test(tab)$p.value, getOR(tab))
  })))
  colnames(deModEnrich) = c("Pvalue", "OR")
  
  numGenes <-  table(net_Lithium$colorsLab)  %>%  as.data.frame() %>% setNames(c("Module" , "Size" )) 
  
  deModEnrich <-  deModEnrich %>% rownames_to_column("Module") %>%
    mutate(across(where(is.numeric), signif, 3)) %>%
    merge(numGenes , by = "Module") %>% 
    merge(tt_df , by = "Module") %>% 
    merge(prop.table_tt_df , by = "Module") %>% 
    setNames(c("Module", "Pvalue",   "OR" ,      "numGenes" , "#Not DEfs" ,  "#DEfs" , "%Not DEfs" ,  "%DEfs" ))
  
  
  write.csv(prop.table_tt , here::here("results" , "wgcna_Enrichment" , "BPseq_model" , paste0("prop.table of  availability of DEf in WGCNA modules" , name , ".csv" ) ), row.names = T)
  
  write.csv(deModEnrich , here::here("results" , "wgcna_Enrichment" , "BPseq_model" , paste0("Enrichment on WGCNA modules stat graphs " , name , ".csv" ) ), row.names = T)
  
  
  jpeg(here::here("graphs" , "BPseq_model" , "wgcna" ,  paste0("Enrichment on WGCNA modules stat graphs " , name , ".jpg" ) ) , res = 300 , width = 5000 , height = 4000)
  
  par( mfrow= c(2,1)  , mar=c(4,4,4,2))
  
  plot(prop.table(tt,1)[,2], -log10(bpdEffect_Lithium$pvalue) ,
       xlab="prop of availability of DEf in modules", ylab="-log10 pvalue of bpdEffect")
  
  plot(-log10(deModEnrich$Pvalue), -log10(bpdEffect_Lithium$pvalue ) ,
       xlab="-log10 pvalue from Chi-square test", ylab="-log10 pvalue of bpdEffect")
  
  
  dev.off()
  
  
  jpeg(here::here("graphs" , "BPseq_model" , "wgcna" ,  paste0("Enrichment DEfs on WGCNA modules stat graphs " , name , ".jpg" ) ) , res = 400 , width = 3000 , height = 4000)
  
  p<-tableGrob(deModEnrich)
  grid.arrange(p)
  
  dev.off()
  
  jpeg(here::here("graphs" , "BPseq_model" , "wgcna" ,  paste0("prop.table of  availability of DEf in WGCNA modules " , name , ".jpg" ) ) , res = 400 , width = 3000 , height = 4000)
  
  p<-tableGrob(prop.table(tt,1))
  grid.arrange(p)
  
  dev.off()
  
}

transcript_plan <-  data.frame( i = 4 , j = c(1 , 2) , name = c("Transcript_Amygdala" , "Transcript_sACC"))

pmap(transcript_plan , Enrichment_fun_transcript)


## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()
