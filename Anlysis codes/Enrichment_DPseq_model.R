### Library #####
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
library(ggplot2)




load(file = here("results" ,"Lithium results DPseqmodel.RDS"))

map_depth(Lithium_DPseqmodel, 2, nrow)

## Define Universe
all_gencode <- map_depth(Lithium_DPseqmodel, 2, "common_gene_id")
map_depth(all_gencode, 2, length)
head(all_gencode$gene$Amygdala)


length(unlist(all_gencode))
# [1] 1503060

## all ENSEMBL
all_ensembl <- unique(ss(unlist(all_gencode),"\\."))
length(all_ensembl)
# [1] 33441

all_entrez <- bitr(all_ensembl, fromType = "ENSEMBL", toType ="ENTREZID", OrgDb="org.Hs.eg.db")
# 30.99% of input gene IDs are fail to map...
nrow(all_entrez)
# [1] 25737
u <- all_entrez$ENTREZ

#### Functions for extracting gene sets ####

get_signif <- function(outFeature, colname = "common_feature_id", cutoff = 0.005, return_unique = FALSE , UpDown ){
  if (UpDown == "none") {
  signif <- outFeature[[colname]][outFeature$P.Value < cutoff]
  if(return_unique) signif <- unique(signif)
  signif <- signif[!is.na(signif)]
  return(signif)
  } else if ( UpDown == "Up") {
    
    signif <- outFeature[colname][c(outFeature$P.Value < cutoff & outFeature$logFC > 0 ), ]

    if(return_unique) signif <- unique(signif)
    signif <- signif[!is.na(signif)]
    return(signif)
    
  } else if (UpDown == "Down") {
    
    
    signif <- outFeature[colname][c(outFeature$P.Value < cutoff & outFeature$logFC < 0 ), ]
    if(return_unique) signif <- unique(signif)
    signif <- signif[!is.na(signif)]
    return(signif)
    
  
  }
}




my_flatten <- function (x, use.names = TRUE, classes = "ANY") {
  #' Source taken from rlist::list.flatten
  len <- sum(rapply(x, function(x) 1L, classes = classes))
  y <- vector("list", len)
  i <- 0L
  items <- rapply(x, function(x) {
    i <<- i + 1L
    y[[i]] <<- x
    TRUE
  }, classes = classes)
  if (use.names && !is.null(nm <- names(items))) 
    names(y) <- nm
  y
}

my_get_entrez <- function(g){
  e <- unique(ss(g, "\\."))
  entrez <- bitr(e, fromType = "ENSEMBL", toType ="ENTREZID", OrgDb="org.Hs.eg.db")
  return(entrez$ENTREZ)
}

#### Get signif genes####
signif_genes_whole_list<- map_depth(Lithium_DPseqmodel, 2, ~get_signif(.x, colname = "common_gene_id", return_unique = TRUE , UpDown = "none" ))

signif_genes_whole_list$All_assigned_genes$Amygdala =  map(signif_genes_whole_list , 1) %>%  list_c() %>% unique()
signif_genes_whole_list$All_assigned_genes$sACC =  map(signif_genes_whole_list , 2) %>%  list_c() %>% unique()


signif_genes_up <- map_depth(Lithium_DPseqmodel, 2, ~get_signif(.x, colname = "common_gene_id", return_unique = TRUE , UpDown = "Up" ))

signif_genes_down<- map_depth(Lithium_DPseqmodel, 2, ~get_signif(.x, colname = "common_gene_id", return_unique = TRUE , UpDown = "Down" ))


### All ###

#### Combine All Features ####
## transpose and flatten so features are combined

signif_genes_af <- transpose(signif_genes_whole_list)


names(signif_genes_af$Amygdala)


signif_genes_af_flat <- my_flatten(map_depth(signif_genes_af, 2, unlist))
map_int(signif_genes_af_flat, length)
 


## Extract unique and convert to entrez
gene_sets_af <- map(signif_genes_af_flat, my_get_entrez)


map_int(gene_sets_af, length)


t(rbind(map_dfr(signif_genes_af_flat, length), map_dfr(gene_sets_af, length)))


str(gene_sets_af)


#### Run Enrichment ####

ont <- c("BP", "CC", "MF", "ALL")
names(ont) <- ont


plan <- tidyr::expand_grid(
  gene_list = gene_sets_af,
  ontology = ont
)

Enrichment_results <-  pmap(plan,  ~enrichGO(.x,
                                                                       universe = u, OrgDb = org.Hs.eg.db,
                                                                       ont = .y , pAdjustMethod = "BH",
                                                                       pvalueCutoff  = .2, qvalueCutoff  = .5,
                                                                       readable= TRUE)   )

names(Enrichment_results) <- c(paste0(names(Enrichment_results), "_", rep(ont)))

gene_set_numbers <-  map_int(Enrichment_results, nrow)

to_remove <-  which(gene_set_numbers == 0)


if (length(to_remove) > 0 ){

Enrichment_results <-  Enrichment_results[-(c(to_remove))]

}


all_plots <-  map2(Enrichment_results, names(Enrichment_results), ~dotplot(.x, title = .y , showCategory = 25))


walk2(names(all_plots), all_plots, ~ggsave(filename = here("graphs" , "DPseq_model" ,  "Enrichment results" , "whole_list" , paste0(.x , ".jpg")), plot = .y, 
                                             height = 14, width = 9))



pdf(here("graphs" , "DPseq_model" ,  "Enrichment results" , "whole_list" ,  "Enrichment results Lithium exposure.pdf") , height = 14 , width = 9)

map2(Enrichment_results, names(Enrichment_results), ~dotplot(.x, title = .y , showCategory = 25))

dev.off()

save(Enrichment_results,file = here("results", "enrichGO_DPseq_model_whole_list.rda"))

### Up ###

#### Combine All Features ####
## transpose and flatten so features are combined

signif_genes_uf <- transpose(signif_genes_up)


names(signif_genes_uf$Amygdala)


signif_genes_uf_flat <- my_flatten(map_depth(signif_genes_uf, 2, unlist))
map_int(signif_genes_uf_flat, length)
 

## Extract unique and convert to entrez
gene_sets_uf <- map(signif_genes_uf_flat, my_get_entrez)


map_int(gene_sets_uf, length)


t(rbind(map_dfr(signif_genes_uf_flat, length), map_dfr(gene_sets_uf, length)))


str(gene_sets_uf)


#### Run Enrichment ####

plan_up <- tidyr::expand_grid(
  gene_list = gene_sets_uf,
  ontology = ont
)

Enrichment_results_up <-  pmap(plan_up,  ~enrichGO(.x,
                                                                       universe = u, OrgDb = org.Hs.eg.db,
                                                                       ont = .y , pAdjustMethod = "BH",
                                                                       pvalueCutoff  = .2, qvalueCutoff  = .5,
                                                                       readable= TRUE)   )

names(Enrichment_results_up) <- c(paste0(names(Enrichment_results_up), "_", rep(ont)))

class(Enrichment_results_up)


gene_set_numbers_up <-  map_int(Enrichment_results_up, nrow)

to_remove_up <-  which(gene_set_numbers_up == 0)

if (length(to_remove_up) > 0 ){

Enrichment_results_up <-  Enrichment_results_up[-(c(to_remove_up))]

}


all_plots_up <-  map2(Enrichment_results_up, names(Enrichment_results_up), ~dotplot(.x, title = .y , showCategory = 25))


walk2(names(all_plots_up), all_plots_up, ~ggsave(filename = here("graphs" , "DPseq_model" ,  "Enrichment results" , "up" , paste0(.x , ".jpg")), plot = .y, 
                                             height = 14, width = 9))




pdf(here("graphs" , "DPseq_model" ,  "Enrichment results" , "up" , "Enrichment results Lithium exposure_up.pdf") , height = 14 , width = 9)

map2(Enrichment_results_up, names(Enrichment_results_up), ~dotplot(.x, title = .y , showCategory = 25))

dev.off()

save(Enrichment_results_up,file = here("results", "enrichGO_DPseq_model_up.rda"))


### down###

#### Combine All Features ####
## transpose and flatten so features are combined

signif_genes_df <- transpose(signif_genes_down)


names(signif_genes_df$Amygdala)


signif_genes_df_flat <- my_flatten(map_depth(signif_genes_df, 2, unlist))
map_int(signif_genes_df_flat, length)
 

## Extract unique and convert to entrez
gene_sets_df <- map(signif_genes_df_flat, my_get_entrez)


map_int(gene_sets_df, length)


t(rbind(map_dfr(signif_genes_df_flat, length), map_dfr(gene_sets_df, length)))


str(gene_sets_df)


#### Run Enrichment ####


plan_down <- tidyr::expand_grid(
  gene_list = gene_sets_df,
  ontology = ont
)

Enrichment_results_down<-  pmap(plan_down ,  ~enrichGO(.x,
                                                                       universe = u, OrgDb = org.Hs.eg.db,
                                                                       ont = .y , pAdjustMethod = "BH",
                                                                       pvalueCutoff  = .2, qvalueCutoff  = .5,
                                                                       readable= TRUE)   )

names(Enrichment_results_down) <- c(paste0(names(Enrichment_results_down), "_", rep(ont)))

class(Enrichment_results_down)


gene_set_numbers_down<-  map_int(Enrichment_results_down, nrow)


to_remove_down<-  which(gene_set_numbers_down == 0)


if (length(to_remove_down) > 0 ){

Enrichment_results_down<-  Enrichment_results_down[-(c(to_remove_down))]

}


all_plots_down <-  map2(Enrichment_results_down, names(Enrichment_results_down), ~dotplot(.x, title = .y , showCategory = 25))


walk2(names(all_plots_down), all_plots_down, ~ggsave(filename = here("graphs" , "DPseq_model" , "Enrichment results" , "down" , paste0(.x , ".jpg")), plot = .y, 
                                             height = 14, width = 9))


pdf(here("graphs" , "DPseq_model" , "Enrichment results" , "down"  , "Enrichment results Lithium exposure_down.pdf") , height = 14 , width = 9)

map2(Enrichment_results_down, names(Enrichment_results_down), ~dotplot(.x, title = .y , showCategory = 25))

dev.off()

save(Enrichment_results_down,file = here("results", "enrichGO_DPseq_model_down.rda"))






## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()

