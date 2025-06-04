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




load(file = here("results" ,"Lithium results DPseqmodel.RDS"))

map_depth(Outlist, 2, nrow)

## Define Universe
all_gencode <- map_depth(Outlist, 2, "common_gene_id")
map_depth(all_gencode, 2, length)
head(all_gencode$gene$Amygdala)


length(unlist(all_gencode))
# [1] 1522730

## all ENSEMBL
all_ensembl <- unique(ss(unlist(all_gencode),"\\."))
length(all_ensembl)
# [1] 33932

all_entrez <- bitr(all_ensembl, fromType = "ENSEMBL", toType ="ENTREZID", OrgDb="org.Hs.eg.db")
# 31.02% of input gene IDs are fail to map...
nrow(all_entrez)
# [1] 23577
u <- all_entrez$ENTREZ

#### Functions for extracting gene sets ####
get_signif <- function(outFeature, colname = "common_feature_id", cutoff = 0.005, return_unique = FALSE){
  signif <- outFeature[[colname]][outFeature$P.Value < cutoff]
  if(return_unique) signif <- unique(signif)
  signif <- signif[!is.na(signif)]
  return(signif)
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
signif_genes <- map_depth(Outlist, 2, ~get_signif(.x, colname = "common_gene_id", return_unique = TRUE))
map_depth(signif_genes, 2, length)


#### Combine All Features ####
## transpose and flatten so features are combined

signif_genes_af <- transpose(signif_genes)


names(signif_genes_af$Amygdala)


signif_genes_af_flat <- my_flatten(map_depth(signif_genes_af, 2, unlist))
map_int(signif_genes_af_flat, length)
# Amygdala.gene Amygdala.exon  Amygdala.jxn   Amygdala.tx     sACC.gene     sACC.exon      sACC.jxn       sACC.tx 
# 249           1014           1675           287             52            526            1019           320 

## Extract unique and convert to entrez
gene_sets_af <- map(signif_genes_af_flat, my_get_entrez)

# 1: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   4.02% of input gene IDs are fail to map...
# 2: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   8.88% of input gene IDs are fail to map...
# 3: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   7.59% of input gene IDs are fail to map...
# 4: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   16.72% of input gene IDs are fail to map...
# 5: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   26.92% of input gene IDs are fail to map...
# 6: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   9.51% of input gene IDs are fail to map...
# 7: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   8.87% of input gene IDs are fail to map...
# 8: In bitr(e, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db") :
#   10.94% of input gene IDs are fail to map...

map_int(gene_sets_af, length)
# Amygdala.gene Amygdala.exon  Amygdala.jxn   Amygdala.tx     sACC.gene     sACC.exon      sACC.jxn       sACC.tx 
# 243            930            1551           241            38           484              937            286 


t(rbind(map_dfr(signif_genes_af_flat, length), map_dfr(gene_sets_af, length)))


#               [,1] [,2]
# Amygdala.gene  249  243
# Amygdala.exon 1014  930
# Amygdala.jxn  1675 1551
# Amygdala.tx    287  241
# sACC.gene       52   38
# sACC.exon      526  484
# sACC.jxn      1019  937
# sACC.tx        320  286


str(gene_sets_af)
# List of 8
# $ Amygdala.gene: chr [1:243] "7133" "712" "714" "713" ...
# $ Amygdala.exon: chr [1:930] "54587" "728661" "8764" "6518" ...
# $ Amygdala.jxn : chr [1:1551] "83858" "8764" "27237" "5293" ...
# $ Amygdala.tx  : chr [1:241] "473" "4681" "100532736" "712" ...
# $ sACC.gene    : chr [1:38] "149563" "100130630" "6566" "107985216" ...
# $ sACC.exon    : chr [1:484] "8718" "473" "57540" "374946" ...
# $ sACC.jxn     : chr [1:937] "100288175" "83858" "5590" "8514" ...
# $ sACC.tx      : chr [1:286] "55052" "6146" "113452" "106479825" ...


#### Run Enrichment ####

names(gene_sets_af)
# 
# [1] "Amygdala.gene" "Amygdala.exon" "Amygdala.jxn"  "Amygdala.tx"   "sACC.gene"     "sACC.exon"    
# [7] "sACC.jxn"      "sACC.tx" 


ont <- c("BP", "CC", "MF", "ALL")
names(ont) <- ont


plan <- tidyr::expand_grid(
  gene_list = gene_sets_af,
  ontology = ont
)

Enrichment_results <-  map2(plan$gene_list , plan$ontology,  ~enrichGO(.x,
                                                                       universe = u, OrgDb = org.Hs.eg.db,
                                                                       ont = .y , pAdjustMethod = "BH",
                                                                       pvalueCutoff  = .2, qvalueCutoff  = .5,
                                                                       readable= TRUE)   )

names(Enrichment_results) <- c(paste0(names(Enrichment_results), "_", rep(ont)))

class(Enrichment_results)


gene_set_numbers <-  map_int(Enrichment_results, nrow)

to_remove <-  which(gene_set_numbers == 0)

Enrichment_results <-  Enrichment_results[-(c(to_remove))]


pdf(here("graphs" , "Enrichment results Lithium exposure.pdf") , height = 12 , width = 7)

map2(Enrichment_results, names(Enrichment_results), ~dotplot(.x, title = .y , showCategory = 25))

dev.off()

save(Enrichment_results,file = here("results", "enrichGO.rda"))


## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()

