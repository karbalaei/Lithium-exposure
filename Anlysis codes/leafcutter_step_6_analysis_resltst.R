

leafviz(infile = "d:/BPseq_data/Lithium_project/Compare models/leafcutter/step_5_BPseq_Amygdala.Rdata")



## Function to load:

#run_leafvir.R
#make_cluster_plot.R
#make_gene_plot.R
#filter_intron_table in server.R

setwd("D:/BPseq_data/Lithium_project/")
library(here)
library(tidyverse)
library(gridExtra)
library(grid)
library(gtable)

here::i_am("Compare models/codes/leafcutter_step_6_analysis_resltst.R")

load(here("Compare models/leafcutter/step_5_BPseq_Amygdala.Rdata"))

DEGs <-   Lithium_DPseqmodel$gene$Amygdala %>% dplyr::filter(P.Value < 0.005) %>% pull(Symbol) %>% unique()

clusters <- clusters %>% dplyr::mutate(gene = gsub("<i>|</i>", "", gene, fixed = FALSE)) %>% 
  dplyr::mutate(DEGs = if_else(gene %in% DEGs , "yes" , "no"))

load(here("Compare models/results/Lithium results DPseqmodel.RDS"))

geneplot <-  make_gene_plot(gene_name = "POLR2J2" ,clusterID <- "clu_33263_-" , cluster_list = clusters , introns = introns , introns_to_plot = introns_to_plot , exons_table = exons_table , snp_pos <-  NA , min_exon_length <- 0.5 )

clusterplot <-  make_cluster_plot(cluster_to_plot = "clu_33263_-",exons_table = exons_table ,cluster_ids = cluster_ids , introns =  introns , meta = meta , counts = counts ,snp_pos <-  NA , file_path = here("Compare models/leafcutter/" , paste0(clu , "_graph.pdf")))

clustertable <-  filter_intron_table(introns, "clu_33263_-", toSave=TRUE)

clusters_list =  unique(introns$clusterID)


# filter_intron_table final!

walk(clusters_list, ~filter_intron_table(introns, .x, toSave=TRUE)) 
  
walk( clusters_list , ~make_cluster_plot(cluster_to_plot = .x ,exons_table = exons_table ,cluster_ids = cluster_ids , introns =  introns , meta = meta , counts = counts ,snp_pos <-  NA , file_path = here("Compare models/leafcutter/" , paste0(.x , "_graph.pdf"))))

walk2( clusters$gene , clusters$clusterID , ~make_gene_plot(gene_name =  .x ,clusterID <- .y , cluster_list = clusters , introns = introns , introns_to_plot = introns_to_plot , exons_table = exons_table , snp_pos <-  NA , min_exon_length <- 0.5 ))



pushViewport(viewport(layout = grid.layout(2, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(geneplot, vp = vplayout(1, 1:2))
print(clusterplot, vp = vplayout(2, 1))
print(clustertable, vp = vplayout(2, 2))

grid.arrange(clusterplot,clustertable ,  nrow = 2)



#assigned_genes <-  Lithium_DPseqmodel$exon$Amygdala %>% dplyr::filter(P.Value < 0.005) %>% pull(Symbol) %>% unique()


# Clusters
# un-italicise gene name


#cluster_out <- file.path(outfolder, paste0(code, "_sig_clusters.tsv") )

write.table(clusters, file = here("Compare models/leafcutter/" , paste0("BPseq_Amygdala" , "_sig_clusters.csv")) , sep = ",", 
            quote = FALSE, row.names = FALSE)

# Introns
intron_df <- introns %>%
  mutate(coord = paste0(chr, ":", start, "-", end)) %>%
  mutate(ensemblID = gsub("_[0-9]$", "", ensemblID) ) %>%
  dplyr::select(clusterID, gene, ensemblID, coord, verdict, deltaPSI = deltapsi)

#intron_out <- file.path(outfolder, paste0(code, "_sig_introns.tsv") )

write.table(intron_df, file = here("Compare models/leafcutter/" , paste0("BPseq_Amygdala" , "_sig_introns.csv")), sep = ",", 
            quote = FALSE, row.names = FALSE)


reports_to_make <- clusters %>% dplyr::select(gene,  clusterID) %>% 
  setNames(c("gene_name",  "cluster_id"))
  
# Create an output directory
if (!dir.exists("output_reports")) dir.create("output_reports")

# --- Run the pipeline using purrr::pwalk ---
# pwalk is used for functions called for their "side effect" (like saving a file).
# It will iterate through each row of the `reports_to_make` data frame.
pwalk(reports_to_make, function(gene_name, cluster_id) {
  
  output_file <- file.path("output_reports", paste0(gene_name, "_", cluster_id, ".pdf"))
  
  # Call the main wrapper function
  create_full_report(
    gene_name = gene_name,
    cluster_id = cluster_id,
    file_path = output_file,
    # Pass all the necessary data objects
    cluster_list = clusters,
    introns = introns,
    introns_to_plot = introns_to_plot,
    exons_table = exons_table,
    counts = counts,
    meta = meta,
    cluster_ids = cluster_ids
  )
})

message("Pipeline finished. Check the 'output_reports' directory.")

