library(here)
library(tidyverse)
library(WGCNA)
library(SummarizedExperiment)
library(jaffelab)
library(ggvenn)
library(gridExtra)
library(readr)
library(ragg)


### load BPseq results

message("load BPseq results")
 
load(here("results" ,"BP results.RDS")) # load statout BPseq analysis


### load Lithim exposure results BPseq DPseq models #####

message("load Lithium  exposure DPseq model")

load(here("results" ,"Lithium results DPseqmodel.RDS"))  # load DPseq model  lithium analysis


#### Load WGCNA results and prepare it ####

message("Load WGCNA results and prepare it")

load(here("data" , "zandiHypde_bipolar_rseGene_n511.rda"), verbose = TRUE)
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, 'Length')
geneIndex = rowMeans(assays(rse_gene)$rpkm) > 0.25  ## both regions
rse_gene = rse_gene[geneIndex,]


load(here("data","constructed_network_signed_bicor.rda"), verbose = TRUE)
# net_list
# net
# fNames

nrow(rse_gene) == length(net$colors)
# get colors
net$colorsLab = labels2colors(net$colors)

colorDat = data.frame(num = net$colors, col = net$colorsLab, 
                      stringsAsFactors=FALSE)
colorDat$Label = paste0("ME", colorDat$num)
colorDat = colorDat[order(colorDat$num),]
colorDat = colorDat[!duplicated(colorDat$num),]
colorDat$numGenes = table(net$colors)[as.character(colorDat$num)]

colorDat


 # Because we are only interested in the Amygdala gene results, the i and j would be 1
Compare_models <-  function(outlist_name,  module , type , i , j ) {
   # Retrieve the data object from its name string
  Outlist <- get(outlist_name)

  message("Starting analysis for Red module")

  # Ensure output directories exist before trying to save
  #dir.create(here("graphs", model, "wgcna"), showWarnings = FALSE, recursive = TRUE)
  #dir.create(here("graphs", model, "t_stat"), showWarnings = FALSE, recursive = TRUE)
  #dir.create(here("results/wgcna_Enrichment/", model), showWarnings = FALSE, recursive = TRUE)

  if ( module == "T" ) {

  message("Performing Red Module analysis")

  module_id <-  which(net$colorsLab == "red")
  
  ensemblID_module <-  rowData(rse_gene)$ensemblID[module_id]
  
  BP_results <-  Outlist_BP[[i]][[j]]
  Li_results <-  Outlist[[i]][[j]]
    
  colnames(BP_results) <-  paste0(colnames(BP_results) , "_BP" )
  colnames(Li_results) <-  paste0(colnames(Li_results) , "_Li" )
  
    message("Joining BP and Lithium results")

  All_results <- BP_results %>% right_join(Li_results , by = c("ensemblID_BP" = "ensemblID_Li")) %>%
       dplyr::mutate("Group" = if_else( (adj.P.Val_BP < 0.05  &  P.Value_Li < 0.005 ), "Both DEf" ,
                                     if_else( (adj.P.Val_BP < 0.05  &  P.Value_Li >= 0.005 ), "BPseq DEf" ,
                                              if_else( (adj.P.Val_BP >= 0.05  &  P.Value_Li < 0.005 ), "Lithium DEf" , "None" ))) ,
                     "Module" = if_else(ensemblID_BP %in%  ensemblID_module , "Yes" , "NO"))
   


    message("Creating plots")
                   
  p1 <-  All_results %>%
    ggplot(aes(t_BP , t_Li , color = Module )) +
    geom_point(alpha = 0.7 ,size = 2) +
    theme_classic()+
    scale_color_brewer(palette="Dark2") +
     labs(title = "Availability in Red module",
       x ="BPseq", y = "Lithium exposure") +
  xlim(-15 , 15) + ylim(-15 , 15)+
  theme(legend.position="none",
        axis.text.x = element_text(face="bold",  
                                   size=14),
        axis.text.y = element_text(face="bold",  
                                   size=14) ,
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=20, face="bold.italic"))
  
  p2 <-  All_results %>%
    ggplot(aes(t_BP , t_Li , color = Group )) +
    geom_point(alpha = 0.7 ,size = 2) +
    theme_classic()+
    scale_color_brewer(palette="Dark2") +
    labs(title ="Amygdala DEGs",
       x ="BPseq", y = "Lithium exposure") +
  xlim(-15 , 15) + ylim(-15 , 15) +
  theme(legend.position=c(0.8, 0.2),
        axis.line.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(face="bold",  
                                   size=14),
        legend.background = element_rect(fill="lightblue",
                                         size=0.5, linetype="solid", 
                                         colour ="darkblue"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14) ,
        axis.title.x =element_blank(),
        plot.title = element_text(size=20, face="bold.italic"))
        
       
    plot_lists <- list(p1, p2)
    bottom <- textGrob("BPseq", gp = gpar(fontsize = 20))
    yleft = textGrob("Lithium exposure", rot=90, gp = gpar(fontsize = 20))


  plot_lists <-  list(p1, p2)

  bottom <- textGrob("BPseq", gp = gpar(fontsize = 20))
   yleft = textGrob("Lithium exposure", rot=90 , gp = gpar(fontsize = 20))
  

  message("Saving files")


  #jpeg(here("graphs" , "DPseq_model", "wgcna" ,    "Red module Compare t stat from BP vs, Lithium expoure.jpeg") , width = 6000, height = 3000 , res = 300)
  
  final_plot_grob <- grid.arrange( grobs=plot_lists , ncol= 2 , bottom = bottom , left = yleft , draw = FALSE)
               
  #dev.off()

   ggsave(
  here("graphs", "DPseq_model", "wgcna", "Red_module_Compare_t_stat.pdf"),
  plot = final_plot_grob,
  width = 12,
  height = 6
   )  
  
     ggsave(
    here("graphs", "DPseq_model", "wgcna", "Red_module_Compare_t_stat.png"),
    plot = final_plot_grob,
    device = ragg::agg_png,
    width = 12,
    height = 6,
    dpi = 300
  )

  #pdf(here("graphs" , "DPseq_model" , "wgcna" ,    "Red module Compare t stat from BP vs, Lithium expoure.pdf") , width = 12, height = 6 )
  
  #grid.arrange( grobs=plot_lists , ncol= 2 , bottom = bottom , left = yleft )
               
  #dev.off()

  while (!is.null(dev.list())) {
    dev.off()
  }

  message("Generated plot Red Module: ", Outlist)
   

 write_rds(All_results , file = here("results/wgcna_Enrichment/" , "DPseq_model" ,    "Red module Compare t stat from BP vs. Lithium expoure.RDS"))
  
   } else {

       message("Performing T-stat comparison analysis")


  BP_results <-  Outlist_BP[[i]][[j]]
  Li_results <-  Outlist[[i]][[j]]  
  
  colnames(BP_results) <-  paste0(colnames(BP_results) , "_BP" )
  colnames(Li_results) <-  paste0(colnames(Li_results) , "_Li" )
  
  message("Joining BP and Lithium results")

  All_results <- BP_results %>% right_join(Li_results , by = c("ensemblID_BP" = "ensemblID_Li")) %>%
       dplyr::mutate("Group" = if_else( (adj.P.Val_BP < 0.05  &  P.Value_Li < 0.005 ), "Both DEf" ,
                                     if_else( (adj.P.Val_BP < 0.05  &  P.Value_Li >= 0.005 ), "BPseq DEf" ,
                                              if_else( (adj.P.Val_BP >= 0.05  &  P.Value_Li < 0.005 ), "Lithium DEf" , "None" )))) 
       
  
     message("Creating plots")

  p1 <-  All_results %>%
    ggplot(aes(t_BP , t_Li , color = Group  )) +
    geom_point(alpha = 0.7 ,size = 2) +
    theme_classic()+
    scale_color_brewer(palette="Dark2") +
     labs(title = "Amygdala",
       x ="BPseq", y = "Lithium exposure") +
  xlim(-15 , 15) + ylim(-15 , 15)+
  theme(legend.position="none",
        axis.text.x = element_text(face="bold",  
                                   size=14),
        axis.text.y = element_text(face="bold",  
                                   size=14) ,
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(size=20, face="bold.italic")
  )

     
  p2 <-  All_results %>%
    ggplot(aes(t_BP , t_Li , color = Group )) +
    geom_point(alpha = 0.7 ,size = 2) +
    theme_classic()+
    scale_color_brewer(palette="Dark2") +
     labs(title ="sACC",
       x ="BPseq", y = "Lithium exposure") +
  xlim(-15 , 15) + ylim(-15 , 15) +
  theme(legend.position=c(0.8, 0.2),
        axis.line.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(face="bold",  
                                   size=14),
        legend.background = element_rect(fill="lightblue",
                                         size=0.5, linetype="solid", 
                                         colour ="darkblue"),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14) ,
        axis.title.x =element_blank(),
        plot.title = element_text(size=20, face="bold.italic")
        )

      
  plot_lists <-  list(p1, p2)

  bottom <- textGrob("BPseq", gp = gpar(fontsize = 20))
   yleft = textGrob("Lithium exposure", rot=90 , gp = gpar(fontsize = 20))
  
      message("Saving files")

  #jpeg(here("graphs" , "DPseq_model" , "t_stat" ,    "Compare t stat from BP vs, Lithium exposure.jpeg") , width = 6000, height = 3000 , res = 300)
  
   final_plot_grob <- grid.arrange( grobs=plot_lists , ncol= 2 , bottom = bottom , left = yleft , draw = FALSE )
               
  
  #dev.off()
   

  #pdf(here("graphs" , "DPseq_model" , "t_stat" ,    "Compare t stat from BP vs, Lithium exposure.pdf") , width = 12, height = 6 )
  
  #grid.arrange( grobs=plot_lists , ncol= 2 , bottom = bottom , left = yleft )
               
  
  #dev.off()


     ggsave(
  here("graphs", "DPseq_model", "t_stat", "Compare t stat from BP vs, Lithium exposure.pdf"),
  plot = final_plot_grob,
  width = 12,
  height = 6
   )  
  
  
    ggsave(
    here("graphs", "DPseq_model", "t_stat", "Compare t stat from BP vs, Lithium exposure.png"),
    plot = final_plot_grob,
    device = ragg::agg_png,
    width = 12,
    height = 6,
    dpi = 300
  )   


   while (!is.null(dev.list())) {
    dev.off()
   }


    message("Generated plot T stat: ", paste( type ,  Outlist))


 write_rds(All_results, file = here("results/wgcna_Enrichment/", "DPseq_model",   paste0(type,  "_",   "Compare t-stat from BP vs. Lithium exposure.RDS")))


  message(paste("Finished analysis for: ", "DPseq_model", " " , type ))
 
 }  
}

plot_params <- tibble::tribble(
  ~module , ~outlist_name,         ~type,      ~i,  ~j,  
  "T",       "Lithium_DPseqmodel",    "Gene",       1 ,  1 ,   
  "F",       "Lithium_DPseqmodel",    "Gene",       1 ,  1 ,  
  "F",       "Lithium_DPseqmodel",    "Gene",       1 ,  2 ,   
  "F",       "Lithium_DPseqmodel",    "Exon",       2 ,  1 ,   
  "F",       "Lithium_DPseqmodel",    "Exon",       2 ,  2 ,   
  "F",       "Lithium_DPseqmodel",    "Junction",   3 ,  1 ,   
  "F",       "Lithium_DPseqmodel",    "Junction",   3 ,  2 ,   
  "F",       "Lithium_DPseqmodel",    "Transcript", 4 ,  1 ,   
  "F",       "Lithium_DPseqmodel",    "Transcript", 4 ,  2   
)

message("Running function")


pwalk(plot_params, Compare_models , .progress = T)

## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()
