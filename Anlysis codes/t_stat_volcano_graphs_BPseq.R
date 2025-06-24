 ### Compare BPseq and DPseq models results #####
rm(list= ls())

library(here)
library(readr)
library(stringr)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(dplyr)


here()


Lithium_gene <- read_csv(here("results" , "Lithium_BPseq_model_gene.csv")) %>% dplyr::rename(Id = ...1)
Lithium_exon <- read_csv(here("results" , "Lithium_BPseq_model_exon.csv")) %>% dplyr::rename(Id = ...1)
Lithium_Jxn <- read_csv(here("results" , "Lithium_BPseq_model_Jxn.csv")) %>% dplyr::rename(Id = ...1)
Lithium_Transcript <- read_csv(here("results" , "Lithium_BPseq_model_Transcript.csv")) %>% dplyr::rename(Id = ...1)


scale_fill_Reza <- function(...){
  ggplot2:::manual_scale(
    'color', 
    values = setNames(c('#e41a1c', '#377eb8', '#984ea3', '#4daf4a'), c("Both Region", "Amygdala","sACC", "None" )), 
    ...
  )
}


Volcano_t_graphs <-  function(df , type , cutoff , xlim1 , xlim2 , ylim1 ) {
  
    df_t <-  df %>% 
    dplyr::mutate("DEp" = if_else( (P.Value_Amygdala < cutoff  &  P.Value_sACC < cutoff ), "Both Region" ,
                                   if_else( (P.Value_Amygdala < cutoff  &  P.Value_sACC >= cutoff ), "Amygdala" ,
                                            if_else( (P.Value_Amygdala >= cutoff  &  P.Value_sACC < cutoff ), "sACC" , "None" ))))
  
  
 
 df_t%>%
    ggplot(aes(t_Amygdala , t_sACC , color = DEp )) +
    geom_point(alpha = 0.8) +
    theme_classic()+
    scale_fill_Reza() +
    #scale_color_brewer(palette="Dark2") +
    labs(title = type , subtitle ="Comapring T stat, highlighting DEp",
         x ="t stat from Amygdala", y = "t stat from sACC") +
    xlim(-15 , 15) +
    ylim(-15 , 15)
  
  ggsave(filename = here("graphs",  "BPseq_model" , "t_stat" ,  paste0("Comparing t stat from Amygdala and sACC regions " , type ,  " .jpeg")) , width = 8 , height = 8 )
  
  # pdf(here("graphs", "t_stat" ,  paste0("Comparing t stat from Amygdala and sACC regions " , type ,  " .pdf")) , width = 8 , height = 8 )
  # 
  # t_graph
  # 
  # dev.off()
  
  
 df_V_Amygdala <-   df%>%
    dplyr::mutate("DEp" = if_else( P.Value_Amygdala < cutoff  , "Amygdala" , "None" ) , 
   Expression = case_when(logFC_Amygdala > 0  & P.Value_Amygdala <cutoff ~ "Up-regulated",
                          logFC_Amygdala < 0  & P.Value_Amygdala <cutoff ~ "Down-regulated",
                          TRUE ~ "Unchanged"))
  
 top_genes_Amygdala <- bind_rows(
   df_V_Amygdala %>% 
    dplyr::filter (Expression == 'Up-regulated') %>% 
     dplyr::arrange(P.Value_Amygdala, desc(abs(logFC_Amygdala))) %>% 
     head(topgenes , n = 5),
   df_V_Amygdala %>% 
     dplyr::filter(Expression == 'Down-regulated') %>% 
     dplyr::arrange(P.Value_Amygdala, desc(abs(logFC_Amygdala))) %>% 
     head(topgenes , n = 5)
 )
 
 
 ggplot(df_V_Amygdala, aes(logFC_Amygdala, -log(P.Value_Amygdala,10))) +
   geom_point(aes(color = Expression), size = 1) +
   xlab(expression("log"[2]*"FC")) + 
   ylab(expression("-log"[10]*"p-value")) +
   scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
   guides(colour = guide_legend(override.aes = list(size=1.5))) +
   geom_label_repel(data = top_genes_Amygdala,
                    mapping = aes(label = common_gene_id_Amygdala),
                    size = 2) + theme_classic()+
   theme(legend.position= "none") +
   xlim(xlim1 , xlim2) + ggtitle( paste0("Amygdala " , type)) +
   ylim(0 , ylim1)
   
 ggsave( here("graphs",  "BPseq_model" , "volcano_plot" ,  paste0("Amygdala region " , type ,  " .jpeg")) , width = 8 , height = 8 )
 
 # pdf(here("graphs", "volcano_plot" ,  paste0("Amygdala region " , type ,  " .pdf")) , width = 8 , height = 8 )
 # 
 # Volcano_Amygdala
 # 
 # dev.off()
 
 
 
 df_V_sACC <-   df%>%
   dplyr::mutate("DEp" = if_else( P.Value_sACC < cutoff  , "sACC" , "None" ) , 
                 Expression = case_when(logFC_sACC > 0  & P.Value_sACC <cutoff ~ "Up-regulated",
                                        logFC_sACC < 0  & P.Value_sACC <cutoff ~ "Down-regulated",
                                        TRUE ~ "Unchanged"))
 
 top_genes_sACC <- bind_rows(
   df_V_sACC %>% 
     dplyr::filter (Expression == 'Up-regulated') %>% 
     dplyr::arrange(P.Value_sACC, desc(abs(logFC_sACC))) %>% 
     head(topgenes , n = 5),
   df_V_sACC %>% 
     dplyr::filter(Expression == 'Down-regulated') %>% 
     dplyr::arrange(P.Value_sACC, desc(abs(logFC_sACC))) %>% 
     head(topgenes , n = 5)
 )
 
 
  ggplot(df_V_sACC, aes(logFC_sACC, -log(P.Value_sACC,10))) +
   geom_point(aes(color = Expression), size = 1) +
   xlab(expression("log"[2]*"FC")) + 
   ylab(expression("-log"[10]*"p-value")) +
   scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
   guides(colour = guide_legend(override.aes = list(size=1.5))) +
   geom_label_repel(data = top_genes_sACC,
                    mapping = aes(label = common_gene_id_sACC),
                    size = 2) + theme_classic()+
   theme(legend.position= "none") +
   xlim(xlim1 , xlim2) + ggtitle( paste0("sACC " , type)) +
   ylim(0 , ylim1)
 
 ggsave( here("graphs",   "BPseq_model" , "volcano_plot" ,  paste0("sACC region " , type ,  " ..jpeg")) , width = 8 , height = 8 )
 
 # pdf(here("graphs", "volcano_plot" ,  paste0("sACC region " , type ,  " .pdf")) , width = 8 , height = 8 )
 # 
 # Volcano_sACC
 # 
 # dev.off()
 
}

Volcano_t_graphs(Lithium_gene , "Gene" , 0.005 , -8 , 8 , 7)
Volcano_t_graphs(Lithium_exon , "Exon" , 0.005 , -10 , 10  , 12)
Volcano_t_graphs(Lithium_Jxn , "Junction" , 0.005 , -10 , 10  , 15)
Volcano_t_graphs(Lithium_Transcript , "Transcript" , 0.005 , -5 , 5  , 6)


## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()
