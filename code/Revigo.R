BiocManager::install("rrvgo")

library(rrvgo)
library(here)


go_analysis <- read.delim(system.file("extdata/example.txt", package="rrvgo"))

simMatrix <- calculateSimMatrix(go_analysis$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")






load("D:/BPseq_data/Old/enrichGO.rda")

to_remove <-  which(endsWith(names(Enrichment_results)   , "ALL"))

Enrichment_results <-  Enrichment_results[-(c(to_remove))]

# 
# 
# Enrichment_results <-  map2(plan$gene_list , plan$ontology,  ~enrichGO(.x,
#                                                                        universe = u, OrgDb = org.Hs.eg.db,
#                                                                        ont = .y , pAdjustMethod = "BH",
#                                                                        pvalueCutoff  = .2, qvalueCutoff  = .5,
#                                                                        readable= TRUE)   )


nams_Enrichment_results <-  names(Enrichment_results) %>% str_split("_") %>%   map(., ~paste(.x, collapse = ",")) %>% 
  unlist


plan <-  read.csv(textConnection(nams_Enrichment_results), header = F , col.names = c("data" , "GO"))




simMatrix_results <-  map2(Enrichment_results , plan$GO,  ~calculateSimMatrix(.x$ID, orgdb="org.Hs.eg.db",
                                                                     ont=.y,method="Rel"  ) , .progress = T)


scores_results <-  map(Enrichment_results , ~setNames(.x$qvalue,  .x$ID ) , .progress = T)

reducedTerms_Rel_results <-  map2(simMatrix_results , scores_results ,  ~reduceSimMatrix(.x, .y ,  threshold=0.7,
                                                                                      orgdb="org.Hs.eg.db") , .progress = T)


here::set_here(".")



map2(Enrichment_results, names(Enrichment_results), ~dotplot(.x, title = .y , showCategory = 25))

jpeg("Revigo_heatmap_Rel.jpg" , res = 600 , width = 7000 , height = 5000)

heatmapPlot(simMatrix_Rel,
            reducedTerms_Rel,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=5)

dev.off()



reducedTerms_Rel <- reduceSimMatrix(simMatrix_Rel,
                                    scores,
                                    threshold=0.7,
                                    orgdb="org.Hs.eg.db")




simMatrix_Rel <- calculateSimMatrix(go_analysis$ID,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")



reducedTerms_Rel <- reduceSimMatrix(simMatrix_Rel,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

jpeg("test_heatmap_Rel.jpg" , res = 600 , width = 7000 , height = 5000)

heatmapPlot(simMatrix_Rel,
            reducedTerms_Rel,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=5)

dev.off()


jpeg("test_scatterPlot_Rel.jpg" , res = 300 , width = 10000 , height = 6000)

scatterPlot(simMatrix_Rel, reducedTerms_Rel)

dev.off()


jpeg("test_treemapPlot_Rel.jpg" , res = 600 , width = 8000 , height = 5000)

treemapPlot(reducedTerms_Rel)

dev.off()



jpeg("test_wordcloudPlot_Rel.jpg" , res = 600 , width = 5000 , height = 5000)

wordcloudPlot(reducedTerms_Rel, min.freq=1, colors="black")
dev.off()




simMatrix_Resnik  <- calculateSimMatrix(go_analysis$ID,
                                    orgdb="org.Hs.eg.db",
                                    ont="BP",
                                    method="Resnik")


reducedTerms_Resnik <- reduceSimMatrix(simMatrix_Resnik,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")

jpeg("test_heatmap_Resnik.jpg" , res = 600 , width = 7000 , height = 5000)

heatmapPlot(simMatrix_Resnik,
            reducedTerms_Resnik,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=5)

dev.off()


jpeg("test_scatterPlot_Resnik.jpg" , res = 300 , width = 10000 , height = 6000)

scatterPlot(simMatrix_Resnik, reducedTerms_Resnik)

dev.off()


jpeg("test_treemapPlot_Resnik.jpg" , res = 600 , width = 8000 , height = 5000)

treemapPlot(reducedTerms_Resnik)

dev.off()



jpeg("test_wordcloudPlot_Resnik.jpg" , res = 600 , width = 5000 , height = 5000)

wordcloudPlot(reducedTerms_Resnik, min.freq=1, colors="black")
dev.off()


