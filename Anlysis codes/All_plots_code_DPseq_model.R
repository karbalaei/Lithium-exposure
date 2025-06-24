### Library #####
library(here)
library("ggVennDiagram")
library(purrr)
library("ggvenn")
library(UpSetR)
library(VennDetail)
library(magrittr)
library(plyr)
library(tidyverse)
library(ggplot2)
library(reshape2)


load(file = here::here("results" ,"Lithium results DPseqmodel.RDS"))


#### Functions for extracting gene sets ####
get_signif_genes <- function(outFeature, colname = "common_feature_id", cutoff = 0.005, return_unique = FALSE){
  signif <- outFeature[[colname]][outFeature$P.Value < cutoff]
  if(return_unique) signif <- unique(signif)
  signif <- signif[!is.na(signif)]
  return(signif)
}


get_signif_feature <- function(outFeature, cutoff = 0.005, return_unique = FALSE){
  signif <- rownames(outFeature)[outFeature$P.Value < cutoff]
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

#### Get signif genes####

signif_genes <- map_depth(Lithium_DPseqmodel, 2, ~get_signif_genes(.x, colname = "common_gene_id", return_unique = TRUE))
map_depth(signif_genes, 2, length)


signif_features <- map_depth(Lithium_DPseqmodel, 2, ~get_signif_feature(.x,  return_unique = TRUE))
map_depth(signif_features, 2, length)


signif_genes_flat <- my_flatten(map_depth(signif_genes, 2, unlist))

map_int(signif_genes_flat, length)



signif_feature_flat <- my_flatten(map_depth(signif_features, 2, unlist))

map_int(signif_feature_flat, length)


Amygdala_Genes_All <-   unique(list_c(signif_genes_flat[c(1 , 3, 5 , 7)]))
sACC_Genes_All <-   unique(list_c(signif_genes_flat[c(2 , 4, 6 , 8)]))

List_All <-  list( Amygdala = Amygdala_Genes_All ,  sACC = sACC_Genes_All)

### Venn diagram Genes ####


# --- Step 1: Define all plot parameters in a data frame ---
output_dir <- here::here("graphs", "DPseq_model", "Venn_diagram")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

plot_params <- tibble::tribble(
  ~plot_type, ~file_suffix,             ~data_indices,         ~extra_args,
  "upset",    "Upset graph All",        1:8,                   list(nsets = 8),
  "upset",    "Upset graph Amygdala",   c(1, 3, 5, 7),         list(nsets = 4),
  "upset",    "Upset graph sACC",       c(2, 4, 6, 8),         list(nsets = 4),
  "ggvenn",   "Venn_diagram_Amygdala phenotypes", c(1, 3, 5, 7), list(set_name_size = 3.5, text_size = 3),
  "ggvenn",   "Venn_diagram_sACC phenotypes",     c(2, 4, 6, 8), list(set_name_size = 3.5, text_size = 3),
  "ggvenn",   "Venn_diagram_gene",      c(1, 2),               list(set_name_size = 4, text_size = 4),
  "ggvenn",   "Venn_diagram_exon",      c(3, 4),               list(set_name_size = 4, text_size = 4),
  "ggvenn",   "Venn_diagram_junction",  c(5, 6),               list(set_name_size = 4, text_size = 4),
  "ggvenn",   "Venn_diagram_transcript",c(7, 8),               list(set_name_size = 4, text_size = 4)
)

# --- Step 2: Create a robust, server-safe plotting function (FINAL VERSION) ---
create_plot_non_interactive <- function(plot_type, file_suffix, data_indices, extra_args) {
  plot_data <- signif_genes_flat[data_indices]
  
  if (plot_type == "ggvenn") {
    ggvenn_args <- c(list(plot_data), extra_args, list(stroke_size = 0.5, fill_color = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")))
    p <- rlang::exec(ggvenn, !!!ggvenn_args)
    ggsave(filename = file.path(output_dir, paste0(file_suffix, ".pdf")), plot = p, width = 8, height = 8)
    ggsave(filename = file.path(output_dir, paste0(file_suffix, ".jpeg")), plot = p, width = 10, height = 8.33, dpi = 600)
    
  } else if (plot_type == "upset") {
    upset_args <- c(list(fromList(plot_data)), extra_args, list(matrix.color = "#1b9e77", main.bar.color = "#7570b3", sets.bar.color = "#e7298a"))
    
    # --- PDF Version ---
    pdf(file.path(output_dir, paste0(file_suffix, ".pdf")), width = 8, height = 8)
    print(do.call('upset', upset_args))
    dev.off()
    
    # --- JPEG Version ---
    jpeg(file.path(output_dir, paste0(file_suffix, ".jpeg")), width = 6000, height = 5000, res = 600)
    print(do.call('upset', upset_args))
    dev.off()
  }
  
  message("Generated plot: ", file_suffix)
}


# --- Step 3: Iterate and create all plots ---
pwalk(plot_params, create_plot_non_interactive)

# --- Handle the special case for 'List_All' ---
message("Generating final plot for List_All...")
ggvenn_all_plot <- ggvenn(
  List_All, stroke_size = 0.5, fill_color = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"),
  set_name_size = 3.5, text_size = 3
)
ggsave(
  filename = file.path(output_dir, "Venn_diagram_Amygdala sACC All.pdf"), 
  plot = ggvenn_all_plot, width = 8, height = 8
)
ggsave(
  filename = file.path(output_dir, "Venn_diagram_Amygdala sACC All.jpeg"), 
  plot = ggvenn_all_plot, width = 10, height = 8.33, dpi = 600
)
message("All plots generated successfully.")

### Venn diagram features ####

# --- Step 1: Define all plot parameters in a data frame ---
output_dir <- here::here("graphs", "DPseq_model", "Venn_diagram")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

plot_params_features <- tibble::tribble(
  ~file_suffix,             ~data_indices,         ~extra_args,
  "Venn_diagram_gene",      c(1, 2),               list(set_name_size = 4, text_size = 4),
  "Venn_diagram_exon",      c(3, 4),               list(set_name_size = 4, text_size = 4),
  "Venn_diagram_junction",  c(5, 6),               list(set_name_size = 4, text_size = 4),
  "Venn_diagram_transcript",c(7, 8),               list(set_name_size = 4, text_size = 4)
)

# --- Step 2: Create a robust, server-safe plotting function (FINAL VERSION) ---
create_plot_non_interactive_feature <- function(file_suffix, data_indices, extra_args) {
  plot_data <- signif_feature_flat[data_indices]
  
    ggvenn_args <- c(list(plot_data), extra_args, list(stroke_size = 0.5, fill_color = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")))
    p <- rlang::exec(ggvenn, !!!ggvenn_args)
    ggsave(filename = file.path(output_dir, paste0(file_suffix, "_feature.pdf")), plot = p, width = 8, height = 8)
    ggsave(filename = file.path(output_dir, paste0(file_suffix, "_feature.jpeg")), plot = p, width = 10, height = 8.33, dpi = 600)
    
  message("Generated plot: ", file_suffix)
}


# --- Step 3: Iterate and create all plots ---
pwalk(plot_params_features, create_plot_non_interactive_feature)

# --- Handle the special case for 'List_All' ---
message("Generating final plot for List_All...")
ggvenn_all_plot <- ggvenn(
  List_All, stroke_size = 0.5, fill_color = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"),
  set_name_size = 3.5, text_size = 3
)
ggsave(
  filename = file.path(output_dir, "Venn_diagram_Amygdala sACC All.pdf"), 
  plot = ggvenn_all_plot, width = 8, height = 8
)
ggsave(
  filename = file.path(output_dir, "Venn_diagram_Amygdala sACC All.jpeg"), 
  plot = ggvenn_all_plot, width = 10, height = 8.33, dpi = 600
)
message("All plots generated successfully.")

##### pdf ####


pdf(here::here("graphs" , "DPseq_model" , "Venn_diagram" , "Upset graph All.pdf") , width = 8 , height = 8)

upset(fromList(signif_genes_flat), nsets = 8 , matrix.color = "#1b9e77" , main.bar.color = "#7570b3" , sets.bar.color	= "#e7298a")

dev.off()



pdf(here::here("graphs" , "DPseq_model" , "Venn_diagram" , "Upset graph Amygdala.pdf") , width = 8 , height = 8)

upset(fromList(signif_genes_flat[c(1 , 3, 5 , 7)]), nsets = 4 , matrix.color = "#1b9e77" , main.bar.color = "#7570b3" , sets.bar.color	= "#e7298a"	)

dev.off()


pdf(here::here("graphs" ,  "DPseq_model" ,  "Venn_diagram" , "Upset graph sACC.pdf") , width = 8 , height = 8)

upset(fromList(signif_genes_flat[c(2 , 4, 6 , 8)]), nsets = 4 , matrix.color = "#1b9e77" , main.bar.color = "#7570b3" , sets.bar.color	= "#e7298a"	)

dev.off()

#v.table <- venn(signif_genes_flat[c(1 , 2)])

pdf(here::here("graphs" ,  "DPseq_model" , "Venn_diagram" , "Venn_diagram_Amygdala phenotypes.pdf") , width = 8 , height = 8)

ggvenn(
  signif_genes_flat[c(1 , 3, 5 , 7)], 
  stroke_size = 0.5 , fill_color = c("#66c2a5" , "#fc8d62" , "#8da0cb", "#e78ac3"),
  set_name_size = 3.5 , text_size = 3
)

dev.off()

pdf(here::here("graphs" ,  "DPseq_model" , "Venn_diagram" , "Venn_diagram_sACC phenotypes.pdf") , width = 8 , height = 8)

ggvenn(
  signif_genes_flat[c(2 , 4, 6 , 8)], 
  stroke_size = 0.5 , fill_color = c("#66c2a5" , "#fc8d62" , "#8da0cb", "#e78ac3"),
  set_name_size = 3.5 , text_size = 3
)

dev.off()


pdf(here::here("graphs" ,  "DPseq_model" , "Venn_diagram" , "Venn_diagram_gene.pdf") , width = 8 , height = 8)

ggvenn(
  signif_genes_flat[c(1 ,2)], 
  stroke_size = 0.5 , fill_color = c("#66c2a5" , "#fc8d62" , "#8da0cb", "#e78ac3"),
  set_name_size = 4 , text_size = 4
)

dev.off()

pdf(here::here("graphs" ,  "DPseq_model" , "Venn_diagram" , "Venn_diagram_exon.pdf") , width = 8 , height = 8)

ggvenn(
  signif_genes_flat[c(3, 4)], 
  stroke_size = 0.5 , fill_color = c("#66c2a5" , "#fc8d62" , "#8da0cb", "#e78ac3"),
  set_name_size = 4 , text_size = 4
)

dev.off()



pdf(here::here("graphs" ,  "DPseq_model" , "Venn_diagram" , "Venn_diagram_junction.pdf") , width = 8 , height = 8)

ggvenn(
  signif_genes_flat[c(5, 6)], 
  stroke_size = 0.5 , fill_color = c("#66c2a5" , "#fc8d62" , "#8da0cb", "#e78ac3"),
  set_name_size = 4 , text_size = 4
)

dev.off()

pdf(here::here("graphs" ,  "DPseq_model" , "Venn_diagram" , "Venn_diagram_transcript.pdf") , width = 8 , height = 8)

ggvenn(
  signif_genes_flat[c(7, 8)], 
  stroke_size = 0.5 , fill_color = c("#66c2a5" , "#fc8d62" , "#8da0cb", "#e78ac3"),
  set_name_size = 4 , text_size = 4
)

dev.off()

pdf(here::here("graphs" ,  "DPseq_model" , "Venn_diagram" , "Venn_diagram_Amygdala sACC All.pdf") , width = 8 , height = 8)

ggvenn(
  List_All , 
  stroke_size = 0.5 , fill_color = c("#66c2a5" , "#fc8d62" , "#8da0cb", "#e78ac3"),
  set_name_size = 3.5 , text_size = 3
)


dev.off()

#### Jpeg #####


jpeg(here::here("graphs" ,  "DPseq_model" , "Venn_diagram" , "Upset graph All.jpeg") , width = 6000 , height = 5000 , res = 600)

upset(fromList(signif_genes_flat), nsets = 8 , matrix.color = "#1b9e77" , main.bar.color = "#7570b3" , sets.bar.color	= "#e7298a"	)

dev.off()



jpeg(here::here("graphs" ,  "DPseq_model" , "Venn_diagram" , "Upset graph Amygdala.jpeg") , width = 6000 , height = 5000 , res = 600)

upset(fromList(signif_genes_flat[c(1 , 3, 5 , 7)]), nsets = 4 , matrix.color = "#1b9e77" , main.bar.color = "#7570b3" , sets.bar.color	= "#e7298a"	)

dev.off()


jpeg(here::here("graphs" ,  "DPseq_model" , "Venn_diagram" , "Upset graph sACC.jpeg") , width = 6000 , height = 5000 , res = 600)

upset(fromList(signif_genes_flat[c(2 , 4, 6 , 8)]), nsets = 4 , matrix.color = "#1b9e77" , main.bar.color = "#7570b3" , sets.bar.color	= "#e7298a"	)

dev.off()



jpeg(here::here("graphs" ,  "DPseq_model" , "Venn_diagram" , "Venn_diagram_Amygdala phenotypes.jpeg") , width = 6000 , height = 5000 , res = 600)

ggvenn(
  signif_genes_flat[c(1 , 3, 5 , 7)], 
  stroke_size = 0.5 , fill_color = c("#66c2a5" , "#fc8d62" , "#8da0cb", "#e78ac3"),
  set_name_size = 3.5 , text_size = 3
)

dev.off()

jpeg(here::here("graphs" ,  "DPseq_model" , "Venn_diagram" , "Venn_diagram_sACC phenotypes.jpeg") , width = 6000 , height = 5000 , res = 600)

ggvenn(
  signif_genes_flat[c(2 , 4, 6 , 8)], 
  stroke_size = 0.5 , fill_color = c("#66c2a5" , "#fc8d62" , "#8da0cb", "#e78ac3"),
  set_name_size = 3.5 , text_size = 3
)

dev.off()


jpeg(here::here("graphs" ,  "DPseq_model" , "Venn_diagram" , "Venn_diagram_gene.jpeg") , width = 6000 , height = 5000 , res = 600)

ggvenn(
  signif_genes_flat[c(1 ,2)], 
  stroke_size = 0.5 , fill_color = c("#66c2a5" , "#fc8d62" , "#8da0cb", "#e78ac3"),
  set_name_size = 4 , text_size = 4
)

dev.off()

jpeg(here::here("graphs" ,  "DPseq_model" , "Venn_diagram" , "Venn_diagram_exon.jpeg") , width = 6000 , height = 5000 , res = 600)

ggvenn(
  signif_genes_flat[c(3, 4)], 
  stroke_size = 0.5 , fill_color = c("#66c2a5" , "#fc8d62" , "#8da0cb", "#e78ac3"),
  set_name_size = 4 , text_size = 4
)

dev.off()



jpeg(here::here("graphs" ,  "DPseq_model" , "Venn_diagram" , "Venn_diagram_junction.jpeg") , width = 6000 , height = 5000 , res = 600)

ggvenn(
  signif_genes_flat[c(5, 6)], 
  stroke_size = 0.5 , fill_color = c("#66c2a5" , "#fc8d62" , "#8da0cb", "#e78ac3"),
  set_name_size = 4 , text_size = 4
)

dev.off()

jpeg(here::here("graphs" ,  "DPseq_model" , "Venn_diagram" , "Venn_diagram_transcript.jpeg") , width = 6000 , height = 5000 , res = 600)

ggvenn(
  signif_genes_flat[c(7, 8)], 
  stroke_size = 0.5 , fill_color = c("#66c2a5" , "#fc8d62" , "#8da0cb", "#e78ac3"),
  set_name_size = 4 , text_size = 4
)

dev.off()

jpeg(here::here("graphs" ,  "DPseq_model" , "Venn_diagram" , "Venn_diagram_Amygdala sACC All.jpeg") , width = 6000 , height = 5000 , res = 600)

ggvenn(
  List_All , 
  stroke_size = 0.5 , fill_color = c("#66c2a5" , "#fc8d62" , "#8da0cb", "#e78ac3"),
  set_name_size = 3.5 , text_size = 3
)


dev.off()

### Bar plots #####

##### features frequencies ####
Gene <- venndetail(signif_genes_flat[c( 1, 2)]) 
Exon <- venndetail(signif_genes_flat[c( 3, 4)])
Junction <-  venndetail(signif_genes_flat[c( 5, 6)])
Transcript <- venndetail(signif_genes_flat[c( 7, 8)])

Gene_df <-   as.data.frame(detail(Gene)) %>% set_colnames("Gene") %>% set_rownames(c("Common" , "sACC" , "Amygdala"))
Exon_df <-   as.data.frame(detail(Exon)) %>% set_colnames("Exon") %>% set_rownames(c("Common" , "sACC" , "Amygdala"))
Junction_df <-   as.data.frame(detail(Junction)) %>% set_colnames("Junction") %>% set_rownames(c("Common" , "sACC" , "Amygdala"))
Transcript_df <-   as.data.frame(detail(Transcript)) %>% set_colnames("Transcript") %>% set_rownames(c("Common" , "sACC" , "Amygdala"))

Type_fre_df <-  t(cbind(Gene_df , Exon_df , Junction_df , Transcript_df)) %>% as.data.frame() %>% rownames_to_column("Type") %>%
  pivot_longer(cols = Common: Amygdala , names_to = "Region" , values_to = "Fre")



p1 <-  Type_fre_df %>%
  mutate(Type = fct_reorder(Type , Fre, .fun='mean' , .desc =F )) %>% 
  ggplot( aes(x=Type, y=Fre , fill = Region)) 

pdf(here::here("graphs" ,  "DPseq_model"  , "Barplot_DEf.pdf") , width = 9 , height = 6)


p1 +   geom_bar(stat='identity' ,colour = "black" ,  position = position_dodge() )+
  geom_text(stat='identity'  , position = position_dodge(width = 0.8)  , aes(label=Fre), vjust=-1, size=3 )+
  scale_fill_manual(values=c( "#66c2a5" , "grey" ,  "#fc8d62")) +
  theme_classic()  + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  xlab("") + ylab("No. significant features")

dev.off()

jpeg(here::here("graphs" ,  "DPseq_model"  , "Barplot_DEf.jpg") , width = 6000 , height = 4000 , res = 600)


p1 +   geom_bar(stat='identity' ,colour = "black" ,  position = position_dodge() )+
  geom_text(stat='identity'  , position = position_dodge(width = 0.8)  , aes(label=Fre), vjust=-1, size=3 )+
  scale_fill_manual(values=c( "#66c2a5" , "grey" ,  "#fc8d62")) +
  theme_classic()  + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  xlab("") + ylab("No. significant features")

dev.off()


#### Gene type ######


Amygdala_gene_type <-  Lithium_DPseqmodel$gene$Amygdala  %>% 
  dplyr::filter(P.Value < 0.005) %>% 
  dplyr::select(c(common_gene_id ,gene_type )) %>% 
  mutate("Region" = "Amygdala")

sACC_gene_type <-  Lithium_DPseqmodel$gene$sACC %>% 
  dplyr::filter(P.Value < 0.005) %>% 
  dplyr::select(c(common_gene_id ,gene_type )) %>% 
  mutate("Region" = "sACC")

Gene_type <-  rbind(Amygdala_gene_type , sACC_gene_type) 


unique(Gene_type$gene_type)

Common_genes <-  intersect(sACC_gene_type$common_gene_id , Amygdala_gene_type$common_gene_id)

Gene_type$Region <-replace (Gene_type$Region, Gene_type$common_gene_id %in% Common_genes, "Common") 

Gene_type <-   Gene_type %>% distinct()

Gene_type$Region <- factor(Gene_type$Region , levels = c("Amygdala"  , "Common" ,  "sACC" ))
Gene_type$Region %>% unique()

Other_RNA <-  c("snRNA","misc_RNA"  ,"sense_overlapping"  , "rRNA"  , "snoRNA" )

pseudogene_RNA <-  c("transcribed_unprocessed_pseudogene" , "transcribed_processed_pseudogene"  , "unprocessed_pseudogene"  ,          
                     "processed_pseudogene" )


replace_RNA_vector <-  rep(c("Other_RNA"  , "pseudogene_RNA" ), c(length(Other_RNA),length(pseudogene_RNA)))


Gene_type$gene_type <-  mapvalues(Gene_type$gene_type ,  c(Other_RNA , pseudogene_RNA) , replace_RNA_vector )


unique(Gene_type$gene_type)

Gene_type$gene_type <-  factor(Gene_type$gene_type , levels = c( "protein_coding" , "antisense" , "pseudogene_RNA" , "lincRNA" , "Other_RNA" ))


### Final graph #####
g1 <-  Gene_type %>% 
  mutate(gene_type = fct_reorder(gene_type , gene_type, .fun='length' , .desc = T )) %>% 
  ggplot(aes(x= gene_type ,fill = Region ))

pdf(here::here("graphs" ,  "DPseq_model"  , "Barplot_gene_type.pdf") , width = 9 , height = 6)

g1 + geom_bar(position = "dodge" , stat='count')+
  geom_text(stat = "count",
            aes(
              label = after_stat(count)
            ),
            position = position_dodge(width = 1),
            color = "black",
            size = 4,
            vjust = -1
  ) +
  scale_fill_manual(values=c( "#66c2a5" , "grey" , "#fc8d62")) +
  theme_classic()  + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  xlab("") + ylab("No. significant features") +scale_x_discrete(breaks=c("protein_coding","pseudogene_RNA","lincRNA" , "Other_RNA" , "antisense"),
                                                                labels=c("Coding", "Pseudogene", "lncRNA" , "Other RNA" , "Antisense"))
dev.off()


jpeg(here::here("graphs" ,  "DPseq_model"  , "Barplot_gene_type.jpg") , width = 6000 , height = 4000 , res = 600)

g1 + geom_bar(position = "dodge" , stat='count')+
  geom_text(stat = "count",
            aes(
              label = after_stat(count)
            ),
            position = position_dodge(width = 1),
            color = "black",
            size = 4,
            vjust = -1
  ) +
  scale_fill_manual(values=c( "#66c2a5" , "grey" , "#fc8d62")) +
  theme_classic()  + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  xlab("") + ylab("No. significant features") +scale_x_discrete(breaks=c("protein_coding","pseudogene_RNA","lincRNA" , "Other_RNA" , "antisense"),
                                                                labels=c("Coding", "Pseudogene", "lncRNA" , "Other RNA" , "Antisense"))
dev.off()



dat = dcast (Gene_type, gene_type ~ Region, fun.aggregate = length)
dat.melt = melt(dat, id.vars = "gene_type", measure.vars = c("Amygdala"  , "Common" ,  "sACC"))
dat.melt

pdf(here::here("graphs" ,  "DPseq_model"  , "Barplot_gene_type2.pdf") , width = 9 , height = 6)

ggplot(dat.melt, aes(x = gene_type,y = value, fill = variable)) + 
  geom_bar(stat = "identity", colour = "black", position = position_dodge(width = .8), width = 0.7) +
  #ylim(0, 14) +
  geom_text(aes(label = value), position = position_dodge(width = .8), vjust = -0.5) +
  scale_fill_manual(values=c( "#66c2a5" , "grey" , "#fc8d62")) +
  theme_classic()  + theme(axis.text.x = element_text(size = 15 , angle = 45, vjust = 0.5, hjust=0.5),
                           axis.title.y = element_text(size = 15 , vjust = 3, hjust=0.5)   ,       axis.text.y = element_text(size = 15 )) +
  xlab("") + ylab("\nNo. significant features\n") +scale_x_discrete(breaks=c("protein_coding","pseudogene_RNA","lincRNA" , "Other_RNA" , "antisense"),
                                                                    labels=c("Coding", "Pseudogene", "lncRNA" , "Other RNA" , "Antisense"))
dev.off()

jpeg(here::here("graphs" ,  "DPseq_model"  , "Barplot_gene_type2.jpg") , width = 6000 , height = 4000 , res = 600)

ggplot(dat.melt, aes(x = gene_type,y = value, fill = variable)) + 
  geom_bar(stat = "identity", colour = "black", position = position_dodge(width = .8), width = 0.7) +
  #ylim(0, 14) +
  geom_text(aes(label = value), position = position_dodge(width = .8), vjust = -0.5) +
  scale_fill_manual(values=c( "#66c2a5" , "grey" , "#fc8d62")) +
  theme_classic()  + theme(axis.text.x = element_text(size = 15 , angle = 45, vjust = 0.5, hjust=0.5),
                           axis.title.y = element_text(size = 15 , vjust = 3, hjust=0.5)   ,       axis.text.y = element_text(size = 15 )) +
  xlab("") + ylab("\nNo. significant features\n") +scale_x_discrete(breaks=c("protein_coding","pseudogene_RNA","lincRNA" , "Other_RNA" , "antisense"),
                                                                    labels=c("Coding", "Pseudogene", "lncRNA" , "Other RNA" , "Antisense"))
dev.off()

