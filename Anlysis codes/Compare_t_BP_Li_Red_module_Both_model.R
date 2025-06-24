#################################################################
### SECTION 1: LOAD LIBRARIES
#################################################################
message("Loading libraries...")
library(here)
library(tidyverse)
library(WGCNA)
library(SummarizedExperiment)
library(jaffelab)
library(gridExtra)
library(ragg) # For high-quality plot saving
library(grid) 
#################################################################
### SECTION 2: LOAD DATA
#################################################################
message("Loading data objects...")

# Load Bipolar Disorder differential expression results
load(here("results", "BP results.RDS"), verbose = TRUE) # Assumes this creates 'Outlist_BP'

# Load Lithium exposure results from two different models
load(here("results", "Lithium results DPseqmodel.RDS"), verbose = TRUE)
load(here("results", "Lithium results BPseqmodel.RDS"), verbose = TRUE)

# Load WGCNA network data
load(here("data", "zandiHypde_bipolar_rseGene_n511.rda"), verbose = TRUE)
load(here("data", "constructed_network_signed_bicor.rda"), verbose = TRUE)

#################################################################
### SECTION 3: PREPARE WGCNA DATA
#################################################################
message("Preparing WGCNA data...")

# Normalize and filter gene expression data
assays(rse_gene)$rpkm <- recount::getRPKM(rse_gene, 'Length')
geneIndex <- rowMeans(assays(rse_gene)$rpkm) > 0.25
rse_gene <- rse_gene[geneIndex, ]

# Prepare module color data
net$colorsLab <- labels2colors(net$colors)
colorDat <- data.frame(num = net$colors, col = net$colorsLab, stringsAsFactors = FALSE)
colorDat$Label <- paste0("ME", colorDat$num)
colorDat <- colorDat[order(colorDat$num), ]
colorDat <- colorDat[!duplicated(colorDat$num), ]
colorDat$numGenes <- table(net$colors)[as.character(colorDat$num)]


#################################################################
### SECTION 4: DEFINE THE MAIN PLOTTING FUNCTION
#################################################################
message("Defining the main analysis function...")

Compare_models <- function(module, model, outlist_name, type, i, j) {

  # Retrieve the main data object for this run based on its name
  Outlist <- get(outlist_name)

  message(paste("--> Starting analysis for:", model, type, "Region index:", j, "Module:", module))

  # Robustly create output directories before trying to save files
  dir.create(here("graphs", model, "wgcna"), showWarnings = FALSE, recursive = TRUE)
  dir.create(here("graphs", model, "t_stat"), showWarnings = FALSE, recursive = TRUE)
  dir.create(here("results/wgcna_Enrichment/", model), showWarnings = FALSE, recursive = TRUE)

  # This 'if' block handles the special case for the "Red module" plot
  if (module == "T") {
    message("Performing Red Module analysis")
    module_id <- which(net$colorsLab == "red")
    ensemblID_module <- rowData(rse_gene)$ensemblID[module_id]

    BP_results <- Outlist_BP[[i]][[j]] %>% rownames_to_column("Id")
    Li_results <- Outlist[[i]][[j]] %>% rownames_to_column("Id")

    colnames(BP_results) <- paste0(colnames(BP_results), "_BP")
    colnames(Li_results) <- paste0(colnames(Li_results), "_Li")

    message("... Joining BP and Lithium results")
    All_results <- BP_results %>%
      right_join(Li_results, by = c("Id_BP" = "Id_Li")) %>%
      dplyr::mutate(
        "Group" = if_else((adj.P.Val_BP < 0.05 & P.Value_Li < 0.005), "Both DEf",
          if_else((adj.P.Val_BP < 0.05 & P.Value_Li >= 0.005), "BPseq DEf",
            if_else((adj.P.Val_BP >= 0.05 & P.Value_Li < 0.005), "Lithium DEf", "None")
          )
        ),
        "Module" = if_else(ensemblID_BP %in% ensemblID_module, "Yes", "NO")
      )

    message("... Creating plots")
    # Your original plot definitions are kept
    p1 <- All_results %>% ggplot(aes(t_BP, t_Li, color = Module)) + geom_point(alpha = 0.7) +
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

    p2 <- All_results %>% ggplot(aes(t_BP, t_Li, color = Group)) + geom_point(alpha = 0.7)  +
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
    yleft <- textGrob("Lithium exposure", rot = 90, gp = gpar(fontsize = 20))

    # Arrange the plots WITHOUT drawing to prevent Rplots.pdf
    final_plot_grob <- grid.arrange(grobs = plot_lists, ncol = 2, bottom = bottom, left = yleft, draw = FALSE)

    # --- Save files using the robust ggsave method ---
    message("... Saving PNG file")
  #  ggsave(
   #  here("graphs", model, "wgcna", "Red_module_Compare_t_stat.pdf"),
   #  plot = final_plot_grob, width = 12, height = 6 , device = cairo_pdf
   # )
    ggsave(
      here("graphs", model, "wgcna", "Red_module_Compare_t_stat.png"),
      plot = final_plot_grob, device = ragg::agg_png, width = 12, height = 6, dpi = 300
    )

    message("... Saving RDS file")
    write_rds(All_results, file = here("results/wgcna_Enrichment/", model, "Red_module_Compare_t_stat.RDS"))

  } else {
    # This 'else' block handles all other plots
    message("... Performing general T-stat comparison analysis")
    BP_results <- Outlist_BP[[i]][[j]] %>% rownames_to_column("Id")

    Li_results <- Outlist[[i]][[j]] %>% rownames_to_column("Id")


    colnames(BP_results) <- paste0(colnames(BP_results), "_BP")
    colnames(Li_results) <- paste0(colnames(Li_results), "_Li")

    All_results <- BP_results %>%
            right_join(Li_results, by = c("Id_BP" = "Id_Li")) %>%
      dplyr::mutate("Group" = if_else((adj.P.Val_BP < 0.05 & P.Value_Li < 0.005), "Both DEf",
        if_else((adj.P.Val_BP < 0.05 & P.Value_Li >= 0.005), "BPseq DEf",
          if_else((adj.P.Val_BP >= 0.05 & P.Value_Li < 0.005), "Lithium DEf", "None")
        )
      ))

    p1 <- All_results %>% ggplot(aes(t_BP, t_Li, color = Group)) + geom_point(alpha = 0.7) +
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




    BP_results <- Outlist_BP[[i]][[j+1]] %>% rownames_to_column("Id")

    Li_results <- Outlist[[i]][[j+1]] %>% rownames_to_column("Id")


    colnames(BP_results) <- paste0(colnames(BP_results), "_BP")
    colnames(Li_results) <- paste0(colnames(Li_results), "_Li")

    All_results <- BP_results %>%
            right_join(Li_results, by = c("Id_BP" = "Id_Li")) %>%
      dplyr::mutate("Group" = if_else((adj.P.Val_BP < 0.05 & P.Value_Li < 0.005), "Both DEf",
        if_else((adj.P.Val_BP < 0.05 & P.Value_Li >= 0.005), "BPseq DEf",
          if_else((adj.P.Val_BP >= 0.05 & P.Value_Li < 0.005), "Lithium DEf", "None")
        )
      ))



    p2 <- All_results %>% ggplot(aes(t_BP, t_Li, color = Group)) + geom_point(alpha = 0.7) + 
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

    plot_lists <- list(p1, p2)
    bottom <- textGrob("BPseq", gp = gpar(fontsize = 20))
    yleft <- textGrob("Lithium exposure", rot = 90, gp = gpar(fontsize = 20))

    final_plot_grob <- grid.arrange(grobs = plot_lists, ncol = 2, bottom = bottom, left = yleft, draw = FALSE)

    # --- Save files using the robust ggsave method ---
    message("PNG file")
   # ggsave(
    #  here("graphs", model, "t_stat", paste0(type, "_Compare_t_stat.pdf")),
    # plot = final_plot_grob, width = 12, height = 6 , device = cairo_pdf
    #)
    ggsave(
      here("graphs", model, "t_stat", paste0(type, "_Compare_t_stat.png")),
      plot = final_plot_grob, device = ragg::agg_png, width = 12, height = 6, dpi = 300
    )

    message("... Saving RDS file")
    write_rds(All_results, file = here("results/wgcna_Enrichment/", model, paste0(type, "_Compare_t_stat.RDS")))
  }

  message(paste("--> Finished analysis for:", model, type, "Region index:", j))

  # Final Cleanup: Aggressively close any lingering graphics devices
  while (!is.null(dev.list())) {
    dev.off()
  }
}


#################################################################
### SECTION 5: EXECUTION
#################################################################

# Define the plan for all analyses using character strings for object names
plot_params <- tibble::tribble(
  ~module, ~model,        ~outlist_name,        ~type,        ~i, ~j,
  "T",     "BPseq_model", "Lithium_BPseqmodel", "Gene",       1,  1,
  "T",     "DPseq_model", "Lithium_DPseqmodel", "Gene",       1,  1,
  "F",     "BPseq_model", "Lithium_BPseqmodel", "Gene",       1,  1,
  "F",     "DPseq_model", "Lithium_DPseqmodel", "Gene",       1,  1,
  "F",     "BPseq_model", "Lithium_BPseqmodel", "Exon",       2,  1,
  "F",     "DPseq_model", "Lithium_DPseqmodel", "Exon",       2,  1,
  "F",     "BPseq_model", "Lithium_BPseqmodel", "Junction",   3,  1,
  "F",     "DPseq_model", "Lithium_DPseqmodel", "Junction",   3,  1,
  "F",     "BPseq_model", "Lithium_BPseqmodel", "Transcript", 4,  1,
  "F",     "DPseq_model", "Lithium_DPseqmodel", "Transcript", 4,  1
)

# Use a robust 'for' loop to iterate through the plan
message("Starting analysis loop...")
for (k in 1:nrow(plot_params)) {
  params <- plot_params[k, ]
  message(paste0("\n--- Starting Iteration ", k, " of ", nrow(plot_params), " ---"))
  tryCatch({
    Compare_models(
      module = params$module, model = params$model, outlist_name = params$outlist_name,
      type = params$type, i = params$i, j = params$j
    )
  }, error = function(e) {
    message(paste("!!!!!! ERROR in iteration", k, ":", e$message))
  })
  message("... Forcing garbage collection to release memory.")
  gc() # Manually trigger garbage collection
}


#################################################################
### SECTION 6: REPRODUCIBILITY
#################################################################
message("--- All analyses complete. ---")
Sys.time()
proc.time()
options(width = 120)
library(sessioninfo)
sessioninfo::session_info()