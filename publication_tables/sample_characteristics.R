library(here)
library(dplyr)

### (output of preprocessing/preprocessing_bpseq_model.R) ###

load(
  here(
    "..","preprocessed_data", 
    "rse_gene_lithium_bpseq_model_processed_data.Rdata"
  )
)

### SAMPLE COUNTS ### 

# Crosstab of # samples lithium x brain region
table(rse_gene_lithium$lithium_group , rse_gene_lithium$BrainRegion)
# Amygdala sACC
# 0       25   27
# 1        9   10

# Crosstab of # samples sex x lithium in amygdala samples only
amyggroup = rse_gene_lithium[, rse_gene_lithium$BrainRegion == "Amygdala"]
table(amyggroup$Sex, amyggroup$lithium_group)
# 0  1
# F  9  4
# M 16  5


# Crosstab of # samples sex x lithium in sACC samples only
saccgroup = rse_gene_lithium[, rse_gene_lithium$BrainRegion == "sACC"]

table(saccgroup$Sex, saccgroup$lithium_group)
# 0  1
# F  9  5
# M 18  5

# Count of N samples by brain region only
table(rse_gene_lithium$BrainRegion)
# Amygdala     sACC 
#       34       37 

# Count of N samples by lithium exposure status only - not unique participants
table(rse_gene_lithium$lithium_group)
#  0  1 
# 52 19 


# Information on unique participants
unique_demogs <- unique(
  colData(rse_gene_lithium)[, c("BrNum", "Sex", "lithium_group")]
)
table(unique_demogs$lithium_group)
#  0  1 
# 29 10 

table(unique_demogs$lithium_group, unique_demogs$Sex)
  #    F  M
  # 0  9 20
  # 1  5  5

### RIN, pH, age of death, PMI ###
as.data.frame(colData(rse_gene_lithium)) %>% 
  dplyr::select(BrainRegion, lithium_group, RIN, pH, AgeDeath, PMI) %>%
  group_by(BrainRegion, lithium_group) %>%
  summarise(
    mean_RIN = mean(RIN, na.rm = FALSE), 
    sd_RIN = sd(RIN, na.rm=FALSE),
    mean_pH = mean(pH),
    sd_pH = sd(pH),
    mean_agedeath = mean(AgeDeath),
    sd_agedeath = sd(AgeDeath),
    mean_pmi = mean(PMI),
    sd_pmi = sd(PMI)
  )


# # A tibble: 4 × 10
# # Groups:   BrainRegion [2]
#   BrainRegion lithium_group mean_RIN sd_RIN mean_pH sd_pH mean_agedeath sd_agedeath mean_pmi sd_pmi
#   <fct>       <fct>            <dbl>  <dbl>   <dbl> <dbl>         <dbl>       <dbl>    <dbl>  <dbl>
# 1 Amygdala    0                 7.55  1.03     4.29  3.01          39.8        12.3     26.3   6.55
# 2 Amygdala    1                 7.17  0.660    4.91  2.79          32.5        12.4     26.3   4.42
# 3 sACC        0                 8.02  0.889    4.21  3.04          40.6        12.3     27.9   7.60
# 4 sACC        1                 7.54  0.717    5.04  2.66          33.7        12.3     25.5   4.93
