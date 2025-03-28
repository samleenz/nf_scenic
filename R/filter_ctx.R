# Fri Mar 28 17:16:27 2025 ------------------------------
library(tidyverse)


# read in CTX results file ------------------------------------------------
pth <- "/vast/scratch/users/lee.sa/nextflow/work/03/434d269903692f9f7c63294ca03e60/ctx_output.tsv"


ctx <- read_tsv(pth, skip = 3, col_names = F)

colnames(ctx) <- c("TF",
"motifID",
"AUC",
"NES",
"Context",
"Annotation",
"MotifSimilarityQvalue",
"OrthologousIdentity",
"RankAtMax",
"TargetGenes"
)


# extract the table of targets per TF ------------------------------------

redun <- function(x){
  purrr::reduce(x, union)
}

df <- ctx %>%
  mutate(gene_symbols = lapply(RankAtMax, function(x) {
    matches <- str_match_all(x, "'([^']+)',")[[1]]
    if (nrow(matches) == 0) return(NA_character_)
    matches[, 2]
  })) |>
  summarise(genes = redun(gene_symbols), .by = TF) 

df


# save long table ---------------------------------------------------------
 



###
### collate steps over all files
###


# load all to single long table
##  need to add an indicator for which run
# identify high confidence genes per regulon
# drop regulons (per run) where they do not have at least 5 high confidence genes
# identify high confidence regulons
# save HCregulons in format for AUCell

# run AUcell woith HCregulons -> this is what we analyse.


