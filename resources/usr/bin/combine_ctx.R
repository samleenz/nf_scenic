#! #!/usr/bin/env Rscript

# Sam Lee
# @samleenz

library(tidyverse)

### custom functions -------------------------------------------------------

redun <- function(x){
  # iteratively apply union to a 
  # list of character vectors
  purrr::reduce(x, union)
}

my_split <- function(x){
  # split a data frame by a column
  split(x$genes, x$TF)
}

write_gmt <- function(gene_sets, file) {
  # write a gmt file (each row is a gene set)
  # format is: gene_set_name \t description (NA) \t gene1 \t gene2 \t ...
  lines <- purrr::imap_chr(gene_sets, ~ paste(c(.y, "NA", .x), collapse = "\t"))
  readr::write_lines(lines, file)
}

# read in the ctx run files -----------------------------------------------

pths <- commandArgs(trailingOnly=TRUE)


cts <- purrr::map(pths, \(pth){
  x <- readr::read_tsv(pth, skip = 3, col_names = F, show_col_types = F)
  colnames(x) <- c("TF",
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
  
  x
})

# extract the table of targets per TF ------------------------------------

dfs <- purrr::map(cts, \(ctx){
  ctx |>
    # this mutate transforms the dict. style column `RankAtMax` to a 
    # list-col of clean gene symbols
    mutate(gene_symbols = lapply(RankAtMax, function(x) {
      matches <- str_match_all(x, "'([^']+)',")[[1]]
      if (nrow(matches) == 0){
        NA_character_
      } else{
        matches[, 2]
      }
    })) |>
    # reframe is similar to summarise, but can return
    # an arbitrary number of rows.
    reframe(genes = redun(gene_symbols), .by = TF)
}) |>
  list_rbind(names_to = "run")


# identify HC genes -------------------------------------------------------

## an HC gene is found in >= 80% of runs for a regulon
min_n <- floor(0.8 * max(dfs$run))


# list of HC genes per regulon
HC_genes <- dfs |>
  summarise(n = n(), .by = c(TF, genes)) |>
  filter(n >= min_n) |>
  my_split()

HC_genes_per_regulon <- do.call("rbind", map(HC_genes, length)) |>
  as.data.frame() |>
  rownames_to_column("TF") |>
  rename(n_hc_genes = V1) 
  


# Identify HC regulons ----------------------------------------------------

# a. regulons found in >= 80% of runs
# b. at least 5 HC genes in the regulon
 
regulon_incidence <- dfs |> 
  distinct(run, TF) |> 
  summarise(n_runs = n(), .by = TF)

HC_regulon_df <- regulon_incidence |>
  left_join(HC_genes_per_regulon, by = "TF") |>
  filter(n_runs >= min_n & n_hc_genes >= 5)

HC_regulons <- HC_genes[HC_regulon_df$TF]

# save HC regulons as a GMT file ------------------------------------------
## see here for GMT file format 
## https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats

write_gmt(HC_regulons, "gene_sets.gmt")

