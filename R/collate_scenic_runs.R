# Fri Mar 21 08:53:44 2025 ------------------------------
setwd("/vast/scratch/users/lee.sa/nf_scenic")
library(tidyverse)
library(SCENIC)
library(SCopeLoomR)
library(AUCell)

# regulon extraction function
get_regulon <- function(regulon, rIM = regulons_incidMat){
  rIM[regulon, rIM[regulon, ] == 1 ] |>
    names()
}

# Read in the TF list (all possible regulons)
all_regs <- read_lines("databases/hs_hgnc_tfs.txt") |>
  paste0("(+)")


proj <- "results/exp6_fullRun"

iters <- list.files(proj)

scenic <- map(iters, \(x){
  fname <- file.path(proj, x, "auc_output.loom")
  
  loom <- open_loom(fname)
  regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
  regulonAUC <- get_regulons_AUC(loom, column.attr.name='RegulonsAUC')
  
  # add all possible regulons to incidence matrix
  rIM <- matrix(0, nrow = length(all_regs), ncol = ncol(regulons_incidMat), dimnames = list(all_regs, colnames(regulons_incidMat)))
  rIM[rownames(regulons_incidMat), ] <- regulons_incidMat
  
  # add all possible regulons to AUC matrix
  aucMat <- matrix(0, nrow = length(all_regs), ncol = ncol(regulonAUC), dimnames = list(all_regs, colnames(regulonAUC)))
  aucMat[rownames(regulonAUC), ] <- SummarizedExperiment::assay(regulonAUC, "AUC")

  
  
  # return list object with info
  list(
    rIM = rIM,
    aucMat = aucMat
  )
},
.progress = TRUE
)

## cumulate the incidence matrices
add_them <- function(x, y){ x + y}

rIM_cumulate <- scenic |>
  map(\(x) x$rIM) |>
  purrr::reduce(add_them)

hist(rIM_cumulate[rIM_cumulate>0])

# number of genes per regulon 
genes_per_regulon <- left_join(
  rIM_cumulate |>
    as_tibble(rownames = "regulon") |>
    pivot_longer(-regulon, names_to = "gene", values_to = "incidence") |>
    filter(incidence > 0) |>
    summarise(unfiltered = n(), .by = regulon),
  rIM_cumulate |>
    as_tibble(rownames = "regulon") |>
    pivot_longer(-regulon, names_to = "gene", values_to = "incidence") |>
    filter(incidence > 30) |>
    summarise(filtered = n(), .by = regulon),
  by = "regulon"
) |>
  mutate(filtered = ifelse(is.na(filtered), 0, filtered))

ggplot(genes_per_regulon, aes(unfiltered, filtered)) +
  geom_point() + 
  theme_bw()


cumulate_multi <- purrr::map(c(0, 10, 20, 30, 40, 45), \(x){
  rIM_cumulate |>
    as_tibble(rownames = "regulon") |>
    pivot_longer(-regulon, names_to = "gene", values_to = "incidence") |>
    filter(incidence > x) |>
    summarise(ngenes = n(), .by = regulon) |>
    mutate(cutoff = x)
}, .progress = TRUE)
  
cumulate_multi |>
  list_rbind() |>
  ggplot(aes(factor(cutoff), ngenes)) +
  geom_boxplot() +
  # scale_x_discrete() +
  theme_bw()

cumulate_multi |>
  list_rbind() |>
  ggplot(aes(cutoff, ngenes, group = regulon)) +
  geom_line() +
  # scale_x_discrete() +
  theme_bw()

rIM_cumulate |>
  as_tibble(rownames = "regulon") |>
  pivot_longer(-regulon, names_to = "gene", values_to = "incidence") |>
  filter(incidence > 0) |>
  ggplot(aes(incidence, group = regulon)) +
  geom_bar() +
  scale_x_binned() +
  theme_bw()

## select non-empty regulons
active_regulons <- genes_per_regulon |>
  filter(filtered > 0)

## get the mean AUC across runs
get_mean <- function(X, len){
  X / len
}

aucMat_mean <-  scenic |>
  purrr::map(\(x) x$aucMat) |>
  purrr::reduce(add_them) |>
  get_mean(len = length(scenic))

map(scenic, \(x){
  dim(x$rIM)
})


# detect which regulons have bimodal distribution
library(mclust)
library(Matrix)
# To binarise the AUC scores we must identify an activation threshold. In the
# pySCENIC library this is done via fitting a two-component gaussian mixture
# model so we will follow this approach. They define the threshold as tr =
# mu(comp_2) - 2SD(comp_2) where comp_2 is the component with the greater mean
# (e.g., the 'active component').
#
# We fit the k=2 GMM with mixtools, first trying to fit the GMM on all data
# points for the regulon. If this does not converge, we then run the GMM on only
# the non-zero regulon scores. This prevents the variance of the smaller
# component shrinking to 0 as a result of high sparsity (e.g., a regulon with
# 40% 0s will fill to converge).
#
# The threshold is then defined as the max between the smallest AUC of the
# regulon and two SDs below the larger components mean. The AUCs can then be
# binarised.

binarise_vec <- function(x, verbose = FALSE) {
  ## if all values are the same, all are set to inactive (0)
  if(length(unique(x)) == 1){
    if(verbose){
      message("All values of x the same, all set to inactive")
    }
    return(rep(0, length(x)))
  }
  
  fit_gmm <- function(x){
    mixtools::normalmixEM(x = x, k = 2, epsilon = 0.01)
  }
  
  gmm_2k <- tryCatch(
    {
      fit_gmm(x)
    },
    error = function(e) {
      # Retry with smallest vals (0s) removed only if the original call fails
      if(verbose){
        message("First attempt failed, retrying with x > 0")
      }
      tryCatch(
        {
          fit_gmm(x[x > min(x)])
        },
        error = function(e) {
          if(verbose){
            message("Second attempt failed, returning NA")
          }
          NA
        }
      )
    }
  )
  
  if(! inherits(gmm_2k, "mixEM")){
    if(verbose){
      message("FItting GMM failed twice, all set to inactive")
    }
    return(rep(0, length(x)))
  }
  
  if(verbose){
    print(gmm_2k)
    
    summary(gmm_2k)
  }
  
  threshold <- max(
    max(gmm_2k$mu) - 2*(gmm_2k$sigma[which.max(gmm_2k$mu)]),
    min(x)
  )
  
  # return 1 if x > t, otherwise 0
  as.numeric(
    x > threshold
  )
}

tst <- binarise_vec(scenic[[2]]$aucMat["FOXO4(+)", ], verbose = TRUE)
result <- binarise_vec(scenic[[1]]$aucMat["HIC1(+)", ], verbose = TRUE)

# Run binarise_vec() while capturing and suppressing any printed output (e.g. EM warnings)
# the output is in `result`
# blahblah <- capture.output({
#   result <- binarise_vec(scenic[[2]]$aucMat["FOXO4(+)", ])
# }, type = "output")
# 
# table(result)

# 
# blahblah <- capture.output({
#   bin_mat <- purrr::map(rownames(scenic[[1]]$aucMat), \(regulon){
#     # print(regulon)
#     binarise_vec(scenic[[1]]$aucMat[regulon, ])
#   })
# }, type = "output")
# 
# bin_mat <- purrr::map(
#   genes_per_regulon$regulon, 
#   \(reg){
#     blahblah <- capture.output({
#       ret <- binarise_vec(scenic[[1]]$aucMat[reg, ], verbose = FALSE)
#     },
#     type = "output"
#     )
#     ret
#     },
#   .progress = TRUE
#   )
# 
# 
# mat <- do.call(rbind, bin_mat) |>
#   Matrix::Matrix(sparse = TRUE)

# binarise all the AUC matrices from the scenic runs and then calculate the incidence matrix
library(future)
library(furrr)

plan(multisession, workers = 8)

auc_incidence <- furrr::future_map(scenic, \(sc){
  aucM <- sc$aucMat
  
  bin_mat <- purrr::map(
    genes_per_regulon$regulon, 
    \(reg){
      blahblah <- capture.output({
        ret <- binarise_vec(aucM[reg, ], verbose = FALSE)
      },
      type = "output"
      )
      ret
    },
    .progress = FALSE
  )
  
  mat <- do.call(rbind, bin_mat) |>
    Matrix::Matrix(sparse = TRUE)
  
  mat
},
.progress = TRUE
) |>
  purrr::reduce(add_them); write_rds(
  auc_incidence,
  here::here("results/R/exp6_fullRun_auc_binary_incidence_mat.rds")
  )

# range(auc_incidence)
auc_incidence[auc_incidence > 0] |>
  as.numeric() |>
  hist(breaks = 50)

auc_cumprop <- table(as.numeric(auc_incidence)) |>
  enframe("x", "count") |>
  mutate(x = as.numeric(x)) |>
  arrange(desc(x)) |>
  mutate(cum = cumsum(count)) |>
  mutate(cumprop = cum / prod(dim(auc_incidence)))

auc_cumprop|>
  ggplot(aes(
    x, cumprop
  )) +
  geom_point() +
  theme_bw() +
  ylim(0,1) +
  scale_x_reverse() +
  labs(
    x = "Runs detected in"
  )

auc_cumprop |>
  filter(cumprop > 0.25) |>
  slice_max(x, n = 1)

# create a mask for active cell X regulons
# the TRUE elements are those to be removed
thresh <- 30
auc_mask <- auc_incidence < thresh


# sum the auc values across runs
aucMat_sum <-  scenic |>
  purrr::map(\(x) x$aucMat) |>
  purrr::reduce(add_them)

# get the non-zero auc values across runs
aucMat_nzero <-  scenic |>
  purrr::map(\(x) {
    a <- matrix(0, nrow = nrow(x$aucMat), ncol = ncol(x$aucMat), dimnames = dimnames(x$aucMat))
    a[x$aucMat > 0] <- 1
    
    a
  }) |>
  purrr::reduce(add_them)

aucMat_mean2 <- aucMat_sum / aucMat_nzero
aucMat_mean2[is.nan(aucMat_mean2)] <- 0

# remove the low confidence scores
auc_mat <- aucMat_mean2
auc_mat[as.matrix(auc_mask)] <- 0

# remove undetected regulons
auc_mat <- auc_mat[rowSums(auc_mat) > 0, ]

# remove regulons that have less than 5 genes consistently detected in at least
# the threshold number of iterations

regulon_lengths <- scenic |>
  purrr::map(\(x) x$rIM) |>
  purrr::map(\(x){
    rowSums(x) |>
      enframe("regulon", "ngenes") |>
       filter(ngenes >= 5)
  }) |>
  list_rbind()

regulon_lengths |>
  summarise(n = n(), .by = regulon) |>
  arrange(desc(n)) |>
  View()

genes_per_regulon |>
  mutate(selected = regulon %in% rownames(auc_mat))  |>
  ggplot(aes(unfiltered, filtered, colour = selected)) +
  geom_point() +
  scico::scale_colour_scico_d() +
  theme_bw()
