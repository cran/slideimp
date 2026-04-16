## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(slideimp)
set.seed(1234)

## -----------------------------------------------------------------------------
# 20 rows, 1000 columns, all columns have at least some NA
sim_obj <- sim_mat(n = 20, p = 1000, perc_col_na = 1)
obj <- sim_obj$input
obj[1:4, 1:4]

## -----------------------------------------------------------------------------
na_loc <- sample_na_loc(obj, n_cols = 200, n_rows = 5, n_reps = 5)
length(na_loc)
na_loc[[1]][1:6, ]

## -----------------------------------------------------------------------------
# This custom function imputes missing values with random normal values and takes
# `mean` and `sd` as params
rnorm_imp <- function(obj, mean, sd) {
  na <- is.na(obj)
  obj[na] <- rnorm(sum(na), mean = mean, sd = sd) # <- impute values with rnorm
  return(obj) # <- return an imputed object with the same dim as obj
}

pca_tune <- tune_imp(
  obj,
  .f = "pca_imp",
  na_loc = na_loc,
  parameters = data.frame(ncp = 10)
)

knn_tune <- tune_imp(
  obj,
  .f = "knn_imp",
  na_loc = na_loc,
  parameters = data.frame(k = 10)
)

rnorm_tune <- tune_imp(
  obj,
  .f = rnorm_imp,
  na_loc = na_loc,
  parameters = data.frame(mean = 0, sd = 1) # must match with arguments of `rnorm_imp`
)

## -----------------------------------------------------------------------------
mean(compute_metrics(pca_tune, metrics = "rmse")$.estimate)
mean(compute_metrics(knn_tune, metrics = "rmse")$.estimate)
mean(compute_metrics(rnorm_tune, metrics = "rmse")$.estimate)

## -----------------------------------------------------------------------------
sim_obj <- sim_mat(n = 20, p = 50, n_col_groups = 2)

# Matrix to be imputed
obj <- sim_obj$input
obj[1:5, 1:4]

# Metadata, i.e., which features belong to which group
meta <- sim_obj$col_group
meta[1:5, ]

# We put feature 1 in `group3`
meta[1, 2] <- "group3"
meta[1:5, ]

## -----------------------------------------------------------------------------
set.seed(1234)
group_imp_df <- prep_groups(colnames(obj), group = meta, min_group_size = 10)
group_imp_df$parameters <- list(list(k = 3), list(k = 4), list(k = 5))
group_imp_df

## -----------------------------------------------------------------------------
knn_results <- group_imp(obj, group = group_imp_df, cores = 4, k = 10)
print(knn_results, p = 4)

## -----------------------------------------------------------------------------
set.seed(1234)
sample_names <- paste0("S", 1:10)
n_sites <- 1000

# Simulate positions with 50–500 bp between each site
distances_between <- sample(50:500, size = n_sites, replace = TRUE)
locations <- cumsum(distances_between) # <- important, location vector

methyl <- data.frame(
  chr = "chr1",
  start = locations,
  end = locations,
  strand = "+"
)

for (i in seq_along(sample_names)) {
  methyl[[paste0("numCs", i)]] <- sample.int(100, size = n_sites, replace = TRUE)
  methyl[[paste0("numTs", i)]] <- sample.int(100, size = n_sites, replace = TRUE)
  methyl[[paste0("coverage", i)]] <- methyl[[paste0("numCs", i)]] + methyl[[paste0("numTs", i)]]
}

methyl[1:5, 1:10]

## -----------------------------------------------------------------------------
numCs_matrix <- as.matrix(methyl[, paste0("numCs", seq_along(sample_names))])
cov_matrix <- as.matrix(methyl[, paste0("coverage", seq_along(sample_names))])
beta_matrix <- numCs_matrix / cov_matrix

colnames(beta_matrix) <- sample_names
rownames(beta_matrix) <- methyl$start

beta_matrix <- t(beta_matrix)
# Set 10% of the data to missing
set.seed(1234)
beta_matrix[sample.int(length(beta_matrix), floor(length(beta_matrix) * 0.1))] <- NA
beta_matrix[1:4, 1:4]

## -----------------------------------------------------------------------------
params <- expand.grid(ncp = c(2, 4), window_size = c(5000, 10000))
params$overlap_size <- 1000
params$min_window_n <- 20 # windows with less than 20 columns are dropped

# Increase n_reps from 2 in actual analyses and use another chromosome (i.e., chr22)
tune_slide_pca <- tune_imp(
  obj = beta_matrix,
  parameters = params,
  .f = "slide_imp",
  n_reps = 2,
  location = locations
)

metrics <- compute_metrics(tune_slide_pca)

aggregate(.estimate ~ .metric + ncp + window_size, data = metrics, FUN = mean)

## -----------------------------------------------------------------------------
slide_imp(
  obj = beta_matrix,
  location = locations,
  window_size = 5000,
  overlap_size = 1000,
  ncp = 2,
  min_window_n = 20,
  dry_run = TRUE # <- dry_run to inspect the windows
)

## -----------------------------------------------------------------------------
slide_imp(
  obj = beta_matrix,
  location = locations,
  window_size = 5000,
  overlap_size = 1000,
  ncp = 2,
  min_window_n = 20,
  dry_run = FALSE,
  .progress = FALSE
)

## -----------------------------------------------------------------------------
slide_imp(
  obj = beta_matrix,
  location = locations,
  window_size = 5000,
  ncp = 2,
  min_window_n = 20,
  subset = c("1323", "33810"),
  flank = TRUE,
  dry_run = TRUE
)

