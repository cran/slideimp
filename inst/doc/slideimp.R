## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(slideimp)

## -----------------------------------------------------------------------------
data(khanmiss1)
obj <- t(khanmiss1)
set.seed(1234)

## -----------------------------------------------------------------------------
n_genes <- 30
n_samples <- 5
n_reps <- 10
na_location <- replicate(
  n = n_reps,
  {
    chosen_genes <- sample.int(ncol(obj), size = n_genes)
    values <- lapply(
      chosen_genes, \(x) {
        chosen_samples <- sample.int(nrow(obj), size = n_samples)
        matrix(c(chosen_samples, rep(x, n_samples)), ncol = 2, dimnames = list(NULL, c("row", "col")))
      }
    )
    do.call(rbind, values)
  },
  simplify = FALSE
)
na_location[[1]][1:10, ]

## -----------------------------------------------------------------------------
tune_knn <- tune_imp(obj, parameters = data.frame(k = 10), rep = na_location, cores = 4)
tune_pca <- tune_imp(obj, parameters = data.frame(ncp = 5), rep = na_location)

## -----------------------------------------------------------------------------
rmse_knn <- mean(compute_metrics(tune_knn, metrics = "rmse")$.estimate)
rmse_pca <- mean(compute_metrics(tune_pca, metrics = "rmse")$.estimate)
c("knn" = rmse_knn, "pca" = rmse_pca)

## -----------------------------------------------------------------------------
# 500 features, 50 samples, 2 chromosomes
sim_obj <- sim_mat(n = 500, m = 50, nchr = 2)
# sim_mat returns matrix with features in rows, so we use t() to put features in columns
obj <- t(sim_obj$input)
# `group_feature` contains metadata for the simulated data
head(sim_obj$group_feature)
n_needed <- 100
# Randomly select CpGs that are needed for clock calculation
subset_cpgs <- sample(colnames(obj), size = n_needed)

## -----------------------------------------------------------------------------
# Construct the group data.frame with one row per group (i.e., chromosome)
group_df <- group_features(obj, sim_obj$group_feature, subset = subset_cpgs)
group_df

## -----------------------------------------------------------------------------
pca_group <- group_df
# For PCA imputation, we use ncp = 5 for the first group and ncp = 10 for the second group
pca_group$parameters <- list(data.frame(ncp = 5), data.frame(ncp = 10))
# For K-NN imputation, we use k = 5 for both groups
knn_group <- group_df
knn_group$parameters <- list(data.frame(k = 5), data.frame(k = 5))

pca_group
knn_group

## -----------------------------------------------------------------------------
system.time(
  knn_results <- group_imp(obj, group = knn_group, cores = 4)
)
system.time(
  pca_results <- group_imp(obj, group = pca_group)
)

## -----------------------------------------------------------------------------
sample_names <- paste0("S", 1:10)
positions <- 500

methyl <- tibble::tibble(
  chr = "chr1",
  start = seq_len(positions),
  end = start,
  strand = "+"
)

set.seed(1234)
for (i in seq_along(sample_names)) {
  methyl[[paste0("numCs", i)]] <- sample.int(100, size = positions, replace = TRUE)
  methyl[[paste0("numTs", i)]] <- sample.int(100, size = positions, replace = TRUE)
  methyl[[paste0("coverage", i)]] <- methyl[[paste0("numCs", i)]] + methyl[[paste0("numTs", i)]]
}

methyl[1:5, 1:10]

## -----------------------------------------------------------------------------
methyl <- methyl[order(methyl$start), ] # <---- Important, sort!
numCs_matrix <- as.matrix(methyl[, paste0("numCs", seq_along(sample_names))])
cov_matrix <- as.matrix(methyl[, paste0("coverage", seq_along(sample_names))])
beta_matrix <- numCs_matrix / cov_matrix

colnames(beta_matrix) <- sample_names
rownames(beta_matrix) <- methyl$start

beta_matrix <- t(beta_matrix)
# Set 10% of the data to missing
set.seed(1234)
beta_matrix[sample.int(length(beta_matrix), floor(length(beta_matrix) * 0.1))] <- NA
beta_matrix[1:5, 1:5]

## -----------------------------------------------------------------------------
# For example, we are tuning `k` value here.
slide_knn_params <- tibble::tibble(k = c(10, 20), n_feat = 100, n_overlap = 10)
# Increase rep from 2 in actual analyses
tune_slide_knn <- tune_imp(
  obj = beta_matrix,
  parameters = slide_knn_params,
  cores = 4,
  rep = 2
)
compute_metrics(tune_slide_knn)

## -----------------------------------------------------------------------------
imputed_beta <- slide_imp(obj = beta_matrix, n_feat = 100, n_overlap = 10, k = 10, cores = 4)
imputed_beta

