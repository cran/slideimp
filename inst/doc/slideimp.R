## -----------------------------------------------------------------------------
#| label: setup
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

# Put feature 1 in `group3`
meta[1, 2] <- "group3"
meta[1:5, ]


## -----------------------------------------------------------------------------
set.seed(1234)
group_imp_df <- prep_groups(colnames(obj), group = meta, min_group_size = 10)
group_imp_df$parameters <- list(list(k = 3), list(k = 4), list(k = 5))
group_imp_df


## -----------------------------------------------------------------------------
knn_results <- group_imp(obj, group = group_imp_df, cores = 2, k = 10)
print(knn_results, p = 4)


## -----------------------------------------------------------------------------
set.seed(1234)
chr1_beta <- sim_mat(n = 30, p = 500)$input
location <- cumsum(sample(50:500, ncol(chr1_beta), replace = TRUE))
chr1_beta[1:4, 1:4]


## -----------------------------------------------------------------------------
targets <- colnames(chr1_beta)[c(100, 250, 400)]
plan <- mipca_flank_plan(
  chr1_beta,
  location,
  window_size = 2000, # bp on each side of the target
  subset = targets,
  min_window_n = 10
)
plan


## -----------------------------------------------------------------------------
win <- mipca_flank_window(plan, targets[1])
tuned <- mipca_tune(win, rank = 1:3, transform = "logit_biwhiten", n_cols = 5, seed = 1)
tuned


## -----------------------------------------------------------------------------
fit <- mipca_flank(
  plan,
  rank = tuned$rank[which.min(tuned$crps)],
  seed = 1,
  transform = "logit_biwhiten",
  .progress = FALSE
)
fit
dim(as.matrix(fit, draw = 1))

