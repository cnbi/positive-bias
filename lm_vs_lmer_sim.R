##########################################################################
# Simulation: naive lm vs. true multilevel model (lmer), for the         #
# valence x intervention interaction effect.                             #
#                                                                        #
# Outcomes: bias and Type I error rate when effect_size = 0, and power   #
# (1 - Type II error) when effect_size != 0, and standard errors         #
##########################################################################

# Libraries
library(lme4)     # for MLM
library(lmerTest) # for p-values
library(tidyverse)
library(nanoparquet) # save design matrix as Parquet file

#------------------------- Simulation design -----------------------------------
n_images <- c(25, 40, 55, 70) # Number of images per experimental condition
n_subj <- c(30, 50, 70) # G*power suggestions?
eff_sizes <- c(0, 0.05, 0.2, 0.5)
icc_subj <- c(0.01, 0.1, 0.2)
icc_img <- c(0.01, 0.1, 0.2)
design_matrix <- expand.grid(
  n_images = n_images,
  n_subjects = n_subj,
  effect_size = eff_sizes,
  icc_subj = icc_subj,
  icc_img = icc_img
)
set.seed(511)
design_matrix <- mutate(design_matrix,
                        row_id = 1:n(),
                        seed = as.integer(sample(2^32 / 2, n())))
write_parquet(design_matrix, "data/simulation/sim_design_matrix.parquet")

# Variances but this needs further discussion
# tau2 <- 0.5
# sigma2 <- 1
# alpha <- 0.05

# Required functions for simulation
source("functions.R")

# Data generation
# n_per_image <- 10     # observations per image
# tau2   <- 0.5         # intercept variance
# sigma2 <- 1           # residual variance
# nsim   <- 100        # MC replications. #TODO: Set to 1000
# alpha  <- 0.05        # alpha level for p-values
# tau2_subj <- 0.7     # intercept variance per subject
#
# sample_size <- 10     # images per condition
# effect_size <- 0      # true interaction effect

# valence: negative, neutral
# intervention: no, yes
# group: high, low anxiety

#------------------------ Run simulation ---------------------------------------
run_simulation <- function(nsim,
                           n_subj,
                           n_images,
                           effect_size,
                           alpha = 0.05,
                           progress = TRUE) {
  res <- vector("list", nsim)
  if (progress)
    pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  
  for (i in seq_len(nsim)) {
    dat <- simulate_data(n_subj, n_images, effect_size)
    res[[i]] <- cbind(sim = i, get_estimates(dat))
    if (progress)
      setTxtProgressBar(pb, i)
  }
  
  if (progress)
    close(pb)
  
  out <- do.call(rbind, res)
  out$reject <- out$p_value < alpha
  out
}


################################################################################
# test stuff #
################################################################################

res <- run_simulation(
  nsim = 100,
  n_subj = 100,
  n_images = 25,
  effect_size = 0,
  alpha = 0.05,
  progress = TRUE
)
summarise_sim(sim_out = res, effect_size = 0)

###
res2 <- run_simulation(
  nsim = nsim,
  n_subj = n_subj,
  n_images = n_images,
  effect_size = effect_size,
  alpha = 0.05,
  progress = TRUE
)
summarise_sim(sim_out = res2, effect_size = effect_size)
