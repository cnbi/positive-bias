########################################################
# Simulation: naive lm vs. fixed-effects lm            
# vs. true multilevel model (lmer), for the            
# valence x intervention interaction effect.           

# Outcomes: bias and Type I error rate when effect_size = 0, and power
# (1 - Type II error) when effect_size != 0.
########################################################

library(lme4)     # for MLM
library(lmerTest) # for p-values

set.seed(123)

n_per_image <- 10     # observations per image
tau2   <- 0.5         # intercept variance
sigma2 <- 1           # residual variance
nsim   <- 1000        # MC replications
alpha  <- 0.05        # alpha level for p-values

sample_size <- 10     # images per condition
effect_size <- 0      # true interaction effect

# simulate one dataset
simulate_data <- function(sample_size, effect_size) {
  cond <- expand.grid(valence = c(0, 1), intervention = c(0, 1))
  img  <- cond[rep(1:4, each = sample_size), ]
  img$image <- 1:nrow(img)
  img$int_interv_val <- interaction(img$valence, img$intervention, drop = TRUE)
  img$b_image <- rnorm(nrow(img), 0, sqrt(tau2))   # random intercept per image N(0,tau2)
  
  dat <- img[rep(1:nrow(img), each = n_per_image), ]
  dat$Y <- effect_size * dat$valence * dat$intervention +
    dat$b_image + rnorm(nrow(dat), 0, sqrt(sigma2))
  
  dat$valence      <- factor(dat$valence, labels = c("neg", "pos"))
  dat$intervention <- factor(dat$intervention, labels = c("no", "yes"))
  dat$image        <- factor(dat$image)
  dat
}

# pull estimates and p-values
# when a term is aliased due to exact collinearity, drop that row from the table 
# for now, this happens in every iteration...
safe_row <- function(coef_table, term) {
  if (term %in% rownames(coef_table)) {
    coef_table[term, c("Estimate", "Pr(>|t|)")]
  } else {
    c(Estimate = NA_real_, `Pr(>|t|)` = NA_real_)
  }
}

# fit all three models on one dataset, extract interaction term
get_estimates <- function(dat) {
  m_naive <- lm(Y ~ valence * intervention, data = dat)
  m_fixed <- lm(Y ~ valence * intervention + image, data = dat)
  m_mlm   <- lmer(Y ~ valence * intervention + (1 | int_interv_val:image),
                  data = dat, REML = TRUE)
  
  term <- "valencepos:interventionyes"
  r_naive <- safe_row(summary(m_naive)$coefficients, term)
  r_fixed <- safe_row(summary(m_fixed)$coefficients, term)
  r_mlm   <- safe_row(summary(m_mlm)$coefficients, term)
  
  data.frame(
    model    = c("naive", "fixed_effects", "multilevel"),
    estimate = c(r_naive["Estimate"], r_fixed["Estimate"], r_mlm["Estimate"]),
    p_value  = c(r_naive["Pr(>|t|)"], r_fixed["Pr(>|t|)"], r_mlm["Pr(>|t|)"])
  )
}

# run simulation
sim_out <- do.call(rbind, lapply(1:nsim, function(i) {
  get_estimates(simulate_data(sample_size, effect_size))
}))
sim_out$reject <- sim_out$p_value < alpha

# summarise: bias and rejection rate (Type I error, or power) per model
summary_tbl <- do.call(rbind, lapply(split(sim_out, sim_out$model), function(d) {
  data.frame(model = d$model[1],
             bias = mean(d$estimate, na.rm = TRUE) - effect_size,
             rejection_rate = mean(d$reject, na.rm = TRUE),
             n_dropped = sum(is.na(d$estimate)))
}))

print(summary_tbl)
