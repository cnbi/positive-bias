########################################################
# Simulation: naive lm vs. true multilevel model (lmer), for the            
# valence x intervention interaction effect.           

# Outcomes: bias and Type I error rate when effect_size = 0, and power
# (1 - Type II error) when effect_size != 0, and standard errors
########################################################

# libraries
library(lme4)     # for MLM
library(lmerTest) # for p-values

# Simulation design
n_images <- c(25, 40, 55, 70) # Number of images per experimental condition
n_subj <- c(30, 50, 70) # G*power suggestions?
effect_sizes <- c(0, 0.2, 0.5, 0.8) # remove 0.8?
# Variances but this needs further discussion
tau2 <- 0.5
sigma2 <- 1
alpha <- c(0.05, 0.01)

# Data generation

set.seed(123)

n_per_image <- 10     # observations per image
tau2   <- 0.5         # intercept variance
sigma2 <- 1           # residual variance
nsim   <- 100        # MC replications. #TODO: Set to 1000
alpha  <- 0.05        # alpha level for p-values
tau2_subj <- 0.7     # intercept variance per subject 

sample_size <- 10     # images per condition
effect_size <- 0      # true interaction effect

# valence: negative, neutral
# intervention: no, yes
# group: high, low anxiety

# simulate one dataset
simulate_data <- function(n_subj, n_images, effect_size) {
    
    cond <- expand.grid(valence = c(0, 1), intervention = c(0, 1))
    
    img  <- cond[rep(1:4, each = n_images), ]
    img$image <- 1:nrow(img)
    
    subj_id <- rep(1:n_subj, each = nrow(img))
    
    # Random effects
    b_subj <- rnorm(n_subj, mean = 0, sd = tau2_subj)
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
# new version with SEs
safe_row <- function(coef_table, term) {
    cols <- c("Estimate", "Std. Error", "Pr(>|t|)")
    if (term %in% rownames(coef_table)) {
        out <- as.numeric(coef_table[term, cols])
    } else {
        out <- rep(NA_real_, length(cols))
    }
    names(out) <- c("estimate", "se", "p_value")
    out
}

# fit all three models on one dataset, extract interaction term
get_estimates <- function(dat) {
    m_naive <- lm(Y ~ valence * intervention, data = dat)
    m_mlm   <- lmer(Y ~ valence * intervention + (1 | intervention:valence:image),
                    data = dat, REML = TRUE)
    
    term <- "valencepos:interventionyes"
    r_naive <- safe_row(summary(m_naive)$coefficients, term)
    r_mlm   <- safe_row(summary(m_mlm)$coefficients, term)
    
    data.frame(
        model    = c("naive", "multilevel"),
        estimate = c(r_naive["estimate"], r_mlm["estimate"]),
        se       = c(r_naive["se"],       r_mlm["se"]),
        p_value  = c(r_naive["p_value"],  r_mlm["p_value"]),
        row.names = NULL
    )
}

# run simulation
# new version with progress bar
run_simulation <- function(nsim, n_subj, n_images, effect_size, alpha = 0.05,
                           progress = TRUE) {
    
    res <- vector("list", nsim)
    if (progress) pb <- txtProgressBar(min = 0, max = nsim, style = 3)
    
    for (i in seq_len(nsim)) {
        dat <- simulate_data(n_subj, n_images, effect_size)
        res[[i]] <- cbind(sim = i, get_estimates(dat))
        if (progress) setTxtProgressBar(pb, i)
    }
    
    if (progress) close(pb)
    
    out <- do.call(rbind, res)
    out$reject <- out$p_value < alpha
    out
}

# summarise: bias, SEs and rejection rate (Type I error, or power) per model
summarise_sim <- function(sim_out, effect_size) {
    do.call(rbind, lapply(split(sim_out, sim_out$model), function(d) {
        data.frame(model = d$model[1],
                   bias = mean(d$estimate, na.rm = TRUE) - effect_size,
                   mean_se = mean(d$se, na.rm = TRUE),
                   emp_sd = sd(d$estimate, na.rm = TRUE),  # "true" SE for comparison
                   rejection_rate = mean(d$reject, na.rm = TRUE),
                   n_dropped = sum(is.na(d$estimate)),
                   row.names = NULL)
    }))
}

################################################################################
# test stuff #
################################################################################

res <- run_simulation(nsim=100, n_subj=100, n_images=25, effect_size=0, alpha = 0.05,
                           progress = TRUE)
summarise_sim(sim_out = res, effect_size = 0)

###
res2 <- run_simulation(nsim=nsim, n_subj=n_subj, n_images=n_images, effect_size=effect_size, alpha = 0.05,
                      progress = TRUE)
summarise_sim(sim_out = res2, effect_size = effect_size)
