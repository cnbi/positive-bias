########################### FUNCTIONS #####################################

# simulate one dataset
simulate_data <- function(n_subj, n_images, effect_size, seed) {
    set.seed(seed)
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
    dat$id           <- subj_id
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

# Fit all three models on one dataset, extract interaction term
get_estimates <- function(dat) {
    m_naive <- lm(Y ~ valence * intervention, data = dat)
    m_mlm   <- lmer(
        Y ~ valence * intervention + (1 + valence |
                                          id) + (1 | image:intervention),
        data = dat,
        REML = TRUE
    )
    # previous model formultaion:
    # m_mlm   <- lmer(Y ~ valence * intervention + (1 | intervention:valence:image),
    #                 data = dat, REML = TRUE)
    term <- "valencepos:interventionyes"
    r_naive <- safe_row(summary(m_naive)$coefficients, term)
    r_mlm   <- safe_row(summary(m_mlm)$coefficients, term)
    
    data.frame(
        model    = c("naive", "multilevel"),
        estimate = c(r_naive["estimate"], r_mlm["estimate"]),
        se       = c(r_naive["se"], r_mlm["se"]),
        p_value  = c(r_naive["p_value"], r_mlm["p_value"]),
        row.names = NULL
    )
}

# Summarise: bias, SEs and rejection rate (Type I error, or power) per model
summarise_sim <- function(sim_out, effect_size) {
    do.call(rbind, lapply(split(sim_out, sim_out$model), function(d) {
        data.frame(
            model = d$model[1],
            bias = mean(d$estimate, na.rm = TRUE) - effect_size,
            mean_se = mean(d$se, na.rm = TRUE),
            emp_sd = sd(d$estimate, na.rm = TRUE),
            # "true" SE for comparison
            rejection_rate = mean(d$reject, na.rm = TRUE),
            n_dropped = sum(is.na(d$estimate)),
            row.names = NULL
        )
    }))
}
