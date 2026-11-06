########################### FUNCTIONS #####################################

# Simulate one dataset
simulate_data <- function(n_subjects, n_images, effect_size, icc_img, icc_subj, prop_slope_var, seed) {
    
    # Experimental conditions (Interverntion: labeling, matching. Valence: netural, negative)
    cond <- expand.grid(valence = c(-1, 1), intervention = c(-1, 1))
    # 1 and -1 because 1 and 0 implies that the ICC is describing the negative labeling
    
    # Stimuli
    img <- cond[rep(1:4, each = n_images), ]
    img$image <- 1:nrow(img)
    img$set <- rep(1:4, each = n_images)
    
    # Participants id
    subj_id <- rep(1:n_subjects, each = nrow(img))
    
    # Variances
    total_var <- 1
    var_subj <- total_var * icc_subj
    var_img <- total_var * icc_img
    slope_var <- total_var * prop_slope_var
    var_e <- 1 - var_subj - var_img - slope_var
    
    # Random
    set.seed(seed)
    img$b_subj <- rnorm(n = n_subjects, mean = 0, sd = sqrt(var_subj))   # random intercept per subject N(0,var_subj)
    img$b_image <- rnorm(n = nrow(img), mean = 0, sd = sqrt(var_img))   # random intercept per image N(0,var_img)
    e <- rnorm(n = n_subjects*nrow(img), mean = 0, sd = sqrt(var_e)) # residual
    
    # TODO: Finish the part of this code
    dat <- img[rep(1:nrow(img), each = ), ]
    dat <- cbind(subj_id, dat)
    
    # TODO: I am not sure what should be the model to generate the data.Should be model 4?
    # Effect size would be the effect of the interaction?
    dat$Y <- -0.005 + (-0.686 * dat$valence) + (0.054 * dat$intervention) + effect_size * dat$valence * dat$intervention +
        dat$b_image + dat$b_subj + 
    
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
