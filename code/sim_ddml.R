library(tidyverse)
library(torch)

source("code/helpers_ddml.R")

#---- synth data ------------------------------------------------------
set.seed(3)
n_bootstrap <- 100 
share_labeled <- 0.4
bootstrap_samplesiz <- 10000
epsilon <-  runif(1) # (+) bias of roads on satellite coverage uncorrelated with forest cover.
alpha <- 0.5 # tuning parameter for loss function in DDML
data_sim <- data.frame(W = rpois(20000, lambda = 1)) %>% 
             mutate(roads = pmax(1-W/4, 0),
                    forest_cov = runif(n()),
                    sattelite_cov = forest_cov + ifelse(W > 0, epsilon, 0))

ols_estimates <- data.frame("beta_ols" = NA, "std_ols" = NA)
ddml_estimates <- data.frame("beta_raw_ptd" = NA,
                             "gamma_bias_ptd" = NA, 
                             "beta_corr_ptd" = NA, 
                             "beta_corr_adv" = NA
                             )
x <- data_sim$roads
W <- data_sim$W
y <- data_sim$sattelite_cov
alpha <- 0.5

#---- Bootstrop estimates ------------------------------------------------------
for (i in 1:n_bootstrap){
  bootstrap_sample <- data_sim %>% 
    slice(sample(1:nrow(data_sim), 
                 bootstrap_samplesiz,
                 replace = TRUE
                 ))
  # OLS
  ols_estimates <- rbind(ols_estimates, 
                         get_ols_estimate(bootstrap_sample))
  # Double Debiased lasso. 
  labelled_sample <- sample(1:nrow(bootstrap_sample), 
                          size = floor(share_labeled * nrow(bootstrap_sample)), 
                          replace = FALSE)
  unlablled <- setdiff(1:nrow(bootstrap_sample), labelled_sample)
  labelled_df <- bootstrap_sample[labelled_sample, ]
  unlabelled_df <- bootstrap_sample[unlablled, ]
  
  result_ptd <- ptd_debias(
    y_true_name    = "forest_cov",
    y_pred_name    = "sattelite_cov",
    x_name         = "roads",
    labelled_data  = labelled_df,
    unlabelled_data = unlabelled_df
  )
  
  
  result_adv <- adversarial_debias(
    y_true_name     = "forest_cov",
    y_pred_name     = "sattelite_cov",
    x_name          = "roads",
    W_names         = c("W"
      #"ndvi", "evi", "ndbi"
      ),   # satellite features
    labelled_data   = labelled_df,
    unlabelled_data = unlabelled_df,
    alpha           = 1.0
  )
  
  
  ddml_estimates <- rbind(ddml_estimates, 
                          data.frame(
                           "beta_raw_ptd" = result_ptd$beta_raw,
                           "gamma_bias_ptd" = result_ptd$gamma_bias,
                           "beta_corr_ptd" = result_ptd$beta_corrected,
                           "beta_corr_adv" = result_adv$beta_D
                         ))
}
ddml_estimates <- ddml_estimates %>% slice(-1) # remove initial NA row
ols_estimates <- ols_estimates %>% slice(-1) # remove initial NA row
# Now extract variance from bootstrap estimates and plot CIs
# var: 
ols_var <- var(ols_estimates$beta_ols, na.rm = TRUE)
ptd_var <- var(ddml_estimates$beta_corr_ptd, na.rm = TRUE)
adv_var <- var(ddml_estimates$beta_corr_adv, na.rm = TRUE)

# CIs
ols_ci <- c(mean(ols_estimates$beta_ols, na.rm = TRUE) - 1.96 * sqrt(ols_var),
             mean(ols_estimates$beta_ols, na.rm = TRUE) + 1.96 * sqrt(ols_var))
ptd_ci <- c(mean(ddml_estimates$beta_corr_ptd, na.rm =
 TRUE) - 1.96 * sqrt(ptd_var),
             mean(ddml_estimates$beta_corr_ptd, na.rm = TRUE) + 1.96 * sqrt(ptd_var))
adv_ci <- c(mean(ddml_estimates$beta_corr_adv, na.rm = TRUE)
 - 1.96 * sqrt(adv_var),
             mean(ddml_estimates$beta_corr_adv, na.rm = TRUE) + 1.96 * sqrt(adv_var))

# plot: 
ci_df <- data.frame(
  method = c("OLS", "PTD-corrected", "Adversarial Debiasing"),
  estimate = c(mean(ols_estimates$beta_ols, na.rm = TRUE),
               mean(ddml_estimates$beta_corr_ptd, na.rm = TRUE),
               mean(ddml_estimates$beta_corr_adv, na.rm = TRUE)),
  ci_low = c(ols_ci[1], ptd_ci[1], adv_ci[1]),
  ci_high = c(ols_ci[2], ptd_ci[2], adv_ci[2])
)
names(ddml_estimates)
# join ddml and ols estimates
estimates <- ddml_estimates %>% 
  select(beta_corr_ptd, beta_corr_adv) %>% 
  cbind(ols_estimates %>% select(beta_ols))
# just plot the distribution of bootstrap estimates
p1 <- ggplot(estimates %>%
         pivot_longer(cols = c(beta_corr_ptd, 
                               beta_corr_adv,
                               beta_ols), 
                                  names_to = "method", 
                                  values_to = "estimate")) +
  geom_density(aes(x = estimate,
                   col = method),
               alpha = 0.5,
               ) +
  theme_minimal() +
  labs(title = "Bootstrap Distribution of PTD-corrected and Adversarial Debiasing Estimates",
       x = "Estimated Coefficient (beta)",
       y = "Density") +
       scale_color_manual(values = c("beta_corr_ptd" = "#B22222", 
                                      "beta_corr_adv" = "#DAA520",
                                      "beta_ols" = "#00008B"),
                            labels = c("PTD-corrected", "Adversarial Debiasing", "OLS")) 

p2 <- ggplot(ci_df, aes(x = method, y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2) +
  theme_minimal() +
  labs(title = "Comparison of OLS, PTD-corrected, and Adversarial Debiasing Estimates",
       y = "Estimated Coefficient (beta)",
       x = "Method")

# save plots
ggsave("output/plots/ptd_adv_comparison_density.png", p1, width = 8, height = 6)
ggsave("output/plots/ptd_adv_comparison_ci.png", p2, width = 8, height = 6)
