# Runs the simulations for the results presented in Figure 6 (Appendix E.4)
library(ggplot2)
devtools::install_github("elliot-young/sandwich.boost")

# Data generating mechanism
generate_data <- function(N, I) {
  n_obs = N*I
  beta = 1
  X = runif(n_obs,-2,2)
  xi = rnorm(n_obs,0,3)
  D = -6*exp(-X) + xi
  U = rnorm(n_obs,0,1)
  epsilon = as.matrix(rep(0,n_obs))
  for (i in seq_len(I)) {
    Sigma0 = matrix(0,N,N)
    f_0 <- 2+tanh(D[(N*(i-1)+1):(N*i)]-3*X[(N*(i-1)+1):(N*i)])
    Sigma0 = matrix(0.2,N,N) + diag(0.8,N)
    Sigma0 = diag(f_0)  %*% Sigma0 %*% diag(f_0)
    epsilon[(N*(i-1)+1):(N*i)] = ((expm::sqrtm(Sigma0)))%*%as.vector(U[(N*(i-1)+1):(N*i)])
  }
  Y = beta*D + tanh(X) + epsilon
  rdf <- data.frame(X=X, D=D, Y=Y, id=rep(1:I, each=N))
  return(rdf)
}

# Run the smaller sample simulations (right plot)
M = 500
Ks = c(2,5,10)
Ss = c(1,5)
small_sims <- list()
for (K in Ks) {
  for (S in Ss) {
    betas <- numeric()
    while (length(betas) < M) {
      set.seed(1+length(betas))
      rdf <- generate_data(N=4, I=500)
      beta_SB <- sandwich.boost::PLR_sandboost(l_formula=Y~s(X,bs="cr"), l_learner="gam",
                                               m_formula=D~s(X,bs="cr"), m_learner="gam",
                                               s_formula=u_s~s(X,bs="cr"), s_learner="gam",
                                               proxyCCF="equicorr", data=rdf, K=K, S=S,
                                               variable_steps = FALSE, m_stop=500,
                                               lambda_s=100, lambda_theta=0.1)$coefficients[1]
      betas <- append(betas, beta_SB)
    }
    small_sims[[paste0("K", K, "S", S)]] <- betas
  }
}

# Plot for smaller sample
betas_sq_small <- lapply(small_sims, function(x) 2000*(x-1)^2)
MSE_small <- unlist(lapply(betas_sq_small, mean))
SD_small <- unlist(lapply(betas_sq_small, function(x) sqrt(var(x)/M)))
MSE_small_lb = MSE_small - qnorm(0.975)*SD_small
MSE_small_ub = MSE_small + qnorm(0.975)*SD_small
small_data <- data.frame(
  Folds = factor(rep(c("2","5","10"), each = 2), levels = c("2","5","10")),
  Crossfit_repeats = factor(rep(c("1","5"), each = 1, times=3), levels = c("1","5")),
  MSEs = MSE_small,
  MSEsLB = MSE_small_lb,
  MSEsUB = MSE_small_ub
)
small_plot <- ggplot(small_data, aes(x = Folds, y = MSEs, color = Crossfit_repeats)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = MSEsLB, ymax = MSEsUB),
                position = position_dodge(width = 0.5), width = 0.25) +
  labs(x = "K = Number of folds for cross fitting", y = expression(N * (hat(beta) - beta[0])^2)) +
  scale_color_discrete(name = "S = Number of \ncross fitting repeats") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(expression(paste("Smaller Sample Size N=", 2000)))


# Run the large sample simulations (right plot)
M = 500
Ks = c(2,5,10)
large_sims <- list()
for (K in Ks) {
    betas <- numeric()
    while (length(betas) < M) {
      set.seed(1+length(betas))
      rdf <- generate_data(N=8, I=2^12)
      beta_SB <- sandwich.boost::PLR_sandboost(l_formula=Y~s(X,bs="cr"), l_learner="gam",
                                               m_formula=D~s(X,bs="cr"), m_learner="gam",
                                               s_formula=u_s~s(X,bs="cr"), s_learner="gam",
                                               proxyCCF="equicorr", data=rdf, K=K, S=1,
                                               variable_steps = FALSE, m_stop=500,
                                               lambda_s=100, lambda_theta=0.1)$coefficients[1]
      betas <- append(betas, beta_SB)
    }
    large_sims[[paste0("K", K)]] <- betas
}

# Plot for large sample
betas_sq_large <- lapply(large_sims, function(x) 2^15*(x-1)^2)
MSE_large <- unlist(lapply(betas_sq_large, mean))
SD_large <- unlist(lapply(betas_sq_large, function(x) sqrt(var(x)/M)))
MSE_large_lb = MSE_large - qnorm(0.975)*SD_large
MSE_large_ub = MSE_large + qnorm(0.975)*SD_large
large_data <- data.frame(
  Folds = factor(rep(c("2","5","10"), each = 1), levels = c("2","5","10")),
  MSEs = MSE_large,
  MSEsLB = MSE_large_lb,
  MSEsUB = MSE_large_ub
)
large_sample_plot <- ggplot(boosting_data_largesample, aes(x = Folds, y = Folds)) +
  geom_point(size=2, color=scales::hue_pal()(1)) +
  geom_errorbar(aes(ymin = MSEsLB, ymax = MSEsUB), width = 0.2, color=scales::hue_pal()(1)) +
  theme_minimal() +
  labs(x = "K = Number of folds for cross fitting", y = expression(N * (hat(beta) - beta[0])^2), xlabel="Folds") +
  theme(axis.text.x = element_text(size = 12)) +
  theme(plot.title = element_text(hjust = 0.5)) +
#  ylim(0.0170, 0.0205) +
  ggtitle(expression(paste("Large Sample Size N=", 2^15, " (Simulation 5.1.4)")))

# Plot both plots together
gridExtra::grid.arrange(large_plot, small_plot, ncol=2, widths = c(1, 1.5))
