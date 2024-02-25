# Runs the numerical simulations for Section 5.1
# Note: As the simulations are over a substantial number of model repeats, this
#  takes a substantial length of time to run (it would be recommended to
#  parallelise).
# Note: We provide the outputs for the simulations of Section 5.1 in the folder
#  simulation.5.1.output

library(mvtnorm)
library(fBasics)
library(magic)
library(expm)
library(jmuOutlier)
library(lme4)
library(caret)
library(foreach)
library(mgcv)
library(geepack)
library(nlme)
# Load relevant package for double machine learning & sandwich boosting
devtools::install_github("elliot-young/sandwich.boost")


# Simulation 5.1.4 - Conditional Variance Misspecification
Simulation.5.1.4.Betas <- Simulation.5.1.4.Variances <- list()
Group_Sizes = 2^(1:8) # Group sizes
M = 500 # Number of model repeats
for (N in Group_Sizes) {
  BETAS.DF <- data.frame()
  Variances.DF <- data.frame()
  # Relevant info for data generation mechanism (Example 5.1.4)
  {
    I = 2^15 / N
    n_obs = N*I
    g_0 = function(x) 1*tanh(x[,1])
    m_0 = function(x) -6*exp(-x[,1])
    f_0 = function(d,x) 2+tanh(d[,1]-3*x[,1])
    beta = 1
  }
  while (nrow(BETAS.DF) < M) {
    # Data generation mechanism (Example 5.1.4)
    set.seed(1+nrow(BETAS.DF))
    {
      X = as.matrix(runif(n_obs,-2,2))
      V = as.matrix(rnorm(n_obs,0,1))
      U = as.matrix(rnorm(n_obs,0,1))
      epsilon = as.matrix(rep(0,n_obs))
      xi = as.matrix(rep(0,n_obs))
      rho_const=0

      sqrtOmega0 = diag(rep(3,N))

      for (i in seq_len(I)) {
        xi[(N*(i-1)+1):(N*i),] = (sqrtOmega0)%*%V[(N*(i-1)+1):(N*i),]
      }

      D = m_0(X) + xi

      for (i in seq_len(I)) {
        a = as.matrix(X[(N*(i-1)+1):(N*i),])
        b = as.matrix(D[(N*(i-1)+1):(N*i),])
        Sigma0 = matrix(0,N,N)
        f_0a <- f_0(b,a)
        Sigma0 = matrix(0.2,N,N) + diag(0.8,N)
        Sigma0 = diag(f_0a)  %*% Sigma0 %*% diag(f_0a)
        epsilon[(N*(i-1)+1):(N*i),] = ((sqrtm(Sigma0)))%*%U[(N*(i-1)+1):(N*i),]
      }

      Y = beta*D + g_0(X) + epsilon

      rdf <- data.frame(X=X, D=D, Y=Y, id=rep(1:I, each=N))
    }

    # Split data into 2 folds
    cv_folds <- groupKFold(rdf$id,k=2)

      beta_hat_k_hetero_marginal_gee = beta_hat_k_ar1_hetero_lik =
      beta_hat_k_ar1_hom_lik = beta_hat_k_hom_marginal_gee =
      beta_hat_num_k_boosting = beta_hat_den_k_boosting =
      beta_hat_num_k_ora = beta_hat_den_k_ora =
      beta_hat_num_k_unw = beta_hat_den_k_unw =
      beta_hat_k_hom_marginal_gee_num = beta_hat_k_hom_marginal_gee_den =
      beta_hat_k_het_GEE_num_k = beta_hat_k_het_GEE_den_k =
      beta_hat_k_hom_gls_num = beta_hat_k_hom_gls_den =
      beta_hat_k_het_gls_num = beta_hat_k_het_gls_den =
      V_num_hat_k_boosting = V_num_hat_k_het_gee = V_num_hat_k_hom_gee =
      V_num_hat_k_hom_gls = V_num_hat_k_unw = V_num_hat_k_het_gls =
      V_s_hom_gls = beta_hat_k_hetero_lik = beta_hat_k_hom_lik =
      numeric(2)


    for (k in seq_len(2)) {
      # Separate data into folds for cross-fitting
      cv_fold <- cv_folds[[k]]
      data.nuisance <- rdf[cv_fold,]
      data.beta <- rdf[-cv_fold,]

      # Fit nuisance functions (common across all estimators)
      fit.l <- sandwich.boost:::regress.l(Y ~ s(X, bs = "cr"), "gam", data.nuisance, data.beta)
      l_residuals.nuisance <- fit.l$l_residuals.nuisance; l_residuals.beta <- fit.l$l_residuals.beta

      fit.m <- sandwich.boost:::regress.m(D ~ s(X, bs = "cr"), "gam", data.nuisance, data.beta)
      m_residuals.nuisance <- fit.m$m_residuals.nuisance; m_residuals.beta <- fit.m$m_residuals.beta

      xi_hat.nuisance <- m_residuals.nuisance
      beta_unweighted.nuisance <- sum(l_residuals.nuisance*xi_hat.nuisance)/sum(xi_hat.nuisance^2)
      epsilon_hat.nuisance <- l_residuals.nuisance - beta_unweighted.nuisance*xi_hat.nuisance

      xi_hat.beta <- m_residuals.beta
      beta_unweighted.beta <- sum(l_residuals.beta*xi_hat.beta)/sum(xi_hat.beta^2)
      epsilon_hat.beta <- l_residuals.beta - beta_unweighted.beta*xi_hat.beta

      # Prepare (s,theta) boosting
      boostdf.nuisance <- cbind(data.nuisance, Y_minus_l_hat=l_residuals.nuisance, epsilon_hat=epsilon_hat.nuisance , xi_hat=xi_hat.nuisance)
      boostdf.beta <- cbind(data.beta, Y_minus_l_hat=l_residuals.beta, epsilon_hat=epsilon_hat.beta , xi_hat=xi_hat.beta)
      # Perform (s,theta) boosting
      s.boost.output <- sandwich.boost:::s.boost.constep(boostdf.nuisance, boostdf.beta, s_formula= u_s ~ s(X, bs = "cr"), s_learner="gam", proxyCCF="equicorr", m_stop=500, lambda_s=100, lambda_theta=0.1)
      s.beta <- s.boost.output$s.beta; theta <- s.boost.output$theta
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, s.beta, theta, proxyCCF="equicorr")
      beta_hat_num_k_boosting[k] <- cross_fit_evals$beta_num
      beta_hat_den_k_boosting[k] <- cross_fit_evals$beta_den
      V_num_hat_k_boosting[k] <- cross_fit_evals$V_num

      # Heteroscedastic GEE
      het.gee.nuisance <- cbind(boostdf.nuisance, squared_epsilon_hat=boostdf.nuisance[,"epsilon_hat"]^2)
      het.gee.beta <- boostdf.beta
      gam_fit <- gam(squared_epsilon_hat ~ s(X,bs="cr"), data=het.gee.nuisance) #cubic splines
      sigma.sq.nuisance <- predict(gam_fit, het.gee.nuisance[,"X",drop=FALSE])
      sigma.sq.nuisance[sigma.sq.nuisance<0.01]=0.01
      s.nuisance <- 1/sqrt(sigma.sq.nuisance)
      sigma.sq.beta <- predict(gam_fit, het.gee.beta[,"X",drop=FALSE])
      sigma.sq.beta[sigma.sq.beta<0.01]=0.01
      s.beta <- 1/sqrt(sigma.sq.beta)
      het.gee.nuisance[,c("D","Y","Y_minus_l_hat","epsilon_hat","xi_hat")] <- het.gee.nuisance[,c("D","Y","Y_minus_l_hat","epsilon_hat","xi_hat")]*s.nuisance
      het.gee.beta[,c("D","Y","Y_minus_l_hat","epsilon_hat","xi_hat")] <- het.gee.beta[,c("D","Y","Y_minus_l_hat","epsilon_hat","xi_hat")]*s.beta

      gee_fit <- geeglm(Y_minus_l_hat ~ xi_hat, data=het.gee.nuisance, corstr="exchangeable", id=id)
      alpha_or_rhoo <- as.numeric(gee_fit$geese[2])
      theta <- alpha_or_rhoo/(1-alpha_or_rhoo)
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, s.beta, theta, proxyCCF="equicorr")
      beta_hat_k_het_GEE_num_k[k] <- cross_fit_evals$beta_num
      beta_hat_k_het_GEE_den_k[k] <- cross_fit_evals$beta_den
      V_num_hat_k_het_gee[k] <- cross_fit_evals$V_num

      # Heteroscedastic ML (gls)
      gls_fit <- gls(Y_minus_l_hat ~ xi_hat, data=boostdf.nuisance, correlation = corCompSymm(form = ~ 1 | id), weights = varComb(varExp(form=~X),varExp(form=~I(X^2)),varExp(form=~I(X^3))), control=glsControl(returnObject=TRUE))
      #gls_fit <- gls(Y_minus_l_hat ~ xi_hat, data=boostdf.nuisance, correlation = corCompSymm(form = ~ 1 | id), weights = varComb(varExp(form=~X),varExp(form=~I(X^2)),varExp(form=~I(X^3))), control=glsControl(returnObject = TRUE, singular.ok = TRUE, minAbsParApVar = 0, sigma=sqrt(mean((boostdf.nuisance$Y_minus_l_hat-boostdf.nuisance$xi_hat)^2))))
      beta_hat_k_hetero_lik[k] <- summary(gls_fit)$coefficients[2]
      all_cov_terms       <- as.numeric(coef(gls_fit$modelStruct, unconstrained=FALSE))
      rhoo <- all_cov_terms[1]; theta = rhoo/(1-rhoo)
      var_a <- all_cov_terms[2]; var_b <- all_cov_terms[3]; var_c <- all_cov_terms[4];
      s.beta <- 1/sqrt(exp( 2*var_a*boostdf.beta$X + 2*var_b*boostdf.beta$X^2 + 2*var_c*boostdf.beta$X^3 ))
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, s.beta, theta, proxyCCF="equicorr")
      beta_hat_k_het_gls_num[k] <- cross_fit_evals$beta_num
      beta_hat_k_het_gls_den[k] <- cross_fit_evals$beta_den
      V_num_hat_k_het_gls[k] <- cross_fit_evals$V_num

      # Homogeneous (hom) ML Estimator
      gls_fit <- gls(Y_minus_l_hat ~ xi_hat, data=boostdf.nuisance, correlation = corCompSymm(form = ~ 1 | id), control=glsControl(sigma=sqrt(mean((boostdf.nuisance$Y_minus_l_hat-boostdf.nuisance$xi_hat)^2))))
      rhoo       <- as.numeric(coef(gls_fit$modelStruct, unconstrained=FALSE))
      beta_num = 0; beta_den = 0; V_num = 0
      theta = rhoo/(1-rhoo)
      s.beta <- rep(1,dim(boostdf.beta)[1])
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, s.beta, theta, proxyCCF="equicorr")
      beta_hat_k_hom_gls_num[k] <- cross_fit_evals$beta_num
      beta_hat_k_hom_gls_den[k] <- cross_fit_evals$beta_den
      V_num_hat_k_hom_gls[k] <- cross_fit_evals$V_num
      #V_s_hom_gls[k] <- gls_fit$varBeta[2,2] # Allows as a sanity check, one to check DML1 versus DML2 for this example


      # Homogeneous (hom) GEE Estimator
      gee_fit <- geeglm(Y_minus_l_hat ~ xi_hat, data=boostdf.nuisance, corstr="exchangeable", id=id)
      alpha_or_rhoo <- as.numeric(gee_fit$geese[2])
      theta <- alpha_or_rhoo/(1-alpha_or_rhoo)
      s.beta <- rep(1,dim(boostdf.beta)[1])
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, s.beta, theta, proxyCCF="equicorr")
      beta_hat_k_hom_marginal_gee_num[k] <- cross_fit_evals$beta_num
      beta_hat_k_hom_marginal_gee_den[k] <- cross_fit_evals$beta_den
      V_num_hat_k_hom_gee[k] <- cross_fit_evals$V_num

      # Unweighted Estimator
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, rep(1,nrow(boostdf.beta)), 0, proxyCCF="equicorr")
      beta_hat_num_k_unw[k] <- cross_fit_evals$beta_num
      beta_hat_den_k_unw[k] <- cross_fit_evals$beta_den
      V_num_hat_k_unw[k] <- cross_fit_evals$V_num

    }

    Beta_boosting = sum(beta_hat_num_k_boosting)/sum(beta_hat_den_k_boosting)
    Beta_unw = sum(beta_hat_num_k_unw)/sum(beta_hat_den_k_unw)
    Beta_het_gee = sum(beta_hat_k_het_GEE_num_k)/sum(beta_hat_k_het_GEE_den_k)
    Beta_hom_gee = sum(beta_hat_k_hom_marginal_gee_num)/sum(beta_hat_k_hom_marginal_gee_den)
    Beta_het_lik = sum(beta_hat_k_het_gls_num)/sum(beta_hat_k_het_gls_den)
    Beta_hom_lik = sum(beta_hat_k_hom_gls_num)/sum(beta_hat_k_hom_gls_den)

    # Collect everthing together
    BETAS.DF <- rbind( BETAS.DF, data.frame(Beta_boosting=Beta_boosting, Beta_het_lik=Beta_het_lik, Beta_hom_lik=Beta_hom_lik, Beta_het_gee=Beta_het_gee, Beta_hom_gee=Beta_hom_gee, Beta_unw=Beta_unw) )
    Variances.DF <- rbind(Variances.DF, data.frame(Boosting=I*N*sum(V_num_hat_k_boosting)/(sum(beta_hat_den_k_boosting))^2 ,HetGEE=I*N*sum(V_num_hat_k_het_gee)/(sum(beta_hat_k_het_GEE_den_k))^2, HomGEE=I*N*sum(V_num_hat_k_hom_gee)/(sum(beta_hat_k_hom_marginal_gee_den))^2, Unweighted=I*N*sum(V_num_hat_k_unw)/(sum(beta_hat_den_k_unw))^2, homGLS=I*N*sum(V_num_hat_k_hom_gls)/(sum(beta_hat_k_hom_gls_den))^2, HetGLS=I*N*sum(V_num_hat_k_het_gls)/(sum(beta_hat_k_het_gls_den))^2 ))

#    Uncomment to display progress
#    cat("\n")
#    cat('Progress: ',dim(BETAS.DF)[1],' of ',M,' complete')
  }
  Simulation.5.1.4.Betas[[which(N==Group_Sizes)]] = BETAS.DF
  Simulation.5.1.4.Variances[[which(N==Group_Sizes)]] =Variances.DF
}


# Simulation 5.1.3 - Mild Conditional Correlation Misspecification
Simulation.5.1.3.Betas <- Simulation.5.1.3.Variances <- list()
Group_Sizes = 2^(1:8) # Group sizes
M = 500 # Number of model repeats
for (N in Group_Sizes) {
  BETAS.DF <- data.frame()
  Variances.DF <- data.frame()
  # Relevant info for data generation mechanism (Example 5.1.4)
  {
    I = 2^15/N
    n_obs = N*I
    Sigma0 = toeplitz(ARMAacf(ar=c(0.3,0.6), ma=c(-0.5), lag.max = N-1))
    sqrtSigma0 = sqrtm(Sigma0)
    Omega_0 = matrix(1/8,N,N) + diag(7/8,N)
    sqrtOmega0 = sqrtm(Omega_0)
    beta = 1
  }
  while (nrow(BETAS.DF) < M) {
    # Data generation mechanism (Example 5.1.4)
    set.seed(1+nrow(BETAS.DF))
    {
      X = as.matrix(runif(n_obs))
      V = as.matrix(rnorm(n_obs))
      U = as.matrix(rnorm(n_obs))
      epsilon = xi = as.matrix(rep(0,n_obs))

      for (i in seq_len(I)) {
        xi[(N*(i-1)+1):(N*i),] = (sqrtOmega0)%*%V[(N*(i-1)+1):(N*i),]
        epsilon[(N*(i-1)+1):(N*i),] = (sqrtSigma0)%*%U[(N*(i-1)+1):(N*i),]
      }
      D = xi
      Y = beta*D + epsilon

      rdf <- data.frame(X=X, D=D, Y=Y, id=rep(1:I, each=N))
    }

    # Split data into 2 folds
    cv_folds <- groupKFold(rdf$id,k=2)

    beta_hat_k_hetero_marginal_gee = beta_hat_k_ar1_hetero_lik =
      beta_hat_k_ar1_hom_lik = beta_hat_k_hom_marginal_gee =
      beta_hat_num_k_boosting = beta_hat_den_k_boosting =
      beta_hat_num_k_ora = beta_hat_den_k_ora =
      beta_hat_num_k_unw = beta_hat_den_k_unw =
      beta_hat_k_hom_marginal_gee_num = beta_hat_k_hom_marginal_gee_den =
      beta_hat_k_het_GEE_num_k = beta_hat_k_het_GEE_den_k =
      beta_hat_k_hom_gls_num = beta_hat_k_hom_gls_den =
      beta_hat_k_het_gls_num = beta_hat_k_het_gls_den =
      V_num_hat_k_boosting = V_num_hat_k_het_gee = V_num_hat_k_hom_gee =
      V_num_hat_k_hom_gls = V_num_hat_k_unw = V_num_hat_k_het_gls =
      V_s_hom_gls = beta_hat_k_hetero_lik = beta_hat_k_hom_lik =
      numeric(2)


    for (k in seq_len(2)) {
      # Separate data into folds for cross-fitting
      cv_fold <- cv_folds[[k]]
      data.nuisance <- rdf[cv_fold,]
      data.beta <- rdf[-cv_fold,]

      # Fit nuisance functions (common across all estimators)
      l_residuals.nuisance <- data.nuisance$Y; l_residuals.beta <- data.beta$Y
      m_residuals.nuisance <- data.nuisance$D; m_residuals.beta <- data.beta$D
      xi_hat.nuisance <- m_residuals.nuisance
      beta_unweighted.nuisance <- sum(l_residuals.nuisance*xi_hat.nuisance)/sum(xi_hat.nuisance^2)
      epsilon_hat.nuisance <- l_residuals.nuisance - beta_unweighted.nuisance*xi_hat.nuisance
      xi_hat.beta <- m_residuals.beta
      beta_unweighted.beta <- sum(l_residuals.beta*xi_hat.beta)/sum(xi_hat.beta^2)
      epsilon_hat.beta <- l_residuals.beta - beta_unweighted.beta*xi_hat.beta

      # Prepare (s,theta) boosting
      boostdf.nuisance <- cbind(data.nuisance, Y_minus_l_hat=l_residuals.nuisance, epsilon_hat=epsilon_hat.nuisance , xi_hat=xi_hat.nuisance)
      boostdf.beta <- cbind(data.beta, Y_minus_l_hat=l_residuals.beta, epsilon_hat=epsilon_hat.beta , xi_hat=xi_hat.beta)
      # Perform (s,theta) boosting
      s.boost.output <- sandwich.boost:::s.boost.constep(boostdf.nuisance, boostdf.beta, s_formula= u_s ~ s(X, bs = "cr"), s_learner="gam", proxyCCF="autoreg", m_stop=500, lambda_s=0, lambda_theta=0.0001)
      s.beta <- s.boost.output$s.beta; theta <- s.boost.output$theta
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, s.beta, theta, proxyCCF="autoreg")
      beta_hat_num_k_boosting[k] <- cross_fit_evals$beta_num
      beta_hat_den_k_boosting[k] <- cross_fit_evals$beta_den
      V_num_hat_k_boosting[k] <- cross_fit_evals$V_num


      # Homogeneous (hom) ML Estimator
      gls_fit <- gls(Y_minus_l_hat ~ xi_hat, data=boostdf.nuisance, correlation = corAR1(form = ~ 1 | id))
      rhoo <- as.numeric(coef(gls_fit$modelStruct, unconstrained=FALSE))
      beta_num = 0; beta_den = 0; V_num = 0
      theta = rhoo
      s.beta <- rep(1,dim(boostdf.beta)[1])
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, s.beta, theta, proxyCCF="autoreg")
      beta_hat_k_hom_gls_num[k] <- cross_fit_evals$beta_num
      beta_hat_k_hom_gls_den[k] <- cross_fit_evals$beta_den
      V_num_hat_k_hom_gls[k] <- cross_fit_evals$V_num
      #V_s_hom_gls[k] <- gls_fit$varBeta[2,2] # Allows as a sanity check, one to check DML1 versus DML2 for this example


      # Homogeneous (hom) GEE Estimator
      gee_fit <- geeglm(Y_minus_l_hat ~ xi_hat, data=boostdf.nuisance, corstr="ar1", id=id)
      alpha_or_rhoo <- as.numeric(gee_fit$geese[2])
      theta <- alpha_or_rhoo
      s.beta <- rep(1,dim(boostdf.beta)[1])
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, s.beta, theta, proxyCCF="autoreg")
      beta_hat_k_hom_marginal_gee_num[k] <- cross_fit_evals$beta_num
      beta_hat_k_hom_marginal_gee_den[k] <- cross_fit_evals$beta_den
      V_num_hat_k_hom_gee[k] <- cross_fit_evals$V_num

      # Unweighted Estimator
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, rep(1,nrow(boostdf.beta)), 0, proxyCCF="autoreg")
      beta_hat_num_k_unw[k] <- cross_fit_evals$beta_num
      beta_hat_den_k_unw[k] <- cross_fit_evals$beta_den
      V_num_hat_k_unw[k] <- cross_fit_evals$V_num

    }

    Beta_boosting = sum(beta_hat_num_k_boosting)/sum(beta_hat_den_k_boosting)
    Beta_unw = sum(beta_hat_num_k_unw)/sum(beta_hat_den_k_unw)
    Beta_hom_gee = sum(beta_hat_k_hom_marginal_gee_num)/sum(beta_hat_k_hom_marginal_gee_den)
    Beta_hom_lik = sum(beta_hat_k_hom_gls_num)/sum(beta_hat_k_hom_gls_den)

    # Collect everthing together
    BETAS.DF <- rbind( BETAS.DF, data.frame(Beta_boosting=Beta_boosting, Beta_hom_lik=Beta_hom_lik, Beta_hom_gee=Beta_hom_gee, Beta_unw=Beta_unw) )
    Variances.DF <- rbind(Variances.DF, data.frame(Boosting=I*N*sum(V_num_hat_k_boosting)/(sum(beta_hat_den_k_boosting))^2, HomGEE=I*N*sum(V_num_hat_k_hom_gee)/(sum(beta_hat_k_hom_marginal_gee_den))^2, Unweighted=I*N*sum(V_num_hat_k_unw)/(sum(beta_hat_den_k_unw))^2, homGLS=I*N*sum(V_num_hat_k_hom_gls)/(sum(beta_hat_k_hom_gls_den))^2 ))

    #    Uncomment to display progress
    #    cat("\n")
    #    cat('Progress: ',dim(BETAS.DF)[1],' of ',M,' complete')
  }
  Simulation.5.1.3.Betas[[which(N==Group_Sizes)]] = BETAS.DF
  Simulation.5.1.3.Variances[[which(N==Group_Sizes)]] =Variances.DF
}


# Simulation 5.1.2 - Conditional Covariance Misspecification
Simulation.5.1.2.Betas <- Simulation.5.1.2.Variances <- list()
Etas = seq(10,100,by=10) # Misspecification parameters
M = 500 # Number of model repeats
for (eta in Etas) {
  BETAS.DF <- data.frame()
  Variances.DF <- data.frame()
  # Relevant info for data generation mechanism (Example 5.1.4)
  {
    N = 4
    I = 10^4
    n_obs = N*I
    g_0 = function(x) tanh(x[,1])
    m_0 = function(x) cos(x[,1])
    beta = 1
    Sigma0 = toeplitz(ARMAacf(ar=c(0.2), ma=c(), lag.max = N-1))
    sqrtSigma0 = sqrtm(Sigma0)
    CovX = matrix(0.9,N,N) + diag(0.1,N)
    sqrtCovX = sqrtm(CovX)
  }
  while (nrow(BETAS.DF) < M) {
    # Data generation mechanism (Example 5.1.4)
    set.seed(1+nrow(BETAS.DF))
    {
      X = as.matrix(rnorm(n_obs,0,1))
      for (i in 1:I) X[((i-1)*N+1):(i*N)] = sqrtCovX %*% X[((i-1)*N+1):(i*N)]
      V = as.matrix(rnorm(n_obs,0,1))
      U = as.matrix(rnorm(n_obs,0,1))
      xi = epsilon = D = as.matrix(rep(0,n_obs))

      for (i in seq_len(I)) {
        a = as.matrix(X[(N*(i-1)+1):(N*i),])
        p <- 1/eta*as.numeric(mean(a)<0) + as.numeric(mean(a)>0)
        B <- rbinom(1,1,p)
        zeta <- B/p
        xi[(N*(i-1)+1):(N*i),] = V[(N*(i-1)+1):(N*i),] * sqrt(zeta)
        D[(N*(i-1)+1):(N*i),] = m_0(a) + xi[(N*(i-1)+1):(N*i),]
        epsilon[(N*(i-1)+1):(N*i),] = (sqrtSigma0)%*%U[(N*(i-1)+1):(N*i),] * sqrt(zeta)
      }

      Y = beta*D + g_0(X) + epsilon

      rdf <- data.frame(X=X, D=D, Y=Y, id=rep(1:I, each=N))
    }

    # Split data into 2 folds
    cv_folds <- groupKFold(rdf$id,k=2)

    beta_hat_k_hetero_marginal_gee = beta_hat_k_ar1_hetero_lik =
      beta_hat_k_ar1_hom_lik = beta_hat_k_hom_marginal_gee =
      beta_hat_num_k_boosting = beta_hat_den_k_boosting =
      beta_hat_num_k_ora = beta_hat_den_k_ora =
      beta_hat_num_k_unw = beta_hat_den_k_unw =
      beta_hat_k_hom_marginal_gee_num = beta_hat_k_hom_marginal_gee_den =
      beta_hat_k_het_GEE_num_k = beta_hat_k_het_GEE_den_k =
      beta_hat_k_hom_gls_num = beta_hat_k_hom_gls_den =
      beta_hat_k_het_gls_num = beta_hat_k_het_gls_den =
      V_num_hat_k_boosting = V_num_hat_k_het_gee = V_num_hat_k_hom_gee =
      V_num_hat_k_hom_gls = V_num_hat_k_unw = V_num_hat_k_het_gls =
      V_s_hom_gls = beta_hat_k_hetero_lik = beta_hat_k_hom_lik =
      numeric(2)


    for (k in seq_len(2)) {
      # Separate data into folds for cross-fitting
      cv_fold <- cv_folds[[k]]
      data.nuisance <- rdf[cv_fold,]
      data.beta <- rdf[-cv_fold,]

      # Fit nuisance functions (common across all estimators)
      fit.l <- sandwich.boost:::regress.l(Y ~ s(X, bs = "cr"), "gam", data.nuisance, data.beta)
      l_residuals.nuisance <- fit.l$l_residuals.nuisance; l_residuals.beta <- fit.l$l_residuals.beta

      fit.m <- sandwich.boost:::regress.m(D ~ s(X, bs = "cr"), "gam", data.nuisance, data.beta)
      m_residuals.nuisance <- fit.m$m_residuals.nuisance; m_residuals.beta <- fit.m$m_residuals.beta

      xi_hat.nuisance <- m_residuals.nuisance
      beta_unweighted.nuisance <- sum(l_residuals.nuisance*xi_hat.nuisance)/sum(xi_hat.nuisance^2)
      epsilon_hat.nuisance <- l_residuals.nuisance - beta_unweighted.nuisance*xi_hat.nuisance

      xi_hat.beta <- m_residuals.beta
      beta_unweighted.beta <- sum(l_residuals.beta*xi_hat.beta)/sum(xi_hat.beta^2)
      epsilon_hat.beta <- l_residuals.beta - beta_unweighted.beta*xi_hat.beta

      # Prepare (s,theta) boosting
      boostdf.nuisance <- cbind(data.nuisance, Y_minus_l_hat=l_residuals.nuisance, epsilon_hat=epsilon_hat.nuisance , xi_hat=xi_hat.nuisance)
      boostdf.beta <- cbind(data.beta, Y_minus_l_hat=l_residuals.beta, epsilon_hat=epsilon_hat.beta , xi_hat=xi_hat.beta)
      # Perform (s,theta) boosting
      s.boost.output <- sandwich.boost:::s.boost.constep(boostdf.nuisance, boostdf.beta, s_formula= u_s ~ s(X, bs = "cr"), s_learner="gam", proxyCCF="autoreg", m_stop=200, lambda_s=10, lambda_theta=0.1)
      s.beta <- s.boost.output$s.beta; theta <- s.boost.output$theta
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, s.beta, theta, proxyCCF="autoreg")
      beta_hat_num_k_boosting[k] <- cross_fit_evals$beta_num
      beta_hat_den_k_boosting[k] <- cross_fit_evals$beta_den
      V_num_hat_k_boosting[k] <- cross_fit_evals$V_num

      # Heteroscedastic GEE
      het.gee.nuisance <- cbind(boostdf.nuisance, squared_epsilon_hat=boostdf.nuisance[,"epsilon_hat"]^2)
      het.gee.beta <- boostdf.beta
      gam_fit <- gam(squared_epsilon_hat ~ s(X,bs="cr"), data=het.gee.nuisance) #cubic splines
      sigma.sq.nuisance <- predict(gam_fit, het.gee.nuisance[,"X",drop=FALSE])
      sigma.sq.nuisance[sigma.sq.nuisance<0.01]=0.01
      s.nuisance <- 1/sqrt(sigma.sq.nuisance)
      sigma.sq.beta <- predict(gam_fit, het.gee.beta[,"X",drop=FALSE])
      sigma.sq.beta[sigma.sq.beta<0.01]=0.01
      s.beta <- 1/sqrt(sigma.sq.beta)
      het.gee.nuisance[,c("D","Y","Y_minus_l_hat","epsilon_hat","xi_hat")] <- het.gee.nuisance[,c("D","Y","Y_minus_l_hat","epsilon_hat","xi_hat")]*s.nuisance
      het.gee.beta[,c("D","Y","Y_minus_l_hat","epsilon_hat","xi_hat")] <- het.gee.beta[,c("D","Y","Y_minus_l_hat","epsilon_hat","xi_hat")]*s.beta

      gee_fit <- geeglm(Y_minus_l_hat ~ xi_hat, data=het.gee.nuisance, corstr="ar1", id=id)
      alpha_or_rhoo <- as.numeric(gee_fit$geese[2])
      theta <- alpha_or_rhoo
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, s.beta, theta, proxyCCF="autoreg")
      beta_hat_k_het_GEE_num_k[k] <- cross_fit_evals$beta_num
      beta_hat_k_het_GEE_den_k[k] <- cross_fit_evals$beta_den
      V_num_hat_k_het_gee[k] <- cross_fit_evals$V_num

      # Heteroscedastic ML (gls)
      gls_fit <- gls(Y_minus_l_hat ~ xi_hat, data=boostdf.nuisance, correlation = corAR1(form = ~ 1 | id), weights = varComb(varExp(form=~X),varExp(form=~I(X^2)),varExp(form=~I(X^3))), control=glsControl(returnObject=TRUE))
      #gls_fit <- gls(Y_minus_l_hat ~ xi_hat, data=boostdf.nuisance, correlation = corCompSymm(form = ~ 1 | id), weights = varComb(varExp(form=~X),varExp(form=~I(X^2)),varExp(form=~I(X^3))), control=glsControl(returnObject = TRUE, singular.ok = TRUE, minAbsParApVar = 0, sigma=sqrt(mean((boostdf.nuisance$Y_minus_l_hat-boostdf.nuisance$xi_hat)^2))))
      beta_hat_k_hetero_lik[k] <- summary(gls_fit)$coefficients[2]
      all_cov_terms       <- as.numeric(coef(gls_fit$modelStruct, unconstrained=FALSE))
      rhoo <- all_cov_terms[1]; theta = rhoo
      var_a <- all_cov_terms[2]; var_b <- all_cov_terms[3]; var_c <- all_cov_terms[4];
      s.beta <- 1/sqrt(exp( 2*var_a*boostdf.beta$X + 2*var_b*boostdf.beta$X^2 + 2*var_c*boostdf.beta$X^3 ))
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, s.beta, theta, proxyCCF="autoreg")
      beta_hat_k_het_gls_num[k] <- cross_fit_evals$beta_num
      beta_hat_k_het_gls_den[k] <- cross_fit_evals$beta_den
      V_num_hat_k_het_gls[k] <- cross_fit_evals$V_num

      # Homogeneous (hom) ML Estimator
      gls_fit <- gls(Y_minus_l_hat ~ xi_hat, data=boostdf.nuisance, correlation = corAR1(form = ~ 1 | id), control=glsControl(sigma=sqrt(mean((boostdf.nuisance$Y_minus_l_hat-boostdf.nuisance$xi_hat)^2))))
      rhoo       <- as.numeric(coef(gls_fit$modelStruct, unconstrained=FALSE))
      beta_num = 0; beta_den = 0; V_num = 0
      theta = rhoo
      s.beta <- rep(1,dim(boostdf.beta)[1])
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, s.beta, theta, proxyCCF="autoreg")
      beta_hat_k_hom_gls_num[k] <- cross_fit_evals$beta_num
      beta_hat_k_hom_gls_den[k] <- cross_fit_evals$beta_den
      V_num_hat_k_hom_gls[k] <- cross_fit_evals$V_num
      #V_s_hom_gls[k] <- gls_fit$varBeta[2,2] # Allows as a sanity check, one to check DML1 versus DML2 for this example


      # Homogeneous (hom) GEE Estimator
      gee_fit <- geeglm(Y_minus_l_hat ~ xi_hat, data=boostdf.nuisance, corstr="ar1", id=id)
      alpha_or_rhoo <- as.numeric(gee_fit$geese[2])
      theta <- alpha_or_rhoo
      s.beta <- rep(1,dim(boostdf.beta)[1])
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, s.beta, theta, proxyCCF="autoreg")
      beta_hat_k_hom_marginal_gee_num[k] <- cross_fit_evals$beta_num
      beta_hat_k_hom_marginal_gee_den[k] <- cross_fit_evals$beta_den
      V_num_hat_k_hom_gee[k] <- cross_fit_evals$V_num

      # Unweighted Estimator
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, rep(1,nrow(boostdf.beta)), 0, proxyCCF="autoreg")
      beta_hat_num_k_unw[k] <- cross_fit_evals$beta_num
      beta_hat_den_k_unw[k] <- cross_fit_evals$beta_den
      V_num_hat_k_unw[k] <- cross_fit_evals$V_num

    }

    Beta_boosting = sum(beta_hat_num_k_boosting)/sum(beta_hat_den_k_boosting)
    Beta_unw = sum(beta_hat_num_k_unw)/sum(beta_hat_den_k_unw)
    Beta_het_gee = sum(beta_hat_k_het_GEE_num_k)/sum(beta_hat_k_het_GEE_den_k)
    Beta_hom_gee = sum(beta_hat_k_hom_marginal_gee_num)/sum(beta_hat_k_hom_marginal_gee_den)
    Beta_het_lik = sum(beta_hat_k_het_gls_num)/sum(beta_hat_k_het_gls_den)
    Beta_hom_lik = sum(beta_hat_k_hom_gls_num)/sum(beta_hat_k_hom_gls_den)

    # Collect everthing together
    BETAS.DF <- rbind( BETAS.DF, data.frame(Beta_boosting=Beta_boosting, Beta_het_lik=Beta_het_lik, Beta_hom_lik=Beta_hom_lik, Beta_het_gee=Beta_het_gee, Beta_hom_gee=Beta_hom_gee, Beta_unw=Beta_unw) )
    Variances.DF <- rbind(Variances.DF, data.frame(Boosting=I*N*sum(V_num_hat_k_boosting)/(sum(beta_hat_den_k_boosting))^2 ,HetGEE=I*N*sum(V_num_hat_k_het_gee)/(sum(beta_hat_k_het_GEE_den_k))^2, HomGEE=I*N*sum(V_num_hat_k_hom_gee)/(sum(beta_hat_k_hom_marginal_gee_den))^2, Unweighted=I*N*sum(V_num_hat_k_unw)/(sum(beta_hat_den_k_unw))^2, homGLS=I*N*sum(V_num_hat_k_hom_gls)/(sum(beta_hat_k_hom_gls_den))^2, HetGLS=I*N*sum(V_num_hat_k_het_gls)/(sum(beta_hat_k_het_gls_den))^2 ))

    #    Uncomment to display progress
    #    cat("\n")
    #    cat('Progress: ',dim(BETAS.DF)[1],' of ',M,' complete')
  }
  Simulation.5.1.2.Betas[[which(eta==Etas)]] = BETAS.DF
  Simulation.5.1.2.Variances[[which(eta==Etas)]] =Variances.DF
}


# Simulation 5.1.1 - Increasing Model Complexity
Simulation.5.1.1.Betas <- Simulation.5.1.1.Variances <- list()
Lambdas = seq(0.1,2.5,by=0.2) # Complexity parameters
M = 500 # Number of model repeats
for (lambda in Lambdas) {
  BETAS.DF <- data.frame()
  Variances.DF <- data.frame()
  # Relevant info for data generation mechanism (Example 5.1.4)
  {
    N = 10
    I = 2000
    n_obs = N*I
    g_0 = function(x) tanh(x[,1])
    m_0 = function(x) cos(x[,1])
    f_0 = function(d,x) 2+cos(lambda*x[,1])
    sqrtOmega0 = sqrtm(matrix(0.1,N,N) + diag(0.9,N))
    Sigma0_core = matrix(0.2,N,N) + diag(0.8,N)
    beta = 1
  }
  while (nrow(BETAS.DF) < M) {
    # Data generation mechanism (Example 5.1.4)
    set.seed(1+nrow(BETAS.DF))
    {
      X = as.matrix(runif(n_obs,-5,5))
      V = as.matrix(rnorm(n_obs,0,1))
      U = as.matrix(rnorm(n_obs,0,1))
      epsilon = xi = as.matrix(rep(0,n_obs))
      for (i in seq_len(I)) {
        xi[(N*(i-1)+1):(N*i),] = (sqrtOmega0)%*%V[(N*(i-1)+1):(N*i),]
      }
      D = m_0(X) + xi
      for (i in seq_len(I)) {
        a = as.matrix(X[(N*(i-1)+1):(N*i),])
        b = as.matrix(D[(N*(i-1)+1):(N*i),])
        Sigma0 = matrix(0,N,N)
        f_0a <- f_0(b,a)
        Sigma0 = diag(f_0a)  %*% Sigma0_core %*% diag(f_0a)
        epsilon[(N*(i-1)+1):(N*i),] = ((sqrtm(Sigma0)))%*%U[(N*(i-1)+1):(N*i),]
      }

      Y = beta*D + g_0(X) + epsilon

      rdf <- data.frame(X=X, D=D, Y=Y, id=rep(1:I, each=N))
    }

    # Split data into 2 folds
    cv_folds <- groupKFold(rdf$id,k=2)

    beta_hat_k_hetero_marginal_gee = beta_hat_k_ar1_hetero_lik =
      beta_hat_k_ar1_hom_lik = beta_hat_k_hom_marginal_gee =
      beta_hat_num_k_boosting = beta_hat_den_k_boosting =
      beta_hat_num_k_ora = beta_hat_den_k_ora =
      beta_hat_num_k_unw = beta_hat_den_k_unw =
      beta_hat_k_hom_marginal_gee_num = beta_hat_k_hom_marginal_gee_den =
      beta_hat_k_het_GEE_num_k = beta_hat_k_het_GEE_den_k =
      beta_hat_k_hom_gls_num = beta_hat_k_hom_gls_den =
      beta_hat_k_het_gls_num = beta_hat_k_het_gls_den =
      V_num_hat_k_boosting = V_num_hat_k_het_gee = V_num_hat_k_hom_gee =
      V_num_hat_k_hom_gls = V_num_hat_k_unw = V_num_hat_k_het_gls =
      V_s_hom_gls = beta_hat_k_hetero_lik = beta_hat_k_hom_lik =
      V_num_hat_k_ora = numeric(2)


    for (k in seq_len(2)) {
      # Separate data into folds for cross-fitting
      cv_fold <- cv_folds[[k]]
      data.nuisance <- rdf[cv_fold,]
      data.beta <- rdf[-cv_fold,]

      # Fit nuisance functions (common across all estimators)
      fit.l <- sandwich.boost:::regress.l(Y ~ s(X, bs = "cr"), "gam", data.nuisance, data.beta)
      l_residuals.nuisance <- fit.l$l_residuals.nuisance; l_residuals.beta <- fit.l$l_residuals.beta

      fit.m <- sandwich.boost:::regress.m(D ~ s(X, bs = "cr"), "gam", data.nuisance, data.beta)
      m_residuals.nuisance <- fit.m$m_residuals.nuisance; m_residuals.beta <- fit.m$m_residuals.beta

      xi_hat.nuisance <- m_residuals.nuisance
      beta_unweighted.nuisance <- sum(l_residuals.nuisance*xi_hat.nuisance)/sum(xi_hat.nuisance^2)
      epsilon_hat.nuisance <- l_residuals.nuisance - beta_unweighted.nuisance*xi_hat.nuisance

      xi_hat.beta <- m_residuals.beta
      beta_unweighted.beta <- sum(l_residuals.beta*xi_hat.beta)/sum(xi_hat.beta^2)
      epsilon_hat.beta <- l_residuals.beta - beta_unweighted.beta*xi_hat.beta

      # Prepare (s,theta) boosting
      boostdf.nuisance <- cbind(data.nuisance, Y_minus_l_hat=l_residuals.nuisance, epsilon_hat=epsilon_hat.nuisance , xi_hat=xi_hat.nuisance)
      boostdf.beta <- cbind(data.beta, Y_minus_l_hat=l_residuals.beta, epsilon_hat=epsilon_hat.beta , xi_hat=xi_hat.beta)
      # Perform (s,theta) boosting
      s.boost.output <- sandwich.boost:::s.boost.varstep(boostdf.nuisance, boostdf.beta, s_formula= u_s ~ s(X, bs = "cr"), s_learner="gam", proxyCCF="equicorr", m_stop=2000, Lambda_s=c(0.001,10), lambda_theta=0.1, mu_s=0.1)
      s.beta <- s.boost.output$s.beta; theta <- s.boost.output$theta
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, s.beta, theta, proxyCCF="equicorr")
      beta_hat_num_k_boosting[k] <- cross_fit_evals$beta_num
      beta_hat_den_k_boosting[k] <- cross_fit_evals$beta_den
      V_num_hat_k_boosting[k] <- cross_fit_evals$V_num

      # Heteroscedastic GEE
      het.gee.nuisance <- cbind(boostdf.nuisance, squared_epsilon_hat=boostdf.nuisance[,"epsilon_hat"]^2)
      het.gee.beta <- boostdf.beta
      gam_fit <- gam(squared_epsilon_hat ~ s(X,bs="cr"), data=het.gee.nuisance) #cubic splines
      sigma.sq.nuisance <- predict(gam_fit, het.gee.nuisance[,"X",drop=FALSE])
      sigma.sq.nuisance[sigma.sq.nuisance<0.01]=0.01
      s.nuisance <- 1/sqrt(sigma.sq.nuisance)
      sigma.sq.beta <- predict(gam_fit, het.gee.beta[,"X",drop=FALSE])
      sigma.sq.beta[sigma.sq.beta<0.01]=0.01
      s.beta <- 1/sqrt(sigma.sq.beta)
      het.gee.nuisance[,c("D","Y","Y_minus_l_hat","epsilon_hat","xi_hat")] <- het.gee.nuisance[,c("D","Y","Y_minus_l_hat","epsilon_hat","xi_hat")]*s.nuisance
      het.gee.beta[,c("D","Y","Y_minus_l_hat","epsilon_hat","xi_hat")] <- het.gee.beta[,c("D","Y","Y_minus_l_hat","epsilon_hat","xi_hat")]*s.beta

      gee_fit <- geeglm(Y_minus_l_hat ~ xi_hat, data=het.gee.nuisance, corstr="exchangeable", id=id)
      alpha_or_rhoo <- as.numeric(gee_fit$geese[2])
      theta <- alpha_or_rhoo/(1-alpha_or_rhoo)
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, s.beta, theta, proxyCCF="equicorr")
      beta_hat_k_het_GEE_num_k[k] <- cross_fit_evals$beta_num
      beta_hat_k_het_GEE_den_k[k] <- cross_fit_evals$beta_den
      V_num_hat_k_het_gee[k] <- cross_fit_evals$V_num

      # Heteroscedastic ML (gls)
      gls_fit <- gls(Y_minus_l_hat ~ xi_hat, data=boostdf.nuisance, correlation = corCompSymm(form = ~ 1 | id), weights = varComb(varExp(value=0,form=~X),varExp(value=0,form=~I(X^2)),varExp(value=0,form=~I(X^3)),varExp(value=0,form=~I(X^4))), control=glsControl(maxIter=500, msMaxIter = 500, returnObject=TRUE, msVerbose=TRUE, singular.ok=TRUE, minAbsParApVar=0, sigma=1))
      all_cov_terms <- as.numeric(coef(gls_fit$modelStruct, unconstrained=FALSE))
      rhoo <- all_cov_terms[1]; theta = rhoo/(1-rhoo)
      var_a <- all_cov_terms[2]; var_b <- all_cov_terms[3]; var_c <- all_cov_terms[4]; var_d <- all_cov_terms[5];# var_e <- all_cov_terms[6];
      s.beta <- 1/sqrt(exp( 2*var_a*boostdf.beta$X + 2*var_b*boostdf.beta$X^2 + 2*var_c*boostdf.beta$X^3 + 2*var_d*boostdf.beta$X^4 ))
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, s.beta, theta, proxyCCF="equicorr")
      beta_hat_k_het_gls_num[k] <- cross_fit_evals$beta_num
      beta_hat_k_het_gls_den[k] <- cross_fit_evals$beta_den
      V_num_hat_k_het_gls[k] <- cross_fit_evals$V_num

      # Homogeneous (hom) ML Estimator
      gls_fit <- gls(Y_minus_l_hat ~ xi_hat, data=boostdf.nuisance, correlation = corCompSymm(form = ~ 1 | id), control=glsControl(sigma=sqrt(mean((boostdf.nuisance$Y_minus_l_hat-boostdf.nuisance$xi_hat)^2))))
      rhoo <- as.numeric(coef(gls_fit$modelStruct, unconstrained=FALSE))
      beta_num = 0; beta_den = 0; V_num = 0
      theta = rhoo/(1-rhoo)
      s.beta <- rep(1,dim(boostdf.beta)[1])
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, s.beta, theta, proxyCCF="equicorr")
      beta_hat_k_hom_gls_num[k] <- cross_fit_evals$beta_num
      beta_hat_k_hom_gls_den[k] <- cross_fit_evals$beta_den
      V_num_hat_k_hom_gls[k] <- cross_fit_evals$V_num
      #V_s_hom_gls[k] <- gls_fit$varBeta[2,2] # Allows as a sanity check, one to check DML1 versus DML2 for this example


      # Homogeneous (hom) GEE Estimator
      gee_fit <- geeglm(Y_minus_l_hat ~ xi_hat, data=boostdf.nuisance, corstr="exchangeable", id=id)
      alpha_or_rhoo <- as.numeric(gee_fit$geese[2])
      theta <- alpha_or_rhoo/(1-alpha_or_rhoo)
      s.beta <- rep(1,dim(boostdf.beta)[1])
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, s.beta, theta, proxyCCF="equicorr")
      beta_hat_k_hom_marginal_gee_num[k] <- cross_fit_evals$beta_num
      beta_hat_k_hom_marginal_gee_den[k] <- cross_fit_evals$beta_den
      V_num_hat_k_hom_gee[k] <- cross_fit_evals$V_num

      # Unweighted Estimator
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, rep(1,nrow(boostdf.beta)), 0, proxyCCF="equicorr")
      beta_hat_num_k_unw[k] <- cross_fit_evals$beta_num
      beta_hat_den_k_unw[k] <- cross_fit_evals$beta_den
      V_num_hat_k_unw[k] <- cross_fit_evals$V_num

      # Oracle Estimator
      cross_fit_evals <- sandwich.boost:::cross_fit_eval(boostdf.beta, 1/f_0(as.matrix(boostdf.beta$D),as.matrix(boostdf.beta$X)), 0.2/(1-0.2), proxyCCF="equicorr")
      beta_hat_num_k_ora[k] <- cross_fit_evals$beta_num
      beta_hat_den_k_ora[k] <- cross_fit_evals$beta_den
      V_num_hat_k_ora[k] <- cross_fit_evals$V_num

    }

    Beta_boosting = sum(beta_hat_num_k_boosting)/sum(beta_hat_den_k_boosting)
    Beta_unw = sum(beta_hat_num_k_unw)/sum(beta_hat_den_k_unw)
    Beta_het_gee = sum(beta_hat_k_het_GEE_num_k)/sum(beta_hat_k_het_GEE_den_k)
    Beta_hom_gee = sum(beta_hat_k_hom_marginal_gee_num)/sum(beta_hat_k_hom_marginal_gee_den)
    Beta_het_lik = sum(beta_hat_k_het_gls_num)/sum(beta_hat_k_het_gls_den)
    Beta_hom_lik = sum(beta_hat_k_hom_gls_num)/sum(beta_hat_k_hom_gls_den)
    Beta_ora = sum(beta_hat_num_k_ora)/sum(beta_hat_den_k_ora)

    # Collect everthing together
    BETAS.DF <- rbind( BETAS.DF, data.frame(Beta_oracle=Beta_ora, Beta_boosting=Beta_boosting, Beta_het_lik=Beta_het_lik, Beta_hom_lik=Beta_hom_lik, Beta_het_gee=Beta_het_gee, Beta_hom_gee=Beta_hom_gee, Beta_unw=Beta_unw) )
    Variances.DF <- rbind(Variances.DF, data.frame(Oracle=I*N*sum(V_num_hat_k_ora)/(sum(beta_hat_den_k_ora))^2, Boosting=I*N*sum(V_num_hat_k_boosting)/(sum(beta_hat_den_k_boosting))^2 ,HetGEE=I*N*sum(V_num_hat_k_het_gee)/(sum(beta_hat_k_het_GEE_den_k))^2, HomGEE=I*N*sum(V_num_hat_k_hom_gee)/(sum(beta_hat_k_hom_marginal_gee_den))^2, Unweighted=I*N*sum(V_num_hat_k_unw)/(sum(beta_hat_den_k_unw))^2, homGLS=I*N*sum(V_num_hat_k_hom_gls)/(sum(beta_hat_k_hom_gls_den))^2, HetGLS=I*N*sum(V_num_hat_k_het_gls)/(sum(beta_hat_k_het_gls_den))^2 ))

    #    Uncomment to display progress
    #    cat("\n")
    #    cat('Progress: ',dim(BETAS.DF)[1],' of ',M,' complete')
  }
  Simulation.5.1.1.Betas[[which(lambda==Lambdas)]] = BETAS.DF
  Simulation.5.1.1.Variances[[which(lambda==Lambdas)]] =Variances.DF
}

