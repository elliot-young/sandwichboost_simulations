# Relevant Code for the real data anlaysis examples given in Section 5.2
#  5.2.1 - Orange Juice Price Elasticity
#  5.2.2 - National Longitudinal Survey of Young Working Women

# We load some relevant packages and functions (note many of the functions here
#  are in the package sandwich.boost, but to also define functions to fit mixed
#  effects models and GEEs in the DML framework we include all the relevant
#  for the simulations in Section 5.2 below).
{
  {
    library(mvtnorm)
    library(fBasics)
    library(magic)
    library(expm)
    library(jmuOutlier)
    library(lme4)
    library(caret)
    library(foreach)
    library(mgcv)
    library(Rcpp)
    library(geepack)

    regress.l <- function(l_formula, l_learner, data.nuisance, data.beta) {
      Y <- all.vars(l_formula)[1]; X <- all.vars(l_formula)[-1]
      if (l_learner == "gam") {
        nuisance_l_fit <- gam(l_formula, data=data.nuisance) #cubic splines
        l_residuals.nuisance <- data.nuisance[,Y] - nuisance_l_fit$fitted
        l_residuals.beta <- data.beta[,Y] - predict(nuisance_l_fit, data.beta[,X,drop=FALSE])
      } else if (l_learner == "lm") {#TODO
        nuisance_l_fit <- lm(l_formula, data.nuisance)#TO DO
        l_residuals.nuisance <- data.nuisance[,Y] - nuisance_l_fit$fitted#TO DO
        l_residuals.beta <- data.beta[,Y] - predict(nuisance_l_fit, data.beta[,X,drop=FALSE])
      } else if (l_learner == "randomforest") {
        stuff <- data.nuisance[, !(names(data.nuisance) %in% c("D", "id"))]
        l.task <- makeRegrTask(data=stuff, target="Y")
        res <- tuneRanger(l.task,measure=list(mse), iters = 70, iters.warmup = 30, time.budget = NULL, num.threads = NULL, num.trees = 500, parameters = list(replace = FALSE, respect.unordered.factors="order"), tune.parameters = c("min.node.size"))
        nuisance_l_fit <- ranger(l_formula, data.nuisance, min.node.size=res$recommended.pars$min.node.size, num.trees=500)
        l_residuals.nuisance <- data.nuisance[,Y] - predict(nuisance_l_fit, data.nuisance[,X,drop=FALSE])$predictions
        l_residuals.beta <- data.beta[,Y] - predict(nuisance_l_fit, data.beta[,X,drop=FALSE])$predictions
      } else if (l_learner == "lasso") {#TODO
        nuisance_l_fit <- cv.glmnet(data.nuisance[,X], data.nuisance[,Y])
        l_residuals.nuisance <- data.nuisance[,Y] - predict(nuisance_l_fit, newx = data.nuisance[,X], s = "lambda.min")
        l_residuals.beta <- data.nuisance[,Y] - predict(nuisance_l_fit, newx = data.beta[,X], s = "lambda.min")
      } else if (l_learner == "ridge") {#TODO
        nuisance_l_fit <- cv.glmnet(data.nuisance[,X], data.nuisance[,Y], alpha=0)
        l_residuals.nuisance <- data.nuisance[,Y] - predict(nuisance_l_fit, newx = data.nuisance[,X], s = "lambda.min")
        l_residuals.beta <- data.nuisance[,Y] - predict(nuisance_l_fit, newx = data.beta[,X], s = "lambda.min")
      } else {
        stop("Error: Invalid l_learner. l_learner must be gam, randomforest, lasso, ridge, or lm.")
      }
      return(list(l_residuals.nuisance=l_residuals.nuisance, l_residuals.beta=l_residuals.beta))
    }

    regress.m <- function(m_formula, m_learner, data.nuisance, data.beta) {
      D <- all.vars(m_formula)[1]
      X <- all.vars(m_formula)[-1]
      if (m_learner == "gam") {
        nuisance_m_fit <- gam(m_formula, data=data.nuisance) #cubic splines
        m_residuals.nuisance <- data.nuisance[,D] - nuisance_m_fit$fitted
        m_residuals.beta <- data.beta[,D] - predict(nuisance_m_fit, data.beta[,X,drop=FALSE])
      } else if (m_learner == "lm") {#TODO
        nuisance_m_fit <- lm(m_formula, data.nuisance)#TO DO
        m_residuals.nuisance <- data.nuisance[,D] - nuisance_m_fit$fitted#TO DO
        m_residuals.beta <- data.beta[,D] - predict(nuisance_m_fit, data.beta[,X,drop=FALSE])
      } else if (m_learner == "randomforest") {
        stuff <- data.nuisance[, !(names(data.nuisance) %in% c("Y", "id"))]
        m.task <- makeRegrTask(data=stuff, target="D")
        res <- tuneRanger(m.task,measure=list(mse), iters = 70, iters.warmup = 30, time.budget = NULL, num.threads = NULL, num.trees = 500, parameters = list(replace = FALSE, respect.unordered.factors="order"), tune.parameters = c("min.node.size"))
        nuisance_m_fit <- ranger(m_formula, data.nuisance, min.node.size=res$recommended.pars$min.node.size, num.trees=500)
        m_residuals.nuisance <- data.nuisance[,D] - predict(nuisance_m_fit, data.nuisance[,X,drop=FALSE])$predictions
        m_residuals.beta <- data.beta[,D] - predict(nuisance_m_fit, data.beta[,X,drop=FALSE])$predictions
      } else if (m_learner == "lasso") {#TODO
        nuisance_m_fit <- cv.glmnet(data.nuisance[,X], data.nuisance[,D])
        m_residuals.nuisance <- data.nuisance[,D] - predict(nuisance_m_fit, newx = data.nuisance[,X], s = "lambda.min")
        m_residuals.beta <- data.nuisance[,D] - predict(nuisance_m_fit, newx = data.beta[,X], s = "lambda.min")
      } else if (m_learner == "ridge") {#TODO
        nuisance_m_fit <- cv.glmnet(data.nuisance[,X], data.nuisance[,D], alpha=0)
        m_residuals.nuisance <- data.nuisance[,D] - predict(nuisance_m_fit, newx = data.nuisance[,D], s = "lambda.min")
        m_residuals.beta <- data.nuisance[,D] - predict(nuisance_m_fit, newx = data.beta[,D], s = "lambda.min")
      } else {
        stop("Error: Invalid m_learner. m_learner must be gam, randomforest, lasso, ridge, or lm.")
      }
      return(list(m_residuals.nuisance=m_residuals.nuisance, m_residuals.beta=m_residuals.beta))
    }

    ngradient = function(rdf, groups, I, n_obs, n_i, s, theta) {
      A1 <- rep(0,I)
      A2 <- rep(0,I)
      B1 <- rep(0,I)
      B2 <- rep(0,I)
      A3 <- numeric(0)
      A4 <- numeric(0)
      for (i in seq_len(I)) {
        s_i <- s[is.element(rdf$id,groups[i])]
        xi_hat_i <- rdf[is.element(rdf$id,groups[i]),"xi_hat"]
        epsilon_hat_i <- rdf[is.element(rdf$id,groups[i]),"epsilon_hat"]
        A1[i] <- sum(s_i^2*xi_hat_i^2)-theta/(1+theta*n_i[i])*(sum(s_i*xi_hat_i))^2
        A2[i] <- sum(s_i^2*xi_hat_i*epsilon_hat_i)-theta/(1+theta*n_i[i])*(sum(s_i*xi_hat_i))*(sum(s_i*epsilon_hat_i))
        B1[i] <- sum(s_i*xi_hat_i)
        B2[i] <- sum(s_i*epsilon_hat_i)
        for (j in seq_len(n_i[i])) {
          A3 <- c(A3, 2*s_i[j]*xi_hat_i[j]^2-2*theta/(1+theta*n_i[i])*B1[i]*xi_hat_i[j])
          A4 <- c(A4, 2*s_i[j]*xi_hat_i[j]*epsilon_hat_i[j]-theta/(1+theta*n_i[i])*(  xi_hat_i[j]*B2[i]   +   epsilon_hat_i[j]*B1[i]   ))
        }
      }

      s_scores <- numeric(0)
      sumA1=sum(A1)
      sumA2sq=sum(A2^2)
      cumsum_n_i <- 0
      for (i in seq_len(I)) {
        for (j in seq_len(n_i[i])) {
          s_scores <- c(s_scores, sumA1*A2[i]*A4[cumsum_n_i+j]-sumA2sq*A3[cumsum_n_i+j])
        }
        cumsum_n_i <- cumsum_n_i + n_i[i]
      }
      s_scores <- n_obs*s_scores*2/(sumA1)^3

      theta_score <- -n_obs*2/(sumA1)^3*( sumA1*sum(A2*B1*B2/(1+theta*n_i)^2) - sumA2sq*sum((B1/(1+theta*n_i))^2) )

      scores <- list(s=s_scores, theta=theta_score)

      return(scores)
    }

    risk = function(rdf, groups, I, n_obs, n_i, s, theta) {
      groups <- unique(rdf[,"id"])
      A1 <- rep(0,I)
      A2 <- rep(0,I)
      for (i in seq_len(I)) {
        s_i <- s[is.element(rdf$id,groups[i])]
        xi_hat_i <- rdf[is.element(rdf$id,groups[i]),"xi_hat"]
        epsilon_hat_i <- rdf[is.element(rdf$id,groups[i]),"epsilon_hat"]
        A1[i] <- sum(s_i^2*xi_hat_i^2)-theta/(1+theta*n_i[i])*(sum(s_i*xi_hat_i))^2
        A2[i] <- sum(s_i^2*xi_hat_i*epsilon_hat_i)-theta/(1+theta*n_i[i])*(sum(s_i*xi_hat_i))*(sum(s_i*epsilon_hat_i))
      }
      sumA1=sum(A1)
      sumA2sq=sum(A2^2)
      Risk = sumA2sq/(sumA1)^2
      return(Risk)
    }

    grouping.metadata <- function(boostdf){
      groups.metadata <- unique(boostdf[,"id"])
      I.metadata <- length(groups.metadata)
      n_i.metadata <- rep(0,I.metadata)
      for (i in seq_len(I.metadata)) {
        n_i.metadata[i] <- nrow(boostdf[boostdf$id==groups.metadata[i],])
      }
      n_obs.metadata <- sum(n_i.metadata)
      return(list(groups.metadata=groups.metadata, I.metadata=I.metadata, n_i.metadata=n_i.metadata, n_obs.metadata=n_obs.metadata))
    }

    regress.s <- function(s_formula, s_learner, boostdf.nuisance.train, boostdf.nuisance.test, boostdf.beta, u_s, U_S, X) {
      if (s_learner == "gam") {
        scores.fit <- gam(s_formula, data=cbind(boostdf.nuisance.train[,X,drop=FALSE],u_s)) #cubic splines
        u_s_hat.train <- scores.fit$fitted
        u_s_hat.test <- predict(scores.fit, boostdf.nuisance.test[,X,drop=FALSE])
        u_s_hat.beta <- predict(scores.fit, boostdf.beta[,X,drop=FALSE])
      } else if (s_learner == "lm") {#TODO
        scores.fit <- lm(s_formula, data=cbind(boostdf.nuisance.train[,X,drop=FALSE],u_s))
        u_s_hat.train <- scores.fit$fitted
        u_s_hat.test <- predict(scores.fit, boostdf.nuisance.test[,X,drop=FALSE])
        u_s_hat.beta <- predict(scores.fit, boostdf.beta[,X,drop=FALSE])
      } else if (s_learner == "randomforest") {
        scores.fit <- ranger(s_formula, data=cbind(boostdf.nuisance.train[,X,drop=FALSE],u_s), max.depth=1, num.trees=10)
        u_s_hat.train <- predict(scores.fit, boostdf.nuisance.train[,X,drop=FALSE])$predictions
        u_s_hat.test <- predict(scores.fit, boostdf.nuisance.test[,X,drop=FALSE])$predictions
        u_s_hat.beta <- predict(scores.fit, boostdf.beta[,X,drop=FALSE])$predictions
      } else {
        stop("Error: Invalid s_learner. s_learner must be gam, lm or randomforest")
      }
      return(list(u_s_hat.train=u_s_hat.train, u_s_hat.test=u_s_hat.test, u_s_hat.beta=u_s_hat.beta))
    }

    # Load the relevant C++ functions (to calculate scores)
    {
      cppFunction('List ngradient_cpp(NumericVector epsilon_hat, NumericVector xi_hat, int I, int n_obs, NumericVector n_i, NumericVector s, double theta) {
  NumericVector A1(I);
  NumericVector A2(I);
  NumericVector B1(I);
  NumericVector B2(I);
  NumericVector A3(n_obs);
  NumericVector A4(n_obs);
  int n_each = 0;
  int n_cumsum_start = 0;
  for(int i = 0; i < I; i++) {
    NumericVector s_i(n_i(i));
    NumericVector xi_hat_i(n_i(i));
    NumericVector epsilon_hat_i(n_i(i));
    for(int gg = 0; gg < n_i(i); gg++) {
      s_i(gg) = s(n_cumsum_start + gg);
      xi_hat_i(gg) = xi_hat(n_cumsum_start + gg);
      epsilon_hat_i(gg) = epsilon_hat(n_cumsum_start + gg);
    }
    A1(i) = sum(pow(s_i,2.0)*pow(xi_hat_i,2.0))-theta/(1+theta*n_i(i))*pow((sum(s_i*xi_hat_i)),2.0);
    A2(i) = sum(pow(s_i,2.0)*xi_hat_i*epsilon_hat_i)-theta/(1+theta*n_i(i))*(sum(s_i*xi_hat_i))*(sum(s_i*epsilon_hat_i));
    B1(i) = sum(s_i*xi_hat_i);
    B2(i) = sum(s_i*epsilon_hat_i);
    for(int j = 0; j < n_i(i); j++) {
      A3(n_each) = 2*s_i(j)*pow(xi_hat_i(j),2.0)-2*theta/(1+theta*n_i(i))*B1(i)*xi_hat_i(j);
      A4(n_each) = 2*s_i(j)*xi_hat_i(j)*epsilon_hat_i(j)-theta/(1+theta*n_i(i))*(  xi_hat_i(j)*B2(i)   +   epsilon_hat_i(j)*B1(i)   );
      n_each += 1;
    }
    n_cumsum_start += n_i(i);
  }

  NumericVector s_scores(n_obs);
  double sumA1 = sum(A1);
  double sumA2sq = sum(pow(A2,2.0));
  double scaling = n_obs*2/pow(sumA1,3.0);
  n_each = 0;
  for(int i = 0; i < I; i++) {
    for(int j = 0; j < n_i(i); j++) {
      s_scores(n_each) = scaling*( sumA1*A2(i)*A4(n_each) - sumA2sq*A3(n_each) );
      n_each += 1;
    }
  }

  double theta_score = -scaling*( sumA1*sum(A2*B1*B2/pow((1+theta*n_i),2.0)) - sumA2sq*sum(pow((B1/(1+theta*n_i)),2.0)) );

  List scores;
  scores["s_scores"] = s_scores;
  scores["theta_score"] = theta_score;
  return scores;
}')
    }
    {
      cppFunction('double risk_cpp(NumericVector epsilon_hat, NumericVector xi_hat, int I, int n_obs, NumericVector n_i, NumericVector s, double theta) {
  NumericVector A1(I);
  NumericVector A2(I);
  int n_cumsum_start = 0;
  for(int i = 0; i < I; i++) {
    NumericVector s_i(n_i(i));
    NumericVector xi_hat_i(n_i(i));
    NumericVector epsilon_hat_i(n_i(i));
    for(int gg = 0; gg < n_i(i); gg++) {
      s_i(gg) = s(n_cumsum_start + gg);
      xi_hat_i(gg) = xi_hat(n_cumsum_start + gg);
      epsilon_hat_i(gg) = epsilon_hat(n_cumsum_start + gg);
    }
    A1(i) = sum(pow(s_i,2.0)*pow(xi_hat_i,2.0))-theta/(1+theta*n_i(i))*pow((sum(s_i*xi_hat_i)),2.0);
    A2(i) = sum(pow(s_i,2.0)*xi_hat_i*epsilon_hat_i)-theta/(1+theta*n_i(i))*(sum(s_i*xi_hat_i))*(sum(s_i*epsilon_hat_i));
    n_cumsum_start += n_i(i);
  }

  NumericVector s_scores(n_obs);
  double sumA1 = sum(A1);
  double sumA2sq = sum(pow(A2,2.0));

  double risk = n_obs*sumA2sq/pow(sumA1,2.0);

  return risk;
}')
    }
    {
      cppFunction('List ngradient_equicorr_cpp(NumericVector epsilon_hat, NumericVector xi_hat, int I, int n_obs, NumericVector n_i, NumericVector s, double theta) {
  NumericVector A1(I);
  NumericVector A2(I);
  NumericVector B1(I);
  NumericVector B2(I);
  NumericVector A3(n_obs);
  NumericVector A4(n_obs);
  int n_each = 0;
  int n_cumsum_start = 0;
  for(int i = 0; i < I; i++) {
    NumericVector s_i(n_i(i));
    NumericVector xi_hat_i(n_i(i));
    NumericVector epsilon_hat_i(n_i(i));
    for(int gg = 0; gg < n_i(i); gg++) {
      s_i(gg) = s(n_cumsum_start + gg);
      xi_hat_i(gg) = xi_hat(n_cumsum_start + gg);
      epsilon_hat_i(gg) = epsilon_hat(n_cumsum_start + gg);
    }
    A1(i) = sum(pow(s_i,2.0)*pow(xi_hat_i,2.0))-theta/(1+theta*n_i(i))*pow((sum(s_i*xi_hat_i)),2.0);
    A2(i) = sum(pow(s_i,2.0)*xi_hat_i*epsilon_hat_i)-theta/(1+theta*n_i(i))*(sum(s_i*xi_hat_i))*(sum(s_i*epsilon_hat_i));
    B1(i) = sum(s_i*xi_hat_i);
    B2(i) = sum(s_i*epsilon_hat_i);
    for(int j = 0; j < n_i(i); j++) {
      A3(n_each) = 2*s_i(j)*pow(xi_hat_i(j),2.0)-2*theta/(1+theta*n_i(i))*B1(i)*xi_hat_i(j);
      A4(n_each) = 2*s_i(j)*xi_hat_i(j)*epsilon_hat_i(j)-theta/(1+theta*n_i(i))*(  xi_hat_i(j)*B2(i)   +   epsilon_hat_i(j)*B1(i)   );
      n_each += 1;
    }
    n_cumsum_start += n_i(i);
  }

  NumericVector s_scores(n_obs);
  double sumA1 = sum(A1);
  double sumA2sq = sum(pow(A2,2.0));
  double scaling = n_obs*2/pow(sumA1,3.0);
  n_each = 0;
  for(int i = 0; i < I; i++) {
    for(int j = 0; j < n_i(i); j++) {
      s_scores(n_each) = scaling*( sumA1*A2(i)*A4(n_each) - sumA2sq*A3(n_each) );
      n_each += 1;
    }
  }

  double theta_score = -scaling*( sumA1*sum(A2*B1*B2/pow((1+theta*n_i),2.0)) - sumA2sq*sum(pow((B1/(1+theta*n_i)),2.0)) );

  List scores;
  scores["s_scores"] = s_scores;
  scores["theta_score"] = theta_score;
  return scores;
}')

      cppFunction('double risk_equicorr_cpp(NumericVector epsilon_hat, NumericVector xi_hat, int I, int n_obs, NumericVector n_i, NumericVector s, double theta) {
  NumericVector A1(I);
  NumericVector A2(I);
  int n_cumsum_start = 0;
  for(int i = 0; i < I; i++) {
    NumericVector s_i(n_i(i));
    NumericVector xi_hat_i(n_i(i));
    NumericVector epsilon_hat_i(n_i(i));
    for(int gg = 0; gg < n_i(i); gg++) {
      s_i(gg) = s(n_cumsum_start + gg);
      xi_hat_i(gg) = xi_hat(n_cumsum_start + gg);
      epsilon_hat_i(gg) = epsilon_hat(n_cumsum_start + gg);
    }
    A1(i) = sum(pow(s_i,2.0)*pow(xi_hat_i,2.0))-theta/(1+theta*n_i(i))*pow((sum(s_i*xi_hat_i)),2.0);
    A2(i) = sum(pow(s_i,2.0)*xi_hat_i*epsilon_hat_i)-theta/(1+theta*n_i(i))*(sum(s_i*xi_hat_i))*(sum(s_i*epsilon_hat_i));
    n_cumsum_start += n_i(i);
  }

  NumericVector s_scores(n_obs);
  double sumA1 = sum(A1);
  double sumA2sq = sum(pow(A2,2.0));

  double risk = n_obs*sumA2sq/pow(sumA1,2.0);

  return risk;
}')


      ##### C++ Scores for AR(1) Case #####

      cppFunction('List ngradient_AR_cpp(NumericVector epsilon_hat, NumericVector xi_hat, int I, int n_obs, NumericVector n_i, NumericVector s, double theta) {
  NumericVector A1(I);
  NumericVector A2(I);
  NumericVector alpha1(I);
  NumericVector beta1(I);
  NumericVector gamma1(I);
  NumericVector alpha2(I);
  NumericVector beta2(I);
  NumericVector gamma2(I);
  NumericVector A3(n_obs);
  NumericVector A4(n_obs);
  int n_each = 0;
  int n_cumsum_start = 0;
  for(int i = 0; i < I; i++) {
    NumericVector s_i(n_i(i));
    NumericVector xi_hat_i(n_i(i));
    NumericVector epsilon_hat_i(n_i(i));
    for(int gg = 0; gg < n_i(i); gg++) {
      s_i(gg) = s(n_cumsum_start + gg);
      xi_hat_i(gg) = xi_hat(n_cumsum_start + gg);
      epsilon_hat_i(gg) = epsilon_hat(n_cumsum_start + gg);
    }
        if (n_i(i) > 2){
        NumericVector s_i_sub1(n_i(i)-1);
        NumericVector s_i_sub2(n_i(i)-1);
        NumericVector xi_hat_i_sub1(n_i(i)-1);
        NumericVector xi_hat_i_sub2(n_i(i)-1);
        NumericVector epsilon_hat_i_sub1(n_i(i)-1);
        NumericVector epsilon_hat_i_sub2(n_i(i)-1);
        for(int gg = 0; gg < (n_i(i)-1); gg++) {
          s_i_sub1(gg) = s_i(gg);
          s_i_sub2(gg) = s_i(gg+1);
          xi_hat_i_sub1(gg) = xi_hat_i(gg);
          xi_hat_i_sub2(gg) = xi_hat_i(gg+1);
          epsilon_hat_i_sub1(gg) = epsilon_hat_i(gg);
          epsilon_hat_i_sub2(gg) = epsilon_hat_i(gg+1);
        }
      alpha1(i) = sum(pow(s_i*xi_hat_i,2.0));
      beta1(i) = alpha1(i) - pow(s_i(0)*xi_hat_i(0),2.0) - pow(s_i(n_i(i)-1)*xi_hat_i(n_i(i)-1),2.0);
      gamma1(i) = 2*sum(s_i_sub1*s_i_sub2*xi_hat_i_sub1*xi_hat_i_sub2);
      alpha2(i) = sum(pow(s_i,2.0)*xi_hat_i*epsilon_hat_i);
      beta2(i) = alpha2(i) - pow(s_i(0),2.0)*xi_hat_i(0)*epsilon_hat_i(0) - pow(s_i(n_i(i)-1),2.0)*xi_hat_i(n_i(i)-1)*epsilon_hat_i(n_i(i)-1);
      gamma2(i) = sum(s_i_sub1*s_i_sub2*(xi_hat_i_sub2*epsilon_hat_i_sub1+xi_hat_i_sub1*epsilon_hat_i_sub2));
      A1(i) = alpha1(i) + pow(theta,2.0)*beta1(i) - theta*gamma1(i);
      A2(i) = alpha2(i) + pow(theta,2.0)*beta2(i) - theta*gamma2(i);
        A3(n_each) = 2*( pow(xi_hat_i(0),2.0)*s_i(0) - theta*xi_hat_i(1)*xi_hat_i(0)*s_i(1) );
        A4(n_each) = 2*epsilon_hat_i(0)*xi_hat_i(0)*s_i(0) - theta*( epsilon_hat_i(1)*xi_hat_i(0)*s_i(1) + epsilon_hat_i(0)*xi_hat_i(1)*s_i(1) );
        n_each += 1;
      for(int j = 1; j < (n_i(i)-1); j++) {
        A3(n_each) = 2*((1+pow(theta,2.0))*pow(xi_hat_i(j),2.0)*s_i(j)-theta*(xi_hat_i(j+1)*xi_hat_i(j)*s_i(j+1)+xi_hat_i(j)*xi_hat_i(j-1)*s_i(j-1)));
        A4(n_each) = 2*(1+pow(theta,2.0))*epsilon_hat_i(j)*xi_hat_i(j)*s_i(j)-theta*(epsilon_hat_i(j+1)*xi_hat_i(j)*s_i(j+1) + epsilon_hat_i(j)*xi_hat_i(j-1)*s_i(j-1)) + epsilon_hat_i(j-1)*xi_hat_i(j)*s_i(j-1) + epsilon_hat_i(j)*xi_hat_i(j+1)*s_i(j+1);
      n_each += 1;
      }
        A3(n_each) = 2*(pow(xi_hat_i(n_i(i)-1),2.0)*s_i(n_i(i)-1) - theta*(xi_hat_i(n_i(i)-1)*xi_hat_i(n_i(i)-2)*s_i(n_i(i)-2)) );
        A4(n_each) = 2*epsilon_hat_i(n_i(i)-1)*xi_hat_i(n_i(i)-1)*s_i(n_i(i)-1) - theta*(epsilon_hat_i(n_i(i)-1)*xi_hat_i(n_i(i)-2)*s_i(n_i(i)-2) + epsilon_hat_i(n_i(i)-2)*xi_hat_i(n_i(i)-1)*s_i(n_i(i)-1) );
        n_each += 1;
    }
    else if (n_i(i) == 2){
      alpha1(i) = sum(pow(s_i*xi_hat_i,2.0));
      gamma1(i) = 2*(s_i(1)*s_i(0)*xi_hat_i(1)*xi_hat_i(0));
      alpha2(i) = sum(pow(s_i,2.0)*xi_hat_i*epsilon_hat_i);
      gamma2(i) = s_i(1)*s_i(0)*(xi_hat_i(0)*epsilon_hat_i(1)+xi_hat_i(1)*epsilon_hat_i(0));
      A1(i) = alpha1(i) - theta*gamma1(i);
      A2(i) = alpha2(i) - theta*gamma2(i);
        A3(n_each) = 2*( pow(xi_hat_i(0),2.0)*s_i(0) - theta*xi_hat_i(1)*xi_hat_i(0)*s_i(1) ) ;
        A4(n_each) = 2*epsilon_hat_i(0)*xi_hat_i(0)*s_i(0) - theta*( epsilon_hat_i(1)*xi_hat_i(0)*s_i(1) + epsilon_hat_i(0)*xi_hat_i(1)*s_i(1) );
        n_each += 1;
        A3(n_each) = 2*( pow(xi_hat_i(1),2.0)*s_i(1) - theta*xi_hat_i(1)*xi_hat_i(0)*s_i(0) );
        A4(n_each) = 2*epsilon_hat_i(1)*xi_hat_i(1)*s_i(1) - theta*( epsilon_hat_i(1)*xi_hat_i(0)*s_i(0) + epsilon_hat_i(0)*xi_hat_i(1)*s_i(0) );
        n_each += 1;
    }
    else if (n_i(i) == 1){
      alpha1(i) = 2*pow(s_i(0)*xi_hat_i(0),2.0);
      alpha2(i) = 2*pow(s_i(0),2.0)*xi_hat_i(0)*epsilon_hat_i(0);
      A1(i) = (1-pow(theta,2.0))*alpha1(i);
      A2(i) = (1-pow(theta,2.0))*alpha2(i);
      A3(n_each) = 2*(1-pow(theta,2.0))*pow(xi_hat_i(0),2.0)*s_i(0);
      A4(n_each) = 2*(1-pow(theta,2.0))*epsilon_hat_i(0)*xi_hat_i(0)*s_i(0);
      n_each += 1;
    }
    n_cumsum_start += n_i(i);
  }

  NumericVector s_scores(n_obs);
  double sumA1 = sum(A1);
  double sumA2sq = sum(pow(A2,2.0));
  double scaling = n_obs*2/pow(sumA1,3.0);
  n_each = 0;
  for(int i = 0; i < I; i++) {
    for(int j = 0; j < n_i(i); j++) {
      s_scores(n_each) = scaling*(sumA1*A2(i)*A4(n_each)-sumA2sq*A3(n_each));
      n_each += 1;
    }
  }

  double theta_score = scaling*( sumA1*sum(A2*(2*theta*beta2-gamma2)) - sumA2sq*sum(2*theta*beta1-gamma1) );

  List scores;
  scores["s_scores"] = s_scores;
  scores["theta_score"] = theta_score;
  return scores;
}')

      cppFunction('double risk_AR_cpp(NumericVector epsilon_hat, NumericVector xi_hat, int I, int n_obs, NumericVector n_i, NumericVector s, double theta) {
  NumericVector A1(I);
  NumericVector A2(I);
  NumericVector alpha1(I);
  NumericVector beta1(I);
  NumericVector gamma1(I);
  NumericVector alpha2(I);
  NumericVector beta2(I);
  NumericVector gamma2(I);
  int n_cumsum_start = 0;
  for(int i = 0; i < I; i++) {
    NumericVector s_i(n_i(i));
    NumericVector xi_hat_i(n_i(i));
    NumericVector epsilon_hat_i(n_i(i));
    for(int gg = 0; gg < n_i(i); gg++) {
      s_i(gg) = s(n_cumsum_start + gg);
      xi_hat_i(gg) = xi_hat(n_cumsum_start + gg);
      epsilon_hat_i(gg) = epsilon_hat(n_cumsum_start + gg);
    }
    if (n_i(i) > 2){
        NumericVector s_i_sub1(n_i(i)-1);
        NumericVector s_i_sub2(n_i(i)-1);
        NumericVector xi_hat_i_sub1(n_i(i)-1);
        NumericVector xi_hat_i_sub2(n_i(i)-1);
        NumericVector epsilon_hat_i_sub1(n_i(i)-1);
        NumericVector epsilon_hat_i_sub2(n_i(i)-1);
        for(int gg = 0; gg < (n_i(i)-1); gg++) {
          s_i_sub1(gg) = s_i(gg);
          s_i_sub2(gg) = s_i(gg+1);
          xi_hat_i_sub1(gg) = xi_hat_i(gg);
          xi_hat_i_sub2(gg) = xi_hat_i(gg+1);
          epsilon_hat_i_sub1(gg) = epsilon_hat_i(gg);
          epsilon_hat_i_sub2(gg) = epsilon_hat_i(gg+1);
        }
      alpha1(i) = sum(pow(s_i,2.0)*pow(xi_hat_i,2.0));
      beta1(i) = alpha1(i) - pow(s_i(0)*xi_hat_i(0),2.0) - pow(s_i(n_i(i)-1)*xi_hat_i(n_i(i)-1),2.0);
      gamma1(i) = 2*sum(s_i_sub1*s_i_sub2*xi_hat_i_sub1*xi_hat_i_sub2);
      alpha2(i) = sum(pow(s_i,2.0)*xi_hat_i*epsilon_hat_i);
      beta2(i) = alpha2(i) - pow(s_i(0),2.0)*xi_hat_i(0)*epsilon_hat_i(0) - pow(s_i(n_i(i)-1),2.0)*xi_hat_i(n_i(i)-1)*epsilon_hat_i(n_i(i)-1);
      gamma2(i) = sum(s_i_sub1*s_i_sub2*(xi_hat_i_sub2*epsilon_hat_i_sub1+xi_hat_i_sub1*epsilon_hat_i_sub2));
      A1(i) = alpha1(i) + pow(theta,2.0)*beta1(i) - theta*gamma1(i);
      A2(i) = alpha2(i) + pow(theta,2.0)*beta2(i) - theta*gamma2(i);
    }
    else if (n_i(i) == 2){
      alpha1(i) = sum(pow(s_i,2.0)*pow(xi_hat_i,2.0));
      gamma1(i) = 2*(s_i(1)*s_i(0)*xi_hat_i(1)*xi_hat_i(0));
      alpha2(i) = sum(pow(s_i,2.0)*xi_hat_i*epsilon_hat_i);
      gamma2(i) = s_i(1)*s_i(0)*(xi_hat_i(0)*epsilon_hat_i(1)+xi_hat_i(1)*epsilon_hat_i(0));
      A1(i) = alpha1(i) - theta*gamma1(i);
      A2(i) = alpha2(i) - theta*gamma2(i);
    }
    else if (n_i(i) == 1){
      alpha1(i) = 2*pow(s_i(0),2.0)*pow(xi_hat_i(0),2.0);
      alpha2(i) = 2*pow(s_i(0),2.0)*xi_hat_i(0)*epsilon_hat_i(0);
      A1(i) = (1-pow(theta,2.0))*alpha1(i);
      A2(i) = (1-pow(theta,2.0))*alpha2(i);
    }
    n_cumsum_start += n_i(i);
  }

  NumericVector s_scores(n_obs);
  double sumA1 = sum(A1);
  double sumA2sq = sum(pow(A2,2.0));

  double risk = n_obs*sumA2sq/pow(sumA1,2.0);

  return risk;
}')

      cppFunction('List ngradient_nested_cpp(NumericVector epsilon_hat, NumericVector xi_hat, NumericVector J, int n_obs, NumericMatrix n_ij, NumericVector s, NumericVector theta) {
  int I = J.length();
  NumericVector phi(I);
  NumericVector phidiff(I);
  NumericVector B1_sum1(I);
  NumericVector B1_sum2(I);
  NumericVector B1_sum3(I);
  NumericVector B1_sum4(I);
  NumericVector B2_sum1(I);
  NumericVector B2_sum2(I);
  NumericVector B2_sum3(I);
  NumericVector B2_sum4(I);
  NumericVector A1(I);
  NumericVector A2(I);
  NumericVector C1(I);
  NumericVector C2(I);
  NumericVector A3(n_obs);
  NumericVector A4(n_obs);
  int n_each = 0;
  int n_cumsum_start = 0;
  for(int i = 0; i < I; i++) {
    NumericVector b1(J(i));
    NumericVector b2(J(i));
    NumericVector c1(J(i));
    NumericVector c2(J(i));
    NumericVector N_ij(J(i));
    for(int gg = 0; gg < J(i); gg++) {
      N_ij(gg) = n_ij(i,gg);
    }
    for (int j = 0; j < J(i); j++) {
      NumericVector s_i(N_ij(j));
      NumericVector xi_hat_i(N_ij(j));
      NumericVector epsilon_hat_i(N_ij(j));
      for(int gg = 0; gg < N_ij(j); gg++) {
        s_i(gg) = s(n_cumsum_start + gg);
        xi_hat_i(gg) = xi_hat(n_cumsum_start + gg);
        epsilon_hat_i(gg) = epsilon_hat(n_cumsum_start + gg);
      }
      b1(j) = sum(s_i*xi_hat_i);
      b2(j) = sum(s_i*epsilon_hat_i);
      c1(j) = sum(pow(s_i*xi_hat_i,2.0));
      c2(j) = sum(pow(s_i,2.0)*epsilon_hat_i*xi_hat_i);
    }
    phi(i) = sum(N_ij/(1+theta(0)*N_ij));
    phidiff(i) = -sum(pow(N_ij/(1+theta(0)*N_ij),2.0));
    B1_sum1(i) = sum(b1/(1+theta(0)*N_ij));
    B1_sum2(i) = sum(pow(b1/(1+theta(0)*N_ij),2.0));
    B1_sum3(i) = sum(N_ij/pow(1+theta(0)*N_ij,2.0)*b1);
    B1_sum4(i) = sum(N_ij/(1+theta(0)*N_ij)*pow(b1,2.0));
    B2_sum1(i) = sum(b2/(1+theta(0)*N_ij));
    B2_sum2(i) = sum(b1*b2/pow(1+theta(0)*N_ij,2.0));
    B2_sum3(i) = sum(N_ij/pow(1+theta(0)*N_ij,2.0)*b2);
    B2_sum4(i) = sum(N_ij/(1+theta(0)*N_ij)*b1*b2);
    C1(i) = sum(c1);
    C2(i) = sum(c2);
    A1(i) = (1+theta(0)+theta(1))*C1(i) - B1_sum4(i) - theta(1)/(1+theta(1)*phi(i))*pow(B1_sum1(i),2.0);
    A2(i) = (1+theta(0)+theta(1))*C2(i) - B2_sum4(i) - theta(1)/(1+theta(1)*phi(i))*B1_sum1(i)*B2_sum1(i);

    for (int j=0; j < J(i); j++) {
      NumericVector s_i(N_ij(j));
      NumericVector xi_hat_i(N_ij(j));
      NumericVector epsilon_hat_i(N_ij(j));
      for(int gg = 0; gg < N_ij(j); gg++) {
        s_i(gg) = s(n_cumsum_start + gg);
        xi_hat_i(gg) = xi_hat(n_cumsum_start + gg);
        epsilon_hat_i(gg) = epsilon_hat(n_cumsum_start + gg);
      }
      n_cumsum_start += N_ij(j);
      for (int k=0; k < N_ij(j); k++) {
        A3(n_each) = 2*(1+theta(0)+theta(1))*s_i(k)*pow(xi_hat_i(k),2.0) - 2*theta(0)/(1+theta(0)*N_ij(j))*b1(j)*xi_hat_i(k) - 2*theta(1)/(1+theta(1)*phi(i))*B1_sum1(i)*xi_hat_i(k)/(1+theta(0)*N_ij(j));
        A4(n_each) = 2*(1+theta(0)+theta(1))*s_i(k)*xi_hat_i(k)*epsilon_hat_i(k) - theta(0)/(1+theta(0)*N_ij(j))*( b1(j)*epsilon_hat_i(k) + b2(j)*xi_hat_i(k) ) - theta(1)/(1+theta(1)*phi(i))*( B2_sum1(i)*xi_hat_i(k) + B1_sum1(i)*epsilon_hat_i(k) )/(1+theta(0)*N_ij(j));
        n_each += 1;
      }
    }
  }

  NumericVector s_scores(n_obs);
  double sumA1 = sum(A1);
  double sumA2sq = sum(pow(A2,2.0));
  double scaling = n_obs*2/pow(sumA1,3.0);
  n_each = 0;
  for(int i = 0; i < I; i++) {
    for(int j = 0; j < J(i); j++) {
      for (int k=0; k < n_ij(i,j); k++) {
        s_scores(n_each) = scaling*( sumA1*A2(i)*A4(n_each) - sumA2sq*A3(n_each) );
        n_each += 1;
      }
    }
  }

  NumericVector theta_score(2);
  theta_score(0) = scaling*( sumA1*sum(A2*(C2-B2_sum2+B1_sum1*B2_sum1*pow(theta(1)/(1+theta(1)*phi),2.0)*phidiff+B1_sum1*B2_sum3*theta(1)/(1+theta(1)*phi)+B2_sum1*B1_sum3)) - sumA2sq*sum(C1-B1_sum2+pow(theta(1)/(1+theta(1)*phi),2.0)*phidiff*pow(B1_sum1,2.0)+2*theta(1)/(1+theta(1)*phi)*B1_sum1*B1_sum3) );
  theta_score(1) = scaling*( sumA1*sum(A2*(C2-B2_sum1*B1_sum1/pow(1+theta(1)*phi,2.0))) - sumA2sq*sum(C1-B1_sum1/pow(1+theta(1)*phi,2.0)) );

  List scores;
  scores["s_scores"] = s_scores;
  scores["theta_score"] = theta_score;
  return scores;

}')

      cppFunction('double risk_nested_cpp(NumericVector epsilon_hat, NumericVector xi_hat, NumericVector J, int n_obs, NumericMatrix n_ij, NumericVector s, NumericVector theta) {
  int I = J.length();
  NumericVector phi(I);
  NumericVector phidiff(I);
  NumericVector B1_sum1(I);
  NumericVector B1_sum4(I);
  NumericVector B2_sum1(I);
  NumericVector B2_sum4(I);
  NumericVector A1(I);
  NumericVector A2(I);
  NumericVector C1(I);
  NumericVector C2(I);
  NumericVector A3(n_obs);
  NumericVector A4(n_obs);
  int n_cumsum_start = 0;
  for(int i = 0; i < I; i++) {
    NumericVector b1(J(i));
    NumericVector b2(J(i));
    NumericVector c1(J(i));
    NumericVector c2(J(i));
    NumericVector N_ij(J(i));
    for(int gg = 0; gg < J(i); gg++) {
      N_ij(gg) = n_ij(i,gg);
    }
    for (int j = 0; j < J(i); j++) {
      NumericVector s_i(N_ij(j));
      NumericVector xi_hat_i(N_ij(j));
      NumericVector epsilon_hat_i(N_ij(j));
      for(int gg = 0; gg < N_ij(j); gg++) {
        s_i(gg) = s(n_cumsum_start + gg);
        xi_hat_i(gg) = xi_hat(n_cumsum_start + gg);
        epsilon_hat_i(gg) = epsilon_hat(n_cumsum_start + gg);
      }
      n_cumsum_start += N_ij(j);
      b1(j) = sum(s_i*xi_hat_i);
      b2(j) = sum(s_i*epsilon_hat_i);
      c1(j) = sum(pow(s_i*xi_hat_i,2.0));
      c2(j) = sum(pow(s_i,2.0)*epsilon_hat_i*xi_hat_i);
    }
    phi(i) = sum(N_ij/(1+theta(0)*N_ij));
    phidiff(i) = -sum(pow(N_ij/(1+theta(0)*N_ij),2.0));
    B1_sum1(i) = sum(b1/(1+theta(0)*N_ij));
    B1_sum4(i) = sum(N_ij/(1+theta(0)*N_ij)*pow(b1,2.0));
    B2_sum1(i) = sum(b2/(1+theta(0)*N_ij));
    B2_sum4(i) = sum(N_ij/(1+theta(0)*N_ij)*b1*b2);
    C1(i) = sum(c1);
    C2(i) = sum(c2);
    A1(i) = (1+theta(0)+theta(1))*C1(i) - B1_sum4(i) - theta(1)/(1+theta(1)*phi(i))*pow(B1_sum1(i),2.0);
    A2(i) = (1+theta(0)+theta(1))*C2(i) - B2_sum4(i) - theta(1)/(1+theta(1)*phi(i))*B1_sum1(i)*B2_sum1(i);
  }

  double sumA1 = sum(A1);
  double sumA2sq = sum(pow(A2,2.0));

  double risk = n_obs*sumA2sq/pow(sumA1,2.0);

  return risk;
}')

    }
    {
      cppFunction('double optimal_step_size_cpp(NumericVector epsilon_hat, NumericVector xi_hat, int I, int n_obs, NumericVector n_i, NumericVector s, NumericVector h, double theta) {
  NumericVector A1(I);
  NumericVector A2(I);
  NumericVector Ab(I);
  NumericVector Ac(I);
  NumericVector Bb(I);
  NumericVector Bc(I);
  int n_cumsum_start = 0;
  for(int i = 0; i < I; i++) {
    NumericVector s_i(n_i(i));
    NumericVector h_i(n_i(i));
    NumericVector xi_hat_i(n_i(i));
    NumericVector epsilon_hat_i(n_i(i));
    for(int gg = 0; gg < n_i(i); gg++) {
      s_i(gg) = s(n_cumsum_start + gg);
      h_i(gg) = h(n_cumsum_start + gg);
      xi_hat_i(gg) = xi_hat(n_cumsum_start + gg);
      epsilon_hat_i(gg) = epsilon_hat(n_cumsum_start + gg);
    }
    A1(i) = sum(pow(s_i,2.0)*pow(xi_hat_i,2.0))-theta/(1+theta*n_i(i))*pow((sum(s_i*xi_hat_i)),2.0);
    A2(i) = sum(pow(s_i,2.0)*xi_hat_i*epsilon_hat_i)-theta/(1+theta*n_i(i))*(sum(s_i*xi_hat_i))*(sum(s_i*epsilon_hat_i));
    Ab(i) = 2*sum(xi_hat_i*epsilon_hat_i*h_i*s_i)-theta/(1+theta*n_i(i))*(sum(s_i*xi_hat_i)*sum(h_i*epsilon_hat_i)+sum(h_i*xi_hat_i)*sum(s_i*epsilon_hat_i));
    Ac(i) = sum(pow(xi_hat_i,2.0)*pow(h_i,2.0))-theta/(1+theta*n_i(i))*pow((sum(h_i*xi_hat_i)),2.0);
    Bb(i) = sum(pow(xi_hat_i,2.0)*h_i*s_i)-theta/(1+theta*n_i(i))*(sum(h_i*xi_hat_i))*(sum(s_i*xi_hat_i));
    Bc(i) = sum(pow(xi_hat_i,2.0)*pow(h_i,2.0))-theta/(1+theta*n_i(i))*pow(sum(h_i*xi_hat_i),2.0);
    n_cumsum_start += n_i(i);
  }

  double a1 = sum(pow(A2,2.0));
  double a2 = -2*sum(A2*Ab);
  double a4 = sum(pow(Ab,2.0));
  double b1 = sum(A1);
  double b2 = -2*sum(Bb);
  double b4 = sum(Bc);

  double eta_optimal = 0.5*(2*a1*b2-a2*b1)/(a4*b1 - 2*a2*b2 - 2*a1*b4 + 3*(a1*pow(b2,2.0)/b1));

  return eta_optimal;
}')


      # Corrected Equicorr Optimal Step Size

      cppFunction('double optimal_step_size_cpp(NumericVector epsilon_hat, NumericVector xi_hat, int I, int n_obs, NumericVector n_i, NumericVector s, NumericVector h, double theta) {
  NumericVector A1(I);
  NumericVector A2(I);
  NumericVector Ab(I);
  NumericVector Ac(I);
  NumericVector Bb(I);
  NumericVector Bc(I);
  int n_cumsum_start = 0;
  for(int i = 0; i < I; i++) {
    NumericVector s_i(n_i(i));
    NumericVector h_i(n_i(i));
    NumericVector xi_hat_i(n_i(i));
    NumericVector epsilon_hat_i(n_i(i));
    for(int gg = 0; gg < n_i(i); gg++) {
      s_i(gg) = s(n_cumsum_start + gg);
      h_i(gg) = h(n_cumsum_start + gg);
      xi_hat_i(gg) = xi_hat(n_cumsum_start + gg);
      epsilon_hat_i(gg) = epsilon_hat(n_cumsum_start + gg);
    }
    A1(i) = sum(pow(s_i,2.0)*pow(xi_hat_i,2.0))-theta/(1+theta*n_i(i))*pow((sum(s_i*xi_hat_i)),2.0);
    A2(i) = sum(pow(s_i,2.0)*xi_hat_i*epsilon_hat_i)-theta/(1+theta*n_i(i))*(sum(s_i*xi_hat_i))*(sum(s_i*epsilon_hat_i));
    Ab(i) = 2*sum(xi_hat_i*epsilon_hat_i*h_i*s_i)-theta/(1+theta*n_i(i))*(sum(s_i*xi_hat_i)*sum(h_i*epsilon_hat_i)+sum(h_i*xi_hat_i)*sum(s_i*epsilon_hat_i));
    Ac(i) = sum(xi_hat_i*epsilon_hat_i*pow(h_i,2.0))-theta/(1+theta*n_i(i))*(sum(h_i*xi_hat_i))*(sum(h_i*epsilon_hat_i));
    Bb(i) = sum(pow(xi_hat_i,2.0)*h_i*s_i)-theta/(1+theta*n_i(i))*(sum(h_i*xi_hat_i))*(sum(s_i*xi_hat_i));
    Bc(i) = sum(pow(xi_hat_i,2.0)*pow(h_i,2.0))-theta/(1+theta*n_i(i))*pow(sum(h_i*xi_hat_i),2.0);
    n_cumsum_start += n_i(i);
  }

  double a1 = sum(pow(A2,2.0));
  double a2 = -2*sum(A2*Ab);
  double a4 = sum(pow(Ab,2.0));
  double b1 = sum(A1);
  double b2 = -2*sum(Bb);
  double b4 = sum(Bc);

  double eta_optimal = 0.5*(2*a1*b2-a2*b1)/(a4*b1 - 2*a2*b2 - 2*a1*b4 + 3*(a1*pow(b2,2.0)/b1));

  return eta_optimal;
}')


      # Optimal Step Size for AR(1) proxyCCF
      cppFunction('double optimal_step_size_AR_cpp(NumericVector epsilon_hat, NumericVector xi_hat, int I, int n_obs, NumericVector n_i, NumericVector s, NumericVector h, double theta) {
  NumericVector A1(I);
  NumericVector A2(I);
  NumericVector Ab(I);
  NumericVector Ac(I);
  NumericVector Bb(I);
  NumericVector Bc(I);
  NumericVector alpha1(I);
  NumericVector beta1(I);
  NumericVector gamma1(I);
  NumericVector alpha2(I);
  NumericVector beta2(I);
  NumericVector gamma2(I);
  NumericVector alpha1h(I);
  NumericVector beta1h(I);
  NumericVector gamma1h(I);
  NumericVector alpha1sh(I);
  NumericVector beta1sh(I);
  NumericVector gamma1sh(I);
  NumericVector alpha2h(I);
  NumericVector beta2h(I);
  NumericVector gamma2h(I);
  NumericVector alpha2sh(I);
  NumericVector beta2sh(I);
  NumericVector gamma2sh(I);
  int n_cumsum_start = 0;
  for(int i = 0; i < I; i++) {
    NumericVector s_i(n_i(i));
    NumericVector h_i(n_i(i));
    NumericVector xi_hat_i(n_i(i));
    NumericVector epsilon_hat_i(n_i(i));
    for(int gg = 0; gg < n_i(i); gg++) {
      s_i(gg) = s(n_cumsum_start + gg);
      h_i(gg) = h(n_cumsum_start + gg);
      xi_hat_i(gg) = xi_hat(n_cumsum_start + gg);
      epsilon_hat_i(gg) = epsilon_hat(n_cumsum_start + gg);
    }
    if (n_i(i) > 2){
        NumericVector s_i_sub1(n_i(i)-1);
        NumericVector s_i_sub2(n_i(i)-1);
        NumericVector h_i_sub1(n_i(i)-1);
        NumericVector h_i_sub2(n_i(i)-1);
        NumericVector xi_hat_i_sub1(n_i(i)-1);
        NumericVector xi_hat_i_sub2(n_i(i)-1);
        NumericVector epsilon_hat_i_sub1(n_i(i)-1);
        NumericVector epsilon_hat_i_sub2(n_i(i)-1);
        for(int gg = 0; gg < (n_i(i)-1); gg++) {
          s_i_sub1(gg) = s_i(gg);
          s_i_sub2(gg) = s_i(gg+1);
          h_i_sub1(gg) = h_i(gg);
          h_i_sub2(gg) = h_i(gg+1);
          xi_hat_i_sub1(gg) = xi_hat_i(gg);
          xi_hat_i_sub2(gg) = xi_hat_i(gg+1);
          epsilon_hat_i_sub1(gg) = epsilon_hat_i(gg);
          epsilon_hat_i_sub2(gg) = epsilon_hat_i(gg+1);
        }
      alpha1(i) = sum(pow(s_i,2.0)*pow(xi_hat_i,2.0));
      beta1(i) = alpha1(i) - pow(s_i(0)*xi_hat_i(0),2.0) - pow(s_i(n_i(i)-1)*xi_hat_i(n_i(i)-1),2.0);
      gamma1(i) = 2*sum(s_i_sub1*s_i_sub2*xi_hat_i_sub1*xi_hat_i_sub2);
      alpha2(i) = sum(pow(s_i,2.0)*xi_hat_i*epsilon_hat_i);
      beta2(i) = alpha2(i) - pow(s_i(0),2.0)*xi_hat_i(0)*epsilon_hat_i(0) - pow(s_i(n_i(i)-1),2.0)*xi_hat_i(n_i(i)-1)*epsilon_hat_i(n_i(i)-1);
      gamma2(i) = sum(s_i_sub1*s_i_sub2*(xi_hat_i_sub2*epsilon_hat_i_sub1+xi_hat_i_sub1*epsilon_hat_i_sub2));
      alpha1h(i) = sum(pow(h_i,2.0)*pow(xi_hat_i,2.0));
      beta1h(i) = alpha1(i) - pow(h_i(0)*xi_hat_i(0),2.0) - pow(h_i(n_i(i)-1)*xi_hat_i(n_i(i)-1),2.0);
      gamma1h(i) = 2*sum(h_i_sub1*h_i_sub2*xi_hat_i_sub1*xi_hat_i_sub2);
      alpha1sh(i) = sum(s_i*h_i*pow(xi_hat_i,2.0));
      beta1sh(i) = alpha1(i) - s_i(0)*h_i(0)*pow(xi_hat_i(0),2.0) - s_i(n_i(i)-1)*h_i(n_i(i)-1)*pow(xi_hat_i(n_i(i)-1),2.0);
      gamma1sh(i) = sum(h_i_sub1*s_i_sub2*xi_hat_i_sub1*xi_hat_i_sub2) + sum(s_i_sub1*h_i_sub2*xi_hat_i_sub1*xi_hat_i_sub2);
      alpha2h(i) = sum(pow(h_i,2.0)*xi_hat_i*epsilon_hat_i);
      beta2h(i) = alpha2(i) - pow(h_i(0),2.0)*xi_hat_i(0)*epsilon_hat_i(0) - pow(h_i(n_i(i)-1),2.0)*xi_hat_i(n_i(i)-1)*epsilon_hat_i(n_i(i)-1);
      gamma2h(i) = sum(h_i_sub1*h_i_sub2*(xi_hat_i_sub2*epsilon_hat_i_sub1+xi_hat_i_sub1*epsilon_hat_i_sub2));
      alpha2sh(i) = sum(s_i*h_i*xi_hat_i*epsilon_hat_i);
      beta2sh(i) = alpha2(i) - s_i(0)*h_i(0)*xi_hat_i(0)*epsilon_hat_i(0) - s_i(n_i(i)-1)*h_i(n_i(i)-1)*xi_hat_i(n_i(i)-1)*epsilon_hat_i(n_i(i)-1);
      gamma2sh(i) = sum((s_i_sub1*h_i_sub2+h_i_sub1*s_i_sub2)*0.5*(xi_hat_i_sub2*epsilon_hat_i_sub1+xi_hat_i_sub1*epsilon_hat_i_sub2));
      A1(i) = alpha1(i) + pow(theta,2.0)*beta1(i) - theta*gamma1(i);
      A2(i) = alpha2(i) + pow(theta,2.0)*beta2(i) - theta*gamma2(i);
      Ab(i) = 2*( alpha2sh(i) + pow(theta,2.0)*beta2sh(i) - theta*gamma2sh(i) );
      Ac(i) = alpha2h(i) + pow(theta,2.0)*beta2h(i) - theta*gamma2h(i);
      Bb(i) = alpha1sh(i) + pow(theta,2.0)*beta1sh(i) - theta*gamma1sh(i);
      Bc(i) = alpha1h(i) + pow(theta,2.0)*beta1h(i) - theta*gamma1h(i);
    }
    else if (n_i(i) == 2){
      alpha1(i) = sum(pow(s_i,2.0)*pow(xi_hat_i,2.0));
      gamma1(i) = 2*(s_i(1)*s_i(0)*xi_hat_i(1)*xi_hat_i(0));
      alpha2(i) = sum(pow(s_i,2.0)*xi_hat_i*epsilon_hat_i);
      gamma2(i) = s_i(1)*s_i(0)*(xi_hat_i(0)*epsilon_hat_i(1)+xi_hat_i(1)*epsilon_hat_i(0));
      alpha1h(i) = sum(pow(h_i,2.0)*pow(xi_hat_i,2.0));
      gamma1h(i) = 2*h_i(0)*h_i(1)*xi_hat_i(0)*xi_hat_i(1);
      alpha1sh(i) = sum(s_i*h_i*pow(xi_hat_i,2.0));
      gamma1sh(i) = h_i(0)*s_i(1)*xi_hat_i(0)*xi_hat_i(1) + s_i(0)*h_i(1)*xi_hat_i(0)*xi_hat_i(1);
      alpha2h(i) = sum(pow(h_i,2.0)*xi_hat_i*epsilon_hat_i);
      gamma2h(i) = h_i(0)*h_i(1)*(xi_hat_i(1)*epsilon_hat_i(0)+xi_hat_i(0)*epsilon_hat_i(1));
      alpha2sh(i) = sum(s_i*h_i*xi_hat_i*epsilon_hat_i);
      gamma2sh(i) = (s_i(0)*h_i(1)+h_i(0)*s_i(1))*0.5*(xi_hat_i(1)*epsilon_hat_i(0)+xi_hat_i(0)*epsilon_hat_i(1));
      A1(i) = alpha1(i) - theta*gamma1(i);
      A2(i) = alpha2(i) - theta*gamma2(i);
      Ab(i) = 2*( alpha2sh(i) - theta*gamma2sh(i) );
      Ac(i) = alpha2h(i) - theta*gamma2h(i);
      Bb(i) = alpha1sh(i) - theta*gamma1sh(i);
      Bc(i) = alpha1h(i) - theta*gamma1h(i);
    }
    else if (n_i(i) == 1){
      alpha1(i) = 2*pow(s_i(0),2.0)*pow(xi_hat_i(0),2.0);
      alpha2(i) = 2*pow(s_i(0),2.0)*xi_hat_i(0)*epsilon_hat_i(0);
      alpha1h(i) = pow(h_i(0),2.0)*pow(xi_hat_i(0),2.0);
      alpha1sh(i) = s_i(0)*h_i(0)*pow(xi_hat_i(0),2.0);
      alpha2h(i) = pow(h_i(0),2.0)*xi_hat_i(0)*epsilon_hat_i(0);
      alpha2sh(i) = s_i(0)*h_i(0)*xi_hat_i(0)*epsilon_hat_i(0);
      A1(i) = (1-pow(theta,2.0))*alpha1(i);
      A2(i) = (1-pow(theta,2.0))*alpha2(i);
      Ab(i) = (1-pow(theta,2.0))*2*alpha2sh(i);
      Ac(i) = (1-pow(theta,2.0))*alpha2h(i);
      Bb(i) = (1-pow(theta,2.0))*alpha1sh(i);
      Bc(i) = (1-pow(theta,2.0))*alpha1h(i);
    }
    n_cumsum_start += n_i(i);
  }

  double a1 = sum(pow(A2,2.0));
  double a2 = -2*sum(A2*Ab);
  double a4 = sum(pow(Ab,2.0));
  double b1 = sum(A1);
  double b2 = -2*sum(Bb);
  double b4 = sum(Bc);

  double eta_optimal = 0.5*(2*a1*b2-a2*b1)/(a4*b1 - 2*a2*b2 - 2*a1*b4 + 3*(a1*pow(b2,2.0)/b1));

  return eta_optimal;
}')


    }
  }

  # Functions that incorporate the DML procedure for each method in Section 5.2.1

  WDML_HPLR_oj <- function(l_formula, l_learner, m_formula, m_learner, s_formula, s_learner, boost_control, data, K, m_stop=100, nu_s=0.1, nu_theta=0.1) {

    cv_folds <- groupKFold(data$id,k=K)
    beta_hat_num_k <- beta_hat_den_k <- V_hat_num_k <- numeric(K)
    for (k in seq_len(K)) {

      # Sample splitting
      cv_fold <- cv_folds[[k]]
      data.nuisance <- data[cv_fold,]
      data.beta <- data[-cv_fold,]

      # Nuisance function fitting
      fit.l <- regress.l(l_formula, l_learner, data.nuisance, data.beta)
      l_residuals.nuisance <- fit.l$l_residuals.nuisance; l_residuals.beta <- fit.l$l_residuals.beta

      fit.m <- regress.m(m_formula, m_learner, data.nuisance, data.beta)
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
      s.boost.output <- s.boost(boostdf.nuisance=boostdf.nuisance, boostdf.beta=boostdf.beta, s_formula=s_formula, s_learner=s_learner, m_stop=m_stop, early_stopping=TRUE, nu_s=nu_s, nu_theta=nu_theta)
      s.beta <- s.boost.output$s.beta; theta <- s.boost.output$theta
      plot(boostdf.beta$time[order(boostdf.beta$time)],max(s.beta[order(boostdf.beta$time)])/s.beta[order(boostdf.beta$time)],type="l",xlab="x",ylab="1/s(x)")
      groups.beta <- unique(boostdf.beta[,"id"])
      I.beta <- length(groups.beta)
      beta_num = 0; beta_den = 0; V_num = 0
      for (i in seq_len(I.beta)) {
        Y_minus_l_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"Y_minus_l_hat"]
        xi_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"xi_hat"]
        epsilon_hat_i = boostdf.beta[is.element(data.beta$id,groups.beta[i]),"epsilon_hat"]
        s_i <- s.beta[is.element(boostdf.beta$id,groups.beta[i])]
        N <- length(s_i)
        beta_num = beta_num + sum(s_i^2*Y_minus_l_hat_i*xi_hat_i) - (theta/(1+N*theta))*sum(s_i*Y_minus_l_hat_i)*sum(s_i*xi_hat_i)
        beta_den = beta_den + sum(s_i^2*xi_hat_i^2) - (theta/(1+N*theta))*(sum(s_i*xi_hat_i))^2
        V_num = V_num + (sum(s_i^2*epsilon_hat_i*xi_hat_i) - (theta/(1+N*theta))*sum(s_i*epsilon_hat_i)*sum(s_i*xi_hat_i))^2
      }
      beta_hat_num_k[k] <- beta_num
      beta_hat_den_k[k] <- beta_den
      V_hat_num_k[k] <- V_num

      #cat("\n")
      #cat('Progress: Fold ',k,' of ',K,' complete')
    }

    beta_hat <- sum(beta_hat_num_k)/sum(beta_hat_den_k)
    V_hat_div_n <- sum(V_hat_num_k)/(sum(beta_hat_den_k))^2#for CIS

    PLR_output <- list(beta_hat = beta_hat, CI_width = sqrt(V_hat_div_n))

    return(PLR_output)
  }

  WDML_HPLR_MEM_intercept = function(l_formula, l_learner, m_formula, m_learner, s_formula, s_learner, boost_control, data, K) {

    cv_folds <- groupKFold(data$id,k=K)
    beta_hat_num_k <- beta_hat_den_k <- V_hat_num_k <- numeric(K)
    rhoo=numeric(2)
    for (k in seq_len(K)) {

      cv_fold <- cv_folds[[k]]
      data.nuisance <- data[cv_fold,]
      data.beta <- data[-cv_fold,]

      fit.l <- regress.l(l_formula, l_learner, data.nuisance, data.beta)
      l_residuals.nuisance <- fit.l$l_residuals.nuisance; l_residuals.beta <- fit.l$l_residuals.beta

      fit.m <- regress.m(m_formula, m_learner, data.nuisance, data.beta)
      m_residuals.nuisance <- fit.m$m_residuals.nuisance; m_residuals.beta <- fit.m$m_residuals.beta

      xi_hat.nuisance <- m_residuals.nuisance
      beta_unweighted.nuisance <- sum(l_residuals.nuisance*xi_hat.nuisance)/sum(xi_hat.nuisance^2)
      epsilon_hat.nuisance <- l_residuals.nuisance - beta_unweighted.nuisance*xi_hat.nuisance

      xi_hat.beta <- m_residuals.beta
      beta_unweighted.beta <- sum(l_residuals.beta*xi_hat.beta)/sum(xi_hat.beta^2)
      epsilon_hat.beta <- l_residuals.beta - beta_unweighted.beta*xi_hat.beta

      #
      boostdf.nuisance <- cbind(data.nuisance, Y_minus_l_hat=l_residuals.nuisance, epsilon_hat=epsilon_hat.nuisance , xi_hat=xi_hat.nuisance)
      boostdf.beta <- cbind(data.beta, Y_minus_l_hat=l_residuals.beta, epsilon_hat=epsilon_hat.beta , xi_hat=xi_hat.beta)

      # lmer MEM
      D <- all.vars(m_formula)[1];   X <- all.vars(m_formula)[-1]
      lmer_response <- l_residuals.nuisance
      data.lmer <- cbind(data.nuisance, lmer_response=lmer_response, cesdd=data.nuisance[,D]-m_residuals.nuisance)##########
      lmer_fit <- lmer(lmer_response ~ cesdd - 1 + (1|id), data=data.lmer)
      lmer_sigma <- sigma(lmer_fit)
      lmer_cov_matrix <- summary(lmer_fit)$varcor$id
      groups.beta <- unique(boostdf.beta[,"id"])

      beta_num = 0; beta_den = 0; V_num = 0
      I.beta <- length(groups.beta)
      for (i in seq_len(I.beta)) {
        Y_minus_l_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"Y_minus_l_hat"]
        xi_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"xi_hat"]
        epsilon_hat_i =  boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"epsilon_hat"]
        N <- length(xi_hat_i)
        Sigma <- matrix(c(rep(1,N)),nrow=N,ncol=1) %*% lmer_cov_matrix %*% t(matrix(c(rep(1,N)),nrow=N,ncol=1)) + lmer_sigma^2*diag(N)
        Sigma_inv <- solve(Sigma)
        beta_num = beta_num + t(Y_minus_l_hat_i)%*%Sigma_inv%*%xi_hat_i
        beta_den = beta_den + t(xi_hat_i)%*%Sigma_inv%*%xi_hat_i
        V_num = V_num + (t(epsilon_hat_i)%*%Sigma_inv%*%xi_hat_i)^2
      }
      beta_hat_num_k[k] <- beta_num
      beta_hat_den_k[k] <- beta_den
      V_hat_num_k[k] <- V_num

      cat("\n")
      cat('Progress: Fold ',k,' of ',K,' complete')
    }

    beta_hat <- sum(beta_hat_num_k)/sum(beta_hat_den_k)
    V_hat_div_n <- sum(V_hat_num_k)/(sum(beta_hat_den_k))^2#for CIS

    PLR_output <- list(beta_hat = beta_hat, CI_width = sqrt(V_hat_div_n), rhoo=rhoo)

    return(PLR_output)
  }

  WDML_HPLR_MEM_int_time = function(l_formula, l_learner, m_formula, m_learner, s_formula, s_learner, boost_control, data, K) {

    cv_folds <- groupKFold(data$id,k=K)
    beta_hat_num_k <- numeric(K); beta_hat_den_k <- numeric(K); V_hat_num_k <- numeric(K)
    rhoo=numeric(2)
    for (k in seq_len(K)) {
      cv_fold <- cv_folds[[k]]
      data.nuisance <- data[cv_fold,]
      data.beta <- data[-cv_fold,]

      fit.l <- regress.l(l_formula, l_learner, data.nuisance, data.beta)
      l_residuals.nuisance <- fit.l$l_residuals.nuisance; l_residuals.beta <- fit.l$l_residuals.beta

      fit.m <- regress.m(m_formula, m_learner, data.nuisance, data.beta)
      m_residuals.nuisance <- fit.m$m_residuals.nuisance; m_residuals.beta <- fit.m$m_residuals.beta

      xi_hat.nuisance <- m_residuals.nuisance
      beta_unweighted.nuisance <- sum(l_residuals.nuisance*xi_hat.nuisance)/sum(xi_hat.nuisance^2)
      epsilon_hat.nuisance <- l_residuals.nuisance - beta_unweighted.nuisance*xi_hat.nuisance

      xi_hat.beta <- m_residuals.beta
      beta_unweighted.beta <- sum(l_residuals.beta*xi_hat.beta)/sum(xi_hat.beta^2)
      epsilon_hat.beta <- l_residuals.beta - beta_unweighted.beta*xi_hat.beta

      #
      boostdf.nuisance <- cbind(data.nuisance, Y_minus_l_hat=l_residuals.nuisance, epsilon_hat=epsilon_hat.nuisance , xi_hat=xi_hat.nuisance)
      boostdf.beta <- cbind(data.beta, Y_minus_l_hat=l_residuals.beta, epsilon_hat=epsilon_hat.beta , xi_hat=xi_hat.beta)

      D <- all.vars(m_formula)[1];   XX <- all.vars(m_formula)[-1]
      lmer_response <- l_residuals.nuisance##########
      data.lmer <- cbind(data.nuisance, lmer_response=lmer_response, cesdd=data.nuisance[,D]-m_residuals.nuisance)##########
      lmer_fit <- lmer(lmer_response ~ cesdd - 1 + (1 + time|id), data=data.lmer)
      lmer_sigma <- sigma(lmer_fit)
      lmer_cov_matrix <- summary(lmer_fit)$varcor$id
      groups.beta <- unique(boostdf.beta[,"id"])

      beta_num = 0; beta_den = 0; V_num = 0
      I.beta <- length(groups.beta)
      for (i in seq_len(I.beta)) {
        Y_minus_l_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"Y_minus_l_hat"]
        xi_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"xi_hat"]
        epsilon_hat_i =  boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"epsilon_hat"]
        N <- length(xi_hat_i)
        Sigma <- matrix(c(rep(1,N),data[is.element(data$id,groups.beta[i]),"time"]),nrow=N,ncol=2) %*% lmer_cov_matrix %*% t(matrix(c(rep(1,N),data[is.element(data$id,groups.beta[i]),"time"]),nrow=N,ncol=2)) + lmer_sigma^2*diag(N)
        Sigma_inv <- solve(Sigma)
        beta_num = beta_num + t(Y_minus_l_hat_i)%*%Sigma_inv%*%xi_hat_i
        beta_den = beta_den + t(xi_hat_i)%*%Sigma_inv%*%xi_hat_i
        V_num = V_num + (t(epsilon_hat_i)%*%Sigma_inv%*%xi_hat_i)^2
      }
      beta_hat_num_k[k] <- beta_num
      beta_hat_den_k[k] <- beta_den
      V_hat_num_k[k] <- V_num

      #cat("\n")
      #cat('Progress: Fold ',k,' of ',K,' complete')
    }

    beta_hat <- sum(beta_hat_num_k)/sum(beta_hat_den_k)
    V_hat_div_n <- sum(V_hat_num_k)/(sum(beta_hat_den_k))^2#for CIS

    PLR_output <- list(beta_hat = beta_hat, CI_width = sqrt(V_hat_div_n), rhoo=rhoo)

    return(PLR_output)
  }

  WDML_HPLR_geepack = function(l_formula, l_learner, m_formula, m_learner, s_formula, s_learner, boost_control, data, K) {

    cv_folds <- groupKFold(data$id,k=K)
    beta_hat_num_k <- numeric(K); beta_hat_den_k <- numeric(K); V_hat_num_k <- numeric(K)
    for (k in seq_len(K)) {
      cv_fold <- cv_folds[[k]]
      data.nuisance <- data[cv_fold,]
      data.beta <- data[-cv_fold,]

      fit.l <- regress.l(l_formula, l_learner, data.nuisance, data.beta)
      l_residuals.nuisance <- fit.l$l_residuals.nuisance; l_residuals.beta <- fit.l$l_residuals.beta

      fit.m <- regress.m(m_formula, m_learner, data.nuisance, data.beta)
      m_residuals.nuisance <- fit.m$m_residuals.nuisance; m_residuals.beta <- fit.m$m_residuals.beta

      xi_hat.nuisance <- m_residuals.nuisance
      beta_unweighted.nuisance <- sum(l_residuals.nuisance*xi_hat.nuisance)/sum(xi_hat.nuisance^2)
      epsilon_hat.nuisance <- l_residuals.nuisance - beta_unweighted.nuisance*xi_hat.nuisance

      xi_hat.beta <- m_residuals.beta
      beta_unweighted.beta <- sum(l_residuals.beta*xi_hat.beta)/sum(xi_hat.beta^2)
      epsilon_hat.beta <- l_residuals.beta - beta_unweighted.beta*xi_hat.beta

      #
      boostdf.nuisance <- cbind(data.nuisance, Y_minus_l_hat=l_residuals.nuisance, epsilon_hat=epsilon_hat.nuisance , xi_hat=xi_hat.nuisance)
      boostdf.beta <- cbind(data.beta, Y_minus_l_hat=l_residuals.beta, epsilon_hat=epsilon_hat.beta , xi_hat=xi_hat.beta)

      # Fit variance function with cubic splines (See Appendix E.1)
      df.nuisance.X1 <- cbind(boostdf.nuisance, squared_epsilon_hat=boostdf.nuisance[,"epsilon_hat"]^2)
      df.beta.X1 <- boostdf.beta
      gam_fit <- gam(squared_epsilon_hat ~ s(X,bs="cr"), data=df.nuisance.X1) #cubic splines
      sigma.sq.nuisance <- predict(gam_fit, df.nuisance.X1[,"X",drop=FALSE])
      sigma.sq.nuisance[sigma.sq.nuisance<0.01]=0.01
      s.nuisance <- 1/sqrt(sigma.sq.nuisance)
      sigma.sq.beta <- predict(gam_fit, df.beta.X1[,"X",drop=FALSE])
      sigma.sq.beta[sigma.sq.beta<0.01]=0.01
      s.beta <- 1/sqrt(sigma.sq.beta)
      df.nuisance.X1[,c("D","Y","Y_minus_l_hat","epsilon_hat","xi_hat")] <- df.nuisance.X1[,c("D","Y","Y_minus_l_hat","epsilon_hat","xi_hat")]*s.nuisance
      df.beta.X1[,c("D","Y","Y_minus_l_hat","epsilon_hat","xi_hat")] <- df.beta.X1[,c("D","Y","Y_minus_l_hat","epsilon_hat","xi_hat")]*s.beta
      # Then fit exchangable GEE structure (see Appendix E.1)
      gee_fit <- geeglm(Y_minus_l_hat ~ xi_hat, data=df.nuisance.X1, corstr="exchangeable", id=id)
      alpha_or_rhoo <- as.numeric(gee_fit$geese[2])
      theta <- alpha_or_rhoo/(1-alpha_or_rhoo)
      # Aggregate
      groups.beta <- unique(boostdf.beta[,"id"])
      I.beta <- length(groups.beta)
      beta_num = 0; beta_den = 0; V_num = 0
      for (i in seq_len(I.beta)) {
        Y_minus_l_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"Y_minus_l_hat"]
        xi_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"xi_hat"]
        epsilon_hat_i = boostdf.beta[is.element(data.beta$id,groups.beta[i]),"epsilon_hat"]
        s_i <- s.beta[is.element(boostdf.beta$id,groups.beta[i])]
        N <- length(s_i)
        beta_num = beta_num + sum(s_i^2*Y_minus_l_hat_i*xi_hat_i) - (theta/(1+N*theta))*sum(s_i*Y_minus_l_hat_i)*sum(s_i*xi_hat_i)
        beta_den = beta_den + sum(s_i^2*xi_hat_i^2) - (theta/(1+N*theta))*(sum(s_i*xi_hat_i))^2
        V_num = V_num + (sum(s_i^2*epsilon_hat_i*xi_hat_i) - (theta/(1+N*theta))*sum(s_i*epsilon_hat_i)*sum(s_i*xi_hat_i))^2
      }
      beta_hat_num_k[k] <- beta_num
      beta_hat_den_k[k] <- beta_den
      V_hat_num_k[k] <- V_num

      beta_hat <- sum(beta_hat_num_k)/sum(beta_hat_den_k)
      V_hat_div_n <- sum(V_hat_num_k)/(sum(beta_hat_den_k))^2#for CIS

      #cat("\n")
      #cat('Progress: Fold ',k,' of ',K,' complete')
    }

    beta_hat <- sum(beta_hat_num_k)/sum(beta_hat_den_k)
    V_hat_div_n <- sum(V_hat_num_k)/(sum(beta_hat_den_k))^2#for CIS

    PLR_output <- list(beta_hat = beta_hat, CI_width = sqrt(V_hat_div_n))

    return(PLR_output)
  }

  WDML_HPLR_geepack_hom = function(l_formula, l_learner, m_formula, m_learner, s_formula, s_learner, boost_control, data, K) {

    cv_folds <- groupKFold(data$id,k=K)
    beta_hat_num_k <- numeric(K); beta_hat_den_k <- numeric(K); V_hat_num_k <- numeric(K)
    for (k in seq_len(K)) {
      cv_fold <- cv_folds[[k]]
      data.nuisance <- data[cv_fold,]
      data.beta <- data[-cv_fold,]

      fit.l <- regress.l(l_formula, l_learner, data.nuisance, data.beta)
      l_residuals.nuisance <- fit.l$l_residuals.nuisance; l_residuals.beta <- fit.l$l_residuals.beta

      fit.m <- regress.m(m_formula, m_learner, data.nuisance, data.beta)
      m_residuals.nuisance <- fit.m$m_residuals.nuisance; m_residuals.beta <- fit.m$m_residuals.beta

      xi_hat.nuisance <- m_residuals.nuisance
      beta_unweighted.nuisance <- sum(l_residuals.nuisance*xi_hat.nuisance)/sum(xi_hat.nuisance^2)
      epsilon_hat.nuisance <- l_residuals.nuisance - beta_unweighted.nuisance*xi_hat.nuisance

      xi_hat.beta <- m_residuals.beta
      beta_unweighted.beta <- sum(l_residuals.beta*xi_hat.beta)/sum(xi_hat.beta^2)
      epsilon_hat.beta <- l_residuals.beta - beta_unweighted.beta*xi_hat.beta

      #
      boostdf.nuisance <- cbind(data.nuisance, Y_minus_l_hat=l_residuals.nuisance, epsilon_hat=epsilon_hat.nuisance , xi_hat=xi_hat.nuisance)
      boostdf.beta <- cbind(data.beta, Y_minus_l_hat=l_residuals.beta, epsilon_hat=epsilon_hat.beta , xi_hat=xi_hat.beta)

      # Fit GEE with equicorrelated proxyCFF
      gee_fit <- geeglm(Y_minus_l_hat ~ xi_hat, data=boostdf.nuisance, corstr="exchangeable", id=id)
      alpha_or_rhoo <- as.numeric(gee_fit$geese[2])
      theta <- alpha_or_rhoo/(1-alpha_or_rhoo)
      # Aggregate
      groups.beta <- unique(boostdf.beta[,"id"])
      I.beta <- length(groups.beta)
      beta_num = 0; beta_den = 0; V_num = 0
      for (i in seq_len(I.beta)) {
        Y_minus_l_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"Y_minus_l_hat"]
        xi_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"xi_hat"]
        epsilon_hat_i = boostdf.beta[is.element(data.beta$id,groups.beta[i]),"epsilon_hat"]
        N <- length(xi_hat_i)
        beta_num = beta_num + sum(Y_minus_l_hat_i*xi_hat_i) - (theta/(1+N*theta))*sum(Y_minus_l_hat_i)*sum(xi_hat_i)
        beta_den = beta_den + sum(xi_hat_i^2) - (theta/(1+N*theta))*(sum(xi_hat_i))^2
        V_num = V_num + (sum(epsilon_hat_i*xi_hat_i) - (theta/(1+N*theta))*sum(epsilon_hat_i)*sum(xi_hat_i))^2
      }
      beta_hat_num_k[k] <- beta_num
      beta_hat_den_k[k] <- beta_den
      V_hat_num_k[k] <- V_num

      beta_hat <- sum(beta_hat_num_k)/sum(beta_hat_den_k)
      V_hat_div_n <- sum(V_hat_num_k)/(sum(beta_hat_den_k))^2#for CIS

      #cat("\n")
      #cat('Progress: Fold ',k,' of ',K,' complete')
    }

    beta_hat <- sum(beta_hat_num_k)/sum(beta_hat_den_k)
    V_hat_div_n <- sum(V_hat_num_k)/(sum(beta_hat_den_k))^2#for CIS

    PLR_output <- list(beta_hat = beta_hat, CI_width = sqrt(V_hat_div_n))

    return(PLR_output)
  }

  # Functions that incorporate the DML procedure for each method in Section 5.2.2

  WDML_HPLR_wages = function(l_formula, l_learner, m_formula, m_learner, s_formula, s_learner, boost_control, data, K, m_stop=100, nu_s=0.1, nu_theta=0.1) {

    cv_folds <- groupKFold(data$id,k=K)
    beta_hat_num_k <- numeric(K); beta_hat_den_k <- numeric(K); V_hat_num_k <- numeric(K)
    for (k in seq_len(K)) {
      cv_fold <- cv_folds[[k]]
      data.nuisance <- data[cv_fold,]
      data.beta <- data[-cv_fold,]

      fit.l <- regress.l(l_formula, l_learner, data.nuisance, data.beta)
      l_residuals.nuisance <- fit.l$l_residuals.nuisance; l_residuals.beta <- fit.l$l_residuals.beta

      fit.m <- regress.m(m_formula, m_learner, data.nuisance, data.beta)
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

      s.boost.output <- s.boost(boostdf.nuisance=boostdf.nuisance, boostdf.beta=boostdf.beta, s_formula=s_formula, s_learner=s_learner, m_stop=m_stop, early_stopping=TRUE, nu_s=nu_s, nu_theta=nu_theta)
      s.beta <- s.boost.output$s.beta; theta <- s.boost.output$theta

      groups.beta <- unique(boostdf.beta[,"id"])
      I.beta <- length(groups.beta)
      beta_num = 0; beta_den = 0; V_num = 0
      for (i in seq_len(I.beta)) {
        Y_minus_l_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"Y_minus_l_hat"]
        xi_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"xi_hat"]
        epsilon_hat_i = boostdf.beta[is.element(data.beta$id,groups.beta[i]),"epsilon_hat"]
        s_i <- s.beta[is.element(boostdf.beta$id,groups.beta[i])]
        N <- length(s_i)
        beta_num = beta_num + sum(s_i^2*Y_minus_l_hat_i*xi_hat_i) - (theta/(1+N*theta))*sum(s_i*Y_minus_l_hat_i)*sum(s_i*xi_hat_i)
        beta_den = beta_den + sum(s_i^2*xi_hat_i^2) - (theta/(1+N*theta))*(sum(s_i*xi_hat_i))^2
        V_num = V_num + (sum(s_i^2*epsilon_hat_i*xi_hat_i) - (theta/(1+N*theta))*sum(s_i*epsilon_hat_i)*sum(s_i*xi_hat_i))^2
      }
      beta_hat_num_k[k] <- beta_num
      beta_hat_den_k[k] <- beta_den
      V_hat_num_k[k] <- V_num

      #cat("\n")
      #cat('Progress: Fold ',k,' of ',K,' complete')
    }

    beta_hat <- sum(beta_hat_num_k)/sum(beta_hat_den_k)
    V_hat_div_n <- sum(V_hat_num_k)/(sum(beta_hat_den_k))^2#for CIS

    #confid <- 0.95#for CIS
    #confinterval <- data.frame(lower=beta_hat-sqrt(V_hat_div_n)*qnorm(1-(1-confid)/2),upper=beta_hat+sqrt(V_hat_div_n)*qnorm(1-(1-confid)/2))#for CIS

    #PLR_output <- list(beta_hat=beta_hat, CI_width=sqrt(V_hat_div_n), CI=confinterval, m.plot=m.plot, l.plot=l.plot)
    PLR_output <- list(beta_hat = beta_hat, CI_width = sqrt(V_hat_div_n))

    return(PLR_output)
  }

  WDML_HPLR_MEM_int_age_tenure_for_wages <- function(l_formula, l_learner, m_formula, m_learner, data, K) {

    cv_folds <- groupKFold(data$id,k=K)
    beta_hat_num_k <- numeric(K); beta_hat_den_k <- numeric(K); V_hat_num_k <- numeric(K)
    rhoo=numeric(2)
    for (k in seq_len(K)) {
      cv_fold <- cv_folds[[k]]
      data.nuisance <- data[cv_fold,]
      data.beta <- data[-cv_fold,]

      fit.l <- regress.l(l_formula, l_learner, data.nuisance, data.beta)
      l_residuals.nuisance <- fit.l$l_residuals.nuisance; l_residuals.beta <- fit.l$l_residuals.beta

      fit.m <- regress.m(m_formula, m_learner, data.nuisance, data.beta)
      m_residuals.nuisance <- fit.m$m_residuals.nuisance; m_residuals.beta <- fit.m$m_residuals.beta

      xi_hat.nuisance <- m_residuals.nuisance
      beta_unweighted.nuisance <- sum(l_residuals.nuisance*xi_hat.nuisance)/sum(xi_hat.nuisance^2)
      epsilon_hat.nuisance <- l_residuals.nuisance - beta_unweighted.nuisance*xi_hat.nuisance

      xi_hat.beta <- m_residuals.beta
      beta_unweighted.beta <- sum(l_residuals.beta*xi_hat.beta)/sum(xi_hat.beta^2)
      epsilon_hat.beta <- l_residuals.beta - beta_unweighted.beta*xi_hat.beta

      #
      boostdf.nuisance <- cbind(data.nuisance, Y_minus_l_hat=l_residuals.nuisance, epsilon_hat=epsilon_hat.nuisance , xi_hat=xi_hat.nuisance)
      boostdf.beta <- cbind(data.beta, Y_minus_l_hat=l_residuals.beta, epsilon_hat=epsilon_hat.beta , xi_hat=xi_hat.beta)

      # Fit lmer mixed effects model
      D <- all.vars(m_formula)[1];   XX <- all.vars(m_formula)[-1]
      lmer_response <- l_residuals.nuisance
      data.lmer <- cbind(data.nuisance, lmer_response=lmer_response, cesdd=data.nuisance[,D]-m_residuals.nuisance)
      lmer_fit <- lmer(lmer_response ~ cesdd - 1 + (1 + time + age|id), data=data.lmer)
      lmer_sigma <- sigma(lmer_fit)
      lmer_cov_matrix <- summary(lmer_fit)$varcor$id
      groups.beta <- unique(boostdf.beta[,"id"])

      beta_num = 0; beta_den = 0; V_num = 0
      I.beta <- length(groups.beta)
      for (i in seq_len(I.beta)) {
        Y_minus_l_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"Y_minus_l_hat"]
        xi_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"xi_hat"]
        epsilon_hat_i =  boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"epsilon_hat"]
        N <- length(xi_hat_i)
        Sigma <- matrix(c(rep(1,N),data[is.element(data$id,groups.beta[i]),"time"],data[is.element(data$id,groups.beta[i]),"age"]),nrow=N,ncol=3) %*% lmer_cov_matrix %*% t(matrix(c(rep(1,N),data[is.element(data$id,groups.beta[i]),"time"],data[is.element(data$id,groups.beta[i]),"age"]),nrow=N,ncol=3)) + lmer_sigma^2*diag(N)
        Sigma_inv <- solve(Sigma)
        beta_num = beta_num + t(Y_minus_l_hat_i)%*%Sigma_inv%*%xi_hat_i
        beta_den = beta_den + t(xi_hat_i)%*%Sigma_inv%*%xi_hat_i
        V_num = V_num + (t(epsilon_hat_i)%*%Sigma_inv%*%xi_hat_i)^2
      }
      beta_hat_num_k[k] <- beta_num
      beta_hat_den_k[k] <- beta_den
      V_hat_num_k[k] <- V_num

      #cat("\n")
      #cat('Progress: Fold ',k,' of ',K,' complete')
    }

    beta_hat <- sum(beta_hat_num_k)/sum(beta_hat_den_k)
    V_hat_div_n <- sum(V_hat_num_k)/(sum(beta_hat_den_k))^2#for CIS

    PLR_output <- list(beta_hat = beta_hat, CI_width = sqrt(V_hat_div_n), rhoo=rhoo)

    return(PLR_output)
  }

  WDML_HPLR_geepack_for_wages <- function(l_formula, l_learner, m_formula, m_learner, data, K) {

    cv_folds <- groupKFold(data$id,k=K)
    beta_hat_num_k <- numeric(K); beta_hat_den_k <- numeric(K); V_hat_num_k <- numeric(K)
    for (k in seq_len(K)) {
      cv_fold <- cv_folds[[k]]
      data.nuisance <- data[cv_fold,]
      data.beta <- data[-cv_fold,]

      fit.l <- regress.l(l_formula, l_learner, data.nuisance, data.beta)
      l_residuals.nuisance <- fit.l$l_residuals.nuisance; l_residuals.beta <- fit.l$l_residuals.beta

      fit.m <- regress.m(m_formula, m_learner, data.nuisance, data.beta)
      m_residuals.nuisance <- fit.m$m_residuals.nuisance; m_residuals.beta <- fit.m$m_residuals.beta

      xi_hat.nuisance <- m_residuals.nuisance
      beta_unweighted.nuisance <- sum(l_residuals.nuisance*xi_hat.nuisance)/sum(xi_hat.nuisance^2)
      epsilon_hat.nuisance <- l_residuals.nuisance - beta_unweighted.nuisance*xi_hat.nuisance

      xi_hat.beta <- m_residuals.beta
      beta_unweighted.beta <- sum(l_residuals.beta*xi_hat.beta)/sum(xi_hat.beta^2)
      epsilon_hat.beta <- l_residuals.beta - beta_unweighted.beta*xi_hat.beta

      #
      boostdf.nuisance <- cbind(data.nuisance, Y_minus_l_hat=l_residuals.nuisance, epsilon_hat=epsilon_hat.nuisance , xi_hat=xi_hat.nuisance)
      boostdf.beta <- cbind(data.beta, Y_minus_l_hat=l_residuals.beta, epsilon_hat=epsilon_hat.beta , xi_hat=xi_hat.beta)

      # Fit variance function with cubic splines (See Appendix E.1)
      df.nuisance.X1 <- cbind(boostdf.nuisance, squared_epsilon_hat=boostdf.nuisance[,"epsilon_hat"]^2)
      boostdf.beta.X1 <- boostdf.beta
      gam_fit <- gam(squared_epsilon_hat ~ s(tenure,bs="cr")+s(age,bs="cr"), data=df.nuisance.X1) #cubic splines
      sigma.sq.nuisance <- predict(gam_fit, df.nuisance.X1[,c("tenure","age"),drop=FALSE])
      sigma.sq.nuisance[sigma.sq.nuisance<0.01]=0.01
      s.nuisance <- 1/sqrt(sigma.sq.nuisance)
      sigma.sq.beta <- predict(gam_fit, boostdf.beta.X1[,c("tenure","age"),drop=FALSE])
      sigma.sq.beta[sigma.sq.beta<0.01]=0.01
      s.beta <- 1/sqrt(sigma.sq.beta)
      df.nuisance.X1[,c("ttl_exp","ln_wage","Y_minus_l_hat","epsilon_hat","xi_hat")] <- df.nuisance.X1[,c("ttl_exp","ln_wage","Y_minus_l_hat","epsilon_hat","xi_hat")]*s.nuisance
      boostdf.beta.X1[,c("ttl_exp","ln_wage","Y_minus_l_hat","epsilon_hat","xi_hat")] <- boostdf.beta.X1[,c("ttl_exp","ln_wage","Y_minus_l_hat","epsilon_hat","xi_hat")]*s.beta
      # Then fit exchangable GEE structure (see Appendix E.1)
      gee_fit <- geeglm(Y_minus_l_hat ~ xi_hat, data=boostdf.nuisance.X1, corstr="exchangeable", id=id)
      alpha_or_rhoo <- as.numeric(gee_fit$geese[2])
      theta <- alpha_or_rhoo/(1-alpha_or_rhoo)
      # Aggregate
      groups.beta <- unique(boostdf.beta[,"id"])
      I.beta <- length(groups.beta)
      beta_num = 0; beta_den = 0; V_num = 0
      for (i in seq_len(I.beta)) {
        Y_minus_l_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"Y_minus_l_hat"]
        xi_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"xi_hat"]
        epsilon_hat_i = boostdf.beta[is.element(data.beta$id,groups.beta[i]),"epsilon_hat"]
        s_i <- s.beta[is.element(boostdf.beta$id,groups.beta[i])]
        N <- length(s_i)
        beta_num = beta_num + sum(s_i^2*Y_minus_l_hat_i*xi_hat_i) - (theta/(1+N*theta))*sum(s_i*Y_minus_l_hat_i)*sum(s_i*xi_hat_i)
        beta_den = beta_den + sum(s_i^2*xi_hat_i^2) - (theta/(1+N*theta))*(sum(s_i*xi_hat_i))^2
        V_num = V_num + (sum(s_i^2*epsilon_hat_i*xi_hat_i) - (theta/(1+N*theta))*sum(s_i*epsilon_hat_i)*sum(s_i*xi_hat_i))^2
      }
      beta_hat_num_k[k] <- beta_num
      beta_hat_den_k[k] <- beta_den
      V_hat_num_k[k] <- V_num

      beta_hat <- sum(beta_hat_num_k)/sum(beta_hat_den_k)
      V_hat_div_n <- sum(V_hat_num_k)/(sum(beta_hat_den_k))^2#for CIS

      #cat("\n")
      #cat('Progress: Fold ',k,' of ',K,' complete')
    }

    beta_hat <- sum(beta_hat_num_k)/sum(beta_hat_den_k)
    V_hat_div_n <- sum(V_hat_num_k)/(sum(beta_hat_den_k))^2#for CIS

    PLR_output <- list(beta_hat = beta_hat, CI_width = sqrt(V_hat_div_n))

    return(PLR_output)
  }

  WDML_HPLR_geepack_hom_wages <- function(l_formula, l_learner, m_formula, m_learner, data, K) {

    cv_folds <- groupKFold(data$id,k=K)
    beta_hat_num_k <- numeric(K); beta_hat_den_k <- numeric(K); V_hat_num_k <- numeric(K)
    for (k in seq_len(K)) {
      cv_fold <- cv_folds[[k]]
      data.nuisance <- data[cv_fold,]
      data.beta <- data[-cv_fold,]

      fit.l <- regress.l(l_formula, l_learner, data.nuisance, data.beta)
      l_residuals.nuisance <- fit.l$l_residuals.nuisance; l_residuals.beta <- fit.l$l_residuals.beta

      fit.m <- regress.m(m_formula, m_learner, data.nuisance, data.beta)
      m_residuals.nuisance <- fit.m$m_residuals.nuisance; m_residuals.beta <- fit.m$m_residuals.beta

      xi_hat.nuisance <- m_residuals.nuisance
      beta_unweighted.nuisance <- sum(l_residuals.nuisance*xi_hat.nuisance)/sum(xi_hat.nuisance^2)
      epsilon_hat.nuisance <- l_residuals.nuisance - beta_unweighted.nuisance*xi_hat.nuisance

      xi_hat.beta <- m_residuals.beta
      beta_unweighted.beta <- sum(l_residuals.beta*xi_hat.beta)/sum(xi_hat.beta^2)
      epsilon_hat.beta <- l_residuals.beta - beta_unweighted.beta*xi_hat.beta

      #
      boostdf.nuisance <- cbind(data.nuisance, Y_minus_l_hat=l_residuals.nuisance, epsilon_hat=epsilon_hat.nuisance , xi_hat=xi_hat.nuisance)
      boostdf.beta <- cbind(data.beta, Y_minus_l_hat=l_residuals.beta, epsilon_hat=epsilon_hat.beta , xi_hat=xi_hat.beta)

      # Fit (homoscedastic) GEE model
      gee_fit <- geeglm(Y_minus_l_hat ~ xi_hat, data=boostdf.nuisance, corstr="exchangeable", id=id)
      alpha_or_rhoo <- as.numeric(gee_fit$geese[2])
      theta <- alpha_or_rhoo/(1-alpha_or_rhoo)

      groups.beta <- unique(boostdf.beta[,"id"])
      I.beta <- length(groups.beta)
      beta_num = 0; beta_den = 0; V_num = 0
      for (i in seq_len(I.beta)) {
        Y_minus_l_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"Y_minus_l_hat"]
        xi_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"xi_hat"]
        epsilon_hat_i = boostdf.beta[is.element(data.beta$id,groups.beta[i]),"epsilon_hat"]
        N <- length(xi_hat_i)
        beta_num = beta_num + sum(Y_minus_l_hat_i*xi_hat_i) - (theta/(1+N*theta))*sum(Y_minus_l_hat_i)*sum(xi_hat_i)
        beta_den = beta_den + sum(xi_hat_i^2) - (theta/(1+N*theta))*(sum(xi_hat_i))^2
        V_num = V_num + (sum(epsilon_hat_i*xi_hat_i) - (theta/(1+N*theta))*sum(epsilon_hat_i)*sum(xi_hat_i))^2
      }
      beta_hat_num_k[k] <- beta_num
      beta_hat_den_k[k] <- beta_den
      V_hat_num_k[k] <- V_num

      beta_hat <- sum(beta_hat_num_k)/sum(beta_hat_den_k)
      V_hat_div_n <- sum(V_hat_num_k)/(sum(beta_hat_den_k))^2#for CIS

      #cat("\n")
      #cat('Progress: Fold ',k,' of ',K,' complete')
    }

    beta_hat <- sum(beta_hat_num_k)/sum(beta_hat_den_k)
    V_hat_div_n <- sum(V_hat_num_k)/(sum(beta_hat_den_k))^2#for CIS

        PLR_output <- list(beta_hat = beta_hat, CI_width = sqrt(V_hat_div_n))

    return(PLR_output)
  }

  s.boost <- function(boostdf.nuisance=boostdf.nuisance, boostdf.beta=boostdf.beta, s_formula=s_formula, s_learner=s_learner, m_stop=NULL, early_stopping=TRUE, nu_s=0.1, nu_theta=0.1) {
    cv.m_stop <- 1000 # selected by cross-validation
    {
      boostdf.nuisance.train <- boostdf.nuisance
      boostdf.nuisance.test <- boostdf.nuisance

      # Grouping Details
      grouping.nuisance.train <- grouping.metadata(boostdf.nuisance.train)
      groups.nuisance.train <- grouping.nuisance.train$groups.metadata
      I.nuisance.train <- grouping.nuisance.train$I.metadata
      n_i.nuisance.train <- grouping.nuisance.train$n_i.metadata
      n_obs.nuisance.train <- grouping.nuisance.train$n_obs.metadata

      grouping.nuisance.test <- grouping.metadata(boostdf.nuisance.test)
      groups.nuisance.test <- grouping.nuisance.test$groups.metadata
      I.nuisance.test <- grouping.nuisance.test$I.metadata
      n_i.nuisance.test <- grouping.nuisance.test$n_i.metadata
      n_obs.nuisance.test <- grouping.nuisance.test$n_obs.metadata

      grouping.beta <- grouping.metadata(boostdf.beta)
      groups.beta <- grouping.beta$groups.metadata
      I.beta <- grouping.beta$I.metadata
      n_i.beta <- grouping.beta$n_i.metadata
      n_obs.beta <- grouping.beta$n_obs.metadata

      # Initalise boosting
      s.train <- rep(1,n_obs.nuisance.train)
      theta <- 0
      s.test <- rep(1,n_obs.nuisance.test)
      s.beta <- rep(1,n_obs.beta)

      U_S <- all.vars(s_formula)[1]
      if (U_S != "u_s") {
        stop("Error: s_formula should have response in the form 'u_s ~ predictors'")
      }
      X <- all.vars(s_formula)[-1]

      for (m in seq_len(cv.m_stop)) {

        # Calculate scores
        scores <- ngradient_cpp(boostdf.nuisance.train[,"epsilon_hat"], boostdf.nuisance.train[,"xi_hat"], I.nuisance.train, n_obs.nuisance.train, n_i.nuisance.train, s.train, theta)
        u_s <- scores$s
        fit.s <- regress.s(s_formula, s_learner, boostdf.nuisance.train, boostdf.nuisance.test, boostdf.beta, u_s, U_S, X)
        u_s_hat.train <- fit.s$u_s_hat.train; u_s_hat.test <- fit.s$u_s_hat.test; u_s_hat.beta <- fit.s$u_s_hat.beta

        # Line Search for optimal step size
        eta_optimal <- optimal_step_size_cpp(boostdf.nuisance.train[,"epsilon_hat"], boostdf.nuisance.train[,"xi_hat"], I.nuisance.train, n_obs.nuisance.train, n_i.nuisance.train, s.train, u_s_hat.train, theta)
        eta_optimal <- max(min(eta_optimal,1), 0.001)
        # s boosting iteration (with shrinkage)
        s.test <- s.test - nu_s * eta_optimal * u_s_hat.test
        s.beta <- s.beta - nu_s * eta_optimal * u_s_hat.beta
        s.train <- s.train - nu_s * eta_optimal * u_s_hat.train
        # theta gradient descent
        u_theta <- scores$theta
        theta <- theta - nu_theta * u_theta
        theta <- min(max(theta, 0),9)

      }
    }
    return(list(s.beta=s.beta, theta=theta))
  }

}


# ----------------------------- SIMULATION 5.2.1 -----------------------------
# ---------------------- Orange Juice Price Elasticity -----------------------

# Import dataset, reformat, and transform relevant variables
library(IndexNumR)
{
  oj <- read.csv("https://msalicedatapublic.z5.web.core.windows.net/datasets/OrangeJuice/oj_large.csv")
  oj[,"week"] <- as.numeric(oj[,"week"])
  oj[,"logmove"] <- as.numeric(oj[,"logmove"])
  oj[,"price"] <- as.numeric(oj[,"price"])
  oj[,"INCOME"] <- as.numeric(oj[,"INCOME"])
  oj[,"HVAL150"] <- as.numeric(oj[,"HVAL150"])
  colnames(oj)[1] <- "id"
  oj_trop_logprice <- oj[oj$brand=="tropicana",]
  colnames(oj_trop_logprice)[3] <- "time"
  colnames(oj_trop_logprice)[6] <- "logprice"
  oj_trop_logprice$logprice <- log(oj_trop_logprice$logprice)
  oj_trop_logpriceTIME <- oj_trop_logprice
  oj_trop_logprice_XDY=oj_trop_logprice
  colnames(oj_trop_logprice_XDY)[c(3,4,6)]=c("X","Y","D")
}
# Diagnostic plots (check nuisance function fits and linearity assumption of the partially linear model)
par(mfrow=c(1,3))
# Nuisance function E[D|X]
plot(oj_trop_logprice$time, oj_trop_logprice$logmove, xlab="time", ylab="logmove")
lines(oj_trop_logprice$time[order(oj_trop_logprice$time)], gam(logmove~s(time,bs="cr"),data=oj_trop_logprice)$fitted[order(oj_trop_logprice$time)], col="red",lwd=2)
# Nuisance function E[Y|X]
plot(oj_trop_logprice$time, oj_trop_logprice$logprice, xlab="time", ylab="logprice")
lines(oj_trop_logprice$time[order(oj_trop_logprice$time)], gam(logprice~s(time,bs="cr"),data=oj_trop_logprice)$fitted[order(oj_trop_logprice$time)], col="red",lwd=2)
# Residual-Residual plot of R_Y versus R_D
plot(gam(logprice~s(time,bs="cr"),data=oj_trop_logprice)$residuals, gam(logmove~s(time,bs="cr"),data=oj_trop_logprice)$residuals)

  # Run Simulation 5.2.1
  set.seed(1)
  A=B=C=D=E=Ab=Bb=Cb=Db=Eb=numeric(0)
  K=5
  par(mfrow=c(1,K))
  while (length(A) < 50) {
    set.seed(1+length(A))
    # Note for the sandwich boosting method (WDML_HPLR_oj) we also plot the s-boosted function fold-by-fold (as in Figure 5a)
    a <- WDML_HPLR_oj(l_formula = logmove ~ s(time,bs="cr"), l_learner="gam", m_formula = logprice ~ s(time,bs="cr"), m_learner="gam", s_formula = u_s ~ s(time,bs="cr"), s_learner="gam", data=oj_trop_logprice, K=K, nu_s=0.1, nu_theta=0.1, m_stop=NULL) # sandwich boosting
    b <- WDML_HPLR_MEM_intercept(l_formula = logmove ~ s(time,bs="cr"), l_learner="gam", m_formula = logprice ~ s(time,bs="cr"), m_learner="gam", data=oj_trop_logprice, K=K) # intercept only mixed effects model (i.e. equicorrelated proxyCCF and maximum likelihood loss function)
    c <- WDML_HPLR_MEM_int_time(l_formula = logmove ~ s(time,bs="cr"), l_learner="gam", m_formula = logprice ~ s(time,bs="cr"), m_learner="gam", data=oj_trop_logprice, K=K) # intercept + time mixed effects model
    d <- WDML_HPLR_geepack(l_formula = Y ~ s(X,bs="cr"), l_learner="gam", m_formula = D ~ s(X,bs="cr"), m_learner="gam", s_formula = u_s ~ X, s_learner="randomforest", data=oj_trop_logprice_XDY, K=K) # heteroscedastic GEE
    e <- WDML_HPLR_geepack_hom(l_formula = Y ~ s(X,bs="cr"), l_learner="gam", m_formula = D ~ s(X,bs="cr"), m_learner="gam", s_formula = u_s ~ X, s_learner="randomforest", data=oj_trop_logprice_XDY, K=K) # homoscedastic GEE

    A = append(A, a$CI_width) # sandwich boosting variances
    B = append(B, b$CI_width) # intercept MEM variances
    C = append(C, c$CI_width) # intercept + time MEM variances
    D = append(D, d$CI_width) # het GEE variances
    E = append(E, e$CI_width) # hom GEE variances
    Ab = append(Ab, a$beta_hat) # sandwich boosting betas
    Bb = append(Bb, b$beta_hat) # intercept MEM betas
    Cb = append(Cb, c$beta_hat) # intercept + time MEM betas
    Db = append(Db, d$beta_hat) # het GEE betas
    Eb = append(Eb, e$beta_hat) # hom GEE betas

    print(paste0("Progress: ",length(A)," of 50 permutations completed"))
  }

  Boosting_Beta = median(Ab)
  IntMEM_Beta = median(Bb)
  IntTimeMEM_Beta = median(Cb)
  HetGEE_Beta = median(Db)
  HomGEE_Beta = median(Eb)

  Boosting_Var = median(A+(Ab-median(Ab))^2)
  IntMEM_Var = median(B+(Bb-median(Bb))^2)
  IntTimeMEM_Var = median(C+(Cb-median(Cb))^2)
  HetGEE_Var = median(D+(Db-median(Db))^2)
  HomGEE_Var = median(E+(Eb-median(Eb))^2)



# ----------------------------- SIMULATION 5.2.2 -----------------------------
# ------------ National Longitudinal Survey of Young Working Women ------------

# Results for Table 3
# Load data
  wagesID=read.csv("~/Downloads/wagesID.csv")
  # NB: Data available to download at: https://www.stata-press.com/data/r10/nlswork.dta
  wagesIDTIME = wagesID
  colnames(wagesIDTIME)[5] = "time"

  A=B=C=D=E=Ab=Bb=Cb=Db=Eb=numeric(0)
  K=5
  while (length(A) < 50) {
    set.seed(1+length(A))
#    a <- WDML_HPLR_wages(l_formula = ln_wage ~ s(age,bs="cr") + s(tenure,bs="cr"), l_learner="gam", m_formula = ttl_exp ~ s(age,bs="cr") + s(tenure,bs="cr"), m_learner="gam", s_formula = u_s ~ s(tenure,bs="cr"), s_learner="gam", data=wagesID, K=K, nu_s=0.1, nu_theta=0.1, m_stop=100)
    b <- WDML_HPLR_MEM_intercept(l_formula = ln_wage ~ s(age,bs="cr") + s(tenure,bs="cr"), l_learner="gam", m_formula = ttl_exp ~ s(age,bs="cr") + s(tenure,bs="cr"), m_learner="gam", data=wagesID, K=K)
    c <- WDML_HPLR_MEM_int_age_tenure_for_wages(l_formula = ln_wage ~ s(time,bs="cr") + s(age,bs="cr"), l_learner="gam", m_formula = ttl_exp ~ s(time,bs="cr") + s(age,bs="cr"), m_learner="gam", data=wagesIDTIME, K=K)
    d <- WDML_HPLR_geepack_for_wages(l_formula = ln_wage ~ s(age,bs="cr") + s(tenure,bs="cr"), l_learner="gam", m_formula = ttl_exp ~ s(age,bs="cr") + s(tenure,bs="cr"), m_learner="gam", data=wagesID, K=K)
    e <- WDML_HPLR_geepack_hom_wages(l_formula = ln_wage ~ s(age,bs="cr") + s(tenure,bs="cr"), l_learner="gam", m_formula = ttl_exp ~ s(age,bs="cr") + s(tenure,bs="cr"), m_learner="gam", data=wagesID, K=K)

    A = append(A, a$CI_width) # sandwich boosting variances
    B = append(B, b$CI_width) # intercept MEM variances
    C = append(C, c$CI_width) # intercept + tenure + age MEM variances
    D = append(D, d$CI_width) # het GEE variances
    E = append(E, e$CI_width) # hom GEE variances
    Ab = append(Ab, a$beta_hat) # sandwich boosting betas
    Bb = append(Bb, b$beta_hat) # intercept MEM betas
    Cb = append(Cb, c$beta_hat) # intercept + tenure + age MEM betas
    Db = append(Db, d$beta_hat) # het GEE betas
    Eb = append(Eb, e$beta_hat) # hom GEE betas

    print(paste0("Progress: ",length(A)," of 50 permutations completed"))
  }

  Boosting_Beta = median(Ab)
  IntMEM_Beta = median(Bb)
  IntTimeMEM_Beta = median(Cb)
  HetGEE_Beta = median(Db)
  HomGEE_Beta = median(Eb)

  Boosting_Var = median(A+(Ab-median(Ab))^2)
  IntMEM_Var = median(B+(Bb-median(Bb))^2)
  IntTimeMEM_Var = median(C+(Cb-median(Cb))^2)
  HetGEE_Var = median(D+(Db-median(Db))^2)
  HomGEE_Var = median(E+(Eb-median(Eb))^2)

# Plot of loss functions (Figure 5b)
  # Generate relevant loss function evaluations
{

  actual_theta <- c(seq(0,0.5,by=0.01), seq(0.5,2,by=0.02), seq(2,4,by=0.1), seq(4,9,by=0.25))
  V_hat_num <- beta_hat_den <- GEE1_theta <- MEM_int_theta <- rep(0,length(actual_theta))

  cv_folds <- groupKFold(wagesID$id,k=5)#Split into K folds (then can paralellise)
  l_formula = ln_wage ~ s(age,bs="cr") + s(tenure,bs="cr"); l_learner="gam"; m_formula = ttl_exp ~ s(age,bs="cr") + s(tenure,bs="cr"); m_learner="gam"; s_formula = u_s ~ s(tenure,bs="cr"); s_learner="gam"; K=5;
  for (k in 1:K) {
    cv_fold <- cv_folds[[k]]
    data.nuisance <- wagesID[cv_fold,]
    data.beta <- wagesID[-cv_fold,]

    fit.l <- regress.l(l_formula, l_learner, data.nuisance, data.beta)
    l_residuals.nuisance <- fit.l$l_residuals.nuisance; l_residuals.beta <- fit.l$l_residuals.beta

    fit.m <- regress.m(m_formula, m_learner, data.nuisance, data.beta)
    m_residuals.nuisance <- fit.m$m_residuals.nuisance; m_residuals.beta <- fit.m$m_residuals.beta

    xi_hat.nuisance <- m_residuals.nuisance
    beta_unweighted.nuisance <- sum(l_residuals.nuisance*xi_hat.nuisance)/sum(xi_hat.nuisance^2)
    epsilon_hat.nuisance <- l_residuals.nuisance - beta_unweighted.nuisance*xi_hat.nuisance

    xi_hat.beta <- m_residuals.beta
    beta_unweighted.beta <- sum(l_residuals.beta*xi_hat.beta)/sum(xi_hat.beta^2)
    epsilon_hat.beta <- l_residuals.beta - beta_unweighted.beta*xi_hat.beta

    boostdf.nuisance <- cbind(data.nuisance, Y_minus_l_hat=l_residuals.nuisance, epsilon_hat=epsilon_hat.nuisance , xi_hat=xi_hat.nuisance)
    boostdf.beta <- cbind(data.beta, Y_minus_l_hat=l_residuals.beta, epsilon_hat=epsilon_hat.beta , xi_hat=xi_hat.beta)

    tracker = 0
    for (THETA in actual_theta) {
      tracker = tracker + 1
      print(THETA)
      s.beta <- rep(1,dim(data.beta)[1]); theta <- THETA
      groups.beta <- unique(boostdf.beta[,"id"])
      I.beta <- length(groups.beta)
      beta_num = 0; beta_den = 0; V_num = 0
      for (i in seq_len(I.beta)) {
        Y_minus_l_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"Y_minus_l_hat"]
        xi_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"xi_hat"]
        epsilon_hat_i = boostdf.beta[is.element(data.beta$id,groups.beta[i]),"epsilon_hat"]
        s_i <- s.beta[is.element(boostdf.beta$id,groups.beta[i])]
        N <- length(s_i)
        beta_den = beta_den + sum(s_i^2*xi_hat_i^2) - (theta/(1+N*theta))*(sum(s_i*xi_hat_i))^2
        V_num = V_num + (sum(s_i^2*epsilon_hat_i*xi_hat_i) - (theta/(1+N*theta))*sum(s_i*epsilon_hat_i)*sum(s_i*xi_hat_i))^2
      }
      beta_hat_den[tracker] <- beta_hat_den[tracker] + beta_den
      V_hat_num[tracker] <- V_hat_num[tracker] + V_num

      # GEE
      sigma_hat_sq <- mean(boostdf.beta$epsilon_hat^2)
      groups.beta <- unique(boostdf.beta[,"id"])
      I.beta <- length(groups.beta)
      GEE1_crit <- 0
      for (i in seq_len(I.beta)) {
        s_i <- s.beta[is.element(boostdf.beta$id,groups.beta[i])]
        N <- length(s_i)
        epsilon_hat_i = boostdf.beta[is.element(data.beta$id,groups.beta[i]),"epsilon_hat"]
        if (N==1) {
          GEE1_crit <- GEE1_crit
        } else {
          epsilon_hat_i = as.vector(epsilon_hat_i)
          EPSILON_HAT_i = ( (epsilon_hat_i) %*% t(epsilon_hat_i) - diag(epsilon_hat_i^2) )/sigma_hat_sq
          THETA_MAT_i = matrix(THETA/(1+THETA), N, N) - diag(THETA/(1+THETA), N)
          GEE1_crit <- GEE1_crit + sum((EPSILON_HAT_i - THETA_MAT_i)^2)
        }
      }
      GEE1_theta[tracker] <- GEE1_theta[tracker] + GEE1_crit

      # MLE
      gls_fit <- gls(Y_minus_l_hat ~ xi_hat, data=boostdf.nuisance, correlation = corCompSymm(form = ~ 1 | id, value=THETA/(1+THETA), fixed="TRUE"))
      MEM_int_theta[tracker] <- MEM_int_theta[tracker] + summary(gls_fit)$logLik

    }
  }

  AV_theta <- (V_hat_num)/(beta_hat_den)^2

}

  # Altenatively, output is given below (uncomment)
{
#  AV_theta <- c(4.531145e-06, 4.298511e-06, 4.110758e-06, 3.956157e-06, 3.826833e-06, 3.717279e-06, 3.623512e-06, 3.542574e-06, 3.472218e-06, 3.410700e-06, 3.356647e-06, 3.308958e-06, 3.266742e-06, 3.229269e-06, 3.195932e-06, 3.166225e-06, 3.139720e-06, 3.116053e-06, 3.094914e-06, 3.076036e-06, 3.059184e-06, 3.044159e-06, 3.030782e-06, 3.018897e-06, 3.008368e-06, 2.999072e-06, 2.990901e-06, 2.983758e-06, 2.977557e-06, 2.972218e-06, 2.967672e-06, 2.963856e-06, 2.960711e-06, 2.958186e-06, 2.956234e-06, 2.954810e-06, 2.953875e-06, 2.953393e-06, 2.953332e-06, 2.953659e-06, 2.954348e-06, 2.955373e-06, 2.956709e-06, 2.958335e-06, 2.960231e-06, 2.962377e-06, 2.964756e-06, 2.967351e-06, 2.970149e-06, 2.973135e-06, 2.976295e-06, 2.976295e-06, 2.983092e-06, 2.990453e-06, 2.998302e-06, 3.006571e-06, 3.015201e-06, 3.024142e-06, 3.033348e-06, 3.042778e-06, 3.052398e-06, 3.062176e-06, 3.072084e-06, 3.082097e-06, 3.092194e-06, 3.102354e-06, 3.112561e-06, 3.122799e-06, 3.133054e-06, 3.143314e-06, 3.153568e-06, 3.163806e-06, 3.174020e-06, 3.184201e-06, 3.194343e-06, 3.204440e-06, 3.214486e-06, 3.224477e-06, 3.234408e-06, 3.244275e-06, 3.254076e-06, 3.263807e-06, 3.273466e-06, 3.283051e-06, 3.292560e-06, 3.301992e-06, 3.311344e-06, 3.320616e-06, 3.329807e-06, 3.338916e-06, 3.347943e-06, 3.356888e-06, 3.365749e-06, 3.374527e-06, 3.383222e-06, 3.391833e-06, 3.400361e-06, 3.408807e-06, 3.417170e-06, 3.425451e-06, 3.433650e-06, 3.441767e-06, 3.449804e-06, 3.457760e-06, 3.465637e-06, 3.473434e-06, 3.481153e-06, 3.488795e-06, 3.496359e-06, 3.503847e-06, 3.511259e-06, 3.518596e-06, 3.525859e-06, 3.533049e-06, 3.540166e-06, 3.547211e-06, 3.554185e-06, 3.561089e-06, 3.567924e-06, 3.574690e-06, 3.581388e-06, 3.588018e-06, 3.594583e-06, 3.601082e-06, 3.607516e-06, 3.613887e-06, 3.620194e-06, 3.620194e-06, 3.650809e-06, 3.679959e-06, 3.707731e-06, 3.734212e-06, 3.759482e-06, 3.783613e-06, 3.806677e-06, 3.828739e-06, 3.849857e-06, 3.870090e-06, 3.889487e-06, 3.908099e-06, 3.925970e-06, 3.943142e-06, 3.959654e-06, 3.975541e-06, 3.990839e-06, 4.005578e-06, 4.019787e-06, 4.033494e-06, 4.033494e-06, 4.065730e-06, 4.095341e-06, 4.122631e-06, 4.147858e-06, 4.171245e-06, 4.192986e-06, 4.213245e-06, 4.232169e-06, 4.249884e-06, 4.266502e-06, 4.282121e-06, 4.296828e-06, 4.310701e-06, 4.323807e-06, 4.336209e-06, 4.347962e-06, 4.359114e-06, 4.369712e-06, 4.379794e-06, 4.389399e-06 )
#  GEE1_theta <- c(373595.8, 371606.4, 369695.6, 367859.8, 366096.0, 364400.9, 362771.6, 361205.4, 359699.7, 358251.8,
#                  356859.4, 355520.2, 354232.1, 352992.9, 351800.7, 350653.7, 349550.0, 348487.9, 347465.8, 346482.1,
#                  345535.4, 344624.1, 343747.1, 342902.9, 342090.3, 341308.2, 340555.4, 339830.7, 339133.3, 338462.0,
#                  337816.0, 337194.2, 336595.8, 336020.1, 335466.1, 334933.1, 334420.4, 333927.2, 333452.9, 332996.7,
#                  332558.2, 332136.6, 331731.3, 331341.9, 330967.8, 330608.4, 330263.3, 329931.9, 329613.9, 329308.7,
#                  329016.0, 329016.0, 328466.3, 327961.6, 327498.9, 327075.8, 326689.5, 326337.8, 326018.5, 325729.7,
#                  325469.3, 325235.7, 325027.2, 324842.2, 324679.4, 324537.4, 324414.9, 324310.7, 324223.8, 324153.1,
#                  324097.7, 324056.6, 324028.9, 324014.0, 324010.9, 324019.1, 324037.8, 324066.4, 324104.4, 324151.1,
#                  324206.0, 324268.8, 324338.7, 324415.6, 324498.8, 324588.1, 324683.1, 324783.4, 324888.7, 324998.6,
#                  325113.0, 325231.4, 325353.7, 325479.6, 325608.9, 325741.3, 325876.6, 326014.7, 326155.3, 326298.2,
#                  326443.4, 326590.6, 326739.7, 326890.5, 327042.9, 327196.8, 327352.1, 327508.6, 327666.2, 327824.9,
#                  327984.5, 328144.9, 328306.0, 328467.9, 328630.3, 328793.2, 328956.5, 329120.2, 329284.3, 329448.5,
#                  329613.0, 329777.5, 329942.2, 330106.9, 330271.6, 330436.2, 330600.7, 330600.7, 331420.4, 332232.5,
#                  333033.9, 333822.1, 334595.6, 335353.0, 336093.6, 336816.7, 337522.3, 338210.0, 338880.1, 339532.6,
#                  340167.8, 340786.1, 341387.7, 341973.0, 342542.6, 343096.7, 343635.9, 344160.6, 344160.6, 345411.8,
#                  346582.0, 347677.6, 348704.7, 349668.9, 350575.3, 351428.6, 352233.1, 352992.5, 353710.4, 354389.9,
#                  355034.0, 355645.1, 356225.8, 356778.1, 357304.0, 357805.4, 358283.8, 358740.8, 359177.8)
#  MEM_int_theta <- c(-62391.36, -60517.08, -58898.20, -57480.85, -56226.42, -55106.31, -54098.69, -53186.56, -52356.37,
#                     -51597.22, -50900.13, -50257.70, -49663.69, -49112.86, -48600.69, -48123.35, -47677.48, -47260.17,
#                     -46868.89, -46501.38, -46155.67, -45830.00, -45522.78, -45232.62, -44958.26, -44698.54, -44452.46,
#                     -44219.06, -43997.50, -43787.02, -43586.89, -43396.48, -43215.19, -43042.47, -42877.83, -42720.79,
#                     -42570.93, -42427.85, -42291.18, -42160.59, -42035.74, -41916.35, -41802.13, -41692.83, -41588.20,
#                     -41488.03, -41392.09, -41300.19, -41212.14, -41127.76, -41046.90, -41046.90, -40895.07, -40755.52,
#                     -40627.21, -40509.25, -40400.79, -40301.11, -40209.53, -40125.46, -40048.33, -39977.66, -39912.99,
#                     -39853.90, -39800.01, -39750.98, -39706.48, -39666.23, -39629.94, -39597.37, -39568.29, -39542.49,
#                     -39519.76, -39499.92, -39482.81, -39468.26, -39456.12, -39446.26, -39438.54, -39432.86, -39429.08,
#                     -39427.12, -39426.87, -39428.24, -39431.14, -39435.49, -39441.22, -39448.24, -39456.50, -39465.94,
#                     -39476.48, -39488.07, -39500.67, -39514.21, -39528.66, -39543.96, -39560.07, -39576.96, -39594.58,
#                     -39612.90, -39631.88, -39651.49, -39671.70, -39692.49, -39713.81, -39735.65, -39757.98, -39780.77,
#                     -39804.01, -39827.67, -39851.72, -39876.16, -39900.96, -39926.10, -39951.57, -39977.35, -40003.42,
#                     -40029.77, -40056.39, -40083.26, -40110.36, -40137.70, -40165.24, -40192.99, -40220.93, -40249.05,
#                     -40277.34, -40277.34, -40421.05, -40567.77, -40716.64, -40866.98, -41018.23, -41169.91, -41321.66,
#                     -41473.16, -41624.18, -41774.49, -41923.93, -42072.36, -42219.67, -42365.78, -42510.61, -42654.10,
#                     -42796.22, -42936.92, -43076.20, -43214.02, -43214.02, -43552.18, -43881.23, -44201.27, -44512.52,
#                     -44815.24, -45109.71, -45396.26, -45675.19, -45946.81, -46211.44, -46469.36, -46720.86, -46966.23,
#                     -47205.71, -47439.55, -47668.00, -47891.28, -48109.60, -48323.16, -48532.15)
}

  # Plot loss functions against theta
  library(RColorBrewer)
  coul <- brewer.pal(3, "Dark2")
  plot(actual_theta[actual_theta<2], AV_theta[actual_theta<2]/min(AV_theta), col=coul[2], lwd=2, type="l", xlab=expression(theta), ylab=expression(paste("Objective function for ", theta)))
  lines(actual_theta[actual_theta<2], (GEE1_theta[actual_theta<2]/min(GEE1_theta)-0.7)/min((GEE1_theta[actual_theta<2]/min(GEE1_theta)-0.7)), col=coul[1], lwd=2)
  lines(actual_theta[actual_theta<2], (-MEM_int_theta[actual_theta<2]/min(-MEM_int_theta)-0.6)/min((-MEM_int_theta[actual_theta<2]/min(-MEM_int_theta)-0.6)), col=coul[3], lwd=2)
  legend("topright", legend=c("ML", "GEE", "Sandwich Loss"), col=c(coul[3], coul[1], coul[2]), lty=1, cex=1, lwd=2)

  # Plot loss functions against rho (=theta/(1+theta))
  plot(actual_theta/(1+actual_theta), AV_theta/min(AV_theta), col=coul[2], lwd=2, type="l", xlab=expression(rho), ylab=expression(paste("Objective function for ", rho)))
  lines(actual_theta/(1+actual_theta), (GEE1_theta/min(GEE1_theta)-0.7)/min((GEE1_theta/min(GEE1_theta)-0.7)), col=coul[1], lwd=2)
  lines(actual_theta/(1+actual_theta), (-MEM_int_theta/min(-MEM_int_theta)-0.6)/min((-MEM_int_theta/min(-MEM_int_theta)-0.6)), col=coul[3], lwd=2)
  legend("topright", legend=c("ML", "GEE", "Sandwich Loss"), col=c(coul[3], coul[1], coul[2]), lty=1, cex=1, lwd=2)
