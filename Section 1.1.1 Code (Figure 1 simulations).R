# Code for Figure 1

# Sort out colours
library(RColorBrewer)
coul <- brewer.pal(3, "Dark2")
t_col <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  invisible(t.col)
}
mycol.red <- t_col(coul[3], perc = 30, name = "lt.red")
mycol.green <- t_col(coul[1], perc = 30, name = "lt.green")

cairo_pdf("correl.pdf", 12, 6)

# 1.1.1: Example 1: Conditional Correlation Misspecification

# Setting (a)
# Pre-run results
{
  actual_theta = c(seq(-0.05,0.99,by=0.01),0.999)
  AV_theta = c(1.226515e-05, 1.201818e-05, 1.177592e-05, 1.153862e-05,
               1.130653e-05, 1.107988e-05, 1.085890e-05, 1.064381e-05,
               1.043480e-05, 1.023207e-05, 1.003578e-05, 9.846095e-06,
               9.663162e-06, 9.487107e-06, 9.318042e-06, 9.156063e-06,
               9.001251e-06, 8.853668e-06, 8.713361e-06, 8.580358e-06,
               8.454674e-06, 8.336303e-06, 8.225226e-06, 8.121406e-06,
               8.024791e-06, 7.935313e-06, 7.852890e-06, 7.777424e-06,
               7.708805e-06, 7.646908e-06, 7.591598e-06, 7.542724e-06,
               7.500129e-06, 7.463642e-06, 7.433083e-06, 7.408265e-06,
               7.388991e-06, 7.375060e-06, 7.366262e-06, 7.362383e-06,
               7.363205e-06, 7.368506e-06, 7.378060e-06, 7.391642e-06,
               7.409023e-06, 7.429976e-06, 7.454272e-06, 7.481686e-06,
               7.511993e-06, 7.544970e-06, 7.580397e-06, 7.618060e-06,
               7.657746e-06, 7.699248e-06, 7.742364e-06, 7.786898e-06,
               7.832658e-06, 7.879460e-06, 7.927125e-06, 7.975482e-06,
               8.024364e-06, 8.073614e-06, 8.123080e-06, 8.172619e-06,
               8.222093e-06, 8.271374e-06, 8.320338e-06, 8.368871e-06,
               8.416865e-06, 8.464217e-06, 8.510836e-06, 8.556632e-06,
               8.601526e-06, 8.645444e-06, 8.688319e-06, 8.730088e-06,
               8.770697e-06, 8.810096e-06, 8.848241e-06, 8.885093e-06,
               8.920620e-06, 8.954792e-06, 8.987586e-06, 9.018984e-06,
               9.048969e-06, 9.077531e-06, 9.104663e-06, 9.130361e-06,
               9.154627e-06, 9.177461e-06, 9.198872e-06, 9.218867e-06,
               9.237459e-06, 9.254660e-06, 9.270488e-06, 9.284960e-06,
               9.298096e-06, 9.309917e-06, 9.320447e-06, 9.329711e-06,
               9.337734e-06, 9.344542e-06, 9.350165e-06, 9.354629e-06,
               9.357965e-06, 9.360026e-06)
  GEE1_theta = c(20840668, 20840334, 20839987, 20839625, 20839247, 20838853,
                 20838442, 20838011, 20837561, 20837091, 20836599, 20836084,
                 20835546, 20834984, 20834396, 20833782, 20833141, 20832471,
                 20831773, 20831044, 20830284, 20829492, 20828666, 20827806,
                 20826911, 20825980, 20825011, 20824003, 20822956, 20821867,
                 20820737, 20819563, 20818344, 20817079, 20815767, 20814406,
                 20812995, 20811533, 20810017, 20808447, 20806820, 20805136,
                 20803392, 20801587, 20799719, 20797786, 20795786, 20793718,
                 20791579, 20789368, 20787082, 20784719, 20782278, 20779756,
                 20777151, 20774461, 20771684, 20768818, 20765861, 20762811,
                 20759667, 20756425, 20753087, 20749649, 20746112, 20742475,
                 20738738, 20734903, 20730972, 20726947, 20722832, 20718634,
                 20714361, 20710023, 20705634, 20701213, 20696781, 20692368,
                 20688011, 20683755, 20679660, 20675799, 20672268, 20669186,
                 20666707, 20665028, 20664404, 20665164, 20667744, 20672716,
                 20680845, 20693168, 20711104, 20736632, 20772559, 20822960,
                 20893895, 20994666, 21140077, 21354776, 21682140, 22203994,
                 23088718, 24721953, 28105341, 35012462)
  MEM_theta = c(-284505.9, -284317.2, -284148.0, -283998.2, -283868.1,
                -283757.5, -283666.7, -283595.7, -283544.5, -283513.1,
                -283501.5, -283509.9, -283538.1, -283586.1, -283654.0,
                -283741.7, -283849.1, -283976.2, -284122.9, -284289.1,
                -284474.7, -284679.7, -284903.9, -285147.1, -285409.4,
                -285690.4, -285990.1, -286308.3, -286644.8, -286999.4,
                -287372.0, -287762.3, -288170.3, -288595.5, -289037.9,
                -289497.2, -289973.2, -290465.6, -290974.2, -291498.8,
                -292039.1, -292594.9, -293165.8, -293751.7, -294352.2,
                -294967.2, -295596.3, -296239.3, -296895.9, -297565.9,
                -298248.9, -298944.6, -299652.9, -300373.5, -301106.0,
                -301850.2, -302605.9, -303372.7, -304150.4, -304938.8,
                -305737.6, -306546.5, -307365.2, -308193.6, -309031.4,
                -309878.3, -310734.1, -311598.6, -312471.5, -313352.7,
                -314241.9, -315138.9, -316043.6, -316955.6, -317874.9,
                -318801.3, -319734.6, -320674.7, -321621.3, -322574.5,
                -323534.1, -324500.0, -325472.1, -326450.4, -327434.8,
                -328425.4, -329422.2, -330425.3, -331434.8, -332450.9,
                -333473.9, -334504.1, -335542.1, -336588.4, -337644.1,
                -338710.1, -339788.2, -340880.6, -341990.4, -343122.4,
                -344283.9, -345487.5, -346756.7, -348144.7, -349821.1,
                -353008.4)
}
# Code to generate data above
{
# {
#   library(expm)
#   library(nlme)
#   cross_fit_eval <- function(boostdf.beta=boostdf.beta, s.beta=s.beta, theta=theta, proxyCCF=proxyCCF) {
#     if (proxyCCF == "equicorr") {
#       groups.beta <- unique(boostdf.beta[,"id"])
#       I.beta <- length(groups.beta)
#       beta_num = 0; beta_den = 0; V_num = 0
#       for (i in seq_len(I.beta)) {
#         Y_minus_l_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"Y_minus_l_hat"]
#         xi_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"xi_hat"]
#         epsilon_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"epsilon_hat"]
#         s_i <- s.beta[is.element(boostdf.beta$id,groups.beta[i])]
#         N <- length(s_i)
#         A1 <- sum(s_i^2*xi_hat_i^2) - (theta/(1+N*theta))*(sum(s_i*xi_hat_i))^2
#         A2_Y_minus_l_hat <- sum(s_i^2*Y_minus_l_hat_i*xi_hat_i) - (theta/(1+N*theta))*sum(s_i*Y_minus_l_hat_i)*sum(s_i*xi_hat_i)
#         A2 <- sum(s_i^2*epsilon_hat_i*xi_hat_i) - (theta/(1+N*theta))*sum(s_i*epsilon_hat_i)*sum(s_i*xi_hat_i)
#         beta_num = beta_num + A2_Y_minus_l_hat
#         beta_den = beta_den + A1
#         V_num = V_num + A2^2
#       }
#     } else if (proxyCCF == "autoreg") {
#       groups.beta <- unique(boostdf.beta[,"id"])
#       I.beta <- length(groups.beta)
#       beta_num = 0; beta_den = 0; V_num = 0
#       for (i in seq_len(I.beta)) {
#         Y_minus_l_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"Y_minus_l_hat"]
#         xi_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"xi_hat"]
#         epsilon_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"epsilon_hat"]
#         s_i <- s.beta[is.element(boostdf.beta$id,groups.beta[i])]
#         N <- length(s_i)
#         if (N>2) {
#           beta1 <- sum((s_i[2:(N-1)]*xi_hat_i[2:(N-1)])^2)
#           alpha1 <- beta1 + (s_i[1]*xi_hat_i[1])^2 + (s_i[N]*xi_hat_i[N])^2
#           gamma1 <- 2*sum(s_i[2:N]*s_i[1:(N-1)]*xi_hat_i[2:N]*xi_hat_i[1:(N-1)])
#           beta2_Y_minus_l_hat <- sum(s_i[2:(N-1)]^2*xi_hat_i[2:(N-1)]*Y_minus_l_hat_i[2:(N-1)])
#           alpha2_Y_minus_l_hat <- beta2_Y_minus_l_hat + s_i[1]^2*xi_hat_i[1]*Y_minus_l_hat_i[1] + s_i[N]^2*xi_hat_i[N]*Y_minus_l_hat_i[N]
#           gamma2_Y_minus_l_hat <- sum(s_i[2:N]*s_i[1:(N-1)]*(xi_hat_i[2:N]*Y_minus_l_hat_i[1:(N-1)] + Y_minus_l_hat_i[2:N]*xi_hat_i[1:(N-1)]))
#           beta2 <- sum(s_i[2:(N-1)]^2*xi_hat_i[2:(N-1)]*epsilon_hat_i[2:(N-1)])
#           alpha2 <- beta2 + s_i[1]^2*xi_hat_i[1]*epsilon_hat_i[1] + s_i[N]^2*xi_hat_i[N]*epsilon_hat_i[N]
#           gamma2 <- sum(s_i[2:N]*s_i[1:(N-1)]*(xi_hat_i[2:N]*epsilon_hat_i[1:(N-1)] + epsilon_hat_i[2:N]*xi_hat_i[1:(N-1)]))
#           A1 <- alpha1 + theta^2*beta1 - theta*gamma1
#           A2_Y_minus_l_hat <- alpha2_Y_minus_l_hat + theta^2*beta2_Y_minus_l_hat - theta*gamma2_Y_minus_l_hat
#           A2 <- alpha2 + theta^2*beta2 - theta*gamma2
#         } else if (N==2) {
#           alpha1 <- (s_i[1]*xi_hat_i[1])^2 + (s_i[N]*xi_hat_i[N])^2
#           gamma1 <- 2*sum(s_i[N]*s_i[N-1]*xi_hat_i[N]*xi_hat_i[N-1])
#           alpha2_Y_minus_l_hat <- s_i[1]^2*xi_hat_i[1]*Y_minus_l_hat_i[1] + s_i[N]^2*xi_hat_i[N]*Y_minus_l_hat_i[N]
#           gamma2_Y_minus_l_hat <- sum(s_i[N]*s_i[N-1]*(xi_hat_i[N]*Y_minus_l_hat_i[N-1] + Y_minus_l_hat_i[N]*xi_hat_i[N-1]))
#           alpha2 <- s_i[1]^2*xi_hat_i[1]*epsilon_hat_i[1] + s_i[N]^2*xi_hat_i[N]*epsilon_hat_i[N]
#           gamma2 <- sum(s_i[N]*s_i[N-1]*(xi_hat_i[N]*epsilon_hat_i[N-1] + epsilon_hat_i[N]*xi_hat_i[N-1]))
#           A1 <- alpha1 - theta*gamma1
#           A2_Y_minus_l_hat <- alpha2_Y_minus_l_hat - theta*gamma2_Y_minus_l_hat
#           A2 <- alpha2 - theta*gamma2
#         } else if (N==1) {
#           alpha1 <- (s_i[1]*xi_hat_i[1])^2 + (s_i[N]*xi_hat_i[N])^2
#           alpha2_Y_minus_l_hat <- s_i[1]^2*xi_hat_i[1]*Y_minus_l_hat_i[1] + s_i[N]^2*xi_hat_i[N]*Y_minus_l_hat_i[N]
#           alpha2 <- s_i[1]^2*xi_hat_i[1]*epsilon_hat_i[1] + s_i[N]^2*xi_hat_i[N]*epsilon_hat_i[N]
#           A1 <- alpha1
#           A2_Y_minus_l_hat <- alpha2_Y_minus_l_hat
#           A2 <- alpha2
#         }
#         beta_num = beta_num + A2_Y_minus_l_hat
#         beta_den = beta_den + A1
#         V_num = V_num + A2^2
#       }
#     } else if (proxyCCF == "nested"){
#       groups.beta <- unique(boostdf.beta[,"id"])
#       I.beta <- length(groups.beta)
#       beta_num = 0; beta_den = 0; V_num = 0
#       for (i in seq_len(I.beta)) {
#         sub.boostdf.beta <- boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),]
#         sub.s.beta <- s.beta[is.element(boostdf.beta$id,groups.beta[i])]
#         sub.groups.beta <- unique(sub.boostdf.beta[,"subid"])
#         J.beta <- length(sub.groups.beta)
#         b1 <- b2 <- b3 <- c1 <- c2 <- c3 <- n_ij <- rep(0,J.beta)
#         for (j in seq_len(J.beta)) {
#           Y_minus_l_hat_ij = sub.boostdf.beta[is.element(sub.boostdf.beta$subid,sub.groups.beta[j]),"Y_minus_l_hat"]
#           xi_hat_ij = sub.boostdf.beta[is.element(sub.boostdf.beta$subid,sub.groups.beta[j]),"xi_hat"]
#           epsilon_hat_ij = sub.boostdf.beta[is.element(sub.boostdf.beta$subid,sub.groups.beta[j]),"epsilon_hat"]
#           s_ij = sub.s.beta[is.element(sub.boostdf.beta$subid,sub.groups.beta[j])]
#           b1[j] = sum(s_ij*xi_hat_ij)
#           b2[j] = sum(s_ij*epsilon_hat_ij)
#           b3[j] = sum(s_ij*Y_minus_l_hat_ij)
#           n_ij[j] = length(s_ij)
#           c1[j] = sum(s_ij^2*xi_hat_ij^2)
#           c2[j] = sum(s_ij^2*xi_hat_ij*epsilon_hat_ij)
#           c3[j] = sum(s_ij^2*xi_hat_ij*Y_minus_l_hat_ij)
#         }
#         phi <- sum(n_ij/(1+theta[1]*n_ij))
#         A1 <- (1+theta[1]+theta[2])*sum(c1) - sum(theta[1]/(1+theta[1]*n_ij)*b1^2) - theta[2]/(1+theta[2]*phi)*(sum(b1/(1+theta[1]*n_ij)))^2
#         A2 <- (1+theta[1]+theta[2])*sum(c2) - sum(theta[1]/(1+theta[1]*n_ij)*b1*b2) - theta[2]/(1+theta[2]*phi)*sum(b1/(1+theta[1]*n_ij))*sum(b2/(1+theta[1]*n_ij))
#         A2_Y_minus_l_hat <- (1+theta[1]+theta[2])*sum(c3) - sum(theta[1]/(1+theta[1]*n_ij)*b1*b3) - theta[2]/(1+theta[2]*phi)*sum(b1/(1+theta[1]*n_ij))*sum(b3/(1+theta[1]*n_ij))
#         beta_num = beta_num + A2_Y_minus_l_hat
#         beta_den = beta_den + A1
#         V_num = V_num + A2^2
#       }
#     }
#     return(list(beta_num=beta_num, beta_den=beta_den, V_num=V_num))
#   }
#
#   set.seed(1)
#   N = 100
#   Sigma0 = toeplitz(ARMAacf(ar=c(0.3,0.6), ma=c(-0.5), lag.max = N-1)) # GEE1/AV = 1.5
#   sqrtSigma0 = sqrtm(Sigma0)
#   Omega0 = matrix(1/8,N,N) + 7/8*diag(rep(1,N))
#   sqrtOmega0 = sqrtm(Omega0)
#
#   {
#     I = 2*10^3
#     n_obs = N*I
#
#     V = as.matrix(rnorm(n_obs,0,1))#V = as.matrix(rlaplace(n_obs,0,1))#V = as.matrix(rnorm(n_obs,0,1))
#     U = as.matrix(rnorm(n_obs,0,1))#U = as.matrix(rlaplace(n_obs,0,1))#U = as.matrix(rnorm(n_obs,0,1))
#     epsilon = as.matrix(rep(0,n_obs))
#     xi = as.matrix(rep(0,n_obs))
#
#     for (i in seq_len(I)) {
#       xi[(N*(i-1)+1):(N*i),] = (sqrtOmega0)%*%V[(N*(i-1)+1):(N*i),]
#     }
#
#     D = xi
#     for (i in seq_len(I)) {
#       epsilon[(N*(i-1)+1):(N*i),] = (sqrtSigma0)%*%U[(N*(i-1)+1):(N*i),]
#     }
#
#     Y = beta*D + epsilon
#
#     rdf <- data.frame(D=D, Y=Y, id=rep(1:I, each=N))
#   }
#   V_hat_num_k1 <- numeric(0)
#   beta_hat_den_k1 <- numeric(0)
#   GEE1_theta_k1 <- numeric(0)
#   MEM_int_theta_k1 <- numeric(0)
#   #actual_theta = c(seq(0,0.99,by=0.01),0.999)
#   actual_theta = c(seq(-0.05,0.99,by=0.01),0.999)
#   { k=1
#
#     for (THETA in actual_theta) {
#       print(THETA)
#       s.beta <- rep(1,dim(rdf)[1]); theta <- THETA
#       groups.beta <- unique(rdf[,"id"])
#       I.beta <- length(groups.beta)
#       rdf = cbind(rdf, Y_minus_l_hat=rdf$Y, epsilon_hat=rdf$Y-rdf$D, xi_hat=rdf$D)
#       data.beta = data.nuisance = boostdf.beta = boostdf.nuisance = rdf
#       cross_fit_evals <- cross_fit_eval(boostdf.beta=rdf, s.beta=s.beta, theta=THETA, proxyCCF="autoreg")
#       beta_hat_den_k1 <- append(beta_hat_den_k1, cross_fit_evals$beta_den)
#       V_hat_num_k1 <- append(V_hat_num_k1, cross_fit_evals$V_num)
#
#
#       groups.beta <- unique(boostdf.beta[,"id"])
#       I.beta <- length(groups.beta)
#       GEE1_crit <- 0
#       sigma_hat_sq = mean((rdf$Y - rdf$D)^2)
#       for (i in seq_len(I.beta)) {
#         s_i <- s.beta[is.element(boostdf.beta$id,groups.beta[i])]
#         N <- length(s_i)
#
#         epsilon_hat_i = boostdf.beta[is.element(data.beta$id,groups.beta[i]),"epsilon_hat"]
#
#         if (N==1) {
#           GEE1_crit <- GEE1_crit
#         } else {
#           epsilon_hat_i = as.vector(epsilon_hat_i)
#           EPSILON_HAT_i = ( (epsilon_hat_i) %*% t(epsilon_hat_i) - diag(epsilon_hat_i^2) )/sigma_hat_sq
#           THETA_MAT_i = toeplitz(ARMAacf(ar=c(THETA), ma=c(0), lag.max = N-1)) - diag(1, N)
#           GEE1_crit <- GEE1_crit + sum((EPSILON_HAT_i - THETA_MAT_i)^2)
#         }
#
#       }
#       GEE1_theta_k1 <- append(GEE1_theta_k1, GEE1_crit)
#
#
#       gls_fit <- gls(Y ~ D, data=rdf, correlation = corAR1(form = ~ 1 | id, value=THETA, fixed="TRUE"))
#       MEM_int_theta_k1 <- append(MEM_int_theta_k1, summary(gls_fit)$logLik)
#
#     }
#
#   }
#   AV_theta = V_hat_num_k1/(beta_hat_den_k1)^2
#   GEE1_theta = GEE1_theta_k1
#   MEM_theta = MEM_int_theta_k1
#
# }
}
# Plots
# NB: the objectives are scaled for clear visuals
{
  par(mfrow=c(1,2))

  plot(actual_theta, AV_theta/min(AV_theta), col=coul[2], lwd=4, type="l",
       xlab=expression(rho), ylab=expression(paste("Objective function for ", rho)),
       main="Setting (a)", ylim=c(1.0-0.02,1.5), xlim=c(0,1), cex=1.2)
  lines(actual_theta, (GEE1_theta/min(GEE1_theta)-0.95)/min((GEE1_theta/min(GEE1_theta)-0.95)), col=coul[1], lwd=2)
  lines(actual_theta, (-MEM_theta/min(-MEM_theta)-0.9)/min(-MEM_theta/min(-MEM_theta)-0.9), col=coul[3], lwd=2)
  lines(actual_theta, AV_theta/min(AV_theta), col=coul[2], lwd=4, type="l")
  legend("topright", legend=c("ML", "GEE", "Asymptotic Variance (SL)"), col=c(coul[3], coul[1], coul[2]), lty=1, cex=1.2, lwd=3)

  RHO_GEE1 = actual_theta[which.min(GEE1_theta)]
  RHO_MEM = actual_theta[which.min(-MEM_theta)]
  RHO_SL = actual_theta[which.min(AV_theta)]

  AV_GEE1_sc = AV_theta[which.min(GEE1_theta)]/min(AV_theta)
  AV_MEM_sc = AV_theta[which.min(-MEM_theta)]/min(AV_theta)

  lines(c(RHO_MEM,RHO_MEM),c(1.000,AV_MEM_sc),lty="dashed", lwd=1, col=coul[3])
  lines(c(RHO_GEE1,RHO_GEE1),c(1.000,AV_GEE1_sc),lty="dashed", lwd=1, col=coul[1])

  AV_all_sc = AV_theta/min(AV_theta)

  #MEM
  AV_all_sc<AV_MEM_sc
  mem_opp_val = 106
  actual_theta[mem_opp_val]
  AV_all_sc[mem_opp_val]
  lines(c(RHO_MEM,actual_theta[mem_opp_val]),c(AV_MEM_sc,AV_MEM_sc),lty="dashed", lwd=1, col=coul[3])

  #GEE
  AV_all_sc<AV_GEE1_sc
  gee1_opp_val = 16
  actual_theta[gee1_opp_val]
  AV_all_sc[gee1_opp_val]
  lines(c(actual_theta[gee1_opp_val],RHO_GEE1),c(AV_GEE1_sc,AV_GEE1_sc),lty="dashed", lwd=1, col=coul[1])
  lines(c(actual_theta[gee1_opp_val],actual_theta[gee1_opp_val]),c(1.000,AV_GEE1_sc),lty="dashed", lwd=1, col=coul[1])

  t_col <- function(color, percent = 50, name = NULL) {
    ## Get RGB values for named color
    rgb.val <- col2rgb(color)
    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 alpha = (100 - percent) * 255 / 100,
                 names = name)
    ## Save the color
    invisible(t.col)
  }
  mycol.red <- t_col(coul[3], perc = 30, name = "lt.red")
  mycol.green <- t_col(coul[1], perc = 30, name = "lt.green")
  lines(c(actual_theta[mem_opp_val],RHO_MEM),c(0.99,0.99), lwd=10, col=mycol.red, lend=1)

  lines(c(actual_theta[gee1_opp_val],RHO_GEE1),c(0.97,0.97), lwd=10, col=mycol.green, lend=1)

  lines(actual_theta, AV_theta/min(AV_theta), col=coul[2], lwd=4, type="l")

}

# Setting (b) (& Table 4)
# Pre-run results
{
  actual_theta = c(-0.999,seq(-0.99,0.99,by=0.01),0.999)
  AV_theta = c(2.314655e-05, 2.313931e-05, 2.312951e-05, 2.311779e-05, 2.310411e-05, 2.308842e-05, 2.307066e-05, 2.305078e-05, 2.302873e-05, 2.300446e-05, 2.297789e-05, 2.294900e-05, 2.291770e-05, 2.288396e-05, 2.284772e-05, 2.280891e-05, 2.276749e-05, 2.272339e-05, 2.267656e-05, 2.262694e-05, 2.257449e-05, 2.251914e-05, 2.246083e-05, 2.239952e-05, 2.233515e-05, 2.226766e-05, 2.219701e-05, 2.212315e-05, 2.204602e-05, 2.196557e-05, 2.188176e-05, 2.179455e-05, 2.170388e-05, 2.160973e-05, 2.151204e-05, 2.141079e-05, 2.130593e-05, 2.119743e-05, 2.108528e-05, 2.096943e-05, 2.084987e-05, 2.072658e-05, 2.059954e-05, 2.046874e-05, 2.033417e-05, 2.019584e-05, 2.005373e-05, 1.990787e-05, 1.975825e-05, 1.960490e-05, 1.944783e-05, 1.928708e-05, 1.912268e-05, 1.895466e-05, 1.878307e-05, 1.860797e-05, 1.842942e-05, 1.824747e-05, 1.806220e-05, 1.787370e-05, 1.768204e-05, 1.748732e-05, 1.728965e-05, 1.708912e-05, 1.688587e-05, 1.668000e-05, 1.647166e-05, 1.626098e-05, 1.604811e-05, 1.583320e-05, 1.561642e-05, 1.539793e-05, 1.517791e-05, 1.495654e-05, 1.473401e-05, 1.451052e-05, 1.428626e-05, 1.406145e-05, 1.383630e-05, 1.361103e-05, 1.338585e-05, 1.316101e-05, 1.293671e-05, 1.271321e-05, 1.249073e-05, 1.226952e-05, 1.204981e-05, 1.183184e-05, 1.161585e-05, 1.140209e-05, 1.119079e-05, 1.098219e-05, 1.077652e-05, 1.057401e-05, 1.037488e-05, 1.017937e-05, 9.987684e-06, 9.800031e-06, 9.616615e-06, 9.437631e-06, 9.263266e-06, 9.093697e-06, 8.929093e-06, 8.769611e-06, 8.615399e-06, 8.466593e-06, 8.323318e-06, 8.185687e-06, 8.053799e-06, 7.927743e-06, 7.807593e-06, 7.693411e-06, 7.585245e-06, 7.483129e-06, 7.387086e-06, 7.297123e-06, 7.213235e-06, 7.135404e-06, 7.063598e-06, 6.997773e-06, 6.937873e-06, 6.883828e-06, 6.835558e-06, 6.792971e-06, 6.755965e-06, 6.724427e-06, 6.698234e-06, 6.677254e-06, 6.661346e-06, 6.650363e-06, 6.644149e-06, 6.642541e-06, 6.645370e-06, 6.652463e-06, 6.663643e-06, 6.678726e-06, 6.697526e-06, 6.719857e-06, 6.745528e-06, 6.774347e-06, 6.806122e-06, 6.840662e-06, 6.877775e-06, 6.917269e-06, 6.958958e-06, 7.002652e-06, 7.048169e-06, 7.095328e-06, 7.143950e-06, 7.193861e-06, 7.244893e-06, 7.296879e-06, 7.349660e-06, 7.403081e-06, 7.456990e-06, 7.511245e-06, 7.565705e-06, 7.620238e-06, 7.674716e-06, 7.729018e-06, 7.783030e-06, 7.836640e-06, 7.889748e-06, 7.942254e-06, 7.994068e-06, 8.045105e-06, 8.095284e-06, 8.144533e-06, 8.192784e-06, 8.239973e-06, 8.286044e-06, 8.330944e-06, 8.374627e-06, 8.417052e-06, 8.458180e-06, 8.497980e-06, 8.536424e-06, 8.573487e-06, 8.609150e-06, 8.643397e-06, 8.676217e-06, 8.707599e-06, 8.737540e-06, 8.766038e-06, 8.793092e-06, 8.818706e-06, 8.842888e-06, 8.865645e-06, 8.886988e-06, 8.906932e-06, 8.925490e-06, 8.942680e-06, 8.958520e-06, 8.973031e-06, 8.986233e-06, 8.998151e-06, 9.008807e-06, 9.018227e-06, 9.026436e-06, 9.033461e-06, 9.038794e-06)
  GEE1_theta = c(13846332, 12865021, 12046304, 11434649, 10972632, 10619909, 10347887,
                 10136106, 9969784, 9838132, 9733195, 9649048, 9581230, 9526352,
                 9481812, 9445593, 9416119, 9392148, 9372691, 9356957, 9344307,
                 9334224, 9326284, 9320138, 9315500, 9312130, 9309829, 9308431,
                 9307795, 9307802, 9308353, 9309362, 9310757, 9312475, 9314464,
                 9316677, 9319077, 9321628, 9324302, 9327075, 9329923, 9332828,
                 9335775, 9338748, 9341736, 9344728, 9347715, 9350688, 9353642,
                 9356569, 9359465, 9362325, 9365146, 9367923, 9370655, 9373339,
                 9375972, 9378553, 9381081, 9383553, 9385970, 9388331, 9390634,
                 9392879, 9395065, 9397193, 9399262, 9401272, 9403223, 9405115,
                 9406947, 9408721, 9410435, 9412090, 9413686, 9415224, 9416703,
                 9418123, 9419486, 9420790, 9422036, 9423224, 9424354, 9425427,
                 9426442, 9427400, 9428300, 9429142, 9429928, 9430656, 9431327,
                 9431940, 9432495, 9432993, 9433433, 9433816, 9434140, 9434405,
                 9434612, 9434760, 9434849, 9434878, 9434847, 9434756, 9434603,
                 9434389, 9434113, 9433774, 9433371, 9432905, 9432373, 9431775,
                 9431111, 9430380, 9429579, 9428709, 9427768, 9426755, 9425669,
                 9424509, 9423272, 9421958, 9420564, 9419090, 9417534, 9415893,
                 9414166, 9412350, 9410444, 9408446, 9406352, 9404160, 9401869,
                 9399474, 9396973, 9394363, 9391641, 9388803, 9385846, 9382766,
                 9379558, 9376220, 9372745, 9369131, 9365371, 9361461, 9357395,
                 9353168, 9348773, 9344205, 9339457, 9334521, 9329391, 9324058,
                 9318515, 9312754, 9306765, 9300539, 9294066, 9287337, 9280340,
                 9273066, 9265503, 9257639, 9249462, 9240961, 9232125, 9222940,
                 9213397, 9203485, 9193194, 9182517, 9171450, 9159992, 9148145,
                 9135920, 9123338, 9110427, 9097237, 9083836, 9070322, 9056830,
                 9043548, 9030733, 9018736, 9008035, 8999281, 8993362, 8991494,
                 8995342, 9007204, 9030257, 9068923, 9129391, 9220367, 9354180,
                 9548393, 9828174, 10229814, 10805930, 11536353)
  MEM_theta = c(-443182.5, -432948.0, -429038.0, -426279.6, -423993.3, -421966.3,
                -420103.8, -418355.3, -416690.9, -415091.1, -413542.8, -412036.7,
                -410566.0, -409125.5, -407711.1, -406319.9, -404949.4, -403597.5,
                -402262.7, -400943.7, -399639.3, -398348.8, -397071.3, -395806.3,
                -394553.2, -393311.8, -392081.6, -390862.4, -389654.0, -388456.3,
                -387269.1, -386092.3, -384926.0, -383770.1, -382624.7, -381489.8,
                -380365.4, -379251.7, -378148.7, -377056.7, -375975.6, -374905.7,
                -373847.2, -372800.2, -371764.9, -370741.5, -369730.2, -368731.3,
                -367744.9, -366771.3, -365810.7, -364863.4, -363929.6, -363009.6,
                -362103.6, -361212.0, -360334.9, -359472.6, -358625.5, -357793.7,
                -356977.7, -356177.5, -355393.6, -354626.3, -353875.6, -353142.1,
                -352425.9, -351727.3, -351046.7, -350384.1, -349740.1, -349114.7,
                -348508.3, -347921.1, -347353.4, -346805.4, -346277.4, -345769.6,
                -345282.3, -344815.5, -344369.7, -343945.0, -343541.5, -343159.6,
                -342799.3, -342460.8, -342144.3, -341850.0, -341578.1, -341328.5,
                -341101.6, -340897.3, -340715.8, -340557.2, -340421.5, -340308.9,
                -340219.4, -340153.0, -340109.7, -340089.7, -340092.8, -340119.2,
                -340168.7, -340241.4, -340337.2, -340456.1, -340598.0, -340762.8,
                -340950.5, -341160.9, -341394.0, -341649.6, -341927.5, -342227.8,
                -342550.1, -342894.4, -343260.5, -343648.1, -344057.2, -344487.5,
                -344938.8, -345410.8, -345903.5, -346416.5, -346949.6, -347502.6,
                -348075.2, -348667.2, -349278.3, -349908.3, -350556.8, -351223.7,
                -351908.7, -352611.4, -353331.6, -354069.0, -354823.4, -355594.5,
                -356381.9, -357185.4, -358004.7, -358839.6, -359689.7, -360554.8,
                -361434.6, -362328.8, -363237.2, -364159.4, -365095.3, -366044.6,
                -367007.0, -367982.2, -368970.1, -369970.3, -370982.7, -372007.0,
                -373043.1, -374090.7, -375149.6, -376219.6, -377300.7, -378392.5,
                -379495.0, -380608.0, -381731.4, -382865.1, -384009.0, -385163.1,
                -386327.3, -387501.6, -388686.0, -389880.6, -391085.4, -392300.5,
                -393526.1, -394762.5, -396009.7, -397268.3, -398538.6, -399821.1,
                -401116.4, -402425.3, -403748.5, -405087.3, -406442.9, -407817.0,
                -409211.4, -410628.7, -412072.1, -413545.4, -415053.8, -416604.2,
                -418206.0, -419872.2, -421622.1, -423486.0, -425514.2, -427801.4,
                -430560.4, -434470.7, -444704.6)
}
# Code to generate data above
{
# {
#   library(expm)
#   library(nlme)
#   cross_fit_eval <- function(boostdf.beta=boostdf.beta, s.beta=s.beta, theta=theta, proxyCCF=proxyCCF) {
#     if (proxyCCF == "equicorr") {
#       groups.beta <- unique(boostdf.beta[,"id"])
#       I.beta <- length(groups.beta)
#       beta_num = 0; beta_den = 0; V_num = 0
#       for (i in seq_len(I.beta)) {
#         Y_minus_l_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"Y_minus_l_hat"]
#         xi_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"xi_hat"]
#         epsilon_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"epsilon_hat"]
#         s_i <- s.beta[is.element(boostdf.beta$id,groups.beta[i])]
#         N <- length(s_i)
#         A1 <- sum(s_i^2*xi_hat_i^2) - (theta/(1+N*theta))*(sum(s_i*xi_hat_i))^2
#         A2_Y_minus_l_hat <- sum(s_i^2*Y_minus_l_hat_i*xi_hat_i) - (theta/(1+N*theta))*sum(s_i*Y_minus_l_hat_i)*sum(s_i*xi_hat_i)
#         A2 <- sum(s_i^2*epsilon_hat_i*xi_hat_i) - (theta/(1+N*theta))*sum(s_i*epsilon_hat_i)*sum(s_i*xi_hat_i)
#         beta_num = beta_num + A2_Y_minus_l_hat
#         beta_den = beta_den + A1
#         V_num = V_num + A2^2
#       }
#     } else if (proxyCCF == "autoreg") {
#       groups.beta <- unique(boostdf.beta[,"id"])
#       I.beta <- length(groups.beta)
#       beta_num = 0; beta_den = 0; V_num = 0
#       for (i in seq_len(I.beta)) {
#         Y_minus_l_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"Y_minus_l_hat"]
#         xi_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"xi_hat"]
#         epsilon_hat_i = boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),"epsilon_hat"]
#         s_i <- s.beta[is.element(boostdf.beta$id,groups.beta[i])]
#         N <- length(s_i)
#         if (N>2) {
#           beta1 <- sum((s_i[2:(N-1)]*xi_hat_i[2:(N-1)])^2)
#           alpha1 <- beta1 + (s_i[1]*xi_hat_i[1])^2 + (s_i[N]*xi_hat_i[N])^2
#           gamma1 <- 2*sum(s_i[2:N]*s_i[1:(N-1)]*xi_hat_i[2:N]*xi_hat_i[1:(N-1)])
#           beta2_Y_minus_l_hat <- sum(s_i[2:(N-1)]^2*xi_hat_i[2:(N-1)]*Y_minus_l_hat_i[2:(N-1)])
#           alpha2_Y_minus_l_hat <- beta2_Y_minus_l_hat + s_i[1]^2*xi_hat_i[1]*Y_minus_l_hat_i[1] + s_i[N]^2*xi_hat_i[N]*Y_minus_l_hat_i[N]
#           gamma2_Y_minus_l_hat <- sum(s_i[2:N]*s_i[1:(N-1)]*(xi_hat_i[2:N]*Y_minus_l_hat_i[1:(N-1)] + Y_minus_l_hat_i[2:N]*xi_hat_i[1:(N-1)]))
#           beta2 <- sum(s_i[2:(N-1)]^2*xi_hat_i[2:(N-1)]*epsilon_hat_i[2:(N-1)])
#           alpha2 <- beta2 + s_i[1]^2*xi_hat_i[1]*epsilon_hat_i[1] + s_i[N]^2*xi_hat_i[N]*epsilon_hat_i[N]
#           gamma2 <- sum(s_i[2:N]*s_i[1:(N-1)]*(xi_hat_i[2:N]*epsilon_hat_i[1:(N-1)] + epsilon_hat_i[2:N]*xi_hat_i[1:(N-1)]))
#           A1 <- alpha1 + theta^2*beta1 - theta*gamma1
#           A2_Y_minus_l_hat <- alpha2_Y_minus_l_hat + theta^2*beta2_Y_minus_l_hat - theta*gamma2_Y_minus_l_hat
#           A2 <- alpha2 + theta^2*beta2 - theta*gamma2
#         } else if (N==2) {
#           alpha1 <- (s_i[1]*xi_hat_i[1])^2 + (s_i[N]*xi_hat_i[N])^2
#           gamma1 <- 2*sum(s_i[N]*s_i[N-1]*xi_hat_i[N]*xi_hat_i[N-1])
#           alpha2_Y_minus_l_hat <- s_i[1]^2*xi_hat_i[1]*Y_minus_l_hat_i[1] + s_i[N]^2*xi_hat_i[N]*Y_minus_l_hat_i[N]
#           gamma2_Y_minus_l_hat <- sum(s_i[N]*s_i[N-1]*(xi_hat_i[N]*Y_minus_l_hat_i[N-1] + Y_minus_l_hat_i[N]*xi_hat_i[N-1]))
#           alpha2 <- s_i[1]^2*xi_hat_i[1]*epsilon_hat_i[1] + s_i[N]^2*xi_hat_i[N]*epsilon_hat_i[N]
#           gamma2 <- sum(s_i[N]*s_i[N-1]*(xi_hat_i[N]*epsilon_hat_i[N-1] + epsilon_hat_i[N]*xi_hat_i[N-1]))
#           A1 <- alpha1 - theta*gamma1
#           A2_Y_minus_l_hat <- alpha2_Y_minus_l_hat - theta*gamma2_Y_minus_l_hat
#           A2 <- alpha2 - theta*gamma2
#         } else if (N==1) {
#           alpha1 <- (s_i[1]*xi_hat_i[1])^2 + (s_i[N]*xi_hat_i[N])^2
#           alpha2_Y_minus_l_hat <- s_i[1]^2*xi_hat_i[1]*Y_minus_l_hat_i[1] + s_i[N]^2*xi_hat_i[N]*Y_minus_l_hat_i[N]
#           alpha2 <- s_i[1]^2*xi_hat_i[1]*epsilon_hat_i[1] + s_i[N]^2*xi_hat_i[N]*epsilon_hat_i[N]
#           A1 <- alpha1
#           A2_Y_minus_l_hat <- alpha2_Y_minus_l_hat
#           A2 <- alpha2
#         }
#         beta_num = beta_num + A2_Y_minus_l_hat
#         beta_den = beta_den + A1
#         V_num = V_num + A2^2
#       }
#     } else if (proxyCCF == "nested"){
#       groups.beta <- unique(boostdf.beta[,"id"])
#       I.beta <- length(groups.beta)
#       beta_num = 0; beta_den = 0; V_num = 0
#       for (i in seq_len(I.beta)) {
#         sub.boostdf.beta <- boostdf.beta[is.element(boostdf.beta$id,groups.beta[i]),]
#         sub.s.beta <- s.beta[is.element(boostdf.beta$id,groups.beta[i])]
#         sub.groups.beta <- unique(sub.boostdf.beta[,"subid"])
#         J.beta <- length(sub.groups.beta)
#         b1 <- b2 <- b3 <- c1 <- c2 <- c3 <- n_ij <- rep(0,J.beta)
#         for (j in seq_len(J.beta)) {
#           Y_minus_l_hat_ij = sub.boostdf.beta[is.element(sub.boostdf.beta$subid,sub.groups.beta[j]),"Y_minus_l_hat"]
#           xi_hat_ij = sub.boostdf.beta[is.element(sub.boostdf.beta$subid,sub.groups.beta[j]),"xi_hat"]
#           epsilon_hat_ij = sub.boostdf.beta[is.element(sub.boostdf.beta$subid,sub.groups.beta[j]),"epsilon_hat"]
#           s_ij = sub.s.beta[is.element(sub.boostdf.beta$subid,sub.groups.beta[j])]
#           b1[j] = sum(s_ij*xi_hat_ij)
#           b2[j] = sum(s_ij*epsilon_hat_ij)
#           b3[j] = sum(s_ij*Y_minus_l_hat_ij)
#           n_ij[j] = length(s_ij)
#           c1[j] = sum(s_ij^2*xi_hat_ij^2)
#           c2[j] = sum(s_ij^2*xi_hat_ij*epsilon_hat_ij)
#           c3[j] = sum(s_ij^2*xi_hat_ij*Y_minus_l_hat_ij)
#         }
#         phi <- sum(n_ij/(1+theta[1]*n_ij))
#         A1 <- (1+theta[1]+theta[2])*sum(c1) - sum(theta[1]/(1+theta[1]*n_ij)*b1^2) - theta[2]/(1+theta[2]*phi)*(sum(b1/(1+theta[1]*n_ij)))^2
#         A2 <- (1+theta[1]+theta[2])*sum(c2) - sum(theta[1]/(1+theta[1]*n_ij)*b1*b2) - theta[2]/(1+theta[2]*phi)*sum(b1/(1+theta[1]*n_ij))*sum(b2/(1+theta[1]*n_ij))
#         A2_Y_minus_l_hat <- (1+theta[1]+theta[2])*sum(c3) - sum(theta[1]/(1+theta[1]*n_ij)*b1*b3) - theta[2]/(1+theta[2]*phi)*sum(b1/(1+theta[1]*n_ij))*sum(b3/(1+theta[1]*n_ij))
#         beta_num = beta_num + A2_Y_minus_l_hat
#         beta_den = beta_den + A1
#         V_num = V_num + A2^2
#       }
#     }
#     return(list(beta_num=beta_num, beta_den=beta_den, V_num=V_num))
#   }
#
#   set.seed(1)
#   N = 30
#   Sigma0 = toeplitz(ARMAacf(ar=c(0.1,0.85), ma=c(-0.4), lag.max = N-1)) # GEE1/AV = 1.5
#   sqrtSigma0 = sqrtm(Sigma0)
#   #Omega0 = matrix(1/8,N,N) + 7/8*diag(rep(1,N))
#   Omega0 = matrix(1/8,N,N) + 7/8*diag(rep(1,N))
#   sqrtOmega0 = sqrtm(Omega0)
#
#   {
#     I = 10^3*8
#     n_obs = N*I
#
#     V = as.matrix(rnorm(n_obs,0,1))#V = as.matrix(rlaplace(n_obs,0,1))#V = as.matrix(rnorm(n_obs,0,1))
#     U = as.matrix(rnorm(n_obs,0,1))#U = as.matrix(rlaplace(n_obs,0,1))#U = as.matrix(rnorm(n_obs,0,1))
#     epsilon = as.matrix(rep(0,n_obs))
#     xi = as.matrix(rep(0,n_obs))
#
#     for (i in seq_len(I)) {
#       xi[(N*(i-1)+1):(N*i),] = (sqrtOmega0)%*%V[(N*(i-1)+1):(N*i),]
#     }
#
#     D = xi
#     for (i in seq_len(I)) {
#       epsilon[(N*(i-1)+1):(N*i),] = (sqrtSigma0)%*%U[(N*(i-1)+1):(N*i),]
#     }
#
#     Y = beta*D + epsilon
#
#     rdf <- data.frame(D=D, Y=Y, id=rep(1:I, each=N))
#   }
#   V_hat_num_k1 <- numeric(0)
#   beta_hat_den_k1 <- numeric(0)
#   GEE1_theta_k1 <- numeric(0)
#   MEM_int_theta_k1 <- numeric(0)
#   #actual_theta = c(seq(0,0.99,by=0.01),0.999)
#   actual_theta = c(-0.999,seq(-0.99,0.99,by=0.01),0.999)
#   #actual_theta = c(-0.999,-0.99,-0.9,-0.6,0,seq(0.2,0.5,by=0.01),0.7,0.9)
#   { k=1
#
#     for (THETA in actual_theta) {
#       print(THETA)
#       s.beta <- rep(1,dim(rdf)[1]); theta <- THETA
#       groups.beta <- unique(rdf[,"id"])
#       I.beta <- length(groups.beta)
#       rdf = cbind(rdf, Y_minus_l_hat=rdf$Y, epsilon_hat=rdf$Y-rdf$D, xi_hat=rdf$D)
#       data.beta = data.nuisance = boostdf.beta = boostdf.nuisance = rdf
#       cross_fit_evals <- cross_fit_eval(boostdf.beta=rdf, s.beta=s.beta, theta=THETA, proxyCCF="autoreg")
#       beta_hat_den_k1 <- append(beta_hat_den_k1, cross_fit_evals$beta_den)
#       V_hat_num_k1 <- append(V_hat_num_k1, cross_fit_evals$V_num)
#
#
#       groups.beta <- unique(boostdf.beta[,"id"])
#       I.beta <- length(groups.beta)
#       GEE1_crit <- 0
#       sigma_hat_sq = mean((rdf$Y - rdf$D)^2)
#       for (i in seq_len(I.beta)) {
#         s_i <- s.beta[is.element(boostdf.beta$id,groups.beta[i])]
#         N <- length(s_i)
#
#         epsilon_hat_i = boostdf.beta[is.element(data.beta$id,groups.beta[i]),"epsilon_hat"]
#
#         if (N==1) {
#           GEE1_crit <- GEE1_crit
#         } else {
#           epsilon_hat_i = as.vector(epsilon_hat_i)
#           EPSILON_HAT_i = ( (epsilon_hat_i) %*% t(epsilon_hat_i) - diag(epsilon_hat_i^2) )/sigma_hat_sq
#           THETA_MAT_i = toeplitz(ARMAacf(ar=c(THETA), ma=c(0), lag.max = N-1)) - diag(1, N)
#           GEE1_crit <- GEE1_crit + sum((EPSILON_HAT_i - THETA_MAT_i)^2)
#         }
#
#       }
#       GEE1_theta_k1 <- append(GEE1_theta_k1, GEE1_crit)
#
#
#       gls_fit <- gls(Y ~ D, data=rdf, correlation = corAR1(form = ~ 1 | id, value=THETA, fixed="TRUE"))
#       MEM_int_theta_k1 <- append(MEM_int_theta_k1, summary(gls_fit)$logLik)
#
#     }
#
#   }
#   AV_theta = V_hat_num_k1/(beta_hat_den_k1)^2
#   GEE1_theta = GEE1_theta_k1
#   MEM_theta = MEM_int_theta_k1
#
# }
}
# Plots
# NB: the objectives are scaled for clear visuals
{
  plot(actual_theta, AV_theta/min(AV_theta), col=coul[2], lwd=4, type="l", xlab=expression(rho),
       ylab=expression(paste("Objective function for ", rho)), ylim=c(1.0-0.1,3.7),
       main="Setting (b)", cex=1.2)
  lines(actual_theta, (GEE1_theta/min(GEE1_theta)-0.95)/min((GEE1_theta/min(GEE1_theta)-0.95)), col=coul[1], lwd=2)
  lines(actual_theta, (-MEM_theta/min(-MEM_theta)-0.9)/min(-MEM_theta/min(-MEM_theta)-0.9), col=coul[3], lwd=2)
  lines(actual_theta, AV_theta/min(AV_theta), col=coul[2], lwd=4, type="l")
  legend("topright", legend=c("ML", "GEE", "Asymptotic variance (SL)"), col=c(coul[3], coul[1], coul[2]), lty=1, cex=1.2, lwd=3)

  RHO_GEE1 = actual_theta[which.min(GEE1_theta)]
  RHO_MEM = actual_theta[which.min(-MEM_theta)]
  RHO_SL = actual_theta[which.min(AV_theta)]

  AV_GEE1_sc = AV_theta[which.min(GEE1_theta)]/min(AV_theta)
  AV_MEM_sc = AV_theta[which.min(-MEM_theta)]/min(AV_theta)

  lines(c(RHO_MEM,RHO_MEM),c(1.000,AV_MEM_sc),lty="dashed", lwd=1, col=coul[3])
  lines(c(RHO_GEE1,RHO_GEE1),c(1.000,AV_GEE1_sc),lty="dashed", lwd=1, col=coul[1])

  AV_all_sc = AV_theta/min(AV_theta)

  #MEM
  AV_all_sc<AV_MEM_sc
  mem_opp_val = length(actual_theta)
  actual_theta[mem_opp_val]
  AV_all_sc[mem_opp_val]
  lines(c(RHO_MEM,actual_theta[mem_opp_val]),c(AV_MEM_sc,AV_MEM_sc),lty="dashed", lwd=1, col=coul[3])

  #MEM
  AV_all_sc<AV_GEE1_sc
  gee1_opp_val = 103
  actual_theta[gee1_opp_val]
  AV_all_sc[gee1_opp_val]
  lines(c(actual_theta[gee1_opp_val],RHO_GEE1),c(AV_GEE1_sc,AV_GEE1_sc),lty="dashed", lwd=1, col=coul[1])
  lines(c(actual_theta[gee1_opp_val],actual_theta[gee1_opp_val]),c(1.000,AV_GEE1_sc),lty="dashed", lwd=1, col=coul[1])

  lines(c(actual_theta[mem_opp_val],RHO_MEM),c(0.95,0.95), lwd=10, col=mycol.red, lend=1)
  lines(c(actual_theta[gee1_opp_val],RHO_GEE1),c(0.85,0.85), lwd=10, col=mycol.green, lend=1)

  lines(actual_theta, AV_theta/min(AV_theta), col=coul[2], lwd=4, type="l")
}

dev.off()
