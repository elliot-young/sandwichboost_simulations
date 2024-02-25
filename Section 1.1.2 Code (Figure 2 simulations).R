# Code for Figure 2
library(RColorBrewer)
coul <- brewer.pal(3, "Dark2")

# 1.1.2: Example 2: Conditional Variance Misspecification
LambdaMAX=10
cairo_pdf("mu_both.pdf", 12, 6)
par(mfrow=c(1,2))
{
  mu <- 0.2

  sech <- function(x) 1/cosh(x)

  a <- 2

  Q <- seq(0, LambdaMAX, by = 0.01)
  Q <- Q[-1]

  AV_vs_unweighted <- rep(0,length(Q))
  MLE_vs_unweighted <- rep(0,length(Q))
  Rsq_vs_unweighted <- rep(0,length(Q))
  ratio <- rep(0,length(Q))
  unweighted <- rep(0,length(Q))

  for (qplace in seq_len(length(Q))) {
    qq <- Q[qplace]
    integral_sig <- function(l,u) {
      int <- 5*(u-l) - 1/qq*(tanh(qq*(u-mu))-tanh(qq*(l-mu))) + 4/qq*(log(cosh(qq*(u-mu))) - log(cosh(qq*(l-mu))))
      return(int)
    }

    gamma <- seq(0,1,0.001)
    lmle <- 2*(log(a+1))*(1-gamma) + integral_sig(0,gamma) + (1/(1+a)^2)*integral_sig(gamma,1)   ##
    gamma_mle <- gamma[which.min(lmle)]
    lAV_mle <- (1/(gamma_mle+(1-gamma_mle)/(1+a)^2)^2)  *  (integral_sig(0,gamma_mle) + (1/(1+a)^4)*integral_sig(gamma_mle,1)) ##
    MLE_vs_unweighted[qplace] <- lAV_mle/integral_sig(0,1)

    lav <- (1/(gamma+(1-gamma)/(1+a)^2)^2)  *  (integral_sig(0,gamma) + (1/(1+a)^4)*integral_sig(gamma,1))
    gamma_AV <- gamma[which.min(lav)]
    lAV_AV <- (1/(gamma_AV+(1-gamma_AV)/(1+a)^2)^2)  *  (integral_sig(0,gamma_AV) + (1/(1+a)^4)*integral_sig(gamma_AV,1)) ##
    AV_vs_unweighted[qplace] <- lAV_AV/integral_sig(0,1)

    ratio[qplace] <- lAV_mle/lAV_AV

    Rsq_int_indef <- function(x,a) {
      int <- (1/qq)*(-8*log(sech(qq*x-qq*mu))+(-8*a^2-16*a+24)*log(cosh(qq*x-qq*mu))+(2*a^2+4*a-23)*tanh(qq*x-qq*mu)+4*(sech(qq*x-qq*mu))^2) - (1/(3*qq))*(tanh(qq*x-qq*mu))^3 + (a^4+4*a^3-4*a^2-16*a+32)*x
      return(int)
    }
    Rsq_int_def <- function(l,u,a) {
      int <- Rsq_int_indef(u,a) - Rsq_int_indef(l,a)
      return(int)
    }
    gamma <- seq(0,1,0.001)
    lRsq <- Rsq_int_def(0,gamma,0) + Rsq_int_def(gamma,1,2)   ##
    gamma_Rsq <- gamma[which.min(lRsq)]
    lAV_Rsq <- (1/(gamma_Rsq+(1-gamma_Rsq)/(1+a)^2)^2)  *  (integral_sig(0,gamma_Rsq) + (1/(1+a)^4)*integral_sig(gamma_Rsq,1)) ##
    Rsq_vs_unweighted[qplace] <- lAV_Rsq/integral_sig(0,1)
    unweighted[qplace] <- integral_sig(0,1)
  }

  dashed_line <- rep(1,length(Q))-1

  dashed_line=dashed_line+1

  par(mar=c(5.1,6.1,4.1,2.1))
  plot(Q,dashed_line,col="white",type="l",lwd=3,lty=2, xlab=expression(lambda),
       ylab="Mean Squared Error", ylim=c(1-0.4,1+0.4),
       main=bquote(paste(mu,' = ',.(mu))), cex.lab=1.2)
  points(Q,MLE_vs_unweighted,col=coul[3],type="l",lwd=3)
  points(Q,Rsq_vs_unweighted,col=coul[1],type="l",lwd=3)
  abline(h=1, lwd=3, lty=2)
  legend("topright", legend=c("Unweighted", "ML", "GEE"),
         col=c("black", coul[3], coul[1]), lty=c(2, 1, 1), lwd=3, cex=1)
}
{

  mu <- 0.5

  sech <- function(x) 1/cosh(x)

  a <- 2

  Q <- seq(0, LambdaMAX, by = 0.01)
  Q <- Q[-1]

  AV_vs_unweighted <- rep(0,length(Q))
  MLE_vs_unweighted <- rep(0,length(Q))
  Rsq_vs_unweighted <- rep(0,length(Q))
  ratio <- rep(0,length(Q))
  unweighted <- rep(0,length(Q))

  for (qplace in seq_len(length(Q))) {
    qq <- Q[qplace]
    integral_sig <- function(l,u) {
      int <- 5*(u-l) - 1/qq*(tanh(qq*(u-mu))-tanh(qq*(l-mu))) + 4/qq*(log(cosh(qq*(u-mu))) - log(cosh(qq*(l-mu))))
      return(int)
    }

    gamma <- seq(0,1,0.001)
    lmle <- 2*(log(a+1))*(1-gamma) + integral_sig(0,gamma) + (1/(1+a)^2)*integral_sig(gamma,1)   ##
    gamma_mle <- gamma[which.min(lmle)]
    lAV_mle <- (1/(gamma_mle+(1-gamma_mle)/(1+a)^2)^2)  *  (integral_sig(0,gamma_mle) + (1/(1+a)^4)*integral_sig(gamma_mle,1)) ##
    MLE_vs_unweighted[qplace] <- lAV_mle/integral_sig(0,1)

    lav <- (1/(gamma+(1-gamma)/(1+a)^2)^2)  *  (integral_sig(0,gamma) + (1/(1+a)^4)*integral_sig(gamma,1))
    gamma_AV <- gamma[which.min(lav)]
    lAV_AV <- (1/(gamma_AV+(1-gamma_AV)/(1+a)^2)^2)  *  (integral_sig(0,gamma_AV) + (1/(1+a)^4)*integral_sig(gamma_AV,1)) ##
    AV_vs_unweighted[qplace] <- lAV_AV/integral_sig(0,1)

    ratio[qplace] <- lAV_mle/lAV_AV



    # Rsq method added
    Rsq_int_indef <- function(x,a) {
      int <- (1/qq)*(-8*log(sech(qq*x-qq*mu))+(-8*a^2-16*a+24)*log(cosh(qq*x-qq*mu))+(2*a^2+4*a-23)*tanh(qq*x-qq*mu)+4*(sech(qq*x-qq*mu))^2) - (1/(3*qq))*(tanh(qq*x-qq*mu))^3 + (a^4+4*a^3-4*a^2-16*a+32)*x
      return(int)
    }
    Rsq_int_def <- function(l,u,a) {
      int <- Rsq_int_indef(u,a) - Rsq_int_indef(l,a)
      return(int)
    }

    gamma <- seq(0,1,0.001)
    lRsq <- Rsq_int_def(0,gamma,0) + Rsq_int_def(gamma,1,2)   ##
    gamma_Rsq <- gamma[which.min(lRsq)]
    lAV_Rsq <- (1/(gamma_Rsq+(1-gamma_Rsq)/(1+a)^2)^2)  *  (integral_sig(0,gamma_Rsq) + (1/(1+a)^4)*integral_sig(gamma_Rsq,1)) ##
    Rsq_vs_unweighted[qplace] <- lAV_Rsq/integral_sig(0,1)
    #Rsq method end
    unweighted[qplace] <- integral_sig(0,1)
  }

  dashed_line <- rep(1,length(Q))-1


  par(mar=c(5.1,6.1,4.1,2.1))
  plot(Q,dashed_line,col="white",type="l",lwd=3,lty=2, xlab=expression(lambda),
       ylab="Mean Squared Error", ylim=c(1-0.6,1+0.8),
       main=bquote(paste(mu,' = ',.(mu))), cex.lab=1.2)
  points(Q,MLE_vs_unweighted,col=coul[3],type="l",lwd=3)
  points(Q,Rsq_vs_unweighted,col=coul[1],type="l",lwd=3)
  abline(h=1, lwd=3, lty=2)
  legend("topright", legend=c("Unweighted", "ML", "GEE"),
         col=c("black", coul[3], coul[1]), lty=c(2, 1, 1,1), lwd=3, cex=1)
}
dev.off()
