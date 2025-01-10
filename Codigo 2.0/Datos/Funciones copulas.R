
#hi=BiCopHinv(u[,1], u[,2], family=5, 2.32, obj = NULL, check.pars = TRUE)
#head(hi$hinv1)
library(copula)
library(MASS)
library(mvtnorm)

#Function for performing (nonlinear) quantile regression
#by means of a fitted bivariate copula
qrCopula <- function(copula, p=.5, marginal, x) {
  #copula<-fitCop2
  mvdcQR<-copula
  
  #copula function
  if (class(mvdcQR@mvdc@copula)[1]=="normalCopula") {
    
    biCopFun <-expression(cdf=integrate(function(t) integrate(function(s)(1/(2*pi*sqrt(1-(rho_param^2))))*exp((2*rho_param*s*t-(s^2)-(t^2))/(2*(1-(rho_param^2)))),lower=-Inf, upper = qnorm(u2) ), lower = -Inf, upper = qnorm(u1)))
    rho_param <- mvdcQR@mvdc@copula@parameters
  
  } else if (class(mvdcQR@mvdc@copula)[1]=="frankCopula" ){
    biCopFun <-expression((-1/alpha)*log(1+((exp(-theta*u1)-1)*(exp(-theta*u2)-1)/(exp(-theta1)-1))))
    alpha <- mvdcQR@mvdc@copula@parameters
  } else if (class(mvdcQR@mvdc@copula)[1]=="tCopula") {
    #exps=(1/(2*pi*sqrt(1-(ro^2))))*((1+((s^2)+(t^2)+2*ro*s*t)/(v*(1-(ro^2))))^(-(v+2)/2))
    
    biCopFun <- expression(cdf=integrate(function(t) integrate(function(s)(1/(2*pi*sqrt(1-(rho^2))))*((1+((s^2)+(t^2)+2*rho*s*t)/(v*(1-(rho^2))))^(-(v+2)/2)),lower=-Inf, upper = qt(u2,df=v) ), lower = -Inf, upper = qt(u1, df=v)) )
    rho.1 <- mvdcQR@mvdc@copula@parameters[1]
    df <- mvdcQR@mvdc@copula@parameters[2]
    
  }else if (class(mvdcQR@mvdc@copula)[1]=="r90GumbelCopula") {
  
    biCopFun <-expression(exp(-(((-log(u1))^alpha)+((-log(u2))^alpha))^(1/alpha)))
    alpha <- mvdcQR@mvdc@copula@parameters
    
  }else if (class(mvdcQR@mvdc@copula)[1]=="surBB1Copula") {
   #u1=x, u   y    u2=y, v
    biCopFun <-expression(cdf=(1+(((((u1^(-theta))-1)^delta)+(((u2^(-theta))-1)^delta))^(1/delta)))^-(1/theta))
    theta <- mvdcQR@mvdc@copula@parameters[1]
    delta <- mvdcQR@mvdc@copula@parameters[2]
    
  }else {biCopFun <- mvdcQR@mvdc@copula@exprdist[1]}
  
  
  library(Deriv)
  #install.packages("Deriv")
  # Calcular la derivada parcial de biCopFun con respecto a u1
  #if (class(mvdcQR@mvdc@copula)[1]!="tCopula") {
   # derivative_u1 <- D(biCopFun$cdf, "u1") 
  #}
  
  
  
  #marginal=2
  i <- marginal #dependent variable, y
  j <- ifelse(i==1, 2, 1) #independent variable, x
  
  
  #construct cdf and quantile function for marginal distributions
  cdf.margin <- copula:::asCall(paste0("p", mvdcQR@mvdc@margins[j]),
                                mvdcQR@mvdc@paramMargins[[j]])
  qdf.margin <- copula:::asCall(paste0("q", mvdcQR@mvdc@margins[i]),
                                mvdcQR@mvdc@paramMargins[[i]])
  #x=xpoints1$tci
  u1 <- eval(cdf.margin, list(x=x))
  #p=.5
  # 
  # #compute u2
 if (class(mvdcQR@mvdc@copula)[1]=="normalCopula"){
   u2= numeric(length(u1))
   for (i in seq_along(u1)) {
     u2[i]<- cCopula(c(u1[i], p), copula=mvdcQR@mvdc@copula, inverse=FALSE)[2]
   }
   #u2 <- pnorm(alpha*qnorm(u1) + sqrt(1-alpha^2)*qnorm(p))
 }else if(class(mvdcQR@mvdc@copula)[1]=="frankCopula"){
   u2= numeric(length(u1))
   for (i in seq_along(u1)) {
     u2[i]<- cCopula(c(u1[i], p), copula=mvdcQR@mvdc@copula, inverse=FALSE)[2]
   }
 }else if (class(mvdcQR@mvdc@copula)[1]=="tCopula"){
   u2= numeric(length(u1))
   for (i in seq_along(u1)) {
     u2[i]<- cCopula(c(u1[i], p), copula=mvdcQR@mvdc@copula, inverse=FALSE)[2]
   }
   
   #u2 <- cCopula(c(u1, p), copula=mvdcQR@mvdc@copula, inverse=TRUE)[2]
 }else if (class(mvdcQR@mvdc@copula)[1]=="surBB1Copula"){
   u2fun <- function (u2, u1, theta,delta, t) {t - eval(derivative_u1,list(u1 = u1, u2 = u2, theta = theta, delta = delta))}
   u2= numeric(length(u1))
   for (i in seq_along(u1)) {
     u2[i]<- uniroot(u2fun, interval=c(1-.99999999999999, .99999999999999),              ###########################################
                   u1=u1[i], delta=delta, theta=theta, t=p, tol=1e-5)$root
     }
   # try(u2 <- uniroot(u2fun, interval=c(1-.99999999999999, .99999999999999),              ###########################################
   #                   u1=u1, delta=delta, theta=theta, t=p, tol=1e-5)$root, silent=FALSE)
 
 }else {
  u2 <- NA
  u2fun <- function (u2, u1, alpha, t) {t - attr(eval(dC.du1),"gradient")}
  try(u2 <- uniroot(u2fun, interval=c(1-.99999999999999, .99999999999999),              ###########################################
                    u1=u1, alpha=alpha, t=p, tol=1e-5)$root, silent=FALSE)}

  x2 <- eval(qdf.margin, list(x=u2))
  #x2
}

