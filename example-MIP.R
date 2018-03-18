#######################################################
## L1-REGUARED MULTIPLE-INFLATION POISSON MODEL (MIP)
## with Local Quadratic Approximation (LQA)
#######################################################
# WRITTEN BY: XIAOGANG SU, PH.D.
# UNIVERSITY OF ALABAMA AT BIRMINGHAM
# EMAIL: XGSU@UAB.EDU
#######################################################

# FIRST DOWNLOAD ALL THREE FILES INTO THE SAME FOLDER. 
# THEN RUN THIS FILE ONLY, WHICH CONTAINS A SIMPLE ILLUSTRATION
# USING A SIMULATED DATA SET.  

source("MIP-MLE.R");
source("MIP-L1.R");

# ======================================================
# FUNCTION THAT GENERATES DATA FROM AN MIP MODEL
# ======================================================

rdat.mip <- function(n=100, gammas=c(-3, -1.5, 3, 2), betas=c(-2, 3, 2), 
              poisson.only=F, clogit.only=F, confound=T){
  x1 <- runif(n); x2 <- runif(n); 
  x3 <- runif(n); x4 <- runif(n);
  jvec <- rep(1, n)
  gamma1 <- gammas[-2];  gamma2 <- gammas[-1]
  G <- cbind(jvec, x1, x3) 
  p1 <- logistic(G%*%gamma1); p2 <- logistic(G%*%gamma2) - p1
  p3 <- 1- logistic(G%*%gamma2)
  if (poisson.only)   {p1<- p2 <- rep(0,n); p3 <- rep(1,n)}
  # print(cbind(p1, p2, p3, p1+p2+p3)) 
 
  if (confound) B <- cbind(jvec, x2, x3) 
  else B <- cbind(jvec, x2, x4) 
  lambda <- exp(B%*%betas) 
  y <- rep(0, n)
  for (i in 1:n){
      if (clogit.only) {
            y[i] <- sample(c(0, 1, 2), size=1, prob=c(p1[i], p2[i], p3[i]))
            lambda <- rep(2, n)}
      else {
            y[i] <- sample(c(0, 1, rpois(1, lambda=lambda[i])), 
                    size=1, prob=c(p1[i], p2[i], p3[i]))}
  }
  
  E.y <- p1*0 + p2*1 + p3*lambda # THE EXPECTED Y
  dat <- data.frame(x1=x1, x2=x2, x3=x3, x4=x4, y=y)
  list(dat=dat, E.y=E.y)
} 
 

# ==========================
# AN ILLUSTRATION OF CODES
# ==========================

set.seed(123)   
n <- 500
rdat <- rdat.mip(n=n, confound=T)
dat <- rdat$dat; 

# ----------------------------------------
# FIRST FIT THE MODEL WITH ALL COVARIATES
# ----------------------------------------


source("MIP-MLE.R");

cols.LM <- 1:4; cols.PM <- 1:4; col.y <- 5; M <- 2  
FIT <- MLE.MIP(dat, cols.LM=cols.LM, cols.PM=cols.PM, col.y=col.y, 
		max.it.EM=3, maxit=100, M=M, epsilon=1e-7, use.gradient=F)
out <- FIT$results
out


# -------------------------------------------
# L1-REGULARIZATION FOR VARIABLE SELECTION
# -------------------------------------------
n <- nrow(dat)
out <- FIT$results
theta.hat <- out$theta.hat    
fit.full <- FIT$fit 
Sigma0 <- fit.full$hessian
results.selection <- LAS.LAR.MI(Sigma0, b0=theta.hat, M=M, p=length(cols.LM), 
	n=n, eps = .Machine$double.eps, max.steps =30) 
results.selection
cbind(theta.hat, beta.unpen=results.selection$beta.unpen, bic=results.selection$beta.bic, aic=results.selection$beta.aic)


Best.step.BIC <- which.min(results.selection$BIC)
Best.step.AIC <- which.min(results.selection$AIC)
betas <- results.selection$beta  		# THESE ARE SLOPES ONLY.
which(betas[Best.step.BIC,] !=0)


# TWO PLOTS FROM L1-REGULARIZATION
# ----------------------------------

# postscript(file="fig1.eps", horizontal=F)
par(mfrow=c(2,1), mar=rep(4,4))
aic <- results.selection$AIC
bic <- results.selection$BIC
plot(x=c(1, nrow(betas)), y=c(min(aic), max(bic)), type="n",
	xlab="step", ylab="information criterion", main="(a) AIC & BIC")
lines(x=1:nrow(betas), y=aic, lty=1, col="green", lwd=1)
lines(x=1:nrow(betas), y=bic, lty=1, col="red", lwd=1)
abline(v=Best.step.BIC, lwd=3, col="blue")
abline(v=Best.step.AIC, lwd=3, col="green")

plot(x=c(1, nrow(betas)+1), y=c(min(betas), max(betas)), type="n",
	xlab="step", ylab="coefficient estimate", main="(b) regularized path")
for (j in 1:ncol(betas)){
	lines(x=1:nrow(betas), y=betas[,j], col="red", lty=1, lwd=1)
}
abline(h=0, col="black", lwd=1)
abline(v=Best.step.BIC, lwd=3, col="blue")
abline(v=Best.step.AIC, lwd=3, col="green")
# dev.off()


# ---------------------------------
# FITTING THE SELECTED MIP MODEL
# ---------------------------------

cols.LM <- c(1, 3); cols.PM <- c(2, 3); col.y <- 5; M <- 2  
FIT1 <- MLE.MIP(dat, cols.LM=cols.LM, cols.PM=cols.PM, col.y=col.y, 
		max.it.EM=3, maxit=100, M=M, epsilon=1e-7, use.gradient=F)
out <- FIT1$results
out

# USE THE report() FUNCTION 
vnames <- colnames(dat)
variables.LM <- vnames[cols.LM]
variables.PM <- vnames[cols.PM]
result <- report(fit.MIP=FIT1, variables.LM=variables.LM, 
		variables.PM=variables.PM, M=M)
result


# -------------------------------
# PREDICTION WITH A TEST SAMPLE
# -------------------------------

test <- rdat.mip(n=n, confound=T)$dat; 
y.hat <- est.MIP.mean(theta=out$theta.hat, M=M, cols.LM=cols.LM, 
		cols.PM = cols.PM, dat=test)
plot(test$y, y.hat, type="p", pch=19, cex=0.5, xlab="observed", 
	ylab="predicted", col="blue")
abline(a=0, b=1, col="red", lwd=2)





#
