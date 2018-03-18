
###########################################################
# ML ESTIAMTION FOR MULTIPLE-INFLATION POISSON (MIP) MODEL
###########################################################
# WRITTEN BY XIAOGANG SU, PH.D.
# UNIVERSITY OF ALABAMA AT BIRMINGHAM
# EMAIL: XGSU@UAB.EDU
###########################################################


options(warn=-1, digits=7);
require(MASS) 
if (!is.element("pscl", installed.packages()[,1])) install.packages("pscl")
library(pscl)


require(MASS)  # Function {MASS}polr() for fitting cumulative logit model    
logistic <- function(x) (tanh(x/2) +1)/2
logit <- function(x) log(x/(1-x))

MLE.MIP <- function(dat, cols.LM, cols.PM, col.y, M, max.it.EM=4,
	maxit=200, epsilon=sqrt(.Machine$double.eps), 
	use.gradient=T){

    require(MASS)  # Function {MASS}polr() for fitting cumulative logit model     
    singular.hess <- F                      #  AN INDICATOR OF WHETHER SINGULAR HESSIAN MATRIX IS RETURNED. 
    n <- nrow(dat); jvec <- rep(1, n)
    y <- dat[, col.y];  
    nvar.PM <- length(cols.PM)              # number of predictors in the Poisson Model
    nvar.LM <- length(cols.LM)              # number of predictors in the cumulative Logit Model
    if (nvar.PM==0) cols.PM <- NA
    if (nvar.LM==0) cols.LM <- NA

    
    # ================================
    #  Initialization for Parameters 
    # ================================
  
    # Initial Values for Gamma
    # ---------------------------
    y1 <- y; y1[y1 >= M] <- M
    # BE CAREFUL WITH SOME SPECIAL CASES -- NO X INVOLVED IN THE (CUMULATIVE) LOGIT MODEL	
    if (nvar.LM==0 || is.na(cols.LM)) {
	dat.Clogit <- data.frame(y=y1)
	G <- as.matrix(jvec)  # Obtain G Matrix for the Cumulative Logit Model
    } else  {
	dat.Clogit <- data.frame(dat[,cols.LM], y=y1)
	G <- as.matrix(cbind(jvec, dat[,cols.LM]))  # Matrix for the Cumulative Logit Model
    }


    fit.Clogit <- polr(factor(y)~.,  data=dat.Clogit, method = "logistic", Hess =T)
    #### Note that it is important to add a minus sign in front of fit.Clogit$coef. ### 
    #### Read help(polr) for details. 
    gamma.vec0 <- as.vector(c(fit.Clogit$zeta, -fit.Clogit$coef))    

    # ------------------------------------------
    # Obtain Initial Value for Beta - Two Steps
    # ------------------------------------------

    if (nvar.PM==0 || is.na(cols.PM)) {
	dat.Pois <-  dat[, col.y]
    	dat.poisson <- dat[dat$y >=M, col.y]
    	B <- as.matrix(jvec)   # Obtain B matrix for Poisson model	
    } else {
    	dat.Pois <-  dat[, c(cols.PM, col.y)]
   	dat.poisson <- dat[dat$y >=M, c(cols.PM, col.y)]
    	B <- as.matrix(cbind(jvec, dat[,cols.PM]))   # matrix for Poisson model
    }

    # poisson model for the data without inflated counts
    fit.PM <- glm(y~., data=dat.poisson, family=poisson(link = "log"))
    beta0 <- as.vector(fit.PM$coef) 


    d <- ncol(B)    
    # Create Delta and delta 
    Delta <- matrix(0, nrow=n, ncol=(M+1))
    for (m in 0:M) {
      nam <- paste("delta", m, sep="")
      if (m < M) {assign(nam, sign(y==m)); Delta[,(m+1)] <- sign(y==m)}
      else  {assign(nam, sign(y>=m)); Delta[,(m+1)] <- sign(y>=m)}
    }
    
    # The loglikelihood function for truncated poisson distribution
    L.trunc.pois <- function(betavec, M, dat){
      y <- dat$y; 
      B <- as.matrix(cbind(rep(1, nrow(dat)), dat[, -ncol(dat)]))
      lambda <- as.vector(exp(B%*%betavec))
      L <- y*log(lambda) - lambda - log(ppois(q=M-1, lambda=lambda, lower.tail = F))  
      return(-sum(L))
    } 
    
    Gradient.trunc.pois <- function(betavec, M, dat){
      y <- dat$y; 
      B <- as.matrix(cbind(rep(1, nrow(dat)), dat[, -ncol(dat)])) 
      lambda <- as.vector(exp(B%*%betavec))
      jcob <- (y-lambda) - dpois(M, lambda=lambda)*M/ppois(q=M-1, lambda=lambda, lower.tail = F) 
      as.vector(t(B)%*%jcob)
    }

    if (use.gradient) {
	fit.trunc.pois <- optim(par=beta0, fn=L.trunc.pois, gr=Gradient.trunc.pois, method = "BFGS", 
		M=M, dat=dat.poisson)
    	beta.vec0 <-fit.trunc.pois$par}
    else {
	fit1.trunc.pois <- optim(par=beta0, fn=L.trunc.pois, gr=NULL, method = "BFGS", 
          M=M, dat=dat.poisson)    
    	beta.vec0 <- fit1.trunc.pois$par
    }
    theta0 <- as.numeric(c(gamma.vec0, beta.vec0)); 
    # print(cbind(theta0, theta.true));   # Initial values for EM are ready.

    
    # =======================================================
    # EM Algorithm to Obtain Initial Values for Quasi-Newton
    # =======================================================
    
    EM <- function(theta0, dat, M, niter=3, tol=0.001){
        n <- nrow(dat)
        p <- length(theta0)  # total number of parameters
        theta <- theta0
        iter <- 0
        dif <- 1e10
        while (iter <= niter && dif > tol) {
            gamma0 <- theta[1:M]
            if ((p-d) >= (M+1)) gamma1 <- theta[(M+1):(p-d)]
		else gamma1 <- NULL
            betavec <- theta[(p-d+1):p]
            gammavec <- theta[1:(p-d)]
            
            # Calculate PI.m and P.m
            PI <- matrix(0, nrow=n, ncol=M) 
            P <- matrix(0, nrow=n, ncol=M+1)
            for (m in 1:M){
             PI[,m] <- logistic(G%*%c(gamma0[m], gamma1))   
             if (m==1) P[,m] <- PI[,m]
             else  P[,m] <- PI[,m] - PI[,(m-1)] 
            }
            P[,(M+1)] <- 1- PI[,M]
        
            # Calculate  Lambda
            lambda <- as.vector(exp(B%*%betavec))
        
            # -------------------
            # E-step: Update Z 
            # -------------------
            Z1 <- (Delta[, -(M+1)]*P[,-(M+1)])/(P[,-(M+1)] + P[, (M+1)]*dpois(y, lambda=lambda))
            Z.M <- as.vector(1- apply(Z1, 1, sum))
            Z <- as.matrix(cbind(Z1, Z.M))
            
            # -------------------
            # M-Step: Update Beta
            # -------------------
            fit.EM.beta <- glm(y~., data=dat.Pois, family=poisson(link = "log"), weights=Z.M)
            betavec <- as.vector(fit.EM.beta$coef) 
            
            # -------------------
            # M-Step: Update Gamma
            # -------------------
            Lc.gamma <- function(gammavec, M=M, Z=Z) {
                  gamma0 <- gammavec[1:M]
			if (length(gammavec) >= (M+1)) gamma1 <- gammavec[(M+1):(length(gammavec))]
			else gamma1 <- NULL
                  
                  PI <- matrix(0, nrow=n, ncol=M) 
                  P <- matrix(0, nrow=n, ncol=M+1)
                  for (m in 1:M){
                       PI[,m] <- logistic(G%*%c(gamma0[m], gamma1))   
                       if (m==1) P[,m] <- PI[,m]
                       else  P[,m] <- PI[,m] - PI[,(m-1)] 
                  }
                  P[,(M+1)] <- 1- PI[,M]
                  -sum(log(P^Z))
                  # Lc <- log(P^Z); Lc[is.na(Lc)] <- 0;    
                  # -sum(Lc)
            }
            fit.EM.gamma <- optim(par=gammavec, fn=Lc.gamma, gr=NULL, 
                  method = "BFGS", M=M, Z=Z)
            gammavec <- fit.EM.gamma$par  
            theta1 <- as.vector(c(gammavec, betavec));
            iter <- iter+1
            dif <- max(abs(theta1 - theta))
            theta <- theta1
            # print(iter); print(theta); print(dif)
        }
        return(theta)
    }

    # theta.0 <- EM(theta0=theta0, dat=dat, M=M, niter=10, tol=0.001)
    theta.0 <- EM(theta0=theta0, dat=dat, M=M, niter=max.it.EM, tol=.Machine$double.eps)
    # HERE THE EM ESTIMATES AFTER A FEW INTERACTIONS WILL BE USED AS INITIAL VALUES FOR MLE
    
    
    # =====================================
    # THE LOGLIKELIHOOD FUNCTION FOR Y 
    # =====================================
    
    L <- function(theta){
        p <- length(theta)  # total number of parameters
        gamma0 <- theta[1:M]
	  if ((p-d) >= (M+1)) 	 gamma1 <- theta[(M+1):(p-d)]
	  else gamma1 <- NULL
        betavec <- theta[(p-d+1):p]
        
        # Calculate PI.m and P.m
        PI <- matrix(0, nrow=n, ncol=M) 
        P <- matrix(0, nrow=n, ncol=M+1)
        for (m in 1:M){
         PI[,m] <- logistic(G%*%c(gamma0[m], gamma1))   
         if (m==1) P[,m] <- PI[,m]
         else  P[,m] <- PI[,m] - PI[,(m-1)] 
        }
        P[,(M+1)] <- 1- PI[,M]
        
        # Calculate  Lambda
        lambda <- as.vector(exp(B%*%betavec))
        # cbind(exp(-lambda)*lambda^y / factorial(y), dpois(y, lambda=lambda))
        
        # Compute the log-likelihood function  
        L <-  log((P[,-(M+1)] + P[,(M+1)]* dpois(y, lambda=lambda))^(Delta[,-(M+1)]))
        # L <-  ifelse(Delta[,-(M+1)],  log(P[,-(M+1)] + P[,(M+1)]* dpois(y, lambda=lambda)), 0)    
        L[is.na(L)] <- 0; # print(L)
        # L <- sum(apply(L, 1, sum) + ifelse(Delta[,(M+1)], log(P[,(M+1)]*dpois(y, lambda=lambda)), 0))
        L <- sum(apply(L, 1, sum) + log((P[,(M+1)]*dpois(y, lambda=lambda))^(Delta[,(M+1)])))
        return(-L)                                                           
    }           
    
    Gradient <- function(theta){
        p <- length(theta)  # total number of parameters
        gamma0 <- theta[1:M]; 
	  if ((p-d) >= (M+1)) gamma1 <- theta[(M+1):(p-d)]
	  else gamma1 <- NULL
        betavec <- theta[(p-d+1):p]
        
        # Calculate PI.m and P.m
        PI <- matrix(0, nrow=n, ncol=M) 
        P <- matrix(0, nrow=n, ncol=M+1)
        for (m in 1:M){
         PI[,m] <- logistic(G%*%c(gamma0[m], gamma1))   
         if (m==1) P[,m] <- PI[,m]
         else  P[,m] <- PI[,m] - PI[,(m-1)] 
        }
        P[,(M+1)] <- 1- PI[,M]
        # Calculate lambda
        lambda <- as.vector(exp(B%*%betavec))
        
        A <- P[,-(M+1)] + P[,(M+1)]* dpois(y, lambda=lambda)     # n X M
        
        d1 <- rep(0, M)
        for (m in 1:(M-1)){
          d1[m] <- sum(PI[,m] *(1-PI[,m]) *(Delta[,m]/A[,m] - Delta[,(m+1)]/A[,(m+1)]))    
        }
        # dim(Delta); dim(A); dim(PI)
        d1[M] <- sum(Delta[,M]* PI[,M] *(1-PI[,M])*(1- as.vector(dpois(M, lambda=lambda)))/A[,M]) - sum(Delta[, (M+1)] * PI[,M])
        
        # Compute Vector d2 
        PI.1 <- PI[,-1]*(1-PI[,-1]) - PI[, -M]*(1-PI[, -M])     
        Temp1 <- ifelse(Delta[,-c(1, (M+1))],  (PI.1 - PI[,M] *(1-PI[,M])*dpois(y, lambda=lambda)) /A[,-1], 0)
        if (!is.null(dim(Temp1))) Temp1 <- apply(Temp1, 1, sum)
        #                
        d2 <- Delta[,1] * (PI[,1] *(1-PI[,1]) - PI[,M] *(1-PI[,M])*exp(-lambda))/A[,1] +    
            +  Temp1 - Delta[,(M+1)]* PI[,M] *(1-PI[,M]) / P[, (M+1)]   
        deriv.gamma1 <- t(G[,-1])%*%d2
        
        # Compute Vector d3
        Temp2 <- ifelse(Delta[, -(M+1)], (dpois(y, lambda=lambda)*(y-lambda))/A, 0)
        # print(apply(Temp2, 1, sum))
        d3 <- apply(Temp2, 1, sum)  + Delta[, (M+1)]*(y-lambda)
        deriv.beta <- t(B)%*% d3
        grad <- -as.vector(c(d1, deriv.gamma1, deriv.beta))
        # print(grad)
        return(grad)
    }
    
    fit.ML <- optim(par=theta.0, fn=L, gr=NULL, method = "BFGS", hessian = T, 
			control=list(maxit=maxit, reltol=epsilon))    # This is the best way!!
    hess <- fit.ML$hessian; 
    if (det(hess) < sqrt(.Machine$double.eps)) {
      singular.hess <- T
      # cat("Hey! You got a singular Hessian matrix!!\n")
      covmat <- ginv(hess) 
    } 
    else covmat <- solve(hess)   # NOTE THAT optim() MINIMIZES by default.
    se <- sqrt(diag(covmat))
    results <- data.frame(theta.hat=fit.ML$par, se=se, z=fit.ML$par/se, p.value=2*(1-pnorm(abs(fit.ML$par/se)))) 
    if (nvar.LM ==0) rownames(results) <- c(paste("gamma0", 0:(M-1), sep=""), paste("beta", 0:nvar.PM, sep="")) 
    else rownames(results) <- c(paste("gamma0", 0:(M-1), sep=""), paste("gamma", 1:nvar.LM, sep=""), paste("beta", 0:nvar.PM, sep="")) 
    list(results=results, fit=fit.ML, singular.hess=singular.hess, covmat=covmat)
}






# =================================================================
# FUNCTION est.MIP.mean() ESTIMATES THE MEAN COUNT IN MIP MODEL
# =================================================================

est.MIP.mean <- function(theta, M, cols.LM, cols.PM, dat){
    p <- length(theta) 
    n <- nrow(dat);  jvec <- rep(1, n)
    d <- length(cols.PM) + 1              		# number of parameters in the Poisson Model
    if (length(cols.LM)==0) cols.LM <- NA
    if (length(cols.PM)==0) cols.PM <- NA

    # Obtain G and B
    if (length(cols.LM) ==0 || is.na(cols.LM)) G <- as.matrix(jvec) 
    else G <- as.matrix(cbind(jvec, dat[,cols.LM]))  	# Matrix for the Cumulative Logit Model
    if (length(cols.PM) ==0 || is.na(cols.PM)) B <- as.matrix(jvec)
    else B <- as.matrix(cbind(jvec, dat[,cols.PM]))   # matrix for Poisson model
    
    gamma0 <- theta[1:M]
    if ((p-d) >= (M+1)) gamma1 <- theta[(M+1):(p-d)]  # DEAL WITH SOME SPCIAL CASES 
    else gamma1 <- NULL 

    betavec <- theta[(p-d+1):p]           
    # Calculate PI.m and P.m
    PI <- matrix(0, nrow=n, ncol=M) 
    P <- matrix(0, nrow=n, ncol=M+1)
    for (m in 1:M){
        PI[,m] <- logistic(G%*%c(gamma0[m], gamma1))   
        if (m==1) P[,m] <- PI[,m]
        else  P[,m] <- PI[,m] - PI[,(m-1)] 
    }
    P[,(M+1)] <- 1- PI[,M]
        
    # Calculate  Lambda
    lambda <- as.vector(exp(B%*%betavec))
    est.mean <- P[, -(M+1)]%*% c(0:(M-1)) + P[, (M+1)]*lambda 
    return(est.mean)
} 


# ============================================
# FUNCTION report() REPORTS THE FITTING RESULT
# ============================================

report <- function(fit.MIP, variables.LM, variables.PM, M){
    		out <- fit.MIP$results
    		# print(out)
    		theta <- out$theta.hat 
    		se <- out$se
		wald.z <- theta/se
		p.value <- 2*pnorm(abs(wald.z), lower.tail=F)
		result <- as.data.frame(cbind(theta, se, wald.z, p.value))
		colnames(result) <- c("estimate", "se", "Wald.z", "P-Value")
		row.names(result) <- c(paste("gamma0.", 0:(M-1), sep=""), 
				paste("gamma.", variables.LM, sep=""), 
				paste("beta.", c(0, variables.PM), sep=""))
		return(result)
}


    