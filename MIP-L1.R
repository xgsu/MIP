

#######################################################
## L1-REGUARED MULTIPLE-INFLATION POISSON MODEL (MIP)
## with Local Quadratic Approximation (LQA)
#######################################################
# WRITTEN BY: XIAOGANG SU, PH.D.
# UNIVERSITY OF ALABAMA AT BIRMINGHAM
# EMAIL: XGSU@UAB.EDU
# PART OF THE CODES ADOPTED FROM WANG AND LENG (2006). 
#######################################################


# Sigma0 ---- the varaince-covaraince matrix; 
# b0 ---- the estiamted parameters from the full model; 
# M ---- the number of inflated count values;
# p ---- number of covariates used in the cumulative logit model
# n ---- the sample size

if (!is.element("lars", installed.packages()[,1])) install.packages("lars")
require(lars)

LAS.LAR.MI <- function(Sigma0, b0,      
                      M,  n,  p, 
                      eps = .Machine$double.eps, max.steps =20) 
{ 
 

  # Sigma0 <- solve(Sigma0); # no need if Hessian matrix is provided directly.
  n1 <- nrow(Sigma0)      # n1 - total number of parameters 
  
  # EXCHANGE THE ELEMENTS TO MAKE SURE INTERCEPTS ARE PLACED TOGETHER
  id.new <- c(1:M, (M+p+1), (M+1):(M+p), (M+p+2):length(b0))
  b.0 <- b0[id.new] 
  b0 <- b.0

  Sigma.0 <- Sigma0[id.new,]
  Sigma.0 <- Sigma.0[,id.new]
  Sigma0 <- Sigma.0;  

  # ----------------------
  ## HANDLE THE INTERCEPT
  #  ---------------------
  a11 <- Sigma0[1:(M+1), 1:(M+1)]                   
  a12 <- Sigma0[(M+2):n1,1:(M+1)]; 
  a22 <- Sigma0[(M+2):n1,(M+2):n1]                  
  Sigma <- a22 - a12%*%solve(a11)%*% t(a12)	    
  b <- b0[(M+2):n1]                                  
  # beta0 <- b%*%a12%*%solve(a11)                    
  beta0 <- solve(a11)%*% t(a12) %*% b		     
  # cat("beta0 = ", beta0, "\n")

  # ----------------------
  # THE LAR PROCEDURE
  # ----------------------
  Sigma <- diag(abs(b))%*%Sigma%*%diag(abs(b)) 	     
  b <- sign(b)                                       
  nm <- dim(Sigma)
  m <- nm[2]                                        # m - number of slopes
  im <- inactive <- seq(m)
  Cvec <- drop(t(b)%*%Sigma)
  ssy <- sum(Cvec*b)
                                                     
  if (missing(max.steps))  max.steps <- 8 * m
  beta <- matrix(0, max.steps + 1, m)
  Gamrat <- NULL
  arc.length <- NULL
  R2 <- 1
  RSS <- ssy
  first.in <- integer(m)
  active <- NULL
  actions <- as.list(seq(max.steps))
  drops <- FALSE
  Sign <- NULL
  R <- NULL
  k <- 0
  ignores <- NULL
  while ((k < max.steps) & (length(active) < m)) {
    print(k)
    action <- NULL
    k <- k + 1
    C <- Cvec[inactive]
    Cmax <- max(abs(C))
    if (!any(drops)) {
      new <- abs(C) >= Cmax - eps
      C <- C[!new]
      new <- inactive[new]
      for (inew in new) {
        R <- updateR(Sigma[inew, inew], R, drop(Sigma[inew, active]),
                     Gram = TRUE,eps=eps)
        if(attr(R, "rank") == length(active)) {
          ##singularity; back out
          nR <- seq(length(active))
          R <- R[nR, nR, drop = FALSE]
          attr(R, "rank") <- length(active)
          ignores <- c(ignores, inew)
          action <- c(action,  - inew)
        }
        else {
          if(first.in[inew] == 0)
            first.in[inew] <- k
          active <- c(active, inew)
          Sign <- c(Sign, sign(Cvec[inew]))
          action <- c(action, inew)
        }
      }
    }
    else action <- -dropid

    ################# 
    print(R)
    ###################
    Gi1 <- backsolve(R, backsolvet(R, Sign))
    dropouts <- NULL
    A <- 1/sqrt(sum(Gi1 * Sign))
    w <- A * Gi1
    if (length(active) >= m) {
      gamhat <- Cmax/A      
    }
    else {        
      a <- drop(w %*% Sigma[active, -c(active,ignores), drop = FALSE])
      gam <- c((Cmax - C)/(A - a), (Cmax + C)/(A + a))
      gamhat <- min(gam[gam > eps], Cmax/A)
    }
    # LASSO Lines  =====
      dropid <- NULL
      b1 <- beta[k, active]
      z1 <- -b1/w
      zmin <- min(z1[z1 > eps], gamhat)
      # cat('zmin ',zmin, ' gamhat ',gamhat,'\n') 
      if (zmin < gamhat) {
        gamhat <- zmin
        drops <- z1 == zmin
      }
      else drops <- FALSE
   # =====================
   
    beta[k + 1, ] <- beta[k, ]
    beta[k + 1, active] <- beta[k + 1, active] + gamhat * w
    Cvec <- Cvec - gamhat * Sigma[, active, drop = FALSE] %*% w   
    Gamrat <- c(Gamrat, gamhat/(Cmax/A))
    arc.length <- c(arc.length, gamhat)
    if (any(drops)) {
       dropid <- seq(drops)[drops]
       for (id in rev(dropid)) {
       R <- downdateR(R,id)
       }
       dropid <- active[drops]
       beta[k + 1, dropid] <- 0
       active <- active[!drops]
       Sign <- Sign[!drops]
    }
    actions[[k]] <- action
    inactive <- im[-c(active)]
  }
  beta <- beta[seq(k + 1), ]                       
  dff <- b-t(beta)                                  
  RSS <- diag(t(dff)%*%Sigma%*%dff)                 
  beta <- t(abs(b0[(M+2):n1])*t(beta))             
  beta0 <-  matrix(rep(beta0, (k+1)), , (k+1)) - t(drop(beta%*%a12%*%solve(a11)))
  # print(beta0) 
  # cat("beta0 = ", beta0, "\n")

  # Compute AIC and BIC
  dof <- apply(abs(beta)>eps,1,sum)
  BIC <- RSS+log(n)*dof
  AIC <- RSS+2*dof

  # OBTAIN THE SELECTION RESULTS - AIC AND BIC
  t1 <- sort(BIC, index.return =T)
  t2 <- sort(AIC, ind=T)
  beta0 <- beta0 +  matrix(rep(b0[1:(M+1)], k+1), , k+1)                 
  beta0 <- t(beta0) 
  beta.bic <- c(beta0[t1$ix[1], ], beta[t1$ix[1],])
  beta.aic <- c(beta0[t2$ix[1], ], beta[t2$ix[1],])
  object <- list(AIC = AIC, BIC = BIC, beta = beta, beta0 = beta0,      # RESULTS FROM LASSO
        beta.unpen=b0, beta.bic=beta.bic, beta.aic = beta.aic)          # THE BEST RESULTS
  object
}

