rm(list=ls())

start=Sys.time()
format(start,"%H:%M:%S")


## Set multiplier
L_t <- seq(11.3,21.9,0.1)

## Set R matrix
r <- 0.1
R <- matrix(0,3,3) 
diag(R) <- r

I <- matrix(0,3,3) 
diag(I) <- 1


## Set In-control parameter
mu0 <- 0 ; p0 <- 0.3 ; c0 <- 3

## Require package
library(mvtnorm)
library(readr)

## original rho
rho123 <- c(-0.3,0,0.3,0.6)

## corresponding rho
rho12 <- c(-0.40,0,0.40,0.80) ; 
rho13 <- c(-0.31,0,0.31,0.62) ; 
rho23 <- c(-0.42,0,0.42,0.78) 

for(d in 1:length(rho123)){
  
  rho <- rho123[d]
  ARL <- numeric() ; time_elapse <- numeric()
  
  for(k in 1:length(L_t)){
    
    start_time <- proc.time()
    RL <- numeric() ; H_T2 <- numeric()
    N <- 10000
    
    for(j in 1:N){
      
      ## Set count for ARL_0 
      count <- 1 
      
      ## Setting n (observations) ; m (subgroup size)
      n <- 1 ;  m <- 1  
      
      CovZ <- matrix(0,3,3) # CovZ will be the covariance matrix of the standard normal vector.
      CovZ[1,2] <- rho12[d]      
      CovZ[1,3] <- rho13[d]
      CovZ[2,3] <- rho23[d]
      CovZ <- CovZ + t(CovZ)  # Get the lower triangle too.
      diag(CovZ) <- 1
      
      muZ <- c(0,0,0)   
      
      sig_x0 <- 1 ; sig_p0 <- sqrt(p0*(1-p0)) ; sig_c0 <- sqrt(c0)
      
      #set.seed(count)
      ## Generate simulated observations
      Z <- array(0,dim = c(n,3,m)) ; U <- array(0,dim = c(n,3,m)) ;
      Y0 <- array(0,dim = c(n,3,m))
      for(ss in 1:m){
        Z[,,ss] <- rmvnorm(n,muZ,CovZ)
        U[,,ss] <- pnorm(Z[,,ss])
        Y0[,,ss] <- matrix( cbind(qnorm(U[,1,ss],mu0,1),qbinom(U[,2,ss],1,p0),
                                  qpois(U[,3,ss],c0)),n,3 )
      }
      
      Y0 <- apply(Y0,c(1,2),mean)
      Y0 <- matrix(Y0,3,1)
      
      
      ## Compute EWMA statistic
      W0 <- matrix(c(mu0,p0,c0),3,1)
      E_W <- R %*% Y0 + (I-R) %*% W0
      
      ## Compute Y   muY and Sigma_Y
      CovY <- matrix(0,3,3) 
      CovY[1,2] <- rho123[d]*sqrt(sig_x0)*sqrt(p0*(1-p0))      
      CovY[1,3] <- rho123[d]*sqrt(sig_x0)*sqrt(c0) 
      CovY[2,3] <- rho123[d]*sqrt(p0*(1-p0))*sqrt(c0) 
      CovY <- CovY + t(CovY)          
      diag(CovY) <- c(sig_x0,p0*(1-p0),c0)
      
      mu_Y <- matrix(c(mu0,p0,c0),3,1) ; S_Y <- matrix(CovY,3,3)
      
      lam <- t(mu_Y) %*% solve(S_Y) %*% mu_Y
      
      ## Compute T2
      S_W <- r*(1-(1-r)^(2*n))/(2-r) * S_Y
      T2 <- t(E_W-mu_Y) %*% solve(S_W) %*% (E_W-mu_Y)
      
      
      n1 <- 1
      repeat{
        
        if(T2 > L_t[k]){
          RL[j] <- count  
          H_T2[j] <- T2
          break
        }
        
        count <- count + 1
        #set.seed(count)
        ## Generate simulated observations
        Z <- array(0,dim = c(n,3,m)) ; U <- array(0,dim = c(n,3,m)) ;
        Y0 <- array(0,dim = c(n,3,m))
        for(ss in 1:m){
          Z[,,ss] <- rmvnorm(n,muZ,CovZ)
          U[,,ss] <- pnorm(Z[,,ss])
          Y0[,,ss] <- matrix( cbind(qnorm(U[,1,ss],mu0,1),qbinom(U[,2,ss],1,p0),
                                    qpois(U[,3,ss],c0)),n,3 )
        }
        
        Y0 <- apply(Y0,c(1,2),mean)
        Y0 <- matrix(Y0,3,1)
        
        n1 <- n1+1
        ## Compute EWMA statistic
        temp_w <- E_W
        E_W <- R %*% Y0 + (I-R) %*% temp_w
        
        ## Compute T2
        S_W <- r*(1-(1-r)^(2*n1))/(2-r) * S_Y
        T2 <- t(E_W-mu_Y) %*% solve(S_W) %*% (E_W-mu_Y)
        
      }
      
    }
    
    
    ARL[k] <- mean(RL)
    
    time <- proc.time() - start_time
    time_elapse[k] <- time[3]
    
    ## save result to df.csv file
    if(k == 1){
      df <- data.frame(rho = rho, Lt = L_t[k], ARL = ARL[k], time_elapse=time_elapse[k])
      write_csv(df,file=paste("result/rho/T2-find-L-rho-2.csv",sep = ""),append = TRUE,
                col_names = TRUE)  
    }else{
      df <- data.frame(rho = rho, Lt = L_t[k], ARL = ARL[k], time_elapse=time_elapse[k])
      write_csv(df,file=paste("result/rho/T2-find-L-rho-2.csv",sep = ""),append = TRUE)
    }
      
  }

}

end=Sys.time()
end-start
