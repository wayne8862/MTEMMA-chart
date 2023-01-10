rm(list=ls())

start=Sys.time()
format(start,"%H:%M:%S")

## Set Lambda
lam <- 0.1

## Set multiple values
Lp <- seq(3.49,3.69,by=0.01) 
#Lp <- c(2.38,2.66,2.83)

## Set In-control parameter
mu0 <- 0 ; p0 <- 0.3 ; c0 <- 3

## Require package
library(mvtnorm)
library(readr)

ARL_p <- numeric() 
time_elapse <- numeric()

for(k in 1:length(Lp)){
  
  RL_p <- numeric() 
  start_time <- proc.time()
  N <- 10000
  for(j in 1:N){
    
    ## Set count for ARL_0 
    count_p <- 1 
    
    ## Setting n (observations) ; m (subgroup size)
    n <- 1 ;  m <- 1  
    
    CorrZ <- matrix(0,3,3) # CorrZ will be the correlation matrix of the standard normal vector.
    CorrZ[1,2] <- 0      
    CorrZ[1,3] <- 0
    CorrZ[2,3] <- 0
    CorrZ <- CorrZ + t(CorrZ)          # Get the lower triangle too.
    diag(CorrZ)<-1
    
    muZ <- c(0,0,0)   
    
    ## Set In-control parameter
    sig_p0 <- sqrt(p0*(1-p0)) 
    
    
    ## Compute standard deviation
    Sp <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n)))*sig_p0/sqrt(m)
 
    
    ## Generate simulated observations
    Z <- array(0,dim = c(n,3,m)) ; U <- array(0,dim = c(n,3,m)) ;
    Y0 <- array(0,dim = c(n,3,m))
    for(ss in 1:m){
      Z[,,ss] <- rmvnorm(n,muZ,CorrZ)
      U[,,ss] <- pnorm(Z[,,ss])
      Y0[,,ss] <- matrix( cbind(qnorm(U[,1,ss],mu0,1),qbinom(U[,2,ss],1,p0),
                                qpois(U[,3,ss],c0)),n,3 )
    }
    
    Y0 <- apply(Y0,c(1,2),mean)
    
    
    ##---- Binomial ----------
    E_p0 <- p0
    
    ## Calculate EWMA charting statistic and standardize
    E_p <- lam*Y0[1,2] + (1-lam)*E_p0
    E_p <- max(p0,E_p)
    E_p_stand <- (E_p - p0)/Sp
    

    n1 <- 1
    repeat{
      
      if(abs(E_p_stand) > Lp[k]){
        RL_p[j] <- count_p
        break
      }
      count_p <- count_p +1
      
      Z <- array(0,dim = c(1,3,m)) ; U <- array(0,dim = c(1,3,m)) ;
      Y <- array(0,dim = c(1,3,m))
      for(ss in 1:m){
        Z[,,ss] <- rmvnorm(1,muZ,CorrZ)
        U[,,ss] <- pnorm(Z[,,ss])
        Y[,,ss] <- matrix( cbind(qnorm(U[,1,ss],mu0,1),qbinom(U[,2,ss],1,p0),
                                 qpois(U[,3,ss],c0)),1,3 )
      }
      
      Y <- apply(Y,c(1,2),mean)
      
      n1 <- n1+1
      Sp <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n1)))*sig_p0/sqrt(m)
      
      ## Calculate EWMA charting statistic
      temp_p <- E_p
      E_p <- lam*Y[1,2]+(1-lam)*temp_p
      E_p <- max(p0,E_p)
      E_p_stand <- (E_p - p0)/Sp
      
    }
    
    RL_p[j] <- count_p 
    
    
  }
  
  ARL_p[k] <- mean(RL_p)

  time <- proc.time() - start_time
  time_elapse[k] <- time[3]
  
  ## save result to df.csv file
  if(k == 1){
    df <- data.frame( Lp = Lp[k], ARL = ARL_p[k], time_elapse=time_elapse[k])
    write_csv(df,file=paste("result/one-side-EWMA/Find-c-P.csv",sep = ""),append = TRUE,
              col_names = TRUE)  
  }else{
    df <- data.frame( Lp = Lp[k], ARL = ARL_p[k], time_elapse=time_elapse[k])
    write_csv(df,file=paste("result/one-side-EWMA/Find-c-P.csv",sep = ""),append = TRUE)
  }
  
}


end=Sys.time()
end-start


