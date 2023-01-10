rm(list=ls())

start=Sys.time()
format(start,"%H:%M:%S")

## Set Lambda
lam <- 0.1


## original rho
rho123 <- c(-0.3,0.3,0.6)

## corresponding rho
rho12 <- c(-0.40,0.40,0.80) ; 
rho13 <- c(-0.31,0.31,0.62) ; 
rho23 <- c(-0.42,0.42,0.78) 


## Set In-control parameter
mu0 <- 0 ; p0 <- 0.3 ; c0 <- 3

## Require package
library(mvtnorm)
library(readr)

for(d in 1:length(rho123)){
  
  rho <- rho123[d]
  
  ARL <- numeric()
  time_elapse <- numeric()
  
  for(k in 1:20){
    
    ## Set multiple values for ARL_0 =200, 400, 600
    Lx1 <- c(2.67,2.92,3.06) ; Lp1 <- c(2.88,3.14,3.29) ; Lc1 <- c(2.85,3.14,3.31)
    L1 <- matrix(rbind(Lx1,Lp1,Lc1),3,3)
    
    Lx2 <- c(2.67,2.92,3.06) ; Lp2 <- c(2.88,3.14,3.29) ; Lc2 <- c(2.85,3.14,3.31)
    L2 <- matrix(rbind(Lx2,Lp2,Lc2),3,3)
    
    Lx3 <- c(2.67,2.92,3.06) ; Lp3 <- c(2.88,3.14,3.29) ; Lc3 <- c(2.85,3.14,3.31)
    L3 <- matrix(rbind(Lx3,Lp3,Lc3),3,3)
    
    L <- array(0,dim = c(3,3,3))
    L[,,1] <- L1
    L[,,2] <- L2
    L[,,3] <- L3
    
    
    L[,3,d] <- L[,3,d] - 0.01*(k+1)
  
    start_time <- proc.time()
    
    RL <- numeric() 
    N <- 10000
    
    for(j in 1:N){
      
      ## Set count for ARL_0 
      count <- 1 
      
      ## Setting n (observations) ; m (subgroup size)
      n <- 1 ;  m <- 1  
      
      CorrZ <- matrix(0,3,3) # CorrZ will be the correlation matrix of the standard normal vector.
      CorrZ[1,2] <- rho12[d]      
      CorrZ[1,3] <- rho13[d]
      CorrZ[2,3] <- rho23[d]
      CorrZ <- CorrZ + t(CorrZ)          # Get the lower triangle too.
      diag(CorrZ)<-1
      
      muZ <- c(0,0,0)   
      
      sig_x0 <- 1 ; sig_p0 <- sqrt(p0*(1-p0)) ; sig_c0 <- sqrt(c0)
      
      ## Compute standard deviation
      Sx <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n)))*sig_x0/sqrt(m)
      Sp <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n)))*sig_p0/sqrt(m)
      Sc <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n)))*sig_c0/sqrt(m)
      
      #set.seed(count)
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
      
      ##----- Normal ------------
      E_x0 <- mu0
      ## Calculate EWMA charting statistic and standardize
      E_x <- lam*Y0[1,1] + (1-lam)*E_x0
      E_x <- max(mu0,E_x)
      E_x_stand <- (E_x - mu0)/Sx
      
      ##---- Binomial ----------
      E_p0 <- p0
      ## Calculate EWMA charting statistic and standardize
      E_p <- lam*Y0[1,2] + (1-lam)*E_p0
      E_p <- max(p0,E_p)
      E_p_stand <- (E_p - p0)/Sp
      
      ##---- Poisson ----------
      E_c0 <- c0
      ## Calculate EWMA charting statistic and standardize
      E_c <- lam*Y0[1,3] + (1-lam)*E_c0
      E_c <- max(c0,E_c)
      E_c_stand <- (E_c - c0)/Sc
      
      E <- c(E_x_stand,E_p_stand,E_c_stand)
      
      E_sort <- sort(E)
      E_ord <- order(E)
      
      n1 <- 1
      repeat{
        
        if(E_sort[3] < L[E_ord[3],3,d]){ 
          
        }else{
          RL[j] <- count  
          break        
        } 
        count <- count + 1
        #set.seed(count)
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
        
        n1 <- n1+1
        ## Compute standard deviation
        Sx <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n1)))*sig_x0/sqrt(m)
        Sp <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n1)))*sig_p0/sqrt(m)
        Sc <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n1)))*sig_c0/sqrt(m)
        
        ##----- Normal ------------
        temp_x <- E_x
        ## Calculate EWMA charting statistic and standardize
        E_x <- lam*Y0[1,1] + (1-lam)*temp_x
        E_x <- max(mu0,E_x)
        E_x_stand <- (E_x - mu0)/Sx
        
        ##---- Binomial ----------
        temp_p <- E_p
        ## Calculate EWMA charting statistic and standardize
        E_p <- lam*Y0[1,2] + (1-lam)*temp_p
        E_p <- max(p0,E_p)
        E_p_stand <- (E_p - p0)/Sp
        
        ##---- Poisson ----------
        temp_c <- E_c
        ## Calculate EWMA charting statistic and standardize
        E_c <- lam*Y0[1,3] + (1-lam)*temp_c
        E_c <- max(c0,E_c)
        E_c_stand <- (E_c - c0)/Sc
        
        E <- c(E_x_stand,E_p_stand,E_c_stand)
        
        E_sort <- sort(E)
        E_ord <- order(E)
        
      }
      
    }
    
    
    time <- proc.time() - start_time
    time_elapse[k] <- time[3]
    
    ARL[k] <- mean(RL)
    
    ## save result to df.csv file
    if(k == 1){
      df <- cbind(rho,k,L[1,1,d],L[2,1,d],L[3,1,d],L[1,2,d],L[2,2,d],L[3,2,d],
                  L[1,3,d],L[2,3,d],L[3,3,d],ARL[k], time_elapse[k])
      df <- data.frame(df)
      names(df) <- c("rho","k","Lx2","Lp2","Lc2","Lx4","Lp4","Lc4","Lx6","Lp6","Lc6",
                     "ARL","time_elapse")
      write_csv(df,file=paste("result/one-side-EWMA/Xpc-ARL0-Holm-rho-find-L-up.csv",sep = ""),append = TRUE,
                col_names = TRUE)
    }else{
      df <- cbind(rho,k,L[1,1,d],L[2,1,d],L[3,1,d],L[1,2,d],L[2,2,d],L[3,2,d],
                  L[1,3,d],L[2,3,d],L[3,3,d],ARL[k], time_elapse[k])
      df <- data.frame(df)
      names(df) <- c("rho","k","Lx2","Lp2","Lc2","Lx4","Lp4","Lc4","Lx6","Lp6","Lc6",
                     "ARL","time_elapse")
      write_csv(df,file=paste("result/one-side-EWMA/Xpc-ARL0-Holm-rho-find-L-up.csv",sep = ""),append = TRUE)
    } 
  
  }

}


end=Sys.time()
end-start
