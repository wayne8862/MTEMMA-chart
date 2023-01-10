rm(list=ls())

start=Sys.time()
format(start,"%H:%M:%S")

## Set Lambda
lam <- 0.1


## original rho
rho123 <- c(-0.3,0,0.3,0.6)

## corresponding rho
rho12 <- c(-0.40,0,0.40,0.80) ; 
rho13 <- c(-0.31,0,0.31,0.62) ; 
rho23 <- c(-0.42,0,0.42,0.78) 


## Set multiple values for individual EWMA chart ARL_0 =400, 800, 1200
Lx1_U <- c(2.67,2.92,3.06) ; Lp1_U <- c(2.88,3.14,3.29) ; Lc1_U <- c(2.85,3.14,3.31)
L1_U <- matrix(rbind(Lx1_U,Lp1_U,Lc1_U),3,3)
L1_U[,3] <- L1_U[,3] - 0.01

Lx2_U <- c(2.67,2.92,3.06) ; Lp2_U <- c(2.88,3.14,3.29) ; Lc2_U <- c(2.85,3.14,3.31)
L2_U <- matrix(rbind(Lx2_U,Lp2_U,Lc2_U),3,3)
#L2_U[,3] <- L2_U[,3] + 0.05

Lx3_U <- c(2.67,2.92,3.06) ; Lp3_U <- c(2.88,3.14,3.29) ; Lc3_U <- c(2.85,3.14,3.31)
L3_U <- matrix(rbind(Lx3_U,Lp3_U,Lc3_U),3,3)
L3_U[,3] <- L3_U[,3] - 0.01

Lx4_U <- c(2.67,2.92,3.06) ; Lp4_U <- c(2.88,3.14,3.29) ; Lc4_U <- c(2.85,3.14,3.31)
L4_U <- matrix(rbind(Lx4_U,Lp4_U,Lc4_U),3,3)
L4_U[,3] <- L4_U[,3] - 0.01*4

L_U <- array(0,dim = c(3,3,4))
L_U[,,1] <- L1_U
L_U[,,2] <- L2_U
L_U[,,3] <- L3_U
L_U[,,4] <- L4_U

## Set multiple values for individual EWMA chart ARL_0 =400, 800, 1200
Lx1_L <- c(-2.66,-2.92,-3.06) ; Lp1_L <- c(-2.215,-2.324,-2.375) ; Lc1_L <- c(-2.46,-2.68,-2.79)
L1_L <- matrix(rbind(Lx1_L,Lp1_L,Lc1_L),3,3)
L1_L[,3] <- L1_L[,3] + 0.01*1

Lx2_L <- c(-2.66,-2.92,-3.06) ; Lp2_L <- c(-2.215,-2.324,-2.375) ; Lc2_L <- c(-2.46,-2.68,-2.79)
L2_L <- matrix(rbind(Lx2_L,Lp2_L,Lc2_L),3,3)
#L2_L[,3] <- L2_L[,3] - 0.01*5

Lx3_L <- c(-2.66,-2.92,-3.06) ; Lp3_L <- c(-2.215,-2.324,-2.375) ; Lc3_L <- c(-2.46,-2.68,-2.79)
L3_L <- matrix(rbind(Lx3_L,Lp3_L,Lc3_L),3,3)
L3_L[,3] <- L3_L[,3] + 0.01*1

Lx4_L <- c(-2.66,-2.92,-3.06) ; Lp4_L <- c(-2.215,-2.324,-2.375) ; Lc4_L <- c(-2.46,-2.68,-2.79)
L4_L <- matrix(rbind(Lx4_L,Lp4_L,Lc4_L),3,3)
L4_L[,3] <- L4_L[,3] + 0.01*4

L_L <- array(0,dim = c(3,3,4))
L_L[,,1] <- L1_L
L_L[,,2] <- L2_L
L_L[,,3] <- L3_L
L_L[,,4] <- L4_L


## Set In-control parameter
mu0 <- 0 ; p0 <- 0.3 ; c0 <- 3

## Require package
library(mvtnorm)
library(readr)

ARL <- numeric() ; RL_se <- numeric()
for(d in c(1,2,3,4)){
  
  rho <- rho123[d]
  
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
    E_x_stand <- (E_x - mu0)/Sx
    
    E_x_max <- max(mu0,E_x)
    E_x_max_stand <- (E_x_max - mu0)/Sx
    
    E_x_min <- min(mu0,E_x)
    E_x_min_stand <- (E_x_min - mu0)/Sx
    
    ##---- Binomial ----------
    E_p0 <- p0
    ## Calculate EWMA charting statistic and standardize
    E_p <- lam*Y0[1,2] + (1-lam)*E_p0
    E_p_stand <- (E_p - p0)/Sp
    
    E_p_max <- max(p0,E_p)
    E_p_max_stand <- (E_p_max - p0)/Sp
    
    E_p_min <- min(p0,E_p)
    E_p_min_stand <- (E_p_min - p0)/Sp
    
    ##---- Poisson ----------
    E_c0 <- c0
    ## Calculate EWMA charting statistic and standardize
    E_c <- lam*Y0[1,3] + (1-lam)*E_c0
    E_c_stand <- (E_c - c0)/Sc
    
    E_c_max <- max(c0,E_c)
    E_c_max_stand <- (E_c_max - c0)/Sc
    
    E_c_min <- min(c0,E_c)
    E_c_min_stand <- (E_c_min - c0)/Sc
    
    E <- c(E_x_stand,E_p_stand,E_c_stand)
    E_max <- c(E_x_max_stand,E_p_max_stand,E_c_max_stand)
    E_min <- c(E_x_min_stand,E_p_min_stand,E_c_min_stand)
    
    E_sort <- sort(E)
    E_ord <- order(E)
    
    E_max_sort <- sort(E_max)
    E_max_ord <- order(E_max)
    
    E_min_sort <- sort(E_min)
    E_min_ord <- order(E_min)
    
    n1 <- 1
    repeat{
      
      # if(abs(E_sort[3]) > L[E_ord[3],3]){ 
      #   RL[j] <- count  
      #   break
      # }else if( abs(E_sort[2]) > L[E_ord[2],2] ){
      #   RL[j] <- count  
      #   break        
      # }else if( abs(E_sort[1]) > L[E_ord[1],1] ){
      #   RL[j] <- count  
      #   break        
      # }
      
      if(E_max_sort[3] > L_U[E_max_ord[3],3,d] || E_min_sort[1] < L_L[E_min_ord[1],3,d]){ 
        # if( E_max_sort[2] < L_U[E_max_ord[2],2,d] || E_min_sort[2] > L_L[E_min_ord[2],2,d]){
        #   if( E_max_sort[1] < L_U[E_max_ord[1],1,d] || E_min_sort[3] > L_L[E_min_ord[3],1,d]){
        #        
        #   }else{
        #     #RL[j] <- count  
        #     #break 
        #   }
        # }else{
        #   #RL[j] <- count  
        #   #break 
        # }
        RL[j] <- count  
        break  
      }else{
        #RL[j] <- count  
        #break        
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
      E_x_stand <- (E_x - mu0)/Sx
      
      temp_x_max <- E_x_max
      E_x_max <- lam*Y0[1,1] + (1-lam)*temp_x_max
      E_x_max <- max(mu0,E_x_max)
      E_x_max_stand <- (E_x_max - mu0)/Sx
      
      temp_x_min <- E_x_min
      E_x_min <- lam*Y0[1,1] + (1-lam)*temp_x_min
      E_x_min <- min(mu0,E_x_min)
      E_x_min_stand <- (E_x_min - mu0)/Sx
      
      ##---- Binomial ----------
      temp_p <- E_p
      ## Calculate EWMA charting statistic and standardize
      E_p <- lam*Y0[1,2] + (1-lam)*temp_p
      E_p_stand <- (E_p - p0)/Sp
      
      temp_p_max <- E_p_max
      E_p_max <- lam*Y0[1,2] + (1-lam)*temp_p_max
      E_p_max <- max(p0,E_p_max)
      E_p_max_stand <- (E_p_max - p0)/Sp
      
      temp_p_min <- E_p_min
      E_p_min <- lam*Y0[1,2] + (1-lam)*temp_p_min
      E_p_min <- min(p0,E_p_min)
      E_p_min_stand <- (E_p_min - p0)/Sp
      
      ##---- Poisson ----------
      temp_c <- E_c
      ## Calculate EWMA charting statistic and standardize
      E_c <- lam*Y0[1,3] + (1-lam)*temp_c
      E_c_stand <- (E_c - c0)/Sc
      
      temp_c_max <- E_c_max
      E_c_max <- lam*Y0[1,3] + (1-lam)*temp_c_max
      E_c_max <- max(c0,E_c_max)
      E_c_max_stand <- (E_c_max - c0)/Sc
      
      temp_c_min <- E_c_min
      E_c_min <- lam*Y0[1,3] + (1-lam)*temp_c_min
      E_c_min <- min(c0,E_c_min)
      E_c_min_stand <- (E_c_min - c0)/Sc
      
      E <- c(E_x_stand,E_p_stand,E_c_stand)
      E_max <- c(E_x_max_stand,E_p_max_stand,E_c_max_stand)
      E_min <- c(E_x_min_stand,E_p_min_stand,E_c_min_stand)
      
      
      E_sort <- sort(E)
      E_ord <- order(E)
      
      E_max_sort <- sort(E_max)
      E_max_ord <- order(E_max)
      
      E_min_sort <- sort(E_min)
      E_min_ord <- order(E_min)
      
    }
    
  }
  
  
  c_df <- data.frame(rho=rho123[d],RL=RL)
  ## save result to df.csv file
  #write_csv(c_df ,file=paste("result/one-side-EWMA/RL/Xpc-ARL0-Holm-all-RL-all-",d,".csv",
  #                             sep = ""),append = TRUE,col_names = TRUE)
  
  ARL[d] <- mean(RL)
  RL_se[d] <- sd(RL)/sqrt(N)
  
  
}

df <- data.frame(rho=rho123,ARL=ARL,se=RL_se)
df
write_csv(df,file=paste("result/one-side-EWMA/Xpc-ARL0-Holm-all.csv",sep = ""),
          append = TRUE,col_names = TRUE)  

end=Sys.time()
end-start
