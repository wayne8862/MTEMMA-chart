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
L1_U <- L1_U - 0.01

Lx2_U <- c(2.67,2.92,3.06) ; Lp2_U <- c(2.88,3.14,3.29) ; Lc2_U <- c(2.85,3.14,3.31)
L2_U <- matrix(rbind(Lx2_U,Lp2_U,Lc2_U),3,3)
#L2_U <- L2_U + 0.04

Lx3_U <- c(2.67,2.92,3.06) ; Lp3_U <- c(2.88,3.14,3.29) ; Lc3_U <- c(2.85,3.14,3.31)
L3_U <- matrix(rbind(Lx3_U,Lp3_U,Lc3_U),3,3)
L3_U <- L3_U + 0.01*1

Lx4_U <- c(2.67,2.92,3.06) ; Lp4_U <- c(2.88,3.14,3.29) ; Lc4_U <- c(2.85,3.14,3.31)
L4_U <- matrix(rbind(Lx4_U,Lp4_U,Lc4_U),3,3)
L4_U <- L4_U - 0.01*4

L_U <- array(0,dim = c(3,3,4))
L_U[,,1] <- L1_U
L_U[,,2] <- L2_U
L_U[,,3] <- L3_U
L_U[,,4] <- L4_U

## Set multiple values for individual EWMA chart ARL_0 =400, 800, 1200
Lx1_L <- c(-2.66,-2.92,-3.06) ; Lp1_L <- c(-2.215,-2.324,-2.375) ; Lc1_L <- c(-2.46,-2.68,-2.79)
L1_L <- matrix(rbind(Lx1_L,Lp1_L,Lc1_L),3,3)
L1_L <- L1_L + 0.01*1

Lx2_L <- c(-2.66,-2.92,-3.06) ; Lp2_L <- c(-2.215,-2.324,-2.375) ; Lc2_L <- c(-2.46,-2.68,-2.79)
L2_L <- matrix(rbind(Lx2_L,Lp2_L,Lc2_L),3,3)
#L2_L <- L2_L - 0.01*4

Lx3_L <- c(-2.66,-2.92,-3.06) ; Lp3_L <- c(-2.215,-2.324,-2.375) ; Lc3_L <- c(-2.46,-2.68,-2.79)
L3_L <- matrix(rbind(Lx3_L,Lp3_L,Lc3_L),3,3)
L3_L <- L3_L + 0.01*1

Lx4_L <- c(-2.66,-2.92,-3.06) ; Lp4_L <- c(-2.215,-2.324,-2.375) ; Lc4_L <- c(-2.46,-2.68,-2.79)
L4_L <- matrix(rbind(Lx4_L,Lp4_L,Lc4_L),3,3)
L4_L <- L4_L + 0.01*4

L_L <- array(0,dim = c(3,3,4))
L_L[,,1] <- L1_L
L_L[,,2] <- L2_L
L_L[,,3] <- L3_L
L_L[,,4] <- L4_L

## Set In-control parameter
mu0 <- 0 ; p0 <- 0.3 ; c0 <- 3
sig0 <- 1
pa <- c(mu0,p0,c0)

## Out-of-control parameter 
s <- c(0.2,0.3,0.4,0.5)

sn <- sig0*s
mu1_v <- round(c(mu0,mu0+sn),2)  

sp <- c(0.05,0.10,0.15,0.20)
p1_v <- round(c(p0,p0+sp),2)   

sc <- c(0.2,0.4,0.6,0.8)
lambda1_v <- round(c(c0,c0+sc),2)  


## Require package
library(mvtnorm)
library(readr)

for(d in c(1,2,3,4)){
  
  rho <- rho123[d]
  
  cou <- 1
  for(s in 1:5){
    for(l in 1:5){
      for (o in 1:5) {
      
        M <- 1
        ARL <- numeric() ; 
        ct <- matrix(0,M,7) ; Rt <- matrix(0,M,7)
        time_elapse <- numeric()

        
        df <- matrix(0,M,16)
        df <- data.frame(df)
        names(df) <- c("c1","c2","c3","c12","c13","c23","c123",
                       "x1","x2","x3","x12","x13","x23","x123",
                       "ARL","time_elapse")
      
        for(k in 1:M){
        
        RL <- numeric() 
        start_time <- proc.time()
        N <- 100000
        c_ind_up <- matrix(0,N,3) ; c_ind_down <- matrix(0,N,3)
        
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
          
          ## Set In-control parameter
          sig_x0 <- 1 ; sig_p0 <- sqrt(p0*(1-p0)) ; sig_c0 <- sqrt(c0)
          
          ## Set Out-of-control parameter
          mu1 <- mu1_v[s] ; p1 <- p1_v[l] ; c1 <- lambda1_v[o]
          sig_x1 <- 1 ; sig_p1 <- sqrt(p1*(1-p1)) ; sig_c1 <- sqrt(c1)
          
          
          ## Compute standard deviation
          Sx <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n)))*sig_x0/sqrt(m)
          Sp <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n)))*sig_p0/sqrt(m)
          Sc <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n)))*sig_c0/sqrt(m)
          
          
          ## Generate simulated observations
          Z <- array(0,dim = c(n,3,m)) ; U <- array(0,dim = c(n,3,m)) ;
          Y1 <- array(0,dim = c(n,3,m))
          for(ss in 1:m){
            Z[,,ss] <- rmvnorm(n,muZ,CorrZ)
            U[,,ss] <- pnorm(Z[,,ss])
            Y1[,,ss] <- matrix( cbind(qnorm(U[,1,ss],mu1,1),qbinom(U[,2,ss],1,p1),
                                      qpois(U[,3,ss],c1)),n,3 )
          }
          
          Y1 <- apply(Y1,c(1,2),mean)
          
          ##----- Normal ------------
          E_x0 <- mu0
          ## Calculate EWMA charting statistic and standardize
          E_x <- lam*Y1[1,1] + (1-lam)*E_x0
          E_x_stand <- (E_x - mu0)/Sx
          
          E_x_max <- max(mu0,E_x)
          E_x_max_stand <- (E_x_max - mu0)/Sx
          
          E_x_min <- min(mu0,E_x)
          E_x_min_stand <- (E_x_min - mu0)/Sx
          
          ##---- Binomial ----------
          E_p0 <- p0
          ## Calculate EWMA charting statistic and standardize
          E_p <- lam*Y1[1,2] + (1-lam)*E_p0
          E_p_stand <- (E_p - p0)/Sp
          
          E_p_max <- max(p0,E_p)
          E_p_max_stand <- (E_p_max - p0)/Sp
          
          E_p_min <- min(p0,E_p)
          E_p_min_stand <- (E_p_min - p0)/Sp
          
          ##---- Poisson ----------
          E_c0 <- c0
          ## Calculate EWMA charting statistic and standardize
          E_c <- lam*Y1[1,3] + (1-lam)*E_c0
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
            
            ## Bon
            if( (E_max_sort[1] > L_U[E_max_ord[1],3,d] || E_max_sort[2] > L_U[E_max_ord[2],3,d] 
                 || E_max_sort[3] > L_U[E_max_ord[3],3,d])  
                 ||(E_min_sort[1] < L_L[E_min_ord[1],3,d] || E_min_sort[2] < L_L[E_min_ord[2],3,d] 
                  || E_min_sort[3] < L_L[E_min_ord[3],3,d])
               ){
              
              RL[j] <- count
              if(E_max_sort[3] > L_U[E_max_ord[3],3,d]){c_ind_up[j,E_max_ord[3]] <- 1}
              if(E_max_sort[2] > L_U[E_max_ord[2],3,d]){c_ind_up[j,E_max_ord[2]] <- 1}
              if(E_max_sort[1] > L_U[E_max_ord[1],3,d]){c_ind_up[j,E_max_ord[1]] <- 1}
              
              if(E_min_sort[1] < L_L[E_min_ord[1],3,d]){c_ind_down[j,E_min_ord[1]] <- 1}
              if(E_min_sort[2] < L_L[E_min_ord[2],3,d]){c_ind_down[j,E_min_ord[2]] <- 1}
              if(E_min_sort[3] < L_L[E_min_ord[3],3,d]){c_ind_down[j,E_min_ord[3]] <- 1}
              
              break
            }
            
            count <- count + 1
            
            ## Generate simulated observations
            Z <- array(0,dim = c(n,3,m)) ; U <- array(0,dim = c(n,3,m)) ;
            Y1 <- array(0,dim = c(n,3,m))
            for(ss in 1:m){
              Z[,,ss] <- rmvnorm(n,muZ,CorrZ)
              U[,,ss] <- pnorm(Z[,,ss])
              Y1[,,ss] <- matrix( cbind(qnorm(U[,1,ss],mu1,1),qbinom(U[,2,ss],1,p1),
                                        qpois(U[,3,ss],c1)),n,3 )
            }
            
            Y1 <- apply(Y1,c(1,2),mean)
            
            n1 <- n1+1
            ## Compute standard deviation
            Sx <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n1)))*sig_x0/sqrt(m)
            Sp <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n1)))*sig_p0/sqrt(m)
            Sc <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n1)))*sig_c0/sqrt(m)
            
            ##----- Normal ------------
            temp_x <- E_x
            ## Calculate EWMA charting statistic and standardize
            E_x <- lam*Y1[1,1] + (1-lam)*temp_x
            E_x_stand <- (E_x - mu0)/Sx
            
            temp_x_max <- E_x_max
            E_x_max <- lam*Y1[1,1] + (1-lam)*temp_x_max
            E_x_max <- max(mu0,E_x_max)
            E_x_max_stand <- (E_x_max - mu0)/Sx
            
            temp_x_min <- E_x_min
            E_x_min <- lam*Y1[1,1] + (1-lam)*temp_x_min
            E_x_min <- min(mu0,E_x_min)
            E_x_min_stand <- (E_x_min - mu0)/Sx
            
            ##---- Binomial ----------
            temp_p <- E_p
            ## Calculate EWMA charting statistic and standardize
            E_p <- lam*Y1[1,2] + (1-lam)*temp_p
            E_p_stand <- (E_p - p0)/Sp
            
            temp_p_max <- E_p_max
            E_p_max <- lam*Y1[1,2] + (1-lam)*temp_p_max
            E_p_max <- max(p0,E_p_max)
            E_p_max_stand <- (E_p_max - p0)/Sp
            
            temp_p_min <- E_p_min
            E_p_min <- lam*Y1[1,2] + (1-lam)*temp_p_min
            E_p_min <- min(p0,E_p_min)
            E_p_min_stand <- (E_p_min - p0)/Sp
            
            ##---- Poisson ----------
            temp_c <- E_c
            ## Calculate EWMA charting statistic and standardize
            E_c <- lam*Y1[1,3] + (1-lam)*temp_c
            E_c_stand <- (E_c - c0)/Sc
            
            temp_c_max <- E_c_max
            E_c_max <- lam*Y1[1,3] + (1-lam)*temp_c_max
            E_c_max <- max(c0,E_c_max)
            E_c_max_stand <- (E_c_max - c0)/Sc
            
            temp_c_min <- E_c_min
            E_c_min <- lam*Y1[1,3] + (1-lam)*temp_c_min
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
        
        
        C <- matrix(cbind(c_ind_up,c_ind_down,RL),N,7)
        C <- data.frame(C)
        names(C) <- c("c1","c2","c3","d1","d2","d3","RL")
        
        ## count
        ct[k,1] <- length(which( ((C$c1==1)|(C$d1==1)) &  
                                   ((C$c2==0)&(C$c3==0)&(C$d2==0)&(C$d3==0)) ))
        ct[k,2] <- length(which( ((C$c2==1)|(C$d2==1)) &  
                                   ((C$c1==0)&(C$c3==0)&(C$d1==0)&(C$d3==0)) ))
        ct[k,3] <- length(which( ((C$c3==1)|(C$d3==1)) &  
                                   ((C$c1==0)&(C$c2==0)&(C$d1==0)&(C$d2==0)) ))
        
        ct[k,4] <- length(which( ((C$c1==1)&(C$c2==1))|((C$d1==1)&(C$d2==1))|
                                 ((C$c1==1)&(C$d2==1))|((C$c2==1)&(C$d1==1))  
                                 & ((C$c3==0)&(C$d3==0)) ))
        ct[k,5] <- length(which( ((C$c1==1)&(C$c3==1))|((C$d1==1)&(C$d3==1))|
                                   ((C$c1==1)&(C$d3==1))|((C$c3==1)&(C$d1==1))  
                                 & ((C$c2==0)&(C$d2==0)) ))
        ct[k,6] <- length(which( ((C$c2==1)&(C$c3==1))|((C$d2==1)&(C$d3==1))|
                                   ((C$c2==1)&(C$d3==1))|((C$c3==1)&(C$d2==1))  
                                 & ((C$c1==0)&(C$d1==0)) ))
        
        ct[k,7] <- length(which( ((C$c1==1)&(C$c2==1)&(C$c3==1))|
                                   ((C$d1==1)&(C$d2==1)&(C$d3==1)) ))
        
        ## RL
        Rt[k,1] <- mean(RL[which( ((C$c1==1)|(C$d1==1)) &  
                                   ((C$c2==0)&(C$c3==0)&(C$d2==0)&(C$d3==0)) )])
        Rt[k,2] <- mean(RL[which( ((C$c2==1)|(C$d2==1)) &  
                                   ((C$c1==0)&(C$c3==0)&(C$d1==0)&(C$d3==0)) )])
        Rt[k,3] <- mean(RL[which( ((C$c3==1)|(C$d3==1)) &  
                                   ((C$c1==0)&(C$c2==0)&(C$d1==0)&(C$d2==0)) )])
        
        Rt[k,4] <- mean(RL[which( ((C$c1==1)&(C$c2==1))|((C$d1==1)&(C$d2==1))|
                                   ((C$c1==1)&(C$d2==1))|((C$c2==1)&(C$d1==1))  
                                 & ((C$c3==0)&(C$d3==0)) )])
        Rt[k,5] <- mean(RL[which( ((C$c1==1)&(C$c3==1))|((C$d1==1)&(C$d3==1))|
                                   ((C$c1==1)&(C$d3==1))|((C$c3==1)&(C$d1==1))  
                                 & ((C$c2==0)&(C$d2==0)) )])
        Rt[k,6] <- mean(RL[which( ((C$c2==1)&(C$c3==1))|((C$d2==1)&(C$d3==1))|
                                   ((C$c2==1)&(C$d3==1))|((C$c3==1)&(C$d2==1))  
                                 & ((C$c1==0)&(C$d1==0)) )])
        
        Rt[k,7] <- mean(RL[which( ((C$c1==1)&(C$c2==1)&(C$c3==1))|
                                  ((C$c1==1)&(C$c2==1)&(C$d3==1))|
                                  ((C$c1==1)&(C$d2==1)&(C$c3==1))|
                                  ((C$c1==1)&(C$d2==1)&(C$d3==1))|
                                  ((C$d1==1)&(C$c2==1)&(C$c3==1))|
                                  ((C$d1==1)&(C$c2==1)&(C$d3==1))|
                                  ((C$d1==1)&(C$d2==1)&(C$c3==1))|  
                                  ((C$d1==1)&(C$d2==1)&(C$d3==1))  )])
        for(f in 1:7){
          if(Rt[k,f] == "NaN"){ Rt[k,f] <- 0 }
        }
        
        ARL[k] <- mean(RL)
        
        
        time <- proc.time() - start_time
        time_elapse[k] <- time[3]
        
        ## save result to df
        df[k,] <- cbind(ct[k,1],ct[k,2],ct[k,3],ct[k,4],ct[k,5],ct[k,6],ct[k,7],
                        Rt[k,1],Rt[k,2],Rt[k,3],Rt[k,4],Rt[k,5],Rt[k,6],Rt[k,7],
                        ARL[k], time_elapse[k])
        
      }
      
      cas <- matrix(0,1,3)
      cas[1,1] <- s-1 ; cas[1,2] <- l-1 ; cas[1,3] <- o-1
      cas <- data.frame(cas)
      names(cas) <- c("Normal","Binomial","Poisson")
      
      df_res <- matrix(colMeans(df),1,16)
      df_res <- data.frame(df_res)
      names(df_res) <- c("c1","c2","c3","c12","c13","c23","c123",
                         "x1","x2","x3","x12","x13","x23","x123",
                         "ARL","time_elapse")
      
      c_df <- cbind(rho,cas,df_res)
      
      ## save result to df.csv file
      if(s == 1 & l == 1 & o == 1){
        write_csv(c_df ,file=paste("result/one-side-EWMA/Xpc-5-ARL1-all-Bon-rho-5.csv",sep = ""),append = TRUE,
                  col_names = TRUE)
      }else{
        write_csv(c_df ,file=paste("result/one-side-EWMA/Xpc-5-ARL1-all-Bon-rho-5.csv",sep = ""),append = TRUE)
      }
      
    }
  }
}

}

end=Sys.time()
end-start
