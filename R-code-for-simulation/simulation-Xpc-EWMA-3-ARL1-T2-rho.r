rm(list=ls())

start=Sys.time()
format(start,"%H:%M:%S")

## Set multiplier
L_t <- c(10.78,10.86,10.81,10.82)

## Set R matrix
r <- 0.1
R <- matrix(0,3,3) 
diag(R) <- r

I <- matrix(0,3,3) 
diag(I) <- 1

## original rho
rho123 <- c(-0.3,0,0.3,0.6)

## corresponding rho
rho12 <- c(-0.40,0,0.40,0.80) ; 
rho13 <- c(-0.31,0,0.31,0.62) ; 
rho23 <- c(-0.42,0,0.42,0.78) 


## Set In-control parameter
mu0 <- 0 ; p0 <- 0.3 ; c0 <- 3
sig0 <- 1

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


for(d in 1:length(rho123)){
  
  rho <- rho123[d]
  
  cou <- 1
  for(s in 1:5){
    for(l in 1:5){
      for (o in 1:5) {

        M <- 1
        ARL <- numeric() ; c.ind_H <- matrix(0,M,7) ; c.ind_c <- matrix(0,M,7)
        time_elapse <- numeric()
        
        df <- matrix(0,M,16)
        df <- data.frame(df)
        names(df) <- c("c1","c2","c3","c12","c13","c23","c123",
                       "x1","x2","x3","x12","x13","x23","x123","ARL","time_elapse")
        
        for(k in 1:M){
          
          RL <- numeric() 
          T_1 <- numeric() ; T_2.1 <- numeric() ; T_3.1_2 <- numeric()
          start_time <- proc.time()
          N <- 100000
          c_ind <- matrix(0,N,3)
          
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
            
            ## Set In-control parameter
            sig_x0 <- 1 ; sig_p0 <- sqrt(p0*(1-p0)) ; sig_c0 <- sqrt(c0)
            
            ## Set Out-of-control parameter
            mu1 <- mu1_v[s] ; p1 <- p1_v[l] ; c1 <- lambda1_v[o]
            sig_x1 <- 1 ; sig_p1 <- sqrt(p1*(1-p1)) ; sig_c1 <- sqrt(c1)
            
            #set.seed(count)
            ## Generate simulated observations
            Z <- array(0,dim = c(n,3,m)) ; U <- array(0,dim = c(n,3,m)) ;
            Y1 <- array(0,dim = c(n,3,m))
            for(ss in 1:m){
              Z[,,ss] <- rmvnorm(n,muZ,CovZ)
              U[,,ss] <- pnorm(Z[,,ss])
              Y1[,,ss] <- matrix( cbind(qnorm(U[,1,ss],mu1,1),qbinom(U[,2,ss],1,p1),
                                        qpois(U[,3,ss],c1)),n,3 )
            }
            
            Y1 <- apply(Y1,c(1,2),mean)
            Y1 <- matrix(Y1,3,1)
            
            
            ## Compute EWMA statistic
            W0 <- matrix(c(mu0,p0,c0),3,1)
            E_W <- R %*% Y1 + (I-R) %*% W0
            
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
              
              if(T2 > L_t[d]){
                RL[j] <- count  
                
                ## Decomposition start
                Y_bar <- matrix(E_W,3,1)
                mu_Z <- matrix(c(mu0,p0,c0),3,1)
                
                ## Compute T_1
                T_1[j] <- (Y_bar[1,1]-mu_Z[1,1])^2 / S_W[1,1]
                
                ## Compute T_2.1
                b2 <- S_W[1,2] / S_W[1,1] 
                X_2.1 <- mu_Z[2,1] + b2 * (Y_bar[1,1]-mu_Z[1,1])
                s_2.1 <- S_W[2,2] - (S_W[1,2])*(S_W[1,1])^(-1)*(S_W[1,2])
                T_2.1[j] <- (Y_bar[2,1]-X_2.1)^2 / s_2.1
                
                
                ## Compute T_3.1,2
                Sxx <- matrix(c(S_W[1,1],S_W[1,2],S_W[2,1],S_W[2,2]),2,2)
                Sx <- matrix(c(S_W[1,3],S_W[2,3]),2,1) 
                b3 <- solve(Sxx) %*% Sx 
                x_2 <- matrix(c(Y_bar[1,1],Y_bar[2,1]),2,1) ; 
                x_2_bar <- matrix(c(mu_Z[1,1],mu_Z[2,1]),2,1)
                X_3.1 <- mu_Z[3,1] + t(b3) %*% (x_2 - x_2_bar )
                s_3.1 <- S_W[3,3] - t(Sx) %*% solve(Sxx) %*% (Sx)
                T_3.1_2[j] <- (Y_bar[3,1]-X_3.1)^2 / s_3.1
                
                La <- qchisq(0.995,1)
                if(T_1[j] > La){
                  c_ind[j,1] <- 1
                }
                if(T_2.1[j] > La){
                  c_ind[j,2] <- 1
                }
                if(T_3.1_2[j] > La){
                  c_ind[j,3] <- 1
                }
                
                break
              }
              
              count <- count + 1
              #set.seed(count)
              ## Generate simulated observations
              Z <- array(0,dim = c(n,3,m)) ; U <- array(0,dim = c(n,3,m)) ;
              Y1 <- array(0,dim = c(n,3,m))
              for(ss in 1:m){
                Z[,,ss] <- rmvnorm(n,muZ,CovZ)
                U[,,ss] <- pnorm(Z[,,ss])
                Y1[,,ss] <- matrix( cbind(qnorm(U[,1,ss],mu1,1),qbinom(U[,2,ss],1,p1),
                                          qpois(U[,3,ss],c1)),n,3 )
              }
              
              Y1 <- apply(Y1,c(1,2),mean)
              Y1 <- matrix(Y1,3,1)
              
              n1 <- n1+1
              ## Compute EWMA statistic
              temp_w <- E_W
              E_W <- R %*% Y1 + (I-R) %*% temp_w
              
              ## Compute T2
              S_W <- r*(1-(1-r)^(2*n1))/(2-r) * S_Y
              T2 <- t(E_W-mu_Y) %*% solve(S_W) %*% (E_W-mu_Y)
              
            }
            
          }
          
          
          ## save 10000 RL result to csv file 
          # c_df <- matrix(cbind(rho,c_ind,RL),N,5)
          # c_df <- data.frame(c_df)
          # names(c_df) <- c("rho","Ix","Ip","Ic","RL")
          # if(k == 1){
          #   write_csv(c_df ,file=paste("result/rho/T2-RL/Xpc-ARL1-T2-all-RL-",s-1,
          #                              "-",l-1,"-",o-1,".csv",sep = ""),
          #             append = TRUE,col_names = TRUE)
          # }else{
          #   write_csv(c_df ,file=paste("result/rho/T2-RL/Xpc-ARL1-T2-all-RL-",s-1,
          #                              "-",l-1,"-",o-1,".csv",sep = ""),
          #             append = TRUE)
          # }
          
          c.ind_H_1 <- which(rowSums(c_ind)==1)
          c.ind_H_2 <- which(rowSums(c_ind)==2)
          c.ind_H_3 <- which(rowSums(c_ind)==3)
        
          ## Count index
          c.ind_c[k,1] <- length(which(c_ind[c.ind_H_1,1]==1))
          c.ind_c[k,2] <- length(which(c_ind[c.ind_H_1,2]==1))
          c.ind_c[k,3] <- length(which(c_ind[c.ind_H_1,3]==1))
          
          c.ind_c[k,4] <- length(which(c_ind[c.ind_H_2,1]==1 & c_ind[c.ind_H_2,2]==1))
          c.ind_c[k,5] <- length(which(c_ind[c.ind_H_2,1]==1 & c_ind[c.ind_H_2,3]==1))
          c.ind_c[k,6] <- length(which(c_ind[c.ind_H_2,2]==1 & c_ind[c.ind_H_2,3]==1))
          
          c.ind_c[k,7] <- length(which(c_ind[c.ind_H_3,1]==1 & c_ind[c.ind_H_3,2]==1 &
                                         c_ind[c.ind_H_3,3]==1))
          
          ## Compute average of RL
          c.ind_H[k,1] <- mean(RL[c.ind_H_1][which(c_ind[c.ind_H_1,1]==1)])
          c.ind_H[k,2] <- mean(RL[c.ind_H_1][which(c_ind[c.ind_H_1,2]==1)])
          c.ind_H[k,3] <- mean(RL[c.ind_H_1][which(c_ind[c.ind_H_1,3]==1)])
        
          c.ind_H[k,4] <- mean(RL[c.ind_H_2][which(c_ind[c.ind_H_2,1]==1 & 
                                                     c_ind[c.ind_H_2,2]==1)])
        
          c.ind_H[k,5] <- mean(RL[c.ind_H_2][which(c_ind[c.ind_H_2,1]==1 & 
                                                     c_ind[c.ind_H_2,3]==1)])
        
          c.ind_H[k,6] <- mean(RL[c.ind_H_2][which(c_ind[c.ind_H_2,2]==1 & 
                                                     c_ind[c.ind_H_2,3]==1)])
        
        
          c.ind_H[k,7] <- mean(RL[c.ind_H_3][which(c_ind[c.ind_H_3,1]==1 & 
                                                     c_ind[c.ind_H_3,2]==1 &
                                                      c_ind[c.ind_H_3,3]==1)])
          for(E in 1:7){
            if(c.ind_H[k,E] == "NaN"){ c.ind_H[k,E] <- 0 }
          }
          ARL[k] <- mean(RL)
        
          
          time <- proc.time() - start_time
          time_elapse[k] <- time[3]
          
          ## save result to df
          df[k,] <- cbind(c.ind_c[k,1],c.ind_c[k,2],c.ind_c[k,3],c.ind_c[k,4],
                          c.ind_c[k,5],c.ind_c[k,6],c.ind_c[k,7],
                          c.ind_H[k,1],c.ind_H[k,2],c.ind_H[k,3],c.ind_H[k,4],
                          c.ind_H[k,5],c.ind_H[k,6],c.ind_H[k,7], ARL[k], time_elapse[k])
          
          
        }
      
      cas <- matrix(0,1,3)
      cas[1,1] <- s-1 ; cas[1,2] <- l-1 ; cas[1,3] <- o-1
      cas <- data.frame(cas)
      names(cas) <- c("Normal","Binomial","Poisson")
      
      df_res <- matrix(colMeans(df),1,16)
      df_res <- data.frame(df_res)
      names(df_res) <- c("c1","c2","c3","c12","c13","c23","c123",
                     "x1","x2","x3","x12","x13","x23","x123","ARL","time_elapse")
      
      c_df <- cbind(rho,cas,df_res)
      
      ## save result to df.csv file
      if(s == 1 & l == 1 & o == 1){
        write_csv(c_df ,file=paste("result/rho/Xpc-ARL1-all-T2-rho-3.csv",sep = ""),append = TRUE,
                  col_names = TRUE)
      }else{
        write_csv(c_df ,file=paste("result/rho/Xpc-ARL1-all-T2-rho-3.csv",sep = ""),append = TRUE)
      }

    }
  }
}

}

end=Sys.time()
end-start
