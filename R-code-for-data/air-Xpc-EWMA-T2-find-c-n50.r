rm(list=ls())

start=Sys.time()
format(start,"%H:%M:%S")

# input data
regTaiwan <- read.csv("DataExample/taiwan-with-PM10quality.csv",header = T)
regTaiwan <- data.frame(regTaiwan)

# ARRANGE THE DATA AS A LIST OF DATA SETS
regionsTaiwan <- as.character(unique(regTaiwan$regnames))
dlistTaiwan <- lapply(regionsTaiwan,function(x) regTaiwan[regTaiwan$regnames==x,])
names(dlistTaiwan) <- regionsTaiwan

df_org <- dlistTaiwan[[1]]

#df <- df_org[(which(df_org$year==2001)[1]:5113),]
df <- df_org

n0 <- 50

## Normal
z1 <- matrix(df$o3[1:(n0*30)],n0,30)
x1 <- rowMeans(z1)
mu0 <- mean(x1)
sig0 <- sd(x1)

## Binomial
z2 <- matrix(df$pm10_quality[1:(n0*30)],n0,30)
x2 <- rowMeans(z2)
p0 <- mean(x2)

## Poisson
z3 <- matrix(df$resp[1:(n0*30)],n0,30)
x3 <- rowMeans(z3)
c0 <- mean(x3)


## correlation
rho12 <- cor(x1,x2)
rho13 <- cor(x1,x3)
rho23 <- cor(x2,x3)

## Set R matrix
r <- 0.1
R <- matrix(0,3,3) 
diag(R) <- r

I <- matrix(0,3,3) 
diag(I) <- 1

## Set multiple values
L_t <- seq(16.28,18.09,by=0.01) 


## Require package
library(mvtnorm)
library(readr)

ARL_t <- numeric() 
time_elapse <- numeric()

for(k in 1:length(L_t)){
  
  RL_x <- numeric() 
  start_time <- proc.time()
  N <- 100000
  for(j in 1:N){
    
    ## Set count for ARL_0 
    count_x <- 1 ; count_p <- 1 ; count_c <- 1
    
    ## Setting n (observations) ; m (subgroup size)
    n <- 1 ;  m <- 1  

    sig_x0 <- sig0 ; sig_p0 <- sqrt(p0*(1-p0)) ; sig_c0 <- sqrt(c0)
    
    #set.seed(count)
    ## Generate simulated observations
    Y0 <- matrix(c(sample(x1,1,replace = T),sample(z2,1,replace = T),
                   sample(z3,1,replace = T)),1,3)
    Y0 <- matrix(Y0,3,1)
    ## Compute EWMA statistic
    W0 <- matrix(c(mu0,p0,c0),3,1)
    E_W <- R %*% Y0 + (I-R) %*% W0
    
    ## Compute Y   muY and Sigma_Y
    CovY <- matrix(0,3,3) 
    CovY[1,2] <- rho12*sig_x0*sqrt(p0*(1-p0))      
    CovY[1,3] <- rho13*sig_x0*sqrt(c0) 
    CovY[2,3] <- rho23*sqrt(p0*(1-p0))*sqrt(c0) 
    CovY <- CovY + t(CovY)          
    diag(CovY) <- c(sig_x0^2,p0*(1-p0),c0)
    
    mu_Y <- matrix(c(mu0,p0,c0),3,1) ; S_Y <- matrix(CovY,3,3)
    ## Compute T2
    S_W <- r*(1-(1-r)^(2*n))/(2-r) * S_Y
    T2 <- t(E_W-mu_Y) %*% solve(S_W) %*% (E_W-mu_Y)    
    
    
    n1 <- 1
    repeat{
      
      if(T2 > L_t[k]){
        RL_x[j] <- count_x
        break
      }
      count_x <- count_x +1
      
      ## Generate simulated observations
      Y0 <- matrix(c(sample(x1,1,replace = T),sample(z2,1,replace = T),
                     sample(z3,1,replace = T)),1,3)
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
  
  ARL_t[k] <- mean(RL_x)

  time <- proc.time() - start_time
  time_elapse[k] <- time[3]
  
  ## save result to df.csv file
  if(k == 1){
    df <- data.frame( Lt = L_t[k], ARL = ARL_t[k], time_elapse=time_elapse[k])
    write_csv(df,file=paste("dataExample/result/T2-Find-c-n50.csv",sep = ""),append = TRUE,
              col_names = TRUE)  
  }else{
    df <- data.frame( Lt = L_t[k], ARL = ARL_t[k], time_elapse=time_elapse[k])
    write_csv(df,file=paste("dataExample/result/T2-Find-c-n50.csv",sep = ""),append = TRUE)
  }

}

end=Sys.time()
end-start


