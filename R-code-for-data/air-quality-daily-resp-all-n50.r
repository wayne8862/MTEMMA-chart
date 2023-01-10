rm(list=ls())

# input data
regTaiwan <- read.csv("DataExample/taiwan-with-PM10quality.csv",header = T)
regTaiwan <- data.frame(regTaiwan)

# ARRANGE THE DATA AS A LIST OF DATA SETS
regionsTaiwan <- as.character(unique(regTaiwan$regnames))
dlistTaiwan <- lapply(regionsTaiwan,function(x) regTaiwan[regTaiwan$regnames==x,])
names(dlistTaiwan) <- regionsTaiwan

## Kaohsing (1)  ; Taipei (2) ; Taichung (3)   
## year 1994 ~ 2007
## N = 5113

## Phase I  36 samples (daily) of 1 observations
## Consider X1 = Ozone ; X2 = PM10_quality ; X3 = cvd

##----------------- x-chart --------------------------------------------
df_org <- dlistTaiwan[[1]]

#df <- df_org[(which(df_org$year==2001)[1]:5113),]
df <- df_org

n0 <- 50

set.seed(4)
## Normal
z1 <- matrix(df$o3[1:(n0*30)],n0,30)
x1 <- sample(z1,n0,replace = T)
mu0 <- mean(x1)
sig0 <- sd(x1)

## Binomial
z2 <- matrix(df$pm10_quality[1:(n0*30)],n0,30)
x2 <- sample(z2,n0,replace = T)
p0 <- mean(x2)

## Poisson
z3 <- matrix(df$resp[1:(n0*30)],n0,30)
x3 <- sample(z3,n0,replace = T)
c0 <- mean(x3)

rho12 <- cor(x1,x2)
rho13 <- cor(x1,x3)
rho23 <- cor(x2,x3)

rho12
rho13
rho23


Box.test(x1, lag=3,type = "Ljung-Box")
Box.test(x2, lag=3,type = "Ljung-Box")
Box.test(x3, lag=3,type = "Ljung-Box")



## Setting n (observations) ; m (subgroup size)
n <- 1 ;  m <- 1 

## Set In-control parameter
sig_x0 <- sd(x1) ; sig_p0 <- sqrt(p0*(1-p0)) ; sig_c0 <- sqrt(c0)

## Set Lambda
lam <- 0.1


E <- matrix(0,n0,3) ; E_U <- matrix(0,n0,3) ; E_L <- matrix(0,n0,3)
Y0 <- matrix(c(x1[n],x2[n],x3[n]),1,3)

## Compute standard deviation
Sx <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n)))*sig_x0/sqrt(m)
Sp <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n)))*sig_p0/sqrt(m)
Sc <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n)))*sig_c0/sqrt(m)

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

E[n,] <- c(E_x_stand,E_p_stand,E_c_stand)
E_U[n,] <- c(E_x_max_stand,E_p_max_stand,E_c_max_stand)
E_L[n,] <- c(E_x_min_stand,E_p_min_stand,E_c_min_stand)
  
n1 <- n
repeat{
  
  n1 <- n1 + 1
  
  if(n1 == 51){ break }
  
  Y1 <- matrix(c(x1[n1],x2[n1],x3[n1]),1,3)
  
  ## Compute standard deviation
  Sx <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n1)))*sig_x0/sqrt(m)
  Sp <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n1)))*sig_p0/sqrt(m)
  Sc <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n1)))*sig_c0/sqrt(m)
  
  ##----- Normal ------------
  temp_x <- E_x
  ## Calculate EWMA charting statistic and standardize
  E_x <- lam*Y1[1,1] + (1-lam)*temp_x
  E_x_stand <- (E_x - mu0)/Sx
  
  E_x_max <- max(mu0,E_x)
  E_x_max_stand <- (E_x_max - mu0)/Sx
  
  E_x_min <- min(mu0,E_x)
  E_x_min_stand <- (E_x_min - mu0)/Sx
  
  
  ##---- Binomial ----------
  temp_p <- E_p
  ## Calculate EWMA charting statistic and standardize
  E_p <- lam*Y1[1,2] + (1-lam)*temp_p
  E_p_stand <- (E_p - p0)/Sp
  
  E_p_max <- max(p0,E_p)
  E_p_max_stand <- (E_p_max - p0)/Sp
  
  E_p_min <- min(p0,E_p)
  E_p_min_stand <- (E_p_min - p0)/Sp
  
  ##---- Poisson ----------
  temp_c <- E_c
  ## Calculate EWMA charting statistic and standardize
  E_c <- lam*Y1[1,3] + (1-lam)*temp_c
  E_c_stand <- (E_c - c0)/Sc
  
  E_c_max <- max(c0,E_c)
  E_c_max_stand <- (E_c_max - c0)/Sc
  
  E_c_min <- min(c0,E_c)
  E_c_min_stand <- (E_c_min - c0)/Sc
  
  E[n1,] <- c(E_x_stand,E_p_stand,E_c_stand)
  E_U[n1,] <- c(E_x_max_stand,E_p_max_stand,E_c_max_stand)
  E_L[n1,] <- c(E_x_min_stand,E_p_min_stand,E_c_min_stand)
}


##===== Phase II ================

set.seed(77)
## Normal
z11 <- matrix(df$o3[(n0*30+1):(n0*60)],n0,30)
x11 <- sample(z11,n0,replace = T)
mu1 <- mean(x11)
sig1 <- sd(x11)

## Binomial
z12 <- matrix(df$pm10_quality[(n0*30+1):(n0*60)],n0,30)
x12 <- sample(z12,n0,replace = T)
p1 <- mean(x12)

## Poisson
z13 <- matrix(df$resp[(n0*30+1):(n0*60)],n0,30)
x13 <- sample(z13,n0,replace = T)
c1 <- mean(x13)


## Setting n (observations) ; m (subgroup size)
n <- 1 ;  m <- 1 

## Set In-control parameter
sig_x0 <- sd(x1) ; sig_p0 <- sqrt(p0*(1-p0)) ; sig_c0 <- sqrt(c0)

## Set Lambda
lam <- 0.1


E1 <- matrix(0,n0,3) ; E1_U <- matrix(0,n0,3) ; E1_L <- matrix(0,n0,3)
Y10 <- matrix(c(x11[n],x12[n],x13[n]),1,3)

## Compute standard deviation
Sx <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n)))*sig_x0/sqrt(m)
Sp <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n)))*sig_p0/sqrt(m)
Sc <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n)))*sig_c0/sqrt(m)

##----- Normal ------------
E_x0 <- mu0
## Calculate EWMA charting statistic and standardize
E_x1 <- lam*Y10[1,1] + (1-lam)*E_x0
E_x1_stand <- (E_x1 - mu0)/Sx

E_x1_max <- max(mu0,E_x1)
E_x1_max_stand <- (E_x1_max - mu0)/Sx

E_x1_min <- min(mu0,E_x1)
E_x1_min_stand <- (E_x1_min - mu0)/Sx


##---- Binomial ----------
E_p0 <- p0
## Calculate EWMA charting statistic and standardize
E_p1 <- lam*Y10[1,2] + (1-lam)*E_p0
E_p1_stand <- (E_p1 - p0)/Sp

E_p1_max <- max(p0,E_p1)
E_p1_max_stand <- (E_p1_max - p0)/Sp

E_p1_min <- min(p0,E_p1)
E_p1_min_stand <- (E_p1_min - p0)/Sp

##---- Poisson ----------
E_c0 <- c0
## Calculate EWMA charting statistic and standardize
E_c1 <- lam*Y10[1,3] + (1-lam)*E_c0
E_c1_stand <- (E_c1 - c0)/Sc

E_c1_max <- max(c0,E_c1)
E_c1_max_stand <- (E_c1_max - c0)/Sc

E_c1_min <- min(c0,E_c1)
E_c1_min_stand <- (E_c1_min - c0)/Sc


E1[n,] <- c(E_x1_stand,E_p1_stand,E_c1_stand)
E1_U[n,] <- c(E_x1_max_stand,E_p1_max_stand,E_c1_max_stand)
E1_L[n,] <- c(E_x1_min_stand,E_p1_min_stand,E_c1_min_stand)

n1 <- n
repeat{
  
  n1 <- n1 + 1
  
  if(n1 == 51){ break }
  
  Y11 <- matrix(c(x11[n1],x12[n1],x13[n1]),1,3)
  
  ## Compute standard deviation
  Sx <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n1)))*sig_x0/sqrt(m)
  Sp <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n1)))*sig_p0/sqrt(m)
  Sc <- sqrt(lam/(2-lam)*(1-(1-lam)^(2*n1)))*sig_c0/sqrt(m)
  
  ##----- Normal ------------
  temp_x <- E_x1
  ## Calculate EWMA charting statistic and standardize
  E_x1 <- lam*Y11[1,1] + (1-lam)*temp_x
  E_x1_stand <- (E_x1 - mu0)/Sx
  
  E_x1_max <- max(mu0,E_x1)
  E_x1_max_stand <- (E_x1_max - mu0)/Sx
  
  E_x1_min <- min(mu0,E_x1)
  E_x1_min_stand <- (E_x1_min - mu0)/Sx
  
  ##---- Binomial ----------
  temp_p <- E_p1
  ## Calculate EWMA charting statistic and standardize
  E_p1 <- lam*Y11[1,2] + (1-lam)*temp_p
  E_p1_stand <- (E_p1 - p0)/Sp
  
  E_p1_max <- max(p0,E_p1)
  E_p1_max_stand <- (E_p1_max - p0)/Sp
  
  E_p1_min <- min(p0,E_p1)
  E_p1_min_stand <- (E_p1_min - p0)/Sp
  
  ##---- Poisson ----------
  temp_c <- E_c1
  ## Calculate EWMA charting statistic and standardize
  E_c1 <- lam*Y11[1,3] + (1-lam)*temp_c
  E_c1_stand <- (E_c1 - c0)/Sc
  
  E_c1_max <- max(c0,E_c1)
  E_c1_max_stand <- (E_c1_max - c0)/Sc
  
  E_c1_min <- min(c0,E_c1)
  E_c1_min_stand <- (E_c1_min - c0)/Sc
  
  
  E1[n1,] <- c(E_x1_stand,E_p1_stand,E_c1_stand)
  E1_U[n1,] <- c(E_x1_max_stand,E_p1_max_stand,E_c1_max_stand)
  E1_L[n1,] <- c(E_x1_min_stand,E_p1_min_stand,E_c1_min_stand)
  
}

cbind(x11,x12,x13)

## Combine two E and E1

B <- rbind(E,E1)
B_U <- rbind(E_U,E1_U)
B_L <- rbind(E_L,E1_L)

## Holm
##------- Plot ---------------------- 
x <- c(1,seq(5,100,by=5))
y <- seq(-3.5,4,by=0.5)
plot(B[,1],type = "b",pch=16, xaxt = "n", yaxt = "n",
     xlim = c(1,100),ylim = c(-3.5,5.5),col="red",lty=1,lwd=1,
     xlab = "Sample Sequence",ylab = "value",cex.lab=1.8,cex.axis=1.3)
points(B[,2],col="blue",type = "b",pch=16,lty=1,lwd=1.3)
points(B[,3],type = "b",pch=16,lty=1,lwd=1.6)
axis(1, labels = x, at = x,cex.axis=1.1)
axis(2, labels = y, at = y)
abline(h=0,lty=1,col="gray",lwd=1.4)
legend(1.2,5.5,legend = c(expression(X[1]),expression(X[2]),expression(X[3])),
       col=c("red","blue",1),lty = c(1,1,1),lwd = c(1,1.3,1.6),cex = 0.8 )
abline(v=50)
text(43,5.15,"Phase I",cex = 3)
text(57,5.15,"Phase II",cex = 3)
abline(h=2.47,lty=1,col=6) # x
abline(h=2.76,lty=1,col=6)
abline(h=2.92,lty=1,col=6)
abline(h=2.42,lty=1,col="steelblue") # p
abline(h=2.68,lty=1,col="steelblue")
abline(h=2.83,lty=1,col="steelblue")
abline(h=2.80,lty=1,col=1) # p
abline(h=3.17,lty=1,col=1)
abline(h=3.37,lty=1,col=1)
text(1.7,2.24,expression(paste(UCL[x(1)]," = 2.47")),cex = 1,col=6)
text(8.7,2.24,expression(paste(UCL[x(2)]," = 2.76")),cex = 1,col=6)
text(15.7,2.24,expression(paste(UCL[x(3)]," = 2.92")),cex = 1,col=6)
text(1.7,3.06,expression(paste(UCL[p(1)]," = 2.42")),cex = 1,col="steelblue")
text(8.7,3.06,expression(paste(UCL[p(2)]," = 2.68")),cex = 1,col="steelblue")
text(15.7,3.06,expression(paste(UCL[p(3)]," = 2.83")),cex = 1,col="steelblue")
text(1.7,3.5,expression(paste(LCL[c(1)]," = 2.80")),cex = 1)
text(8.7,3.5,expression(paste(LCL[c(2)]," = 3.17")),cex = 1)
text(15.7,3.5,expression(paste(LCL[c(3)]," = 3.37")),cex = 1)
abline(h=-2.61,lty=1,col=6) # x
abline(h=-2.48,lty=1,col=6)
abline(h=-2.24,lty=1,col=6)
abline(h=-2.56,lty=1,col="steelblue") # p
abline(h=-2.44,lty=1,col="steelblue")
abline(h=-2.25,lty=1,col="steelblue")
abline(h=-2.73,lty=1,col=1) # p
abline(h=-2.60,lty=1,col=1)
abline(h=-2.35,lty=1,col=1)
text(1.7,-2.85,expression(paste(UCL[x(1)]," = -2.24")),cex = 1,col=6)
text(8.7,-2.85,expression(paste(UCL[x(2)]," = -2.48")),cex = 1,col=6)
text(15.7,-2.85,expression(paste(UCL[x(3)]," = -2.61")),cex = 1,col=6)
text(1.7,-3.20,expression(paste(UCL[p(1)]," = -2.25")),cex = 1,col="steelblue")
text(8.7,-3.20,expression(paste(UCL[p(2)]," = -2.44")),cex = 1,col="steelblue")
text(15.7,-3.20,expression(paste(UCL[p(3)]," = -2.56")),cex = 1,col="steelblue")
text(1.7,-2.10,expression(paste(LCL[c(1)]," = -2.35")),cex = 1)
text(8.7,-2.10,expression(paste(LCL[c(2)]," = -2.60")),cex = 1)
text(15.7,-2.10,expression(paste(LCL[c(3)]," = -2.73")),cex = 1)

rej_c <- which(B[,3] > 3.37)
points(rej_c,B[rej_c,3],pch=8,cex=1.2)
text(rej_c,B[rej_c,3]+0.2,rej_c,cex=1.1)


B_C <- rep(0,100)
##------- Plot (2) ---------------------- 
x <- c(1,seq(5,100,by=5))
y <- seq(-3.5,4,by=0.5)
plot(B_U[,1],type = "b",pch=16, xaxt = "n", yaxt = "n",
     xlim = c(1,100),ylim = c(-3.5,5.5),col="red",lty=1,lwd=1,
     xlab = "Sample Sequence",ylab = "value",cex.lab=1.8,cex.axis=1.3)
points(B_U[,2],col="blue",type = "b",pch=16,lty=1,lwd=1.3)
points(B_U[,3],type = "b",pch=16,lty=1,lwd=1.6)
points(B_L[,1],col="red",type = "b",pch=16,lty=1,lwd=1)
points(B_L[,2],col="blue",type = "b",pch=16,lty=1,lwd=1.3)
points(B_L[,3],type = "b",pch=16,lty=1,lwd=1.6)
lines(B_C,type = "l")
axis(1, labels = x, at = x,cex.axis=1.1)
axis(2, labels = y, at = y)
#abline(h=0,lty=1,col="gray",lwd=1.4)
legend(1.2,5.5,legend = c(expression(X[1]),expression(X[2]),expression(X[3])),
       col=c("red","blue",1),lty = c(1,1,1),lwd = c(1,1.3,1.6),cex = 0.8 )
abline(v=50)
text(43,5.15,"Phase I",cex = 3)
text(57,5.15,"Phase II",cex = 3)
abline(h=2.47,lty=1,col=6) # x
abline(h=2.76,lty=1,col=6)
abline(h=2.92,lty=1,col=6)
abline(h=2.42,lty=1,col="steelblue") # p
abline(h=2.68,lty=1,col="steelblue")
abline(h=2.83,lty=1,col="steelblue")
abline(h=2.80,lty=1,col=1) # p
abline(h=3.17,lty=1,col=1)
abline(h=3.37,lty=1,col=1)
text(1.7,2.24,expression(paste(UCL[x(1)]," = 2.47")),cex = 1,col=6)
text(8.7,2.24,expression(paste(UCL[x(2)]," = 2.76")),cex = 1,col=6)
text(15.7,2.24,expression(paste(UCL[x(3)]," = 2.92")),cex = 1,col=6)
text(1.7,3.06,expression(paste(UCL[p(1)]," = 2.42")),cex = 1,col="steelblue")
text(8.7,3.06,expression(paste(UCL[p(2)]," = 2.68")),cex = 1,col="steelblue")
text(15.7,3.06,expression(paste(UCL[p(3)]," = 2.83")),cex = 1,col="steelblue")
text(1.7,3.5,expression(paste(LCL[c(1)]," = 2.80")),cex = 1)
text(8.7,3.5,expression(paste(LCL[c(2)]," = 3.17")),cex = 1)
text(15.7,3.5,expression(paste(LCL[c(3)]," = 3.37")),cex = 1)
abline(h=-2.61,lty=1,col=6) # x
abline(h=-2.48,lty=1,col=6)
abline(h=-2.24,lty=1,col=6)
abline(h=-2.56,lty=1,col="steelblue") # p
abline(h=-2.44,lty=1,col="steelblue")
abline(h=-2.25,lty=1,col="steelblue")
abline(h=-2.73,lty=1,col=1) # p
abline(h=-2.60,lty=1,col=1)
abline(h=-2.35,lty=1,col=1)
text(1.7,-2.85,expression(paste(UCL[x(1)]," = -2.24")),cex = 1,col=6)
text(8.7,-2.85,expression(paste(UCL[x(2)]," = -2.48")),cex = 1,col=6)
text(15.7,-2.85,expression(paste(UCL[x(3)]," = -2.61")),cex = 1,col=6)
text(1.7,-3.20,expression(paste(UCL[p(1)]," = -2.25")),cex = 1,col="steelblue")
text(8.7,-3.20,expression(paste(UCL[p(2)]," = -2.44")),cex = 1,col="steelblue")
text(15.7,-3.20,expression(paste(UCL[p(3)]," = -2.56")),cex = 1,col="steelblue")
text(1.7,-2.10,expression(paste(LCL[c(1)]," = -2.35")),cex = 1)
text(8.7,-2.10,expression(paste(LCL[c(2)]," = -2.60")),cex = 1)
text(15.7,-2.10,expression(paste(LCL[c(3)]," = -2.73")),cex = 1)

rej_c <- which(B[,3] > 3.37)
points(rej_c,B[rej_c,3],pch=8,cex=1.2)
text(rej_c,B[rej_c,3]+0.2,rej_c,cex=1.1)



## Bon
##------- Plot ---------------------- 
x <- c(1,seq(5,100,by=5))
y <- seq(-3.5,4,by=0.5)
plot(B[,1],type = "b",pch=16, xaxt = "n", yaxt = "n",
     xlim = c(1,100),ylim = c(-3.5,5.5),col="red",lty=1,lwd=1,
     xlab = "Sample Sequence",ylab = "value",cex.lab=1.8,cex.axis=1.3)
points(B[,2],col="blue",type = "b",pch=16,lty=1,lwd=1.3)
points(B[,3],type = "b",pch=16,lty=1,lwd=1.6)
axis(1, labels = x, at = x,cex.axis=1.1)
axis(2, labels = y, at = y)
abline(h=0,lty=1,col="gray",lwd=1.4)
legend(1.2,5.5,legend = c(expression(X[1]),expression(X[2]),expression(X[3])),
       col=c("red","blue",1),lty = c(1,1,1),lwd = c(1,1.3,1.6),cex = 0.8 )
abline(v=50)
text(43,5.15,"Phase I",cex = 3)
text(57,5.15,"Phase II",cex = 3)
abline(h=2.92,lty=1,col=6)# x
abline(h=2.83,lty=1,col="steelblue")# p
abline(h=3.37,lty=1,col=1)# c
text(1.7,3.05,expression(paste(UCL[x(3)]," = 2.92")),cex = 1,col=6)
text(1.7,2.70,expression(paste(UCL[p(3)]," = 2.83")),cex = 1,col="steelblue")
text(1.7,3.50,expression(paste(LCL[c(3)]," = 3.37")),cex = 1)
abline(h=-2.61,lty=1,col=6) # x
abline(h=-2.56,lty=1,col="steelblue") # p
abline(h=-2.70,lty=1,col=1) # p
text(1.7,-2.85,expression(paste(UCL[x(3)]," = -2.61")),cex = 1,col=6)
text(1.7,-2.40,expression(paste(UCL[p(3)]," = -2.56")),cex = 1,col="steelblue")
text(1.7,-2.10,expression(paste(LCL[c(3)]," = -2.73")),cex = 1)

rej_c2 <- which(B[,3] > 3.37)
points(rej_c2,B[rej_c2,3],pch=8,cex=1.2)
text(rej_c2,B[rej_c2,3]+0.2,rej_c2,cex=1.1)


##------- Plot (2)---------------------- 
x <- c(1,seq(5,100,by=5))
y <- seq(-3.5,4,by=0.5)
plot(B_U[,1],type = "b",pch=16, xaxt = "n", yaxt = "n",
     xlim = c(1,100),ylim = c(-3.5,5.5),col="red",lty=1,lwd=1,
     xlab = "Sample Sequence",ylab = "value",cex.lab=1.8,cex.axis=1.3)
points(B_U[,2],col="blue",type = "b",pch=16,lty=1,lwd=1.3)
points(B_U[,3],type = "b",pch=16,lty=1,lwd=1.6)
points(B_L[,1],col="red",type = "b",pch=16,lty=1,lwd=1)
points(B_L[,2],col="blue",type = "b",pch=16,lty=1,lwd=1.3)
points(B_L[,3],type = "b",pch=16,lty=1,lwd=1.6)
lines(B_C,type = "l")
axis(1, labels = x, at = x,cex.axis=1.1)
axis(2, labels = y, at = y)
#abline(h=0,lty=1,col="gray",lwd=1.4)
legend(1.2,5.5,legend = c(expression(X[1]),expression(X[2]),expression(X[3])),
       col=c("red","blue",1),lty = c(1,1,1),lwd = c(1,1.3,1.6),cex = 0.8 )
abline(v=50)
text(43,5.15,"Phase I",cex = 3)
text(57,5.15,"Phase II",cex = 3)
abline(h=2.92,lty=1,col=6)# x
abline(h=2.83,lty=1,col="steelblue")# p
abline(h=3.37,lty=1,col=1)# c
text(1.7,3.05,expression(paste(UCL[x(3)]," = 2.92")),cex = 1,col=6)
text(1.7,2.70,expression(paste(UCL[p(3)]," = 2.83")),cex = 1,col="steelblue")
text(1.7,3.50,expression(paste(LCL[c(3)]," = 3.37")),cex = 1)
abline(h=-2.61,lty=1,col=6) # x
abline(h=-2.56,lty=1,col="steelblue") # p
abline(h=-2.70,lty=1,col=1) # p
text(1.7,-2.85,expression(paste(UCL[x(3)]," = -2.61")),cex = 1,col=6)
text(1.7,-2.40,expression(paste(UCL[p(3)]," = -2.56")),cex = 1,col="steelblue")
text(1.7,-2.10,expression(paste(LCL[c(3)]," = -2.73")),cex = 1)

rej_c2 <- which(B[,3] > 3.37)
points(rej_c2,B[rej_c2,3],pch=8,cex=1.2)
text(rej_c2,B[rej_c2,3]+0.2,rej_c2,cex=1.1)


