

###################################################################
##### This is a collection of functions written for compute   #####
##### point and variance estimate in estimated likelihood     #####
##### method by Pepe 91, in a two-stage study with biased     #####
##### sampling, July 2006                                     #####
###################################################################


scoreMat <- function(X1,X,y1,x1,y,pzw,beta,Nyx){
   fit1  <- drop(exp(X1 %*% beta)/(1+ exp(X1 %*% beta)))
   pfit1 <- ifelse(y1==1,fit1,1-fit1)
   pyxz<- pfit1*pzw
   pyx <- rep(0,length(y1))
   pyx[x1==0 & y1==0] <- sum(pyxz[x1==0 & y1==0])
   pyx[x1==1 & y1==0] <- sum(pyxz[x1==1 & y1==0])
   pyx[x1==0 & y1==1] <- sum(pyxz[x1==0 & y1==1])
   pyx[x1==1 & y1==1] <- sum(pyxz[x1==1 & y1==1])  
   S  <- (y1-fit1)*X1
   B  <- S*pfit1*pzw
   V  <- pfit1*(1-pfit1)
   der.B <- t(S[y1==1 &x1==1,])%*% ((V*X1*pzw/pyx)[y1==1&x1==1,]) - t((X1*V*pfit1*pzw/pyx)[y1==1&x1==1,]) %*% (X1[y1==1&x1==1,])
   der.A <- apply((X1*V*pzw)[y1==1 & x1==1,], 2, sum)
   B1    <- apply(B[y1==1 & x1==1,],2,sum)
   inf11  <- der.B -outer(der.A,B1/((unique(pyx[y1==1 & x1==1]))^2))

   der.B <- t(S[y1==1 &x1==0,])%*% ((V*X1*pzw/pyx)[y1==1&x1==0,]) - t((X1*V*pfit1*pzw/pyx)[y1==1&x1==0,]) %*% (X1[y1==1&x1==0,])
   der.A <- apply((X1*V*pzw)[y1==1 & x1==0,], 2, sum)
   B1    <- apply(B[y1==1 & x1==0,],2,sum)
   inf10  <- der.B -outer(der.A,B1/((unique(pyx[y1==1 & x1==0]))^2))

   der.B <- (-t(S[y1==0 & x1==1,])%*% ((V*X1*pzw/pyx)[y1==0 & x1==1,]) - t((X1*V*pfit1*pzw/pyx)[y1==0 & x1==1,]) %*% (X1[y1==0 & x1==1,])) 
   der.A <- (-apply((X1*V*pzw)[y1==0 & x1==1,], 2, sum))
   B0    <- apply(B[y1==0 & x1==1,],2,sum)
   inf01 <- der.B -outer(der.A,B0/((unique(pyx[y1==0 & x1==1]))^2))

   der.B <- (-t(S[y1==0 & x1==0,])%*% ((V*X1*pzw/pyx)[y1==0 & x1==0,]) - t((X1*V*pfit1*pzw/pyx)[y1==0 & x1==0,]) %*% (X1[y1==0 & x1==0,])) 
   der.A <- (-apply((X1*V*pzw)[y1==0 & x1==0,], 2, sum))
   B0    <- apply(B[y1==0 & x1==0,],2,sum)
   inf00 <- der.B -outer(der.A,B0/((unique(pyx[y1==0 & x1==0]))^2))

   fit  <- drop(exp(X %*% beta)/(1+ exp(X %*% beta)))
   v    <- fit*(1-fit)
   vinf <- t(X)%*% (v*X)
   cheese <- (-inf11*unique(Nyx[y1==1 & x1==1]) - inf10*unique(Nyx[y1==1 & x1==0]) -inf01*unique(Nyx[y1==0 & x1==1]) - inf00*unique(Nyx[y1==0 & x1==0])  + vinf)
   
   score  <- apply(rbind((y-fit)*X,Nyx*B/pyx),2,sum)   
   return(list(cheese=cheese,score=score))
}


##### The extra variance drawn from estimated density when using independence assumption #####


VarDensity <- function(beta,X1,y1,x1,pzw,Nyx,z.from.y,n1,n0,n, N1, N0, N){
   fit1  <- drop(exp(X1 %*% beta)/(1+ exp(X1 %*% beta)))
   pfit1 <- ifelse(y1==1,fit1,1-fit1)
   D1 <- ifelse(y1==1,pfit1*(1-pfit1),-pfit1*(1-pfit1))*X1
   DZW1 <- D1*pzw
   D2 <- matrix(0,length(y1),ncol(X1))
   D2[x1==0 & y1==0,] <- matrix(rep(apply(DZW1[x1==0 & y1==0,],2,sum),n),nrow=n,byrow=T)
   D2[x1==1 & y1==0,] <- matrix(rep(apply(DZW1[x1==1 & y1==0,],2,sum),n),nrow=n,byrow=T)
   D2[x1==0 & y1==1,] <- matrix(rep(apply(DZW1[x1==0 & y1==1,],2,sum),n),nrow=n,byrow=T)
   D2[x1==1 & y1==1,] <- matrix(rep(apply(DZW1[x1==1 & y1==1,],2,sum),n),nrow=n,byrow=T)
   pyx <- rep(0,length(y1))
   pyxz<- pfit1*pzw
   pyx[x1==0 & y1==0] <- sum(pyxz[x1==0 & y1==0])
   pyx[x1==1 & y1==0] <- sum(pyxz[x1==1 & y1==0])
   pyx[x1==0 & y1==1] <- sum(pyxz[x1==0 & y1==1])
   pyx[x1==1 & y1==1] <- sum(pyxz[x1==1 & y1==1])
   weight <- ifelse(z.from.y==1,1/n1,1/n0)   
   W <- (D1/pyx - D2/(pyx)^2 *pfit1)*Nyx * weight
   Wbar <- matrix(0,n,ncol(X1))
   num <- rep(1:n,4)
   for (i in 1:n) Wbar[i,] <- apply(W[num==i,],2,sum)         
   from.y <- z.from.y[1:n]
   I <- var(Wbar[from.y==1,])*(N1/N)^2*n1 + var(Wbar[from.y==0,])*(N0/N)^2*n0
   I
}

## Ting (2011-03-22): VarDensity2 is not used anywhere at all
####### The extra variance drawn from estimated density when not using independence assumption #####
##
##VarDensity2 <- function(beta,X1,y1,x1,pzx,Nyx,from.y,from.x,n1,n0,n){
##   fit1  <- drop(exp(X1 %*% beta)/(1+ exp(X1 %*% beta)))
##   pfit1 <- ifelse(y1==1,fit1,1-fit1)
##   D1 <- ifelse(y1==1,pfit1*(1-pfit1),-pfit1*(1-pfit1))*X1
##   DZW1 <- D1*pzx
##   D2 <- matrix(0,length(y1),4)
##   D2[x1==0 & y1==0,] <- matrix(rep(apply(DZW1[x1==0 & y1==0,],2,sum),nz0),nrow=nz0,byrow=T)
##   D2[x1==1 & y1==0,] <- matrix(rep(apply(DZW1[x1==1 & y1==0,],2,sum),nz1),nrow=nz1,byrow=T)
##   D2[x1==0 & y1==1,] <- matrix(rep(apply(DZW1[x1==0 & y1==1,],2,sum),nz0),nrow=nz0,byrow=T)
##   D2[x1==1 & y1==1,] <- matrix(rep(apply(DZW1[x1==1 & y1==1,],2,sum),nz1),nrow=nz1,byrow=T)
##   pyx <- rep(0,length(y1))
##   pyxz<- pfit1*pzx
##   pyx[x1==0 & y1==0] <- sum(pyxz[x1==0 & y1==0])
##   pyx[x1==1 & y1==0] <- sum(pyxz[x1==1 & y1==0])
##   pyx[x1==0 & y1==1] <- sum(pyxz[x1==0 & y1==1])
##   pyx[x1==1 & y1==1] <- sum(pyxz[x1==1 & y1==1])
##   weight <- rep(0,length(y1))
##   for (i in 1:length(y1)) weight[i] <- 1/sum(x==from.x[i] & y==from.y[i])
##   W <- (D1/pyx - D2/(pyx)^2 *pfit1)*Nyx * weight
##   Wbar <- matrix(0,n,4)
##   num <- rep(1:n,2)
##   for (i in 1:n) Wbar[i,] <- apply(W[num==i,],2,sum)         
##   from.y <- c(y[x==1],y[x==0])
##   from.x <- c(x[x==1],x[x==0])   
##   pyx[x1==0 & y1==0] <- NXY1/(NXY1+NXY3)  
##   pyx[x1==1 & y1==0] <- NXY2/(NXY2+NXY4)
##   pyx[x1==0 & y1==1] <- NXY3/(NXY1+NXY3)
##   pyx[x1==1 & y1==1] <- NXY4/(NXY2+NXY4)
##   I1 <- var(Wbar[from.y==1 & from.x==1,])*(unique(pyx[y1==1&x1==1]))^2*n11 + var(Wbar[from.y==1 &from.x==0,])*(unique(pyx[y1==1&x1==0]))^2*n10
##   I2 <- var(Wbar[from.y==0 & from.x==1,])*(unique(pyx[y1==0&x1==1]))^2*n01 + var(Wbar[from.y==0 &from.x==0,])*(unique(pyx[y1==0&x1==0]))^2*n00                                                                             
##   I <- I1 + I2
##   I
##}   




#### EM-like algorithm ####

## Ting (2011-03-22): estlik.EML is not used anywhere at all
##estlik.EML <- function(x,z,y,NXY1,NXY2,NXY3,NXY4,N1,N0) {
##   n1 <- sum(y==1)
##   n0 <- sum(y==0)
##   z.from.y <- ifelse(y==1,1,0)
##   z1 <- rep(z,4)
##   z.from.y<- rep(z.from.y,4)
##   n <- length(y)
##   x1 <- c(rep(1,n),rep(0,n),rep(1,n),rep(0,n))
##   y1 <- c(rep(1,2*n),rep(0,2*n))
##   X1 <- cbind(1,x1,z1,x1*z1)
##   X  <- cbind(1,x,z,x*z) 
##   pzw <- ifelse(z.from.y==0, 1/n0 *(N0/N), 1/n1 *(N1/N))
##   wgt0 <- ifelse(y==1,N1/n1,N0/n0)
##   sfit <- glm(y~x+z + x*z, family=binomial,weights=wgt0)
##   betastart <- sfit$coef   
##   fit1 <- drop(exp(X1 %*% betastart)/(1+ exp(X1 %*% betastart)))
##   pfit1 <- ifelse(y1==1,fit1,1-fit1)
##   Nyx  <- rep(0,length(y1))
##   pyx  <- rep(0,length(y1))
##   for (i in 1:length(y1)) {
##
##    if (x1[i]==0 & y1[i]==0) {
##                     Nyx[i] <- NXY1 - sum(x==0 & y==0)
##                     pyx[i] <- NXY1/(NXY1+NXY3)
##    }
##    if (x1[i]==1 & y1[i]==0) { 
##                     Nyx[i] <- NXY2 - sum(x==1 & y==0)
##                     pyx[i] <- NXY2/(NXY2+NXY4)
##
##    } 
##    if (x1[i]==0 & y1[i]==1) {
##                     Nyx[i] <- NXY3 - sum(x==0 & y==1)
##                     pyx[i] <- NXY3/(NXY1+NXY3)
##
##    }
##    if (x1[i]==1 & y1[i]==1) {
##                     Nyx[i] <- NXY4 - sum(x==1 & y==1)
##                     pyx[i] <- NXY4/(NXY2 + NXY4)
##    } 
##   }
##   
##   diff <- 1
##   tol <- 1e-7
##   it  <- 0
##   start <- betastart
##   beta <- start
##   xx <- c(x,x1)
##   zz <- c(z,z1)
##   XX <- cbind(1,xx,zz,xx*zz)
##   yy <- c(y,y1)
##   while ( diff>=tol & it < 500 ) {
##      it <- it +1
##      wgt1  <- Nyx*pfit1*pzw/pyx
##      wgt   <- c(rep(1,n),wgt1)
##      out  <- iwls(XX,yy,wgt,start,family=binomial(),maxit=300,epsilon=1e-8)
##      pfit1 <- ifelse(yy==1,out$fitted,1-out$fitted)[-(1:n)]
##      betanew <- out$beta  
##      diff <- max(abs(beta-betanew))
##      beta <- betanew
##      pyxz<- pfit1*pzw
##      pyx[x1==0 & y1==0] <- sum(pyxz[x1==0 & y1==0])
##      pyx[x1==1 & y1==0] <- sum(pyxz[x1==1 & y1==0])
##      pyx[x1==0 & y1==1] <- sum(pyxz[x1==0 & y1==1])
##      pyx[x1==1 & y1==1] <- sum(pyxz[x1==1 & y1==1])              
##      if (it==500)  beta <- rep(999,4)
##   }
##   out <- scoreMat(X1,X,y1,x1,y,pzw,beta,Nyx)
##   cheese <- out$cheese   
##   I3 <- VarDensity(beta,X1,y1,x1,pzw,Nyx,z.from.y,n1,n0,n)   
##   CovMat <- solve(cheese) %*% (cheese+I3) %*% solve(cheese)
##   return(list(beta=beta,CovMat=CovMat))
##}


###### EM-like algorithm for estimated likelihood without using independence #####

## Ting (2011-03-22): estlik.EMLnoind is not used anywhere at all
##estlik.EMLnoind <- function(x,z,y,NXY1,NXY2,NXY3,NXY4,N1,N0) {
##   n1 <- sum(y==1)
##   n0 <- sum(y==0)
##   n <- length(y)
##   wgt0 <- ifelse(y==1,N1/n1,N0/n0)
##   sfit <- glm(y~x+z + x*z, family=binomial,weights=wgt0)
##   betastart <- sfit$coef
##   n11 <- sum(y==1 & x==1)
##   n10 <- sum(y==1 & x==0)
##   n01 <- sum(y==0 & x==1)
##   n00 <- sum(y==0 & x==0)  
##   z0  <- z[x==0]
##   z1  <- z[x==1]      
##   nz0 <- length(z0)
##   nz1 <- length(z1)
##   z1  <- c(z1,z0,z1,z0)
##   from.y <- c(y[x==1],y[x==0],y[x==1],y[x==0])
##   z.from.y<-from.y
##   from.x <- c(x[x==1],x[x==0],x[x==1],x[x==0])
##   x1 <- c(rep(1,nz1),rep(0,nz0),rep(1,nz1),rep(0,nz0))
##   y1 <- c(rep(1,n),rep(0,n))
##   X1 <- cbind(1,x1,z1,x1*z1)
##   X <- cbind(1,x,z,x*z)
##   pzx <- rep(0,2*n)
##   Nyx <- pzx
##   pyx <- pzx
##   Nyx[x1==0 & y1==0] <- NXY1 - sum(x==0 & y==0)
##   pyx[x1==0 & y1==0] <- NXY1/(NXY1+NXY3)  
##   Nyx[x1==1 & y1==0] <- NXY2 - sum(x==1 & y==0)
##   pyx[x1==1 & y1==0] <- NXY2/(NXY2+NXY4)
##   Nyx[x1==0 & y1==1] <- NXY3 - sum(x==0 & y==1)
##   pyx[x1==0 & y1==1] <- NXY3/(NXY1+NXY3)
##   Nyx[x1==1 & y1==1] <- NXY4 - sum(x==1 & y==1)
##   pyx[x1==1 & y1==1] <- NXY4/(NXY2 + NXY4)
##   weight <- rep(0,length(y1))
##   for (i in 1:length(y1))  weight[i] <- 1/sum(x==from.x[i] & y==from.y[i])    
##   pzx[from.y==1 & from.x==1] <- weight[from.y==1 & from.x==1] * unique(pyx[y1==1 & x1==1])
##   pzx[from.y==1 & from.x==0] <- weight[from.y==1 & from.x==0] * unique(pyx[y1==1 & x1==0])
##   pzx[from.y==0 & from.x==1] <- weight[from.y==0 & from.x==1] * unique(pyx[y1==0 & x1==1])
##   pzx[from.y==0 & from.x==0] <- weight[from.y==0 & from.x==0] * unique(pyx[y1==0 & x1==0])
##
##   fit1 <- drop(exp(X1 %*% betastart)/(1+ exp(X1 %*% betastart)))
##   pfit1 <- ifelse(y1==1,fit1,1-fit1)
##      
##   diff <- 1
##   tol <- 1e-7
##   it  <- 0
##   start <- betastart
##   beta <- start
##   xx <- c(x,x1)
##   zz <- c(z,z1)
##   XX <- cbind(1,xx,zz,xx*zz)
##   yy <- c(y,y1)
##   while ( diff>=tol & it < 500 ) {
##      it <- it +1
##      wgt1  <- Nyx*pfit1*pzx/pyx
##      wgt   <- c(rep(1,n),wgt1)
##      out  <- iwls(XX,yy,wgt,start,family=binomial(),maxit=300,epsilon=1e-8)
##      pfit1 <- ifelse(yy==1,out$fitted,1-out$fitted)[-(1:n)]
##      betanew <- out$beta  
##      diff <- max(abs(beta-betanew))
##      beta <- betanew
##      pyxz<- pfit1*pzx
##      pyx[x1==0 & y1==0] <- sum(pyxz[x1==0 & y1==0])
##      pyx[x1==1 & y1==0] <- sum(pyxz[x1==1 & y1==0])
##      pyx[x1==0 & y1==1] <- sum(pyxz[x1==0 & y1==1])
##      pyx[x1==1 & y1==1] <- sum(pyxz[x1==1 & y1==1])              
##      if (it==500)  beta <- rep(999,4)
##   }
##   out <- scoreMat(X1,X,y1,x1,y,pzx,beta,Nyx)
##   cheese <- out$cheese   
##
##   fit1  <- drop(exp(X1 %*% beta)/(1+ exp(X1 %*% beta)))
##   pfit1 <- ifelse(y1==1,fit1,1-fit1)
##   D1 <- ifelse(y1==1,pfit1*(1-pfit1),-pfit1*(1-pfit1))*X1
##   DZW1 <- D1*pzx
##   D2 <- matrix(0,length(y1),4)
##   D2[x1==0 & y1==0,] <- matrix(rep(apply(DZW1[x1==0 & y1==0,],2,sum),nz0),nrow=nz0,byrow=T)
##   D2[x1==1 & y1==0,] <- matrix(rep(apply(DZW1[x1==1 & y1==0,],2,sum),nz1),nrow=nz1,byrow=T)
##   D2[x1==0 & y1==1,] <- matrix(rep(apply(DZW1[x1==0 & y1==1,],2,sum),nz0),nrow=nz0,byrow=T)
##   D2[x1==1 & y1==1,] <- matrix(rep(apply(DZW1[x1==1 & y1==1,],2,sum),nz1),nrow=nz1,byrow=T)
##   pyx <- rep(0,length(y1))
##   pyxz<- pfit1*pzx
##   pyx[x1==0 & y1==0] <- sum(pyxz[x1==0 & y1==0])
##   pyx[x1==1 & y1==0] <- sum(pyxz[x1==1 & y1==0])
##   pyx[x1==0 & y1==1] <- sum(pyxz[x1==0 & y1==1])
##   pyx[x1==1 & y1==1] <- sum(pyxz[x1==1 & y1==1])
##
##   W <- (D1/pyx - D2/(pyx)^2 *pfit1)*Nyx * weight
##
##   Wbar <- matrix(0,n,4)
##   num <- rep(1:n,2)
##   for (i in 1:n) Wbar[i,] <- apply(W[num==i,],2,sum)         
##   from.y <- c(y[x==1],y[x==0])
##   from.x <- c(x[x==1],x[x==0])   
##
##   m1 <- apply( Wbar[from.y==1 & from.x==1,],2,mean)*NXY4/(NXY2+NXY4) + apply( Wbar[from.y==1 & from.x==0,],2,mean)*NXY3/(NXY1+NXY3)
##   nxy <- sum(from.y==1 & from.x==1)
##   d1 <- NXY4/(NXY2+NXY4)*(Wbar[from.y==1 & from.x==1,] - matrix(m1,nrow=nxy,ncol=4,byrow=T))
##   nxy <- sum(from.y==1 & from.x==0)
##   d2 <- NXY3/(NXY1+NXY3)*(Wbar[from.y==1 & from.x==0,] - matrix(m1,nrow=nxy,ncol=4,byrow=T))
##   I1 <- t(rbind(d1,d2)) %*% (rbind(d1,d2))
##
##   m1 <- apply( Wbar[from.y==0 & from.x==1,],2,mean)*NXY2/(NXY2+NXY4) + apply( Wbar[from.y==0 & from.x==0,],2,mean)*NXY1/(NXY1+NXY3)
##   nxy <- sum(from.y==0 & from.x==1)
##   d3 <- NXY2/(NXY2+NXY4)*(Wbar[from.y==0 & from.x==1,] - matrix(m1,nrow=nxy,ncol=4,byrow=T))
##   nxy <- sum(from.y==0 & from.x==0)
##   d4 <- NXY1/(NXY1+NXY3)*(Wbar[from.y==0 & from.x==0,] - matrix(m1,nrow=nxy,ncol=4,byrow=T))
##   I2 <- t(rbind(d3,d4)) %*% (rbind(d3,d4))
##   I3 <- I1 + I2
##   
##   CovMat <- solve(cheese) %*% (cheese+I3) %*% solve(cheese)
##   return(list(beta=beta,CovMat=CovMat))
##}
