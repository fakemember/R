letval <- function(x, k = 4) {
  LV <- c("M", "F", "E", "D", "C", "B", "A", "Z", "Y", "X", "W")
  out <- array(NA, c(k, 6))
  lx <- rx <- sort(x)
  dimnames(out) <-
    list(LV[1:k],c("LOWER","UPPER","DEPTH","MID","SPREAD","TAIL"))
  for(i in 1:k) {
    out[i, 1:2] <- c(median(lx), median(rx))
    nn <- (length(lx) + 1)/2
    lx <- lx[1:nn]
    rx <- rev(rev(rx)[1:nn])
    out[i, 3] <- nn
  }
  out[, 4] <- (out[, 1] + out[, 2])/2
  out[, 5] <- out[, 2] - out[, 1]
  out[, 6] <- c(0, out[-1, 5]/2/qnorm(1 - 1/2^(2:k)))
  out
}
letval(stack.loss,7)

letval<-function(x,k=4){
  LV= c("M", "F", "E", "D", "C", "B", "A", "Z", "Y", "X", "W")
  out=array(NA, c(k, 6))
  lx<-rx<-sort(x)
  dimnames(out)=list(LV[1:k],c("LOWER","UPPER","DEPTH","MID","SPREAD","TAIL"))
  for(i in 1:k){
    out[i,1:2]=c(median(lx), median(rx))
    nn <- (length(lx) + 1)/2
    lx <- lx[1:nn]
    rx <- rev(rev(rx)[1:nn])
    out[i, 3] <- nn

  }
  out[, 4] <- (out[, 1] + out[, 2])/2
  out[, 5] <- out[, 2] - out[, 1]
  out[, 6] <- c(0, out[-1, 5]/2/qnorm(1 - 1/2^(2:k)))
  out
}



data(InsectSprays)
boxplot(count~spray,data=InsectSprays,col="lightgray")

boxplot(log(1+count)~spray,data=InsectSprays,col="lightgray")
boxplot((count)^0.25~spray,data=InsectSprays,col="lightgray")



letval2 <- function(x, k = 4)
{
  LV <- c("M", "F", "E", "D", "C", "B", "A", "Z", "Y", "X", "W")
  out <- array(NA, c(k, 6))
  lx <- rx <- sort(x)
  dimnames(out) <- list(LV[1:k],c("LOWER","UPPER","DEPTH","MID","SPREAD","TAIL"))
  for(i in 1:k) {
    out[i, 1:2] <- c(median(lx), median(rx))
    nn <- (length(lx) + 1)/2
    lx <- lx[1:nn]
    rx <- rev(rev(rx)[1:nn])
    out[i, 3] <- nn
  }
  out[, 4] <- (out[, 1] + out[, 2])/2
  out[, 5] <- out[, 2] - out[, 1]
  out[, 6] <- c(0, out[-1, 5]/2/qnorm(1 - 1/2^(2:k)))
  out
}
sym <- function(x, k = 7) {
  z <- letval2(x, k)
  M <- z[1,1]
  u <- ((z[-1,2]-M)^2 + (z[-1,1]-M)^2)/(4*M)
  v <- z[-1, 4] - M
  plot(u, v)
  b <- lsfit(u, v, int = F)$coef
  abline(0, b)
  hist(x^(1 - b))
  1 - b
}


par(mfrow=c(2,2))
data(rivers)
hist(rivers)
sym(rivers,4)
hist(log(rivers))




SprVsLevel = function(u1,u2,...) {
  x = list(u1,u2,...)
  z <- sapply(x, function(z){u=letval(z);c(u[1,1],u[2,5])})
  z <- t(z)
  b <- lsfit(log(z[, 1]), log(z[, 2]))$coef
  boxplot(x)
  plot(log(z))
  abline(b)
  lz <- sapply(x, function(u, uu) u^uu, uu=1-b[2])
  if( is.matrix(lz))
    boxplot(data.frame(lz)) else boxplot(lz)
  1 - b[2]
}
x1<-c(0.250,0.375,0.250,0.250,0.250,0.250,1.250,3.250,2.000,1.375,34.750,
      0.625,37.125,12.000,18.125,5.000,10.000,5.250,4.875,0.375,7.000,1.000,0.250,
      9.375,7.125,24.750)
y1<-c(85.000,68.875,97.000,83.500,65.625,52.500,51.000,77.950,85.125,66.625,
      87.625,94.750,58.875,49.650,72.250,48.625,69.250,82.750,52.125,
      89.500,65.125,80.375,47.875,86.125,79.375,37.875,76.000,49.250,86.125,88.000)
par(mfrow=c(1,1))

SprVsLevel = function(u1,...) {
  if(is.list(u1) ) x = u1 else x = list(u1,...)
  z <- sapply(x, function(z){u=letval2(z);c(u[1,1],u[2,5])})
  z <- t(z)
  b <- lsfit(log(z[, 1]), log(z[, 2]))$coef
  boxplot(x)
  plot(log(z))
  abline(b)
  lz <- sapply(x, function(u, uu) u^uu, uu=1-b[2])
  if( is.matrix(lz))
    boxplot(data.frame(lz)) else boxplot(lz)
  1 - b[2]
}


SprVsLevel(x1,y1)







salary=c(150
         ,25
         ,34
         ,23
         , 36
         , 21
         ,19
         ,143
         , 78
         , 210
         ,35
         , 49
         , 40
         ,   49
         ,  29
         ,30
         ,54
         , 37
         , 26
         ,  70
         ,    28
         ,15
         ,136
         ,   3
         ,   10
         ,  22
         ,21
         ,24
         ,27
         , 36
         ,6
         ,  6
         , 15
         ,16
         ,  7
         , 37
         ,66
         , 24
         ,24
         ,25)



##smooth spline
spl <- function(x,y,p=10){
  xx = unique(quantile(x,c(0:(p-1))/(p-1)))
  p = length(xx)
  x = c(seq(from=min(x),to=max(x),length=200),x)
  z = cbind(x,x^2)
  for(i in 1:(p-1)) z <- cbind(z, pmax(x-xx[i],0)^3)
  yy <- data.frame(y=y ,z[-(1:200),]) 
  ff= predict(xlm<-lm(y~.,data=yy),data.frame(z[(1:200),])) 
  plot(x[-(1:200)],y)
  lines(x[1:200],ff)
  list(x=x,y=predict(xlm))
}

x=rnorm(50)
y=rt(50,df=2)
xx = unique(quantile(x,c(0:(10-1))/(10-1)))
xx
x = c(seq(from=min(x),to=max(x),length=200),x)
spl(x,y)








letval2 <- function(x, k = 4)
{
  LV <- c("M", "F", "E", "D", "C", "B", "A", "Z", "Y", "X", "W")
  out <- array(NA, c(k, 6))
  lx <- rx <- sort(x)
  dimnames(out) <- list(LV[1:k],c("LOWER","UPPER","DEPTH","MID","SPREAD","TAIL"))
  for(i in 1:k) {
    out[i, 1:2] <- c(median(lx), median(rx))
    nn <- (length(lx) + 1)/2
    lx <- lx[1:nn]
    rx <- rev(rev(rx)[1:nn])
    out[i, 3] <- nn
  }
  out[, 4] <- (out[, 1] + out[, 2])/2
  out[, 5] <- out[, 2] - out[, 1]
  out[, 6] <- c(0, out[-1, 5]/2/qnorm(1 - 1/2^(2:k)))
  out
}
sym <- function(x, k = 7) {
  z <- letval2(x, k)
  M <- z[1,1]
  u <- ((z[-1,2]-M)^2 + (z[-1,1]-M)^2)/(4*M)
  v <- z[-1, 4] - M
  plot(u, v)
  b <- lsfit(u, v, int = F)$coef
  abline(0, b)
  hist(x^(1 - b))
  1 - b
}
data(rivers)
hist(rivers)
sym(rivers,4)
hist(log(rivers))




####emax model
m1 = nls(ch8~e0 +(emax * dose)/(dose+ed50), start= c(ed50=68, e0=2.4,emax=16),trace=T,data=edd)

m1 = nls(ch8~e0 +(emax * dose)/(dose+ed50), start= c(ed50=68, e0=2.4,emax=16),control=nls.control(maxiter=100),trace=TRUE,na.action=na.omit,data=edd)



##anova
aov



##robust regression
library(MASS)
stack.lm = lm(stack.loss~stack.x)
stack.m = rlm(stack.loss~stack.x,method="M")
stack.mm = rlm(stack.loss~stack.x,method="MM")
stack.lms = lmsreg(stack.x,stack.loss)
res= cbind(stack.lm$resid,stack.m$resid,stack.mm$resid,stack.lms$resid)
par(mfrow=c(3,2))
plot(res[,1:2],pch=20,col=2); abline(0,1)
plot(res[,c(1,3)],pch=20,col=2); abline(0,1)
plot(res[,c(1,4)],pch=20,col=2); abline(0,1)
plot(res[,2:3],pch=20,col=2); abline(0,1)
plot(res[,c(2,4)],pch=20,col=2); abline(0,1)
plot(res[,3:4],pch=20,col=2); abline(0,1)








RU=c(150
     , 25
     ,34
     ,23
     , 36
     ,21
     ,19
     ,143
     ,78
     , 210
     , 35
     ,49
     ,40
     , 49
     ,29
     ,30
     ,54
     ,37
     , 26
     , 70)



CU=c(28
     ,15
     , 136
     , 3
     ,10
     ,22
     ,21
     ,24
     ,27
     ,36
     ,6
     , 6
     , 15
     , 16
     , 7
     ,37
     , 66
     ,24
     , 24
     ,25)
SprVsLevel(RU,CU)



sresult=SprVsLevel(c1,c2)





c1=RU;c2=CU
sresult=SprVsLevel(c1,c2)
while (sresult!=1) {c1=c1^sresult
c2=c2^sresult
sresult=SprVsLevel(c1,c2) 
}
c1=RU^sresult
c2=CU^sresult



sym(c1,4)
sym(c2,4)
