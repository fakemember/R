ss=read.csv('sss.csv')
ss=ss[,-c(1,2)]
ss
boxplot(salary~gender*diploma*school,data = ss,main='comparison of eight groups on salary',
        ylab='salary')


RMN=ss[(ss[,1]=='RU')&(ss[,2]=='M')&(ss[,3]=='N'),]
RMG=ss[(ss[,1]=='RU')&(ss[,2]=='M')&(ss[,3]=='G'),]
RFN=ss[(ss[,1]=='RU')&(ss[,2]=='F')&(ss[,3]=='N'),]
RFG=ss[(ss[,1]=='RU')&(ss[,2]=='F')&(ss[,3]=='G'),]
CMN=ss[(ss[,1]=='CN')&(ss[,2]=='M')&(ss[,3]=='N'),]
CMG=ss[(ss[,1]=='CN')&(ss[,2]=='M')&(ss[,3]=='G'),]
CFN=ss[(ss[,1]=='CN')&(ss[,2]=='F')&(ss[,3]=='N'),]
CFG=ss[(ss[,1]=='CN')&(ss[,2]=='F')&(ss[,3]=='G'),]



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


letval(RMNS)
letval(RMGS)
letval(RFNS)
letval(RFGS)
letval(CMNS)
letval(CMGS)
letval(CFNS)
letval(CFGS)

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


RMNS=RMN$salary;RMGS=RMG$salary;RFNS=RFN$salary;RFGS=RFG$salary;CMNS=CMN$salary;CMGS=CMG$salary;CFNS=CFN$salary;CFGS=CFG$salary
sresult=SprVsLevel(RMNS,RMGS,RFNS,RFGS,CMNS,CMGS,CFNS,CFGS)
while (sresult!=1) {
  RMNS=RMNS^sresult
  RMGS=RMGS^sresult
  RFNS=RFNS^sresult
  RFGS=RFGS^sresult
  CMNS=CMNS^sresult
  CMGS=CMGS^sresult
  CFNS=CFNS^sresult
  CFGS=CFGS^sresult
  sresult=SprVsLevel(RFGS,CFGS,RFNS,CFNS,RMNS,CMNS,CMGS,RMGS)}



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

sym(c(RMNS,RMGS,RFNS,RFGS,CMNS,CMGS,CFNS,CFGS),k=7)



