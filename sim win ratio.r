
library(sas7bdat)
library(dplyr)
satt=read.sas7bdat('satt.sas7bdat')


store=data.frame()
for (num in 4690:10000){
  
  dat=satt[satt[,1]== num,c(2,4:8)]
  #### Treatment Group
  trt = dat[dat[,1]==1,]
  ## remove the variable TRT
  trt = trt[, 2:6]
  names(trt) = c('evt1_trt',  'evt2_trt',  'evt3_trt', 'evt4_trt','pid_trt')
  
  #### Control Group
  con = dat[dat[,1]==2,]
  con = con[, 2:6]
  names(con) = c('evt1_con', 'evt2_con', 'evt3_con', 'evt4_con','pid_con')
  
  
  #### Join the two groups
  trtcon = merge(trt, con)
  #### Control patient won, per Endpoint 4
  a =  (trtcon$evt4_trt == 1 & trtcon$evt4_con < 1 ) 
  
  #### Treatment patient won, per Endpoint 4
  b =  (trtcon$evt4_trt < 1 & trtcon$evt4_con == 1 ) 
  
  #### Control patient won, per Endpoint 3
  c =  (trtcon$evt3_trt == 1 & trtcon$evt3_con < 1 ) 
  
  #### Treatment patient won, per Endpoint 3
  d =  (trtcon$evt3_trt < 1 & trtcon$evt3_con == 1 ) 
  
  
  #### Control patient won, per Endpoint 2
  e = ifelse(trtcon$evt2_trt  < trtcon$evt2_con -5,1,0) 
  
  
  #### Treatment patient won, per Endpoint 2
  f =  ifelse(trtcon$evt2_con  < trtcon$evt2_trt -5,1,0)
  
  
  #### Control patient won, per Endpoint 2
  g =  ifelse(trtcon$evt1_trt-5 > trtcon$evt1_con,1,0 )
  
  #### Treatment patient won, per Endpoint 2
  h =  ifelse(trtcon$evt1_con-5 > trtcon$evt1_trt,1,0)
  
  #### prioritize
  for (i in 1:length(a) ) {
    if (a[i] == 1 | b[i] == 1) { 
      c[i] = 0
      d[i] = 0
    }
    
    if (a[i] == 1 | b[i] == 1 | c[i] == 1 | d[i] == 1) { 
      e[i] = 0
      f[i] = 0
    }
    
    if (a[i] == 1 | b[i] == 1 | c[i] == 1 | d[i] == 1|  e[i]==1  |  f[i]==1) { 
      g[i] = 0
      h[i] = 0
    }
    
  }
  
  #a=as.data.frame(a)
  #a$stratum=1
  #b=as.data.frame(b)
  #b$stratum=1
  #c=as.data.frame(c)
  #c$stratum=1
  #d=as.data.frame(d)
  #d$stratum=1
  #e=as.data.frame(e)
  #e$stratum=1
  #f=as.data.frame(f)
  #f$stratum=1
  #g=as.data.frame(g)
  #g$stratum=1
  #h=as.data.frame(h)
  #h$stratum=1
  #############################################################################################
  #### Stratum-specific win ratios
  #############################################################################################
  
  na = sum(a)
  nb = sum(b)
  nc = sum(c)
  nd = sum(d)
  ne = sum(e)
  nf = sum(f)
  ng = sum(g)
  nh = sum(h)
  
  
  ### number of wins per stratum
  win_con = sum(na,nc,ne,ng)
  win_trt= sum(nb,nd,nf,nh)
  
  ### stratum-specific win ratio
  WR = win_trt/win_con
  
  ### number of patients per stratum
  N_trt = nrow(trt)
  N_con = nrow(con)
  
  N_trt_con = cbind( N_trt, N_con )
  colnames(N_trt_con)=c( 'N2_trt', 'N2_con')
  
  #### theta K0/L0
  theta_KL_0 = (win_trt + win_con)/(2*nrow(trtcon))
  
  #############################################################################################
  #### variances and covariances per stratum
  #############################################################################################
  
  ### Kernel function K and L
  K = b + d + f+ h
  L = a + c + e+ g
  
  KL = cbind( trtcon, K, L)
  
  sum_k_trt = aggregate(K, by=list( KL$pid_trt), FUN = sum)
  sum_k_con = aggregate(K, by=list( KL$pid_con), FUN = sum)
  
  sum_L_trt = aggregate(L, by=list( KL$pid_trt), FUN = sum)
  sum_L_con = aggregate(L, by=list( KL$pid_con), FUN = sum)
  
  
  
  names(sum_k_trt) = c( 'pid_trt', 'sum_k_trt')
  names(sum_k_con) = c( 'pid_con', 'sum_k_con')
  
  names(sum_L_trt) = c(  'pid_trt', 'sum_L_trt')
  names(sum_L_con) = c( 'pid_con', 'sum_L_con')
  
  KL = merge(KL, sum_k_trt, by=c(  'pid_trt'))
  KL = merge(KL, sum_k_con, by=c( 'pid_con'))
  KL = merge(KL, sum_L_trt, by=c(  'pid_trt'))
  KL = merge(KL, sum_L_con, by=c(  'pid_con'))
  
  
  
  KL = merge(KL, N_trt_con)
  KL$theta_KL_0=theta_KL_0
  
  
  sig2_trt1 = N_trt*N_con*sum((KL$K-theta_KL_0)*(KL$sum_k_trt - KL$K - (N_con - 1)*theta_KL_0 ))/(N_con-1)
  sig2_trt2 = N_con*N_trt*sum((KL$K-theta_KL_0)*(KL$sum_k_con - KL$K - (N_trt - 1)*theta_KL_0 ))/(N_trt-1)
  
  
  sig2_con1 = N_con*N_trt*sum((KL$L-theta_KL_0)*(KL$sum_L_con - KL$L - (N_trt - 1)*theta_KL_0 ))/(N_trt-1)
  sig2_con2 = N_trt*N_con*sum((KL$L-theta_KL_0)*(KL$sum_L_trt - KL$L - (N_con - 1)*theta_KL_0 ))/(N_con-1)
  
  sig_trt_con1 =  N_con*N_trt*sum((KL$K-theta_KL_0)*(KL$sum_L_trt - KL$L - (N_con - 1)*theta_KL_0 ))/(N_con-1)
  sig_trt_con2 = N_trt*N_con*sum((KL$K-theta_KL_0)*(KL$sum_L_con - KL$L - (N_trt - 1)*theta_KL_0 ))/(N_trt-1)
  
  sig2_trt = sig2_trt1/N_trt + sig2_trt2/N_con
  sig2_con = sig2_con1/N_con + sig2_con2/N_trt
  sig_trt_con = sig_trt_con1/N_trt + sig_trt_con2/N_con
  
  
  #############################################################################################
  #### Stratified win ratio
  #############################################################################################
  
  ## Total sample size per stratum
  N = N_trt + N_con
  
  ### Stratified win ratio wit Mantel-Haenszel-type wights
  WR = win_trt/win_con
  
  #############################################################################################
  #### Homogeneity test per Cochran (1954)
  #############################################################################################
  
  ### Variance estimate of log (win ratio) per stratum
  var = 4*(sig2_trt + sig2_con - 2*sig_trt_con)/((win_trt + win_con)*(win_trt + win_con))
  
  wr_z= abs(log(win_trt/win_con)/sqrt(var));
  
  wr_p=1-pnorm(wr_z)
  
  store[num,1]=wr_z
  store[num,2]=wr_p
}

