library(foreach)
library(doSNOW)
library(flexmix)
library(copula)
library(snipEM)
library(MultiLCIRT)
library(MASS)
library(compiler)

load("Yit.RData")
load("Xit.RData")
load("RI.RData")

n=dim(Yit)[1]
Ti=dim(Yit)[2]
d<-3
Cat<-4
nui<-list(NULL,NULL,NULL)
kmax=4

cl = makeCluster(5, type = "SOCK")
registerDoSNOW(cl)

for(k in 2:kmax){
  
  # Frank
  rf = foreach(rp=1:5)%dopar%{
    source("frankcopula.r")
    return(EMC_f(Yit,Xit,RI[[rp]][[k-1]],1e-2,k,TRUE,FALSE))
  }
  
  # Clayton
  rc = foreach(rp=1:5)%dopar%{
    source("claytoncopula.r")
    return(EMC_c(Yit,Xit,RI[[rp]][[k-1]],1e-2,k,TRUE,FALSE))
  }
  
  # Joe
  rj = foreach(rp=1:5)%dopar%{
    source("joecopula.r")
    return(EMC_j(Yit,Xit,RI[[rp]][[k-1]],1e-2,k,TRUE,FALSE))
  }
  
  # Gumbel
  rg = foreach(rp=1:5)%dopar%{
    source("gumbelcopula.r")
    return(EMC_g(Yit,Xit,RI[[rp]][[k-1]],1e-2,k,TRUE,FALSE))
  }
  
  MChoice[[k-1]]=list(rf=rf,rc=rc,rj=rj,rg=rg)
  
}

{
  k=5
  # Frank
  rf = foreach(rp=1:5)%dopar%{
    source("frankcopula.r")
    return(EMC_f(Yit,Xit,RI[[rp]][[k-1]],1e-2,k,TRUE,FALSE))
  }
  
  # Clayton
  rc = foreach(rp=1:5)%dopar%{
    source("claytoncopula.r")
    return(EMC_c(Yit,Xit,RI[[rp]][[k-1]],1e-2,k,TRUE,FALSE))
  }
  
  # Joe
  rj = foreach(rp=1:5)%dopar%{
    source("joecopula.r")
    return(EMC_j(Yit,Xit,RI[[rp]][[k-1]],1e-2,k,TRUE,FALSE))
  }
  
  # Gumbel
  rg = foreach(rp=1:5)%dopar%{
    source("gumbelcopula.r")
    return(EMC_g(Yit,Xit,RI[[rp]][[k-1]],1e-2,k,TRUE,FALSE))
  }
  
  MChoice5=list(rf=rf,rc=rc,rj=rj,rg=rg)
  
}

stopCluster(cl)

# save(MChoice5, file = "Mchoice5.RData")

load("Mchoice.RData")
load("Mchoice5.RData")
load("outk1.RData")

liks = matrix(NA,5,4)

liks[1,] = c(max(sapply(1:5, function(x) outk1$rf1[[x]]$LU)),
             max(sapply(1:5, function(x) outk1$rc1[[x]]$LU)),
             max(sapply(1:5, function(x) outk1$rj1[[x]]$LU)),
             max(sapply(1:5, function(x) outk1$rg1[[x]]$LU)))

for(k in 2:4){
  
  liks[k,] = c(max(sapply(1:5, function(x) MChoice[[k-1]]$rf[[x]]$LU$lik)),
               max(sapply(1:5, function(x) MChoice[[k-1]]$rc[[x]]$LU$lik)),
               max(sapply(1:5, function(x) MChoice[[k-1]]$rj[[x]]$LU$lik)),
               max(sapply(1:5, function(x) MChoice[[k-1]]$rg[[x]]$LU$lik)))
}

liks[5,] = c(max(sapply(1:5, function(X) MChoice5$rf[[X]]$LU$lik)),
             max(sapply(1:5, function(X) MChoice5$rc[[X]]$LU$lik)),
             max(sapply(1:5, function(X) MChoice5$rj[[X]]$LU$lik)),
             max(sapply(1:5, function(X) MChoice5$rg[[X]]$LU$lik)))

BICemc=function(Y,X,k,lik,cpar,nui){
  
  if(k==1){
    n=dim(Y)[1]
    d=dim(Y)[3]
    p=dim(X)[3]
    Cat = length(table(Y[,,3]))-1
    np = ((d-1)*k)+(Cat*k)+(p*d)+(k*(k-1))+k-1+cpar+sum(!is.na(nui))
  }
  if(k >1){
    n=dim(Y)[1]
    d=dim(Y)[3]
    p=dim(X)[4]
    Cat = length(table(Y[,,3]))-1
    np = ((d-1)*k)+(Cat*k)+(p*d)+(k*(k-1))+k-1+cpar+sum(!is.na(nui))
  }
  -2*lik+log(n)*np
}

BICS = matrix(NA,5,4)

for(k in 1:5){
  for(m in 1:4){
    BICS[k,m] = BICemc(Yit,Xit,k,liks[k,m],1,NA)
  }
}
