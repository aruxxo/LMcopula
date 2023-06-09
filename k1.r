library(copula)
library(snipEM)
library(MultiLCIRT)
library(MASS)
library(flexmix)
library(compiler)
library(foreach)
library(doSNOW)

g=c(function(x) exp(x)/(1+exp(x)),
    function(x) exp(x), 
    function(x){ 
      pbs = inv_glob(x,type = "g")$p
      return(pbs)
    })
qdens=c(function(qu,par,nui) qbinom(qu,1,par),
        function(qu,par,nui) qpois(qu,par),
        function(qu,par,nui) (which(cumsum(par)>=qu)[1])-1)
pdens=c(function(x,par,nui) pbinom(x,1,par),function(x,par,nui) ppois(x,par),function(x,par,nui) ifelse(x>=0,cumsum(par)[x+1],0))


getInits1 = function(Y,X,d){
  n=dim(Y)[1]
  p=dim(X)[4]
  d=dim(Y)[3]
  initials=list()
  initials$al=rep(NA,(d+Cat-2))
  initials$be=matrix(NA,d,p)
  pY=as.vector(Y[,,1])
  pX=matrix(X[,,1,],n*Ti,p)
  indata=data.frame(pY,pX)
  m1=glm(cbind(pY,1-pY)~X1+X2+X3+X4,family = "binomial",data = indata)
  initials$al[1]=m1$coefficients[1]
  initials$be[1,]=m1$coefficients[2:(p+1)]
  pY=as.vector(Y[,,2])
  pX=matrix(X[,,2,],n*Ti,p)
  indata=data.frame(pY,pX)
  m2=glm(pY~X1+X2+X3+X4,family = "poisson",data = indata)
  initials$al[2]=m2$coefficients[1]
  initials$be[2,]=m2$coefficients[2:(p+1)]
  pY=as.vector(as.factor(Y[,,3]))
  pX=matrix(X[,,3,],n*Ti,p)
  indata=data.frame(pY,pX)
  m3 = polr(pY~1+X1+X2+X3+X4,data = indata)
  initials$al[-c(1:2)]=-m3$zeta
  initials$be[3,]=m3$coefficients
  return(initials)
}

Cmass <- function(Y,X,al,be,d,p,cop){
  
  n=dim(Y)[1]
  Ti=dim(Y)[2]
  d=dim(Y)[3]
  Cat = length(table(Y[,,d]))
  ms=ms_1=linpred<-array(NA,c(n,Ti,d))
  logits = array(NA,c(n,Ti,Cat-1))
  Px = array(NA,c(n,Ti,Cat))
  Dens = matrix(NA,n,Ti)
  
  for(ti in 1:Ti){
      for(h in 1:(d-1)){
        linpred[,ti,h]  = g[[h]](al[h]+X[,ti,h,]%*%be[h,])
      }
      logits[,ti,] = X[,ti,d,]%*%be[d,]
      logits[,ti,] = logits[,ti,]+t(matrix(rep(al[-c(1:(d-1))],n),Cat-1,n))
      Px[,ti,] =  t(matrix(apply(logits[,ti,],1,g[[3]]),Cat,n))  
    }
  
  
  for(ti in 1:Ti){
      for(h in 1:(d-1)){
        ms[,ti,h]=pdens[[h]](Y[,ti,h],linpred[,ti,h])
        ms_1[,ti,h]=pdens[[h]](Y[,ti,h]-1,linpred[,ti,h])
      }
      ms[,ti,d] = sapply(1:n, function(i) pdens[[d]](Y[i,ti,d],Px[i,ti,]))
      m_1 = as.numeric(sapply(1:n, function(i) pdens[[d]](Y[i,ti,d]-1,Px[i,ti,])))
      m_1[is.na(m_1)]=0
      ms_1[,ti,d] = m_1
    }
  
  
  for(ti in 1:Ti){
      ds = pCopula(ms[,ti,],cop)-pCopula(cbind(ms_1[,ti,1],ms[,ti,2],ms[,ti,3]),cop) - pCopula(cbind(ms[,ti,1],ms_1[,ti,2],ms[,ti,3]),cop)-pCopula(cbind(ms[,ti,1],ms[,ti,2],ms_1[,ti,3]),cop)+
        pCopula(cbind(ms_1[,ti,1],ms_1[,ti,2],ms[,ti,3]),cop)+pCopula(cbind(ms[,ti,1],ms_1[,ti,2],ms_1[,ti,3]),cop) + pCopula(cbind(ms_1[,ti,1],ms[,ti,2],ms_1[,ti,3]),cop)-pCopula(ms_1[,ti,],cop)
      ds[ds < .Machine$double.eps]=.Machine$double.eps
      Dens[,ti]=log(ds)
    }
  
  return(list(Dens=Dens,linpred=linpred,ms=ms,ms_1=ms_1))
}

tgt <- function(theta,Y,X,al,be,p,d,cop,r,ms,ms_1,linpred,logits,Px,Dens,Cat){
  res = 0 
  if(r==d){
    al[r:(r-1+Cat-1)] = theta[1:((Cat-1))]  
    be[r,] = theta[((Cat-1)+1):((Cat-1)+p)]

  al[-c(1:(d-1))]=sort(al[-c(1:(d-1))],decreasing=TRUE)

  }else{
    al[r]=theta[1]
    be[r,]=theta[2:(p+1)]
  }
  
  n=dim(Y)[1]
  Ti=dim(Y)[2]
  d=dim(Y)[3]
  
  for(ti in 1:Ti){
      for(h in 1:(d-1)){
        linpred[,ti,h]  = g[[h]](al[h]+X[,ti,h,]%*%be[h,])
      }
      logits[,ti,] = X[,ti,d,]%*%be[d,]
      logits[,ti,] = logits[,ti,]+t(matrix(rep(al[-c(1:(d-1))],n),Cat-1,n))
      Px[,ti,] =  t(matrix(apply(logits[,ti,],1,g[[3]]),Cat,n))  
    
  }
  
  
  for(ti in 1:Ti){
      for(h in 1:(d-1)){
        ms[,ti,h]=pdens[[h]](Y[,ti,h],linpred[,ti,h])
        ms_1[,ti,h]=pdens[[h]](Y[,ti,h]-1,linpred[,ti,h])
      }
      ms[,ti,d] = sapply(1:n, function(i) pdens[[d]](Y[i,ti,d],Px[i,ti,]))
      ms_1[,ti,d] = as.numeric(sapply(1:n, function(i) pdens[[d]](Y[i,ti,d]-1,Px[i,ti,])))
      ds = pCopula(ms[,ti,],cop)-pCopula(cbind(ms_1[,ti,1],ms[,ti,2],ms[,ti,3]),cop) - pCopula(cbind(ms[,ti,1],ms_1[,ti,2],ms[,ti,3]),cop)-pCopula(cbind(ms[,ti,1],ms[,ti,2],ms_1[,ti,3]),cop)+
        pCopula(cbind(ms_1[,ti,1],ms_1[,ti,2],ms[,ti,3]),cop)+pCopula(cbind(ms[,ti,1],ms_1[,ti,2],ms_1[,ti,3]),cop) + pCopula(cbind(ms_1[,ti,1],ms[,ti,2],ms_1[,ti,3]),cop)-pCopula(ms_1[,ti,],cop)
      ds[ds < .Machine$double.eps]=.Machine$double.eps
      res=res-sum(log(ds))
    }
  res
  
  
}

load("Yit.RData")
load("Xit.RData")
n=dim(Yit)[1]
d=dim(Yit)[length(dim(Yit))]
Cat=length(table(Yit[,,d]))
Ti=dim(Yit)[2]
nui=list(NULL,NULL,NULL)
inits1 = getInits1(Yit,Xit,d)
inits1$tcop = max(apply(Yit,2:3,var))

RI1 = list()
RI1[[1]]=inits1
for(rp in 2:5){
set.seed(rp)
al1 = inits1$al + runif(length(inits1$al),-0.75,0.75)
be1 = inits1$be + runif(prod(dim(inits1$be)),-0.5,0.5)
tcop1 = inits1$tcop + runif(1, -0.5,0.5)
RI1[[rp]]=list(al=al1,be=be1,tcop=tcop1)
}

cl = makeCluster(5,type = "SOCK")
registerDoSNOW(cl)


{
# Frank

tgtC=function(theta,ms,ms_1,Ti) {
  res=0
  copU=frankCopula(theta,dim = d)
  
  for(ti in 1:Ti){
      ds = pCopula(ms[,ti,],copU)-pCopula(cbind(ms_1[,ti,1],ms[,ti,2],ms[,ti,3]),copU) - pCopula(cbind(ms[,ti,1],ms_1[,ti,2],ms[,ti,3]),copU)-pCopula(cbind(ms[,ti,1],ms[,ti,2],ms_1[,ti,3]),copU)+
        pCopula(cbind(ms_1[,ti,1],ms_1[,ti,2],ms[,ti,3]),copU)+pCopula(cbind(ms[,ti,1],ms_1[,ti,2],ms_1[,ti,3]),copU) + pCopula(cbind(ms_1[,ti,1],ms[,ti,2],ms_1[,ti,3]),copU)-pCopula(ms_1[,ti,],copU)
      ds[ds < .Machine$double.eps]=.Machine$double.eps
      res=res-sum(log(ds))
    }
  res
}

Mcopula_f <- function(Y,X,inits,eps,up.theta,verbose){
  
  n = dim(Y)[1]
  Ti = dim(Y)[2]
  d = dim(Y)[3]
  p = dim(X)[4]
  al = inits$al
  be = inits$be
  tcop = 0
  if(up.theta){
    tcop=inits$tcop
  }
  cop = frankCopula(tcop,dim = d)
  cm = Cmass(Y,X,al,be,d,p,cop)
  LC = sum(cm$Dens)
  LU = LC
  lcurr = lupd <- LC
  diff=Inf
  it=0
  ms=ms_1=linpred<-array(NA,c(n,Ti,d))
  logits = array(NA,c(n,Ti,Cat-1))
  Px = array(NA,c(n,Ti,Cat))
  if(verbose){
    cat("------------|-------------|-------------|\n")
    cat("     it     |     llik    |     diff    |\n")
    cat("------------|-------------|-------------|\n")
  }
  
  while(diff > eps){  
    it = it+1
    Dens = cm$Dens
    lcurr <- lupd
    LC <- LU
    for(r in sample(d,d)){
      if(r==d){
        op = optim(c(al[r:(r-1+Cat-1)],be[r,]), function(x) tgt(x,Y,X,al,be,p,d,cop,r,ms,ms_1,linpred,logits,Px,Dens,Cat),method = "BFGS",control = list(maxit=1))  
        al[r:(r-1+Cat-1)] = op$par[1:((Cat-1))]
        al[r:(r-1+Cat-1)] = sort(al[r:(r-1+Cat-1)], decreasing = T)
        be[r,] = op$par[((Cat-1)+1):((Cat-1)+p)] 
      }else{
        op = optim(c(al[r],be[r,]), function(x) tgt(x,Y,X,al,be,p,d,cop,r,ms,ms_1,linpred,logits,Px,Dens,Cat),method = "BFGS",control = list(maxit=1))  
        al[r]=op$par[1]
        be[r,]=op$par[2:(p+1)]}
    }
    
    if(up.theta){
      upM = Cmass(Y,X,al,be,d,p,cop)
      uptcop = optimize(tgtC,c(0,15),upM$ms,upM$ms_1,Ti)
      tcop=uptcop$minimum
      cop = frankCopula(tcop,dim = d)
    }
      
    cm = Cmass(Y,X,al,be,d,p,cop)
    lupd = sum(cm$Dens)
    diff = lupd - lcurr
    cat(sprintf("%11g",c(it,lupd,diff)),sep=" | ","\n")
  }
  
  return(list(LU=lupd,tcop=tcop,al=al,be=be,nui=nui,cm=cm))
}

rf1<- foreach(rp=1:5,.packages = c("copula","MultiLCIRT"))%dopar%{
return(Mcopula_f(Yit,Xit,RI1[[rp]],1e-2,TRUE,FALSE))
}


rfb<- foreach(rp=1:49,.packages = c("copula","MultiLCIRT"))%dopar%{
w=sample(n,n,replace=TRUE)
return(Mcopula_f(Yit[w,,],Xit[w,,,],rf1,1e-2,TRUE,FALSE))
}


# Clayton

tgtC=function(theta,ms,ms_1,Ti) {
  res=0
  copU=claytonCopula(theta,dim = d)
  
  for(ti in 1:Ti){
    ds = pCopula(ms[,ti,],copU)-pCopula(cbind(ms_1[,ti,1],ms[,ti,2],ms[,ti,3]),copU) - pCopula(cbind(ms[,ti,1],ms_1[,ti,2],ms[,ti,3]),copU)-pCopula(cbind(ms[,ti,1],ms[,ti,2],ms_1[,ti,3]),copU)+
      pCopula(cbind(ms_1[,ti,1],ms_1[,ti,2],ms[,ti,3]),copU)+pCopula(cbind(ms[,ti,1],ms_1[,ti,2],ms_1[,ti,3]),copU) + pCopula(cbind(ms_1[,ti,1],ms[,ti,2],ms_1[,ti,3]),copU)-pCopula(ms_1[,ti,],copU)
    ds[ds < .Machine$double.eps]=.Machine$double.eps
    res=res-sum(log(ds))
  }
  res
}

Mcopula_c <- function(Y,X,inits,eps,up.theta,verbose){
  
  n = dim(Y)[1]
  Ti = dim(Y)[2]
  d = dim(Y)[3]
  p = dim(X)[4]
  al = inits$al
  be = inits$be
  tcop = 0
  if(up.theta){
    tcop=inits$tcop
  }
  cop = claytonCopula(tcop,dim = d)
  cm = Cmass(Y,X,al,be,d,p,cop)
  LC = sum(cm$Dens)
  LU = LC
  lcurr = lupd <- LC
  diff=Inf
  it=0
  ms=ms_1=linpred<-array(NA,c(n,Ti,d))
  logits = array(NA,c(n,Ti,Cat-1))
  Px = array(NA,c(n,Ti,Cat))
  if(verbose){
    cat("------------|-------------|-------------|\n")
    cat("     it     |     llik    |     diff    |\n")
    cat("------------|-------------|-------------|\n")
  }
  
  while(diff > eps){  
    it = it+1
    Dens = cm$Dens
    lcurr <- lupd
    LC <- LU
    for(r in sample(d,d)){
      if(r==d){
        op = optim(c(al[r:(r-1+Cat-1)],be[r,]), function(x) tgt(x,Y,X,al,be,p,d,cop,r,ms,ms_1,linpred,logits,Px,Dens,Cat),method = "BFGS",control = list(maxit=1))  
        al[r:(r-1+Cat-1)] = op$par[1:((Cat-1))]
        al[r:(r-1+Cat-1)] = sort(al[r:(r-1+Cat-1)], decreasing = T)
        be[r,] = op$par[((Cat-1)+1):((Cat-1)+p)] 
      }else{
        op = optim(c(al[r],be[r,]), function(x) tgt(x,Y,X,al,be,p,d,cop,r,ms,ms_1,linpred,logits,Px,Dens,Cat),method = "BFGS",control = list(maxit=1))  
        al[r]=op$par[1]
        be[r,]=op$par[2:(p+1)]}
    }
    
    if(up.theta){
      upM = Cmass(Y,X,al,be,d,p,cop)
      uptcop = optimize(tgtC,c(0,15),upM$ms,upM$ms_1,Ti)
      tcop=uptcop$minimum
      cop = claytonCopula(tcop,dim = d)
    }
    
    cm = Cmass(Y,X,al,be,d,p,cop)
    lupd = sum(cm$Dens)
    diff = lupd - lcurr
    cat(sprintf("%11g",c(it,lupd,diff)),sep=" | ","\n")
  }
  
  return(list(LU=lupd,tcop=tcop,al=al,be=be,nui=nui,cm=cm))
}

rc1<- foreach(rp=1:5,.packages = c("copula","MultiLCIRT"))%dopar%{
  return(Mcopula_c(Yit,Xit,RI1[[rp]],1e-2,TRUE,FALSE))
}


# Joe

tgtC=function(theta,ms,ms_1,Ti) {
  res=0
  copU=joeCopula(theta,dim = d)
  
  for(ti in 1:Ti){
    ds = pCopula(ms[,ti,],copU)-pCopula(cbind(ms_1[,ti,1],ms[,ti,2],ms[,ti,3]),copU) - pCopula(cbind(ms[,ti,1],ms_1[,ti,2],ms[,ti,3]),copU)-pCopula(cbind(ms[,ti,1],ms[,ti,2],ms_1[,ti,3]),copU)+
      pCopula(cbind(ms_1[,ti,1],ms_1[,ti,2],ms[,ti,3]),copU)+pCopula(cbind(ms[,ti,1],ms_1[,ti,2],ms_1[,ti,3]),copU) + pCopula(cbind(ms_1[,ti,1],ms[,ti,2],ms_1[,ti,3]),copU)-pCopula(ms_1[,ti,],copU)
    ds[ds < .Machine$double.eps]=.Machine$double.eps
    res=res-sum(log(ds))
  }
  res
}

Mcopula_j <- function(Y,X,inits,eps,up.theta,verbose){
  
  n = dim(Y)[1]
  Ti = dim(Y)[2]
  d = dim(Y)[3]
  p = dim(X)[4]
  al = inits$al
  be = inits$be
  tcop = 0
  if(up.theta){
    tcop=inits$tcop
  }
  cop = joeCopula(tcop,dim = d)
  cm = Cmass(Y,X,al,be,d,p,cop)
  LC = sum(cm$Dens)
  LU = LC
  lcurr = lupd <- LC
  diff=Inf
  it=0
  ms=ms_1=linpred<-array(NA,c(n,Ti,d))
  logits = array(NA,c(n,Ti,Cat-1))
  Px = array(NA,c(n,Ti,Cat))
  if(verbose){
    cat("------------|-------------|-------------|\n")
    cat("     it     |     llik    |     diff    |\n")
    cat("------------|-------------|-------------|\n")
  }
  
  while(diff > eps){  
    it = it+1
    Dens = cm$Dens
    lcurr <- lupd
    LC <- LU
    for(r in sample(d,d)){
      if(r==d){
        op = optim(c(al[r:(r-1+Cat-1)],be[r,]), function(x) tgt(x,Y,X,al,be,p,d,cop,r,ms,ms_1,linpred,logits,Px,Dens,Cat),method = "BFGS",control = list(maxit=1))  
        al[r:(r-1+Cat-1)] = op$par[1:((Cat-1))]
        al[r:(r-1+Cat-1)] = sort(al[r:(r-1+Cat-1)], decreasing = T)
        be[r,] = op$par[((Cat-1)+1):((Cat-1)+p)] 
      }else{
        op = optim(c(al[r],be[r,]), function(x) tgt(x,Y,X,al,be,p,d,cop,r,ms,ms_1,linpred,logits,Px,Dens,Cat),method = "BFGS",control = list(maxit=1))  
        al[r]=op$par[1]
        be[r,]=op$par[2:(p+1)]}
    }
    
    if(up.theta){
      upM = Cmass(Y,X,al,be,d,p,cop)
      uptcop = optimize(tgtC,c(1,60),upM$ms,upM$ms_1,Ti)
      tcop=uptcop$minimum
      cop = joeCopula(tcop,dim = d)
    }
    
    cm = Cmass(Y,X,al,be,d,p,cop)
    lupd = sum(cm$Dens)
    diff = lupd - lcurr
    cat(sprintf("%11g",c(it,lupd,diff)),sep=" | ","\n")
  }
  
  return(list(LU=lupd,tcop=tcop,al=al,be=be,nui=nui,cm=cm))
}

rj1<- foreach(rp=1:5,.packages = c("copula","MultiLCIRT"))%dopar%{
  return(Mcopula_j(Yit,Xit,RI1[[rp]],1e-2,TRUE,FALSE))
}

# Gumbel

tgtC=function(theta,ms,ms_1,Ti) {
  res=0
  copU=gumbelCopula(theta,dim = d)
  
  for(ti in 1:Ti){
    ds = pCopula(ms[,ti,],copU)-pCopula(cbind(ms_1[,ti,1],ms[,ti,2],ms[,ti,3]),copU) - pCopula(cbind(ms[,ti,1],ms_1[,ti,2],ms[,ti,3]),copU)-pCopula(cbind(ms[,ti,1],ms[,ti,2],ms_1[,ti,3]),copU)+
      pCopula(cbind(ms_1[,ti,1],ms_1[,ti,2],ms[,ti,3]),copU)+pCopula(cbind(ms[,ti,1],ms_1[,ti,2],ms_1[,ti,3]),copU) + pCopula(cbind(ms_1[,ti,1],ms[,ti,2],ms_1[,ti,3]),copU)-pCopula(ms_1[,ti,],copU)
    ds[ds < .Machine$double.eps]=.Machine$double.eps
    res=res-sum(log(ds))
  }
  res
}

Mcopula_g <- function(Y,X,inits,eps,up.theta,verbose){
  
  n = dim(Y)[1]
  Ti = dim(Y)[2]
  d = dim(Y)[3]
  p = dim(X)[4]
  al = inits$al
  be = inits$be
  tcop = 0
  if(up.theta){
    tcop=inits$tcop
  }
  cop = gumbelCopula(tcop,dim = d)
  cm = Cmass(Y,X,al,be,d,p,cop)
  LC = sum(cm$Dens)
  LU = LC
  lcurr = lupd <- LC
  diff=Inf
  it=0
  ms=ms_1=linpred<-array(NA,c(n,Ti,d))
  logits = array(NA,c(n,Ti,Cat-1))
  Px = array(NA,c(n,Ti,Cat))
  if(verbose){
    cat("------------|-------------|-------------|\n")
    cat("     it     |     llik    |     diff    |\n")
    cat("------------|-------------|-------------|\n")
  }
  
  while(diff > eps){  
    it = it+1
    Dens = cm$Dens
    lcurr <- lupd
    LC <- LU
    for(r in sample(d,d)){
      if(r==d){
        op = optim(c(al[r:(r-1+Cat-1)],be[r,]), function(x) tgt(x,Y,X,al,be,p,d,cop,r,ms,ms_1,linpred,logits,Px,Dens,Cat),method = "BFGS",control = list(maxit=1))  
        al[r:(r-1+Cat-1)] = op$par[1:((Cat-1))]
        al[r:(r-1+Cat-1)] = sort(al[r:(r-1+Cat-1)], decreasing = T)
        be[r,] = op$par[((Cat-1)+1):((Cat-1)+p)] 
      }else{
        op = optim(c(al[r],be[r,]), function(x) tgt(x,Y,X,al,be,p,d,cop,r,ms,ms_1,linpred,logits,Px,Dens,Cat),method = "BFGS",control = list(maxit=1))  
        al[r]=op$par[1]
        be[r,]=op$par[2:(p+1)]}
    }
    
    if(up.theta){
      upM = Cmass(Y,X,al,be,d,p,cop)
      uptcop = optimize(tgtC,c(1,60),upM$ms,upM$ms_1,Ti)
      tcop=uptcop$minimum
      cop = gumbelCopula(tcop,dim = d)
    }
    
    cm = Cmass(Y,X,al,be,d,p,cop)
    lupd = sum(cm$Dens)
    diff = lupd - lcurr
    cat(sprintf("%11g",c(it,lupd,diff)),sep=" | ","\n")
  }
  
  return(list(LU=lupd,tcop=tcop,al=al,be=be,nui=nui,cm=cm))
}

rg1<- foreach(rp=1:5,.packages = c("copula","MultiLCIRT"))%dopar%{
  return(Mcopula_g(Yit,Xit,RI1[[rp]],1e-2,TRUE,FALSE))
}
}

outk1=list(rf1=rf1,rc1=rc1,rj1=rj1,rg1=rg1)

save(outk1,file = "outk1.RData")

