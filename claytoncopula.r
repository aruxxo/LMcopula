library(copula)
library(snipEM)
library(MultiLCIRT)
library(MASS)
library(flexmix)
library(compiler)
Cat<-4
d<-3
nui<-list(NULL,NULL,NULL)

g=c(function(x) exp(x)/(1+exp(x)),
    function(x) exp(x), 
    function(x){ 
      pbs = inv_glob(x,type = "g")$p
      return(pbs)
    })
qdens=c(function(qu,par,nui) qbinom(qu,1,par),
        function(qu,par,nui) qpois(qu,par),
        function(qu,par,nui) (which(cumsum(par)>=qu)[1])-1)
pdens=c(function(x,par,nui) pbinom(x,1,par),function(x,par,nui) ppois(x,par),function(x,par,nui) cumsum(par)[x+1])


complik=function(Y,X,al,be,pis,PI,k,Ti,Dens) {
  
  n = dim(Y)[1]
  d = dim(Y)[3]
  fwd=bwd<-array(0,c(n,Ti,k))
  
  for(j in 1:k) {fwd[,1,j]=Dens[,1,j]+log(pis[j])}
  fwd[,1,][!is.finite(fwd[,1,])]=-600
  
  for(ti in 2:Ti) {
    for(j in 1:k) {
      fwd[,ti,j] <- Dens[,ti,j]+apply(sweep(fwd[,ti-1,],2,log(PI[,j]),"+"),1,sumlog)}
    fwd[,ti,][!is.finite(fwd[,ti,])]=-600}
  
  for(ti in (Ti-1):1){
    for(j in 1:k){
      bwd[,ti,j] <- apply(sweep(bwd[,ti+1,]+Dens[,ti+1,],2,log(PI[j,]),"+"),1,sumlog)}
    bwd[,ti,][!is.finite(bwd[,ti,])]=-600
  }
  
  lik = sum(apply(fwd[,Ti,],1,sumlog))
  return(list(lik=lik,fwd=fwd,bwd=bwd))
}

tgt <- function(theta,Y,X,W,al,be,k,p,d,cop,r,ms,ms_1,linpred,logits,Px,Dens,Cat){
  res = 0 
  if(r==d){
    al[r:(r-1+Cat-1),] = theta[1:((Cat-1)*k)]  
    be[r,] = theta[((Cat-1)*k+1):((Cat-1)*k+p)]
    for(iter in 1:k){
      al[-c(1:(d-1)),iter]=sort(al[-c(1:(d-1)),iter],decreasing=TRUE)
    }
  }else{
    al[r,]=theta[1:k]
    be[r,]=theta[(1:p)+k]
  }
  
  n=dim(Y)[1]
  Ti=dim(Y)[2]
  d=dim(Y)[3]
  
  for(ti in 1:Ti){
    for(j in 1:k){  
      for(h in 1:(d-1)){
        linpred[,ti,h,j]  = g[[h]](al[h,j]+X[,ti,h,]%*%be[h,])
      }
      logits[,ti,,j] = X[,ti,d,]%*%be[d,]
    }
    for(j in 1:k){
      logits[,ti,,j] = logits[,ti,,j]+t(matrix(rep(al[-c(1:(d-1)),j],n),Cat-1,n))
      Px[,ti,,j] =  t(matrix(apply(logits[,ti,,j],1,g[[3]]),Cat,n))  
    }
  }
  
  
  for(ti in 1:Ti){
    for(j in 1:k){
      for(h in 1:(d-1)){
        ms[,ti,h,j]=pdens[[h]](Y[,ti,h],linpred[,ti,h,j])
        ms_1[,ti,h,j]=pdens[[h]](Y[,ti,h]-1,linpred[,ti,h,j])
      }
      ms[,ti,d,j] = sapply(1:n, function(i) pdens[[d]](Y[i,ti,d],Px[i,ti,,j]))
      m_1 = as.numeric(sapply(1:n, function(i) pdens[[d]](Y[i,ti,d]-1,Px[i,ti,,j])))
      m_1[is.na(m_1)]=0
      ms_1[,ti,d,j] = m_1
      ds = pCopula(ms[,ti,,j],cop)-pCopula(cbind(ms_1[,ti,1,j],ms[,ti,2,j],ms[,ti,3,j]),cop) - pCopula(cbind(ms[,ti,1,j],ms_1[,ti,2,j],ms[,ti,3,j]),cop)-pCopula(cbind(ms[,ti,1,j],ms[,ti,2,j],ms_1[,ti,3,j]),cop)+
        pCopula(cbind(ms_1[,ti,1,j],ms_1[,ti,2,j],ms[,ti,3,j]),cop)+pCopula(cbind(ms[,ti,1,j],ms_1[,ti,2,j],ms_1[,ti,3,j]),cop) + pCopula(cbind(ms_1[,ti,1,j],ms[,ti,2,j],ms_1[,ti,3,j]),cop)-pCopula(ms_1[,ti,,j],cop)
      ds[ds < .Machine$double.eps]=.Machine$double.eps
      res=res-sum(log(ds)*W[,ti,j])
    }}
  res
  
  
}

tgtC=function(theta,ms,ms_1,W,k,Ti){
  res=0
  copU=claytonCopula(theta,dim = d)
  
  for(ti in 1:Ti){
    for(j in 1:k){
      ds = pCopula(ms[,ti,,j],copU)-pCopula(cbind(ms_1[,ti,1,j],ms[,ti,2,j],ms[,ti,3,j]),copU) - pCopula(cbind(ms[,ti,1,j],ms_1[,ti,2,j],ms[,ti,3,j]),copU)-pCopula(cbind(ms[,ti,1,j],ms[,ti,2,j],ms_1[,ti,3,j]),copU)+
        pCopula(cbind(ms_1[,ti,1,j],ms_1[,ti,2,j],ms[,ti,3,j]),copU)+pCopula(cbind(ms[,ti,1,j],ms_1[,ti,2,j],ms_1[,ti,3,j]),copU) + pCopula(cbind(ms_1[,ti,1,j],ms[,ti,2,j],ms_1[,ti,3,j]),copU)-pCopula(ms_1[,ti,,j],copU)
      ds[ds < .Machine$double.eps]=.Machine$double.eps
      res=res-sum(log(ds)*W[,ti,j])
    }}
  res
}

Cmass <- function(Y,X,al,be,k,d,p,cop){
  
  n=dim(Y)[1]
  Ti=dim(Y)[2]
  d=dim(Y)[3]
  Cat = length(table(Y[,,d]))
  ms=ms_1=linpred<-array(NA,c(n,Ti,d,k))
  logits = array(NA,c(n,Ti,Cat-1,k))
  Px = array(NA,c(n,Ti,Cat,k))
  Dens = array(NA,c(n,Ti,k))
  
  for(ti in 1:Ti){
    for(j in 1:k){  
      for(h in 1:(d-1)){
        linpred[,ti,h,j]  = g[[h]](al[h,j]+X[,ti,h,]%*%be[h,])
      }
      logits[,ti,,j] = X[,ti,d,]%*%be[d,]
    }
    for(j in 1:k){
      logits[,ti,,j] = logits[,ti,,j]+t(matrix(rep(al[-c(1:(d-1)),j],n),Cat-1,n))
      Px[,ti,,j] =  t(matrix(apply(logits[,ti,,j],1,g[[3]]),Cat,n))  
    }
  }
  
  
  for(ti in 1:Ti){
    for(j in 1:k){
      for(h in 1:(d-1)){
        ms[,ti,h,j]=pdens[[h]](Y[,ti,h],linpred[,ti,h,j])
        ms_1[,ti,h,j]=pdens[[h]](Y[,ti,h]-1,linpred[,ti,h,j])
      }
      ms[,ti,d,j] = sapply(1:n, function(i) pdens[[d]](Y[i,ti,d],Px[i,ti,,j]))
      m_1 = as.numeric(sapply(1:n, function(i) pdens[[d]](Y[i,ti,d]-1,Px[i,ti,,j])))
      m_1[is.na(m_1)]=0
      ms_1[,ti,d,j] = m_1
    }}
  
  
  for(ti in 1:Ti){
    for(j in 1:k){
      ds = pCopula(ms[,ti,,j],cop)-pCopula(cbind(ms_1[,ti,1,j],ms[,ti,2,j],ms[,ti,3,j]),cop) - pCopula(cbind(ms[,ti,1,j],ms_1[,ti,2,j],ms[,ti,3,j]),cop)-pCopula(cbind(ms[,ti,1,j],ms[,ti,2,j],ms_1[,ti,3,j]),cop)+
        pCopula(cbind(ms_1[,ti,1,j],ms_1[,ti,2,j],ms[,ti,3,j]),cop)+pCopula(cbind(ms[,ti,1,j],ms_1[,ti,2,j],ms_1[,ti,3,j]),cop) + pCopula(cbind(ms_1[,ti,1,j],ms[,ti,2,j],ms_1[,ti,3,j]),cop)-pCopula(ms_1[,ti,,j],cop)
      ds[ds < .Machine$double.eps]=.Machine$double.eps
      Dens[,ti,j]=log(ds)
    }}
  
  return(list(Dens=Dens,linpred=linpred,ms=ms,ms_1=ms_1))
}

EMC_c <- function(Y,X,inits,eps,k,up.theta,verbose){
  
  n = dim(Y)[1]
  Ti = dim(Y)[2]
  d <- dim(Y)[3]
  Cat <- length(table(Y[,,d]))
  p = dim(X)[4]
  al = inits$al
  be = inits$be
  pis = inits$pis
  PI = inits$PI
  tcop = 1.1
  if(up.theta){
    tcop=inits$tcop
  }
  cop = claytonCopula(tcop,dim = d)
  cm = Cmass(Y,X,al,be,k,d,p,cop)
  LC = complik(Y,X,al,be,pis,PI,k,Ti,cm$Dens)
  LU = LC
  lcurr = lupd <- LC$lik
  W <- array(NA,c(n,Ti,k))
  Z <- array(NA,c(n,Ti-1,k,k))
  diff=Inf
  it=0
  if(verbose){
    cat("------------|-------------|-------------|-------------|-------------|\n")
    cat("     it     |     llik    |     diff    |   al[1,1]   |    tcop     |\n")
    cat("------------|-------------|-------------|-------------|-------------|\n")
  }
  
  
  ms=ms_1=linpred<-array(NA,c(n,Ti,d,k))
  logits = array(NA,c(n,Ti,Cat-1,k))
  Px = array(NA,c(n,Ti,Cat,k))
  Dens = array(NA,c(n,Ti,k))
  
  while(diff > eps){
    it = it+1
    lcurr = lupd
    LC = LU
    Dens = cm$Dens
    fwd = LC$fwd
    bwd = LC$bwd
    iL = apply(fwd[,Ti,],1,sumlog)
    for(ti in 1:Ti){
      W[,ti,] = exp(sweep(fwd[,ti,]+bwd[,ti,],1,iL))
    }
    pis <- apply(W[,1,],2,sum)
    pis <- pis/sum(pis)  
    for(ti in 1:(Ti-1)){
      for(j in 1:k){
        for(v in 1:k){
          Z[,ti,j,v] <- exp(bwd[,ti+1,v]+Dens[,ti+1,v]+log(PI[j,v])+fwd[,ti,j]-iL)
        }}}  
    
    PI = apply(Z,3:4,sum)
    PI <- PI/rowSums(PI)
    PI[PI<1e-15]=1e-15
    PI <- PI/rowSums(PI)
    
    for(r in sample(d,d)){
      if(r==d){
        op = optim(c(al[r:(r-1+Cat-1),],be[r,]), function(x) tgt(x,Y,X,W,al,be,k,p,d,cop,r,ms,ms_1,linpred,logits,Px,Dens,Cat),method = "BFGS",control = list(maxit=1))  
        al[r:(r-1+Cat-1),] = op$par[1:((Cat-1)*k)]
        al[r:(r-1+Cat-1),] = apply(al[r:(r-1+Cat-1),],2,sort,decreasing=T)
        be[r,] = op$par[((Cat-1)*k+1):((Cat-1)*k+p)] 
      }else{
        op = optim(c(al[r,],be[r,]), function(x) tgt(x,Y,X,W,al,be,k,p,d,cop,r,ms,ms_1,linpred,logits,Px,Dens,Cat),method = "BFGS",control = list(maxit=1))  
        al[r,]=op$par[1:k]
        be[r,]=op$par[k+1:p]}
    }
    
    if(up.theta){
      upM = Cmass(Y,X,al,be,k,d,p,cop)
      uptcop = optimize(tgtC,c(1e-6,60-1e-5),upM$ms,upM$ms_1,W,k,Ti)
      tcop=uptcop$minimum
      cop = claytonCopula(tcop,dim = d)
    }
    
    o=order(al[1,])
    al[1,] = sort(al[1,])
    al[-1,] = al[-1,o]
    for(ti in 1:(Ti-1)){
      W[,ti,] = W[,ti,o]
      Z[,ti,,] = Z[,ti,o,o]
    }
    W[,Ti,] = W[,Ti,o]
    pis=pis[o]
    PI = PI[o,o]
    
    
    cm = Cmass(Y,X,al,be,k,d,p,cop)  
    LU = complik(Y,X,al,be,pis,PI,k,Ti,cm$Dens)  
    lupd = LU$lik
    diff = lupd-lcurr
    if(verbose){cat(sprintf("%11g",c(it,lupd,diff,al[1,1],tcop)),sep=" | ","\n")}
  }
  
  return(list(LU=LU,pis=pis,PI=PI,tcop=tcop,al=al,be=be,nui=nui,W=W,cm=cm))
  
}

EMC_c=cmpfun(EMC_c)