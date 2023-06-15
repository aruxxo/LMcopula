library(copula)
library(snipEM)
library(MultiLCIRT)
library(MASS)
library(flexmix)
library(compiler)

g = c(function(x) {x=x},
      function(x) exp(x)/(1+exp(x)),
      function(x) exp(x))
qdens = c(function(qu,par,nui) qnorm(qu,par,nui),
          function(qu,par,nui) qbinom(qu,1,par),
          function(qu,par,nui) qpois(qu,par))

pdens = c(function(x,par,nui) pnorm(x,par,nui),
          function(x,par,nui) pbinom(x,1,par),
          function(x,par,nui) ppois(x,par))

ddens <- function(x,par,nui){return(dnorm(x,par,nui,log = T))}

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

tgt_c <- function(theta,Y,X,W,al,be,nui,k,p,d,tcop,r,ms,ms_1,linpred){
  res = 0
  if(r==1){
    al[r,]=theta[1:k]
    be[r,]=theta[(1:p)+k]
    nui[[r]] = theta[length(theta)]
  }else{
    al[r,]=theta[1:k]
    be[r,]=theta[(1:p)+k]  
  }
  
  n=dim(Y)[1]
  Ti=dim(Y)[2]
  d=dim(Y)[3]
  
  for(ti in 1:Ti){
    for(j in 1:k){  
      for(h in 1:d){
        linpred[,ti,h,j]  = g[[h]](al[h,j]+X[,ti,h,]%*%be[h,])}
      for(h in 1:d){
        ms[,ti,h,j]=pdens[[h]](Y[,ti,h],linpred[,ti,h,j],exp(nui[[h]]))
        #  ms[,ti,h,j]=pdens[[h]](Y[,ti,h],linpred[,ti,h,j],nui[[h]])
        if(h != 1){
          ms_1[,ti,h,j]=pdens[[h]](Y[,ti,h]-1,linpred[,ti,h,j],nui[[h]])}
      }}}
  
  for(ti in 1:Ti){
    for(j in 1:k){
      fd = ddens(Y[,ti,1],linpred[,ti,1,j],exp(nui[[1]]))
      #  fd = ddens(Y[,ti,1],linpred[,ti,1,j],nui[[1]])
      ds = apply(ms[,ti,,j],1,ctilde,tcop)-apply(cbind(ms[,ti,1,j],ms_1[,ti,2,j],ms[,ti,3,j]),1,ctilde,tcop)-apply(cbind(ms[,ti,1,j],ms[,ti,2,j],ms_1[,ti,3,j]),1,ctilde,tcop)+apply(cbind(ms[,ti,1,j],ms_1[,ti,2,j],ms_1[,ti,3,j]),1,ctilde,tcop)
      ds[ds < .Machine$double.eps]=.Machine$double.eps
      res=res-sum((log(ds)+fd)*W[,ti,j])
    }
  }
  res
}


tgtC_c=function(theta,ms,ms_1,W,k,Ti,nui,linpred) {
  res=0
  copU=frankCopula(theta,dim = d)
  tcopU = theta
  
  for(ti in 1:Ti){
    for(j in 1:k){
      fd = ddens(Y[,ti,1],linpred[,ti,1,j],exp(nui[[1]]))
      ds = apply(ms[,ti,,j],1,ctilde,tcopU)-apply(cbind(ms[,ti,1,j],ms_1[,ti,2,j],ms[,ti,3,j]),1,ctilde,tcopU)-apply(cbind(ms[,ti,1,j],ms[,ti,2,j],ms_1[,ti,3,j]),1,ctilde,tcopU)+apply(cbind(ms[,ti,1,j],ms_1[,ti,2,j],ms_1[,ti,3,j]),1,ctilde,tcopU)
      ds[ds < .Machine$double.eps]=.Machine$double.eps
      res=res-sum((log(ds)+fd)*W[,ti,j])
    }}
  res
}

ctilde <- function(mg,xi){
  if(xi !=0){
  dF = (exp(-xi*mg[1])*prod(exp(-xi*mg[-1])-1))/(((exp(-xi)-1)^(length(mg)-1)) * (1+prod(exp(-xi*mg)-1)/(exp(-xi)-1)^(length(mg)-1)))
  }else{
  dF = prod(mg)  
  }
  return(dF)
}

Cmass_c <- function(Y,X,al,be,k,d,p,tcop,nui){
  
  n=dim(Y)[1]
  Ti=dim(Y)[2]
  d=dim(Y)[3]
  ms=ms_1=linpred<-array(NA,c(n,Ti,d,k))
  Dens = array(NA,c(n,Ti,k))
  
  for(ti in 1:Ti){
    for(j in 1:k){  
      for(h in 1:d){
        linpred[,ti,h,j]  = g[[h]](al[h,j]+X[,ti,h,]%*%be[h,])
      }
      for(h in 1:d){
        ms[,ti,h,j]=pdens[[h]](Y[,ti,h],linpred[,ti,h,j],exp(nui[[h]]))
        if(h != 1){
          ms_1[,ti,h,j]=pdens[[h]](Y[,ti,h]-1,linpred[,ti,h,j],nui[[h]])}
      }
    }
  }
  
  for(ti in 1:Ti){
    for(j in 1:k){
      fd = ddens(Y[,ti,1],linpred[,ti,1,j],exp(nui[[1]]))
      ds = apply(ms[,ti,,j],1,ctilde,tcop)-apply(cbind(ms[,ti,1,j],ms_1[,ti,2,j],ms[,ti,3,j]),1,ctilde,tcop)-apply(cbind(ms[,ti,1,j],ms[,ti,2,j],ms_1[,ti,3,j]),1,ctilde,tcop)+apply(cbind(ms[,ti,1,j],ms_1[,ti,2,j],ms_1[,ti,3,j]),1,ctilde,tcop)
      ds[ds < .Machine$double.eps]=.Machine$double.eps
      Dens[,ti,j]=log(ds)+fd
    }}
  
  return(list(Dens=Dens,linpred=linpred,ms=ms,ms_1=ms_1))
}

EMC_fc <- function(Y,X,inits,eps,k,up.theta,verbose){
  
  n = dim(Y)[1]
  Ti = dim(Y)[2]
  d = dim(Y)[3]
  p = dim(X)[4]
  al = inits$al
  be = inits$be
  nui = inits$nui
  pis = inits$pis
  PI = inits$PI
  tcop=inits$tcop
  cop = frankCopula(tcop,dim = d)
  cm = Cmass_c(Y,X,al,be,k,d,p,tcop,nui)
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
      if(r==1){
        # op = optim(c(al[r,],be[r,],log(nui[[r]])), function(x) tgt_c(x,Y,X,W,al,be,nui,k,p,d,cop,r,ms,ms_1,linpred),method = "BFGS",control = list(maxit=1))  
        op = optim(c(al[r,],be[r,],nui[[r]]), function(x) tgt_c(x,Y,X,W,al,be,nui,k,p,d,tcop,r,ms,ms_1,linpred),method = "BFGS",control = list(maxit=1))  
        al[r,] = op$par[1:k]
        be[r,] = op$par[k+1:p]
        nui[[r]] <- op$par[length(op$par)]
      }else{
        op = optim(c(al[r,],be[r,]), function(x) tgt_c(x,Y,X,W,al,be,nui,k,p,d,tcop,r,ms,ms_1,linpred),method = "BFGS",control = list(maxit=1))  
        al[r,]=op$par[1:k]
        be[r,]=op$par[k+1:p]}
    }
    
    if(up.theta){
      upM = Cmass_c(Y,X,al,be,k,d,p,tcop,nui)
      uptcop = optimize(tgtC_c,c(0,15),ms=upM$ms,ms_1=upM$ms_1,W=W,k=k,Ti=Ti,nui=nui,linpred=upM$linpred)
      tcop=uptcop$minimum
      cop = frankCopula(tcop,dim = d)
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
    
    
    cm = Cmass_c(Y,X,al,be,k,d,p,tcop,nui)  
    LU = complik(Y,X,al,be,pis,PI,k,Ti,cm$Dens)  
    lupd = LU$lik
    diff = lupd-lcurr
    if(verbose){cat(sprintf("%11g",c(it,lupd,diff,al[1,1],tcop)),sep=" | ","\n")}
  }
  
  return(list(LU=LU,pis=pis,PI=PI,tcop=tcop,al=al,be=be,nui=nui,W=W,cm=cm))
  
}

EMC_fc=cmpfun(EMC_fc)

simu_c=function(n,Ti,d,p,al,be,cop,pis,PI,g,qdens,nui) {
  
  X=array(rnorm(n*Ti*d*p,sd=0.1),c(n,Ti,d,p))
  Y=array(NA,c(n,Ti,d)) 
  k=length(pis)
  U=matrix(NA,n,Ti)
  U[,1]=sample(k,n,replace=T,prob=pis)
  for(t in 2:Ti) {
    for(j in 1:k) {
      if(any(U[,t-1]==j)) {
        w=which(U[,t-1]==j)
        U[w,t]=sample(k,length(w),replace=T,prob=PI[j,])} 
    }}
  
  jc = array(NA,c(n,Ti,d))
  linpred = array(NA,c(n,Ti,d))
  for(ti in 1:Ti) {
    jc[,ti,] <- rCopula(n, cop)
    for(h in 1:(d)){
      for(i in 1:n){
        linpred[i,ti,h] <- g[[h]](al[h,U[i,ti]] + X[i,ti,h,]%*%be[h,])
        Y[i,ti,h] <- qdens[[h]](jc[i,ti,h],linpred[i,ti,h],exp(nui[[h]]))
      }
    }
  }
  
  return(list(Y=Y,X=X,U=U,al=al,be=be,cop=cop,pis=pis,PI=PI,qdens=qdens,g=g,nui=nui,jc=jc,linpred=linpred))}

getInits_c=function(Y,X,d,k){
  n=dim(Y)[1]
  Ti=dim(Y)[2]
  d=dim(Y)[3]
  p = dim(X)[4]
  initials = list()
  initials$al=matrix(NA,d,k)
  initials$be=matrix(NA,d,p)
  initials$nui <- list(NA,NA,NA)
  mapnew=rep(NA,n*Ti)
  pY=as.vector(Y[,,1])
  pX=matrix(X[,,1,],n*Ti,p)
  indata=data.frame(pY,pX,idy=rep(1:n,Ti))
  m1 = flexmix(pY ~ 1+X1+X2|idy, data = indata, k=k, model = FLXMRglmfix(family = "gaussian",fixed = ~X1+X2,varFix = T))
  clmap=m1@cluster
  odp = order(parameters(m1)[1,])
  initials$al[1,]=parameters(m1)[1,odp]
  initials$be[1,]=parameters(m1)[2:(p+1),odp[1]]
  initials$nui[[1]] <- log(parameters(m1)[nrow(parameters(m1)),1])
  for(cl in 1:k){
    mapnew[which(clmap==odp[cl])]=cl  
  }
  pY=as.vector(Y[,,2])
  pX=matrix(X[,,2,],n*Ti,p)
  indata=data.frame(pY,pX,idy=rep(1:n,Ti))
  m2 = flexmix(cbind(pY,1-pY)~1|idy, data = indata,cluster = mapnew,k=k,model = FLXMRglmfix(family = "binomial",fixed = ~pX))
  initials$al[2,]=parameters(m2)[p+1,]
  initials$be[2,]=parameters(m2)[1:p,1]
  pY = as.vector(Y[,,3])
  pX = matrix(X[,,3,],n*Ti,p)
  indata = data.frame(pY,pX,idy=rep(1:n,Ti))
  m3 = flexmix(pY~1|idy,data = indata,k=k,cluster = mapnew,model = FLXMRglmfix(family = "poisson",fixed = ~pX))
  initials$al[3,] <- parameters(m3)[p+1,]
  initials$be[3,] <- parameters(m3)[1:p,1]
  initials$pis=rep(1/k,k)
  initials$PI=matrix((1-0.85)/(k-1),k,k)
  diag(initials$PI)=0.85
  idcorr = combn(1:d,(d-1))
  rhos=rep(NA,d)
  initials$tcop = max(apply(Y,2:3,sd))
  initials$cop = frankCopula(initials$tcop,dim = d)
  # cm = Cmass_c(Y,X,initials$al,initials$be,k,d,p,initials$cop,initials$nui)
  # LC = complik(Y,X,initials$al,initials$be,pis,PI,k,Ti,cm$Dens)
  # fwd = LC$fwd
  # bwd = LC$bwd
  # iL = apply(fwd[,Ti,],1,sumlog)
  # W=array(NA,c(n,Ti,k))
  # for(ti in 1:Ti){
  #   W[,ti,] = exp(sweep(fwd[,ti,]+bwd[,ti,],1,iL))
  # }
  # cpI = optimize(tgtC_c,c(0+1e-6,15-1e-5),cm$ms,cm$ms_1,W,k,Ti,initials$nui,cm$linpred)
  # initials$tcop = cpI$minimum
  # initials$cop = frankCopula(initials$tcop,dim = d)
  return(initials)
}


