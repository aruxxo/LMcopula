source("frankcopula.r")

d=3
p=2
k=2
Ti=4
n=800
pis=c(0.5,0.5)
PI=matrix((1-0.8)/(k-1),k,k)
diag(PI)=0.8
be=matrix(1,d,p)
Cat=4
al= matrix(NA,d+(Cat-2),k)
al[1:(d-1),] = t(matrix(rep(c(-1.5,1.5),d-1),k,d-1))
al[-c(1:(d-1)),] = c(-0.5,-1.5,-2.5,2.5,1.5,0.5)
nui=list(NA,NA,NA)
tcop = 2
cop = frankCopula(tcop,dim = d)
pm = list(n=n,Ti=Ti,k=k,al=al,be=be,pis=pis,PI=PI,tcop=tcop,cop=cop)

set.seed(123)
dat=simu(n,Ti,d,p,al,be,cop,pis,PI,g,qdens,nui)
  
inits=list()
inits <- getInits(dat$Y,dat$X,d,k)

out = EMC_f(dat$Y,dat$X,inits,1e-2,k,TRUE,TRUE)

out$al
out$be
out$pis
out$PI
out$tcop

