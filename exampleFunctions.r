####################################################################################
### This script contains functions that estimate the Minimally Dysregulated (MD) set 
### of features. The main functions are MDset() and MDset2(), the former performs
### computations serially, the latter performs computations in parallel using the 
### doParallel package.
####################################################################################

# compute the jj block of V
Vj.funt = function(sigma,j)
{
  return(outer(sigma[j,],sigma[j,],"*")+sigma[j,j]*sigma)
}

# compute the jh block of V
Vjh.funt = function(sigma,j,h)
{
  return(outer(sigma[h,],sigma[j,],"*")+sigma[j,h]*sigma) 
}

# compute DG[,(j-1)*p+j], gradient wrt sigma[j,j]
DGjj.funt = function(sigma,j)
{
  R = diag(1/sqrt(diag(sigma))) %*% sigma %*% diag(1/sqrt(diag(sigma)))
  p = nrow(sigma)
  temp = matrix(0,nrow=p,ncol=p)
  temp[j,] = -R[j,]/(2*sigma[j,j])
  temp[,j] = -R[j,]/(2*sigma[j,j])
  temp[j,j] = 0
  return(as.vector(temp))
}

# compute DG[,(j-1)*p+i], gradient wrt sigma[i,j]
DGij.funt = function(sigma,i,j)
{
  R = diag(1/sqrt(diag(sigma))) %*% sigma %*% diag(1/sqrt(diag(sigma)))
  p = nrow(sigma)
  temp = matrix(0,nrow=p,ncol=p)  
  temp[i,j] = 1/sqrt(sigma[i,i]*sigma[j,j])
  return(as.vector(temp)) 
}

# a function that computes Gamma
Gamma.funt = function(sigma)
{
  p = nrow(sigma)
  
  # compute DG
  DG = matrix(NA,nrow=p^2,ncol=p^2)
  for (j in 1:p) for (i in 1:p) 
  {
    k = (j-1)*p+i
    if(i==j){DG[,k] = DGjj.funt(sigma,j)}
    if(i!=j){DG[,k] = DGij.funt(sigma,i,j)}
  }
  
  # compute V
  V = matrix(NA,nrow=p^2,ncol=p^2)
  for (j in 1:p) for (h in 1:p) 
  {
    if(j==h){V[(j-1)*p+(1:p),(h-1)*p+(1:p)] = Vj.funt(sigma,j)} 
    if(j!=h){V[(j-1)*p+(1:p),(h-1)*p+(1:p)] = Vjh.funt(sigma,j,h)} 
  }
  
  return(DG%*%V%*%t(DG))
}

Gamma.funt3 = function(sigma)
{
  p = nrow(sigma)
  
  # compute DG
  DG = matrix(NA,nrow=p^2,ncol=p^2)
  for (j in 1:p) for (i in 1:p) 
  {
    k = (j-1)*p+i
    if(i==j){DG[,k] = DGjj.funt(sigma,j)}
    if(i!=j){DG[,k] = DGij.funt(sigma,i,j)}
  }
  
  # compute block diagonal of Gamma
  preGamma = matrix(NA,nrow=p^2,ncol=p^2)
  Gamma = matrix(NA,nrow=p^2,ncol=p^2)
  
  for (j in 1:p) for (h in 1:p)
  {
    Vj..h = if(j!=h) {Vjh.funt(sigma,j,h)} else {Vj.funt(sigma,j)}
    #       a = sapply(1:p,function(x){1/sqrt(sigma[j,j]*sigma[x,x])}); a[j]=0 # not necessarily faster...
    #       b = -1/(2*sigma[j,j]^(3/2)) * sapply(1:p,function(x){sigma[x,j]/sqrt(sigma[x,x])}); b[j]=0
    #       temp = a*Vj..h + t(sapply(1:p,function(x){b[x]*Vj..h}))
    temp = DG[(j-1)*p+1:p,(j-1)*p+1:p] %*% Vj..h 
    for(k in (1:p)[-j])
    {
      DG_jkkk = -sigma[k,j]/(2*sqrt(sigma[j,j]*sigma[k,k])*sigma[k,k])
      V_kk.h = 2 * sigma[k,h] * sigma[k,] 
      temp[k,] = temp[k,] + DG_jkkk * V_kk.h      
    }
    
    preGamma[(j-1)*p+1:p,(h-1)*p+1:p] = temp
  }   
  for (j in 1:p)
  {
    temp =  preGamma[(j-1)*p+1:p,(j-1)*p+1:p] %*% t(DG)[(j-1)*p+1:p,(j-1)*p+1:p] 
    for(k in (1:p)[-j])
    {
      PG_j.kk = preGamma[(j-1)*p+1:p,(k-1)*p+1:p][,k]
      tDG_kkkj = -sigma[k,j]/(2*sqrt(sigma[j,j]*sigma[k,k])*sigma[k,k])    
      temp[,k] = temp[,k] + tDG_kkkj * PG_j.kk     
    }
    Gamma[(j-1)*p+1:p,(j-1)*p+1:p] = temp
  }
  
  return(Gamma)
}

# compute the correlation test statistic for j=1,...,p
cor.test.stat = function(R.hat1,R.hat2,Gamma.hat1,Gamma.hat2,n1,n2)
{
  p = ncol(R.hat1)
  test.stat = NULL
  for (j in 1:p)
  {
    Rj.hat1 = R.hat1[j,-j]
    Rj.hat2 = R.hat2[j,-j]
    Gammaj.hat1 = Gamma.hat1[(j-1)*p+(1:p)[-j],(j-1)*p+(1:p)[-j]]
    Gammaj.hat2 = Gamma.hat2[(j-1)*p+(1:p)[-j],(j-1)*p+(1:p)[-j]]
    
    Tj = n1*n2/(n1+n2) * (Rj.hat1 - Rj.hat2) %*% solve(n2/(n1+n2)*Gammaj.hat1 + n1/(n1+n2)*Gammaj.hat2) %*%
      (Rj.hat1 - Rj.hat2)
    
    test.stat = c(test.stat,Tj)
  }
  
  return(data.frame("covariate"=1:p,test.stat))
}

max.step = function(M1,M2,V1,V2,n1,n2)
{
  test.step = cor.test.stat(M1,M2,V1,V2,n1,n2)  
  return(c(which.max(test.step$test.stat),max(test.step$test.stat)))
}

# determine whether to reject or not reject for one test statistic
decision.funt = function(test.stat,level,deg)
{
  if(is.na(test.stat)) {return(NA)} else if(test.stat>qchisq(1-level,df=deg)) 
  {return("reject")} else {return("do not reject")}
}

MDset = function(data1,data2,alpha=0.05)
{
  n1 = nrow(data1)
  n2 = nrow(data2)
  s1 = cov(data1)
  s2 = cov(data2)
  M1 = diag(1/sqrt(diag(s1))) %*% s1 %*% diag(1/sqrt(diag(s1)))
  M2 = diag(1/sqrt(diag(s2))) %*% s2 %*% diag(1/sqrt(diag(s2)))
  V1 = Gamma.funt3(s1)
  V2 = Gamma.funt3(s2) 
  
  p = nrow(M1)
  covariates = 1:p
  steps = p-1
  
  test.stat = data.frame(matrix(NA,nrow=steps,ncol=p+2)); names(test.stat) = c("j.max","test.stat",1:p)
  
  # compute the first test stat 
  temp = max.step(M1,M2,V1,V2,n1,n2)
  test.stat[1,] = c(temp[1],temp[2],rep(F,p))
  
  # compute the ith test stat
  i = 1
  while(decision.funt(test.stat[i,2],level=alpha/(p-i+1),deg=p-i)=="reject")
  {
    # compute test statistics after removing i variables 
    combos = rbind(combn(p,i),NA,NA)
    for (j in 1:ncol(combos))
    {
      v.remove = combos[1:i,j]
      M1.temp = M1[-v.remove,-v.remove]
      M2.temp = M2[-v.remove,-v.remove]
      V1.temp = V1[-c(sapply(v.remove,function(x){(x-1)*p+1:p}),sapply(v.remove,function(x){p*(0:(p-1))+x})),
                   -c(sapply(v.remove,function(x){(x-1)*p+1:p}),sapply(v.remove,function(x){p*(0:(p-1))+x}))]
      V2.temp = V2[-c(sapply(v.remove,function(x){(x-1)*p+1:p}),sapply(v.remove,function(x){p*(0:(p-1))+x})),
                   -c(sapply(v.remove,function(x){(x-1)*p+1:p}),sapply(v.remove,function(x){p*(0:(p-1))+x}))]
      temp = max.step(M1.temp,M2.temp,V1.temp,V2.temp,n1,n2)
      combos[-(1:i),j] = c(covariates[-v.remove][temp[1]],temp[2])
    }
    
    # identify the min among all combos
    v.remove = combos[1:i,which.min(combos[i+2,])]
    test.stat[i+1,1] = combos[i+1,which.min(combos[i+2,])]
    test.stat[i+1,2] = min(combos[i+2,])
    test.stat[i+1,-(1:2)] = (1:p)%in%v.remove
    
    i=i+1
  }
  
  # return the MD set
  mdset = (1:p)[test.stat[i,-(1:2)]==1]  
  if(length(mdset)==0) return(0) else return(mdset)   
}

library(doParallel)
library(dplyr)
MDset2 = function(data1,data2,alpha=0.05)
{
  n1 = nrow(data1)
  n2 = nrow(data2)
  s1 = cov(data1)
  s2 = cov(data2)
  M1 = diag(1/sqrt(diag(s1))) %*% s1 %*% diag(1/sqrt(diag(s1)))
  M2 = diag(1/sqrt(diag(s2))) %*% s2 %*% diag(1/sqrt(diag(s2)))
  V1 = Gamma.funt3(s1)
  V2 = Gamma.funt3(s2) 
  
  p = nrow(M1)
  covariates = 1:p
  steps = p-1
  
  test.stat = data.frame(matrix(NA,nrow=steps,ncol=p+2)); names(test.stat) = c("j.max","test.stat",1:p)
  
  # compute the first test stat 
  temp = max.step(M1,M2,V1,V2,n1,n2)
  test.stat[1,] = c(temp[1],temp[2],rep(F,p))

  # set up for parallel computations
  nworkers <- detectCores()
  cl <- makeCluster(nworkers)
  registerDoParallel(cl)
  
  # compute the ith test stat
  i = 1
  while(decision.funt(test.stat[i,2],level=alpha/(p-i+1),deg=p-i)=="reject")
  {
    # compute test statistics after removing i variables
    combos = rbind(combn(p,i),NA,NA)
    results = foreach(j = 1:ncol(combos),  
                      .combine      = "bind_rows",
                      .multicombine = TRUE,
                      .export=c("max.step","cor.test.stat")) %dopar% { 
      v.remove = combos[1:i,j]
      M1.temp = M1[-v.remove,-v.remove]
      M2.temp = M2[-v.remove,-v.remove]
      V1.temp = V1[-c(sapply(v.remove,function(x){(x-1)*p+1:p}),sapply(v.remove,function(x){p*(0:(p-1))+x})),
                 -c(sapply(v.remove,function(x){(x-1)*p+1:p}),sapply(v.remove,function(x){p*(0:(p-1))+x}))]
      V2.temp = V2[-c(sapply(v.remove,function(x){(x-1)*p+1:p}),sapply(v.remove,function(x){p*(0:(p-1))+x})),
                 -c(sapply(v.remove,function(x){(x-1)*p+1:p}),sapply(v.remove,function(x){p*(0:(p-1))+x}))]
      temp = max.step(M1.temp,M2.temp,V1.temp,V2.temp,n1,n2)
      data.frame("feature"=covariates[-v.remove][temp[1]],"teststat"=temp[2]) }
    combos[-(1:i),] = t(results)
    
    # identify the min among all combos
    v.remove = combos[1:i,which.min(combos[i+2,])]
    test.stat[i+1,1] = combos[i+1,which.min(combos[i+2,])]
    test.stat[i+1,2] = min(combos[i+2,])
    test.stat[i+1,-(1:2)] = (1:p)%in%v.remove
    
    i=i+1
  }
  stopCluster(cl)

  # return the MD set
  mdset = (1:p)[test.stat[i,-(1:2)]==1]  
  if(length(mdset)==0) return(0) else return(mdset)  
}