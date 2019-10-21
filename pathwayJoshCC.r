#
# Cleaned on July 11th.
# Need to conduct simple sub-selection test with p-values from quadratic test
# Return p-values, set of selected genes and gene p-values and inidivudal assocaitions

calcWithJosh =  function(YY,E,M,Cov,KList,level){
  lc = 2+length(Cov[1,])
pvalQYVw1 = NULL
pvalQEVw1 = NULL

Y = YY
Xc = scale(M,center = TRUE,scale = FALSE)
Y = scale(Y,center = TRUE,scale = FALSE)

rr = lm(Y~E)
Y = rr$residuals
AA = t(Xc)%*%Xc

sigmay = as.numeric(var(Y))
SigmaY = sigmay*AA*0
sigmae = as.numeric(var(E))
SigmaE = sigmae*AA*0

SY = t(Xc)%*%Y*0
SE = t(Xc)%*%E*0



Y =YY ##<--on original Y


LLL = length(KList)
for (jjjj in 1:LLL){
  
  kkk = which((KList[[jjjj]])!=0)
  XX = M[,kkk]

  aa = glm(Y~E+Cov+XX,family='binomial')
  betaY = summary(aa)$coefficients[-c(1:lc),1]
  YY = cbind(rep(1,length(E)),E,Cov,XX)
  MY = vcov(aa)[-c(1:lc),-c(1:lc)]
  vY = MY
  
  ZZ = NULL
  ME = NULL
  betaE = betaY*0
  for (jj in 1:length(XX[1,])){
    
    XXc = XX[,jj]
    XXc[Y==1] = XX[Y==1,jj] - mean(XX[Y==1,jj])
    XXc[Y==0] = XX[Y==0,jj] - mean(XX[Y==0,jj])
    Ec = E
    Ec[Y==1] = E[Y==1] - mean(E[Y==1])
    Ec[Y==0] = E[Y==0] - mean(E[Y==0])
    bb = lm(XXc~Ec+Cov)
    betaE[jj] = summary(bb)$coefficients[2,1]
    ZZ = cbind(ZZ,bb$residuals)
    ME = c(ME,summary(bb)$coefficients[2,2])
    
  }
  
  vE =diag((ME))%*%cor(ZZ)%*%diag((ME))
  
  SigmaY[kkk,kkk] = vY
  SigmaE[kkk,kkk] = vE
  
  SY[kkk,1] = betaY
  SE[kkk,1] = betaE
} 


pvalueYw1= NULL
pvalueEw1= NULL

for (j in 1:LLL){

kkk = which((KList[[j]])!=0)


SigmaYS = SigmaY[kkk,kkk]
SigmaES = SigmaE[kkk,kkk]
ZE = SE[kkk]
ZY = SY[kkk]

lambda = eigen(SigmaES)$values
lambda= Re(lambda)
lambda[lambda<0]=0
qE = sum(ZE^2)
pvalQE = liu(qE, lambda)
pvalueEw1 = c(pvalueEw1,pvalQE)

lambda = eigen(SigmaYS)$values
lambda= Re(lambda)
lambda[lambda<0]=0
qY = sum(ZY^2)
pvalQY = liu(qY, lambda)
pvalueYw1 = c(pvalueYw1,pvalQY)
}

ky = length(which(pvalueYw1<level))
ke = length(which(pvalueEw1<level))

kkkGOpt  =   which(pvalueYw1<=min(0.025/ke,0.025) & pvalueEw1<=min(0.025/ky,0.025))

pvalQEV= pmin(pvalueEw1*max(ky,1),1)
pvalQYV= pmin(pvalueYw1*max(ke,1),1)

pppGOpt =pmin(2*pmax(pvalQEV,pvalQYV),1)

return (list(pppNaiveQ = pppGOpt,selectedNaive = kkkGOpt, pvalQYwithW1=pvalueYw1,pvalQEwithW1 = pvalueEw1))
}
