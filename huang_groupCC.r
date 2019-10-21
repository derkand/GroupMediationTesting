#
# Group testing by Huang 2016 and 2018 Biometrics and Annals of Applied Statistics
#
library('CompQuadForm')
library('mvtnorm')

source('huang2.r')
#
# Z scores
# w weights
# Sigma
# 
# Return conditional p-value
#



calcHuangCC =  function(YY,E,M,Cov,KList,level){
  lc = 2+length(Cov[1,])
  Y = YY
  
  pSumL = NULL
  pQ = NULL
  pMaxL = NULL
  
  Xc = scale(M,center = TRUE,scale = FALSE)
  Y = scale(Y,center = TRUE,scale = FALSE)
  rr = lm(Y~E)
  Y = rr$residuals
  L = length(Xc[1,])
  AA = t(Xc)%*%Xc
  
  sigmay = as.numeric(var(Y))
  SigmaY = sigmay*AA*0
  sigmae = as.numeric(var(E))
  SigmaE = sigmae*AA*0
  
  SY = t(Xc)%*%Y*0
  SE = t(Xc)%*%E*0
  
  SYR = SY
  SER = SE
  SigmaYR = SigmaY
  SigmaER = SigmaE
  
  Y = YY ### <-  on original Y
  
  LLL = length(KList)
  for (jjjj in 1:LLL){
    
    kkk = which((KList[[jjjj]])!=0)
    
    XX = M[,kkk]
    
    aa = glm(Y~E+Cov+XX,family='binomial')
    betaY = summary(aa)$coefficients[-c(1:lc),1]
    vY = vcov(aa)[-c(1:lc),-c(1:lc)]
    
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
  
  
  
  LLL = length(KList)
  for (jjjj in 1:LLL){
    
    kkk = which((KList[[jjjj]])!=0)
    
    XX = M[,kkk]
    XXR = NULL
    
    for (jj in 1:length(XX[1,])){
      
      XXc = XX[,jj]
      XXc[Y==1] = XX[Y==1,jj] - mean(XX[Y==1,jj])
      XXc[Y==0] = XX[Y==0,jj] - mean(XX[Y==0,jj])
      Ec = E
      Ec[Y==1] = E[Y==1] - mean(E[Y==1])
      Ec[Y==0] = E[Y==0] - mean(E[Y==0])
      
      r = residuals(lm(XXc~Ec+Cov))
      XXR = cbind(XXR,as.vector(r))
    }
    
    Q = eigen(t(XXR)%*%XXR)$vectors
    #Find Rotation Q
    
    XX = M[,kkk]%*%Q
    
    aa = glm(Y~E+Cov+XX,family='binomial')
    betaY = summary(aa)$coefficients[-c(1:lc),1]
    YY = cbind(rep(1,length(E)),E,XX)
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
    
    SigmaYR[kkk,kkk] = vY
    SigmaER[kkk,kkk] = vE
    
    SYR[kkk,1] = betaY
    SER[kkk,1] = betaE
    
  }
  
  
  
  
  LL = length(KList)
  
  for (j in 1:LL){
    #cat(j,'\n','\n')
    kkk = which((KList[[j]])!=0)
    
    SYS = SY[kkk]
    SES = SE[kkk]
    SigmaYS = SigmaY[kkk,kkk]
    SigmaES = SigmaE[kkk,kkk]
    
    SYSR = SYR[kkk]
    SESR = SER[kkk]
    SigmaYSR = SigmaYR[kkk,kkk]
    SigmaESR = SigmaER[kkk,kkk]
    
    
    xx = gen_boot_samples(SYS,SigmaYS,SES,SigmaES,perm=10000)
    pQ =c(pQ,mean(xx>sum(SYS^2*SES^2)))
    
    pL = get_pvalueSumL(SYS,SigmaYS,SES,SigmaES)
    pSumL = c(pSumL,pL)
    pM = get_pvalueMaxL(SYSR,diag(SigmaYSR),SESR,diag(SigmaESR))
    pMaxL = c(pMaxL,pM)
  }
  
  
  return (list(pvalSumL = pSumL,pvalMaxL=pMaxL,pvalQ=pQ))
}






