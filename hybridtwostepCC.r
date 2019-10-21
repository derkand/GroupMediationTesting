#
# post selection procedure 
#

## 1) Test each familiy 
##    Select component that pass both tests  (alpha<level and beta<0.1)
## 2) Create set of values restricted by selection 
## 3) Get conditional p-values
## 4) Get test statistics for testing mediation 


library('CompQuadForm')
library('MASS')


source('kernelHybrid.r')
source('subetMediation.r')
#source('calc_Summary.r')


source('intertesection.r')
source('calc_cond_pvalue.r')



#
# Calculate conditional p-value based on the predefined families defined by Kernels (p by p)
# level - level of the first selection
#  

calcpvalueTS =  function(YY,E,M,Cov,KList,level){
  lc = 2+length(Cov[1,])
Y=YY
Xc = scale(M,center = TRUE,scale = FALSE)
Y = scale(Y,center = TRUE,scale = FALSE)
E  = scale(E,center = TRUE,scale = FALSE)
rr = lm(Y~E)
Y = rr$residuals
L = length(Xc[1,])
AA = t(Xc)%*%Xc

sigmay = as.numeric(var(Y))
SigmaY = sigmay*AA
sigmae = as.numeric(var(E))
SigmaE = sigmae*AA

SY = t(Xc)%*%Y
SE = t(Xc)%*%E
Y=YY

LLL = length(KList)


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


pvalueYN = rep(NA,L)
pvalueEN = rep(NA,L)

ZYsq = rep(NA,L)
ZEsq = rep(NA,L)

for (i in 1:L){

pvalueYN[i] =2*(1-pnorm(abs(SY[i]/sqrt(SigmaY[i,i]))))
pvalueEN[i] =2*(1-pnorm(abs(SE[i]/sqrt(SigmaE[i,i]))))

ZYsq[i] =abs(SY[i]/sqrt(SigmaY[i,i]))^2
ZEsq[i] =abs(SE[i]/sqrt(SigmaE[i,i]))^2
}


ObjSingle = subsetMediatonWithoutW(pvalueYN,pvalueEN,KList,level)
pvalMin = ObjSingle$pvalMin



KernelsT = KernelsruleHybrid(KList,SE,SY,Xc,SigmaE,SigmaY,level)
KernelsSE = KernelsT$KernelsSE
KernelsSY = KernelsT$KernelsSY

pvalueYC = rep(1,L)
pvalueEC = rep(1,L)


if (length(KernelsSE)==0 & length(KernelsSY)==0){
  a=NULL
list(ObjMin =ObjSingle,pvalMin = pvalMin,ObjTS=a,pvalTS =rep(1,LLL))
}

LL = length(KernelsSE)
cat('Assoc with E',LL,'\n')
if (LL>0){
for (jjc in 1:LL){
KernelS = list()
SET = which((KernelsSE[[jjc]])!=0)
kkk = SET

XXX = diag(SY[SET]^2)
KernelS[[1]] = XXX
for (i in 1:length(SET)){
SetsE = create_set_Ith(i,SE[kkk],KernelS,SigmaE[kkk,kkk],level)
pvalueEE = get_pvalue_i(i,SE[kkk],SigmaE[kkk,kkk],SetsE)
pvalueEC[SET[i]] =pvalueEE$C
}
}
}

LL = length(KernelsSY)
cat('Assoc with Y',LL,'\n')
cat(LL,'\n')
if (LL>0){
cat(LL,'\n')
#save(KernelsSY,file='save')
for (jjc in 1:LL){
KernelS = list()
SET = which((KernelsSY[[jjc]])!=0)
kkk = SET

XXX = diag(SE[SET]^2)
KernelS[[1]] = XXX
for (i in 1:length(SET)){
SetsY = create_set_Ith(i,SY[kkk],KernelS,SigmaY[kkk,kkk],level)
pvalueYY = get_pvalue_i(i,SY[kkk],SigmaY[kkk,kkk],SetsY)
pvalueYC[SET[i]] =pvalueYY$C
}
}
}


ObjTS = subsetMediatonWithW(pvalueYC,pvalueEC,pvalueYN,pvalueEN,ZYsq,ZEsq,KList,level)
pvalTS = ObjTS$pvalTS 

return (list(ObjMin =ObjSingle,pvalMin = pvalMin,ObjTS =ObjTS,pvalTS = pvalTS))
}










