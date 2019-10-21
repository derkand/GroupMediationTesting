subsetMediatonWithoutW = function(PY,PE,KList,ttt){
  PY[is.na(PY)]=1
  PE[is.na(PE)]=1
  ky = length(which(PY<ttt))
  ke = length(which(PE<ttt))
  
  kkk =   which(PY<=min(0.025/ke,0.025) & PE<=min(0.025/ky,0.025))
  
  pvalQEV= pmin(PE*max(ky,1),1)
  pvalQYV= pmin(PY*max(ke,1),1)
  
  pppGOpt =2*pmax(pvalQEV,pvalQYV)
  
  
  L = length(KList)
  
  pathWayG = NULL
  pvalG = NULL
  for (i in 1:L){
    kkk = which((KList[[i]])!=0)
    pval = pppGOpt[kkk]
    pathWayG = c(pathWayG,i)
    pvalG  = c(pvalG,min(min(pval),1))
  }
  
  #cat('>>',pvalG,'\n')
  return (list(pvalAdjE = pvalQEV,pvalAdjY = pvalQYV,indexSG = pathWayG,pvalMin = pvalG))
}




subsetMediatonWithW = function(PY,PE,PYN,PEN,ZYsq,ZEsq,KList,ttt){
  
  PY1 = PY
  PE1 = PE
  
  PYN[(PE)==1]=1
  PEN[(PY)==1]=1
  
  ZYsq[PYN>ttt]=0
  ZEsq[PEN>ttt]=0
  
  ke = sum(ZEsq)/ZEsq
  ky = sum(ZYsq)/ZYsq
  
  if (sum(ZEsq)==0){ke=rep(Inf,length(ke))}
  if (sum(ZYsq)==0){ky=rep(Inf,length(ky))}
  
  SelectedM =  which(PY<=pmin(ttt/ke,0.025) & PE<=pmin(ttt/ky,0.025))
  
  PY = pmin(PY*ke,1)
  PE = pmin(PE*ky,1)
  
  pvalSingelM = 2*pmax(PE,PY)
  
  
  L = length(KList)
  
  pathWayG = NULL
  pvalG = NULL
  
  for (i in 1:L){
    kkk = which((KList[[i]])!=0)
    pval = pvalSingelM[kkk]
    pathWayG = c(pathWayG,i)
    pvalG  = c(pvalG,min(min(pval),1))
  }
  
  
  return (list(selectedM = SelectedM,pEAdj = PE,pYAdj = PY,selectedG = pathWayG,pvalTS = pvalG))
}
