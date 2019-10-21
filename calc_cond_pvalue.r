#:
# Calculate conditional p-value
#

CalcCondPvalue = function(t,Set,se){
L = length(Set)
whole = 0
top = 0
t = as.numeric(t)

for (i in 1:L){
top = top + (pnorm(min(Set[[i]][2]/se,-t))-pnorm(min(Set[[i]][1]/se,-t)))
whole = whole + (pnorm(Set[[i]][2]/se)-pnorm(Set[[i]][1]/se)) 
}

for (i in 1:L){
top = top + (pnorm(max(Set[[i]][2]/se,t))-pnorm(max(Set[[i]][1]/se,t)))
#whole = whole + (pnorm(Set[[i]][2])-pnorm(Set[[i]][1])) 
}

if (whole==0){
save(Set,se,file='Set')
return (-Inf)
}

return (top/whole)
}

#
# Create restricted space for ith biomarker
# Extend to correlated pathways!
#

create_set_Ith = function(I,S,SelectedKernels,Sigma,level){
Set = list()
L = length(SelectedKernels)
LL = length(Sigma[1,])
he =  matrix(0,nrow=LL,ncol=1)
he[I,1] = 1/sqrt(Sigma[I,I])

#
# creating sets
#
kk = 1
for (i in 1:L){
K = SelectedKernels[[i]]
#cat('>>',dim(K),'\n')

kkk = sort(unique(c(I,which(colSums(abs(K))!=0),which(rowSums(abs(K))!=0))))
#cat(length(kkk),'\n')
heS = he[kkk,]
KlistS = K[kkk,kkk]
SS = S[kkk]
SigmaS = Sigma[kkk,kkk]

LL = length(SigmaS[1,])
threshold = calcCrit(KlistS,SigmaS,level)
#cat('threshold',threshold,'\n')

if (sum(K)<0){threshold = -threshold}

cc = 1/as.numeric(t(heS)%*%SigmaS%*%heS)*SigmaS%*%heS 
zz = (diag(1,LL) - cc%*%t(heS))%*%SS
a = t(cc)%*%KlistS%*%cc
b = 2*t(zz)%*%KlistS%*%cc
c =  t(zz)%*%KlistS%*%zz -threshold
#cat(a,b,c,'\n')
if (a==0 & b==0){
}else{
Set[[kk]] = get_sets(a,b,c)
kk = kk + 1
}
}

if (length(Set)==0){
Set[[1]] = c(0,0)
return (Set)
}

#
# intersection sets
#

A = Set[[1]] 
if (L>1){
for (i in 2:L){
A =  intersectAB(A,Set[[i]])
}
}
return (A)
}  


#
# Adjust p-value on the selection
#

get_pvalue_i = function(i,S,Sigma,Set){

Si = S[i]
se = Sigma[i,i]
t = abs(Si/sqrt(se))

pvalueN = 2*(1-pnorm(abs(t)))

pvalueC =  CalcCondPvalue(t,Set,1)
return (list(C=pvalueC,N=pvalueN))
}
