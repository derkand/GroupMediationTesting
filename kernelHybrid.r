#
# Helper functions Kernel Test
#

#
# Get p-value for TQ = S'KS
# S - vector of beta coeffients
# Sigma = Cov(S)
# K - kernel
# Method use Liu's method
#
#H. Liu, Y. Tang, H.H. Zhang, A new chi-square approximation to the distribution of non-negative
# definite quadratic forms ..., 2009
#

get_pvalueQ = function(TQ,K,Sigma){
dd=eigen(Sigma)$values
dd[dd<0]=0

D = diag(dd)
P = eigen(Sigma)$vectors
MM = sqrt(D)%*%t(P)%*%K%*%P%*%sqrt(D)
lambda = eigen(MM)$values
lambda = Re(lambda)
lambda[lambda<0] = 0
pval = liu(TQ,lambda=lambda)

return (pval)
}

#
# Get critical values
#

calcCrit = function(K,Sigma,level){
dd=eigen(Sigma)$values


dd[dd<0]=0
D = diag(dd)
P = eigen(Sigma)$vectors
MM = sqrt(D)%*%t(P)%*%K%*%P%*%sqrt(D)
lambda = eigen(MM)$values
lambda = Re(lambda)

lambda[lambda<0] = 0



	
c1 = sum(lambda^1)
c2 = sum(lambda^2)
c3 = sum(lambda^3)
c4 = sum(lambda^4)

s1 = c3/c2^(3/2)
s2 = c4/c2^2
    #save(c1,c2,c3,c4,pj,J,file='temp')
if ((s1^2)>s2){
a = 1/(s1-sqrt(s1^2-s2))
delta = s1*a^3 - a^2
l = a^2 - 2*delta
	#cat(l,'AL',EV,'\n')
}else{
	l =c2^3/c3^2
	delta=0
	#cat(l,'A',EV,'\n')
}

ct = qchisq(1-level,df=l,nc=delta)

	mux = l+delta
	muq = c1
	sigmax = sqrt(2*(l+2*delta)) 
	sigmaq = sqrt(2*c2)

	ctt = sigmaq/sigmax*(ct - mux) + muq
return (ctt)
}



#
# Families will be tested by quadratic test 
# let Klist  - list of kernels (families)
# TQ = S' Klist[[1]] S
# SE = vector of p-estimated effects of E on M, M~E
# SY = vector of p-estimated effects of M on Y, Y~M+E 
# sigmae = covariance (SE) - block diagonal matrix
# sigmay = covariance (SY) - block diagonal matrix
# level - screening p-value 



KernelsruleHybrid = function(Klist,SE,SY,Xc,sigmae,sigmay,level){
selected_Kernel_KernelsSE = list()
selected_Kernel_KernelsSY = list()
L = length(Klist)
VarE = sigmae
VarY = sigmay
kkSE = 1
kkSY = 1
for (i in 1:L){

kkk =which((Klist[[i]])!=0)
KlistS1 = diag(Klist[[i]][kkk])
SYS = SY[kkk]
SES = SE[kkk]

TY = t(SYS)%*%KlistS1%*%(SYS)
KlistS2 = diag(SYS^2)
TE = t(SES)%*%KlistS2%*%(SES)

pQY = get_pvalueQ(TY,KlistS1,VarY[kkk,kkk])
pQE = get_pvalueQ(TE,KlistS2,VarE[kkk,kkk])


if (pQY<=0.1 & pQE<=level){
selected_Kernel_KernelsSE[[kkSE]] = Klist[[i]]
kkSE = kkSE +1 
}


KlistS2 = diag(SES^2)
TY = t(SYS)%*%KlistS2%*%(SYS)
TE = t(SES)%*%KlistS1%*%(SES)

pQY = get_pvalueQ(TY,KlistS2,VarY[kkk,kkk])
pQE = get_pvalueQ(TE,KlistS1,VarE[kkk,kkk])

if (pQY<=level & pQE<=0.1){
selected_Kernel_KernelsSY[[kkSY]] = Klist[[i]]
kkSY = kkSY +1 
}

}

if (length(selected_Kernel_KernelsSY)==0){
selected_Kernel_KernelsSY = NULL
}

if (length(selected_Kernel_KernelsSE)==0){
selected_Kernel_KernelsSE = NULL
}

return (list(KernelsSE = selected_Kernel_KernelsSE,KernelsSY = selected_Kernel_KernelsSY))
}


