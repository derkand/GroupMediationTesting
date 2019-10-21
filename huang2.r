#
# Huang helper functions
#


#
# Generate bootsrap samples
#

gen_boot_samples = function(mu1,se1,mu2,se2,perm=10000){
Sigma1 = se1
Sigma2 = se2
b1 = mvrnorm(perm,mu1,Sigma1)
b1  = b1 
b2 = mvrnorm(perm,mu2,Sigma2)
b2  = b2
xx = colMeans(b1*b2)
AA = matrix(xx,nrow=length(xx),ncol=perm) 
xx  = rowSums((b1*b2 -t(AA))^2) 
return(xx)
}

#
# Get p-value using naive-linear test
#

get_pvalueSumL = function(mu1,se1,mu2,se2){
Sigma1 = se1
Sigma2 = se2

Sigma = mu1%*%Sigma2%*%mu1 + mu2%*%Sigma1%*%mu2 

T = sum(mu1%*%mu2)/sqrt(as.numeric(Sigma))
xx = 2*(1-pnorm(abs(T)))
return(xx)
}

#
# Get p-value using max-linear test
#

get_pvalueMaxL= function(mu1,se1,mu2,se2){
Sigma1 = se1
Sigma2 = se2

Sigma = NULL
V = NULL
for (i in 1:length(se1)){
Sigma = c(Sigma,max(mu2[i]^2*se1[i],mu1[i]^2*se2[i]))
}

T = sum(mu1%*%mu2)/sqrt(sum(Sigma))
xx = 2*(1-pnorm(abs(T)))

return(xx)
}









