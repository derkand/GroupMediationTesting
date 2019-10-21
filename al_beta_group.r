#
# This is approximation to quadratic test
# No theoretical justification, remove from Huang test
#
#



library('CompQuadForm')

get_pval_ALBE = function(al,beta,SigmaL,SigmaB){

SigmaX1  = diag(beta)%*%SigmaL%*%diag(beta)
SigmaX2  = diag(al)%*%SigmaB%*%diag(al)

Sigma = SigmaX1 + SigmaX2

lambda = eigen(Sigma)$values
qM = sum((al*beta)^2)

pvalQE = liu(qM, lambda)

return (pvalQE)
}

