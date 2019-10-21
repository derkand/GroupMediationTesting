#
# Interesection of two sets
# Checked on july 10
#
#

intersectAB  = function(A,B){
LA = length(A)
LB = length(B)
C = list()
k=1 
for (i in 1:LA){
for (j in 1:LB){

c1 = max(A[[i]][1],B[[j]][1])
c2 = min(A[[i]][2],B[[j]][2])
C[[k]] = c(c1,c2)
if (c2<c1){
C[[k]] = c(0,0)
}
k=k+1
}
}

return(C)
}


#
# Get sets from quadratic equation
#

get_sets = function(a,b,c){
Set = list()



D = b^2 - 4*a*c
#cat(D,'\n')
if (D<=0 & a<=0){
Set[[1]] = c(0,0)
return (Set)
}
if (D<0 & a>0){
Set[[1]] = c(-Inf,+Inf)
return (Set)
}

a1 = (-b + sqrt(D))/(2*a)
a2 = (-b - sqrt(D))/(2*a)

b2 = max(a1,a2)
b1 = min(a1,a2)
a2 = b2
a1 = b1 

if (D>0 & a>0){
Set[[1]] = c(-Inf,a1)
Set[[2]] = c(a2,+Inf)
return (Set)
}

if (D>0 & a<0){
Set[[1]] = c(a1,a2)
return (Set)
}

}
