library(Rcpp)

sourceCpp("pairmean.cpp")

my.crps = function(xobs,xpred){
  Ef1 = mean(abs(xpred - xobs) )
  Ef2 = pairmean(xpred)/2
  return(Ef1 - Ef2)
}

my.ES = function(x1obs,x1pred,m1,s1,x2obs,x2pred,m2,s2){
  xpred = cbind((x1pred - m1) / s1 , (x2pred - m2) / s2)
  xobs = c ( (x1obs - m1) / s1, (x2obs - m2) / s2 )
  
  Ef1 = mean(apply(xpred,1,function(x){ sqrt(sum((x - xobs)^2))}))
  Ef2 = pairnorm(xpred)
  return(Ef1 - Ef2)
}

# my.ES2 = function(xobs,xpred){
#   mean(apply(xpred,1,function(x){ sqrt(sum((x - xobs)^2))})) - 
#     mean(as.matrix(rdist(xpred)))/2
# }

