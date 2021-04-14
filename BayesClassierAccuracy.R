t1 = Sys.time()
###decision function
p1 = 0.4
p2 = 0.6
decision_new = function(X, p, Mean, Var, rho){
  x1 = X[1]
  x2 = X[2]
  m1 = Mean[1]
  m2 = Mean[2]
  Var1 = Var[1]
  Var2 = Var[2]
  d = log(p) - (1/2)*log(Var1*Var2*(1-rho^2))- (1/2)*(1/(1-rho^2))*
    (((x1-m1)^2)/Var1-2*rho*((x1-m1)/sqrt(Var1))*((x2-m2)/sqrt(Var2))+((x2-m2)^2)/Var2)
  return(d)
}
library(mvtnorm)
#population mean vector, covariance matrix
m1 = matrix(c(80, 120), nrow =  2, ncol =  1)
m2 = matrix(c(140, 150), 2, 1)
cov1 = matrix(c(1225, -525, -525, 400), 2, 2, byrow = T)
cov2 = matrix(c(900, 390, 390, 400), 2, 2, byrow = T)
#population
set.seed(947)
pop1 = rmvnorm(80000, m1, cov1)
v1 = c(cov1[1, 1], cov1[2, 2])
r1 = cov1[1, 2] / sqrt(cov1[1, 1]*cov1[2, 2])
pop2 = rmvnorm(120000, m2, cov2)
v2 = c(cov2[1, 1], cov2[2, 2])
r2 = cov2[1, 2] / sqrt(cov2[1, 1]*cov2[2, 2])
classify_pop = function(X){
  decision_new(X, p1, m1, v1, r1) - decision_new(X, p2, m2, v2, r2)
}
clpop1 = apply(pop1, 1, classify_pop)
clpop2 = apply(pop2, 1, classify_pop)
papop1 = sum(as.numeric(clpop1 > 0)) / length(clpop1)
papop2 = sum(as.numeric(clpop2 < 0)) / length(clpop2)
papop = c(papop1, papop2)
#sample mean vector, covariance matrix
n = 100
mbar = function(X){
  apply(X, 2, mean)
}
varbar = function(X){
  apply(X, 2, var)
}
rhobar = function(X){
  cor(X[,1], X[,2])
}
c1 = replicate(n, list(rmvnorm(80, m1, cov1)))
c2 = replicate(n, list(rmvnorm(120, m2, cov2)))
#####clean#####
rm(m1, m2, cov1, cov2)
rm(v1, v2, r1, r2)
#bootstrap sample
bootstrap = function(X){
  n = length(X)
  sam = vector()
  sam = data.frame(X[1:(n/2)], X[((n/2+1):n)])
  #for(i in 1:(n/2)){
  #  m_temp = data.frame(X[i], X[((n/2)+i)])
  #  sam = rbind(sam, m_temp)
  #}
  bs = sam[sample(nrow(sam), (n/2), replace = T),]
  return(bs)
}
nbootstrap = function(X){
  replicate(n, list(bootstrap(X)))
} 
bootc1 = lapply(c1, nbootstrap)
bootc2 = lapply(c2, nbootstrap)
#bootstrap sample parameters
boombar = function(X){
  lapply(X, mbar)
} 
boovbar = function(X){
  lapply(X, varbar)
}
boorho = function(X){
  lapply(X, rhobar)
}
c1boombar = lapply(bootc1, boombar)
c2boombar = lapply(bootc2, boombar)
c1boovbar = lapply(bootc1, boovbar)
c2boovbar = lapply(bootc2, boovbar)
c1boorho = lapply(bootc1, boorho)
c2boorho = lapply(bootc2, boorho)

interval_first = vector(mode = "list")
for (i in 1:n){
  classify1 = vector(mode = "list")
  classify2 = vector(mode = "list")
  for(j in 1:n){
    classify_new = function(X){
      decision_new(X, p1, c1boombar[[i]][[j]], c1boovbar[[i]][[j]], c1boorho[[i]][[j]]) - decision_new(X, p2, c2boombar[[i]][[j]], c2boovbar[[i]][[j]], c2boorho[[i]][[j]])
    }
    a = apply(bootc1[[i]][[j]], 1, classify_new)
    classify1 = append(classify1, list(a))
    b = apply(bootc2[[i]][[j]], 1, classify_new)
    classify2 = append(classify2, list(b))
  }
  pa_1 = function(X){
    sum(as.numeric(unlist(X) > 0)) / length(unlist(X))
  }
  pa_2 = function(X){
    sum(as.numeric(unlist(X) < 0)) / length(unlist(X))
  }
  pa_interval = function(X, Y){
    pa1 = sapply(X, pa_1)
    pa2 = sapply(Y, pa_2)
    interval1 = quantile(unlist(pa1), c(0.025, 0.975))
    interval2 = quantile(unlist(pa2), c(0.025, 0.975))
    list(interval1, interval2)
  }
  i_temp = pa_interval(classify1, classify2)
  interval_first = append(interval_first, list(i_temp))
}
#####clean#####
rm(bootc1, bootc2)
rm(c1boombar, c2boombar, c1boovbar, c2boovbar, c1boorho, c2boorho)
rm(classify1, classify2)

#bootstrap sample
bootc1 = lapply(c1, nbootstrap)
bootc2 = lapply(c2, nbootstrap)
#####clean#####
rm(c1, c2)
#bootstrap sample parameters
c1boombar = lapply(bootc1, boombar)
c2boombar = lapply(bootc2, boombar)
c1boovbar = lapply(bootc1, boovbar)
c2boovbar = lapply(bootc2, boovbar)
c1boorho = lapply(bootc1, boorho)
c2boorho = lapply(bootc2, boorho)

interval_second = vector(mode = "list")
for (i in 1:n){
  classify1 = vector(mode = "list")
  classify2 = vector(mode = "list")
  for(j in 1:n){
    classify_new = function(X){
      decision_new(X, p1, c1boombar[[i]][[j]], c1boovbar[[i]][[j]], c1boorho[[i]][[j]]) - decision_new(X, p2, c2boombar[[i]][[j]], c2boovbar[[i]][[j]], c2boorho[[i]][[j]])
    }
    a = apply(bootc1[[i]][[j]], 1, classify_new)
    classify1 = append(classify1, list(a))
    b = apply(bootc2[[i]][[j]], 1, classify_new)
    classify2 = append(classify2, list(b))
  }
  pa_1 = function(X){
    sum(as.numeric(unlist(X) > 0)) / length(unlist(X))
  }
  pa_2 = function(X){
    sum(as.numeric(unlist(X) < 0)) / length(unlist(X))
  }
  pa_interval = function(X, Y){
    pa1 = sapply(X, pa_1)
    pa2 = sapply(Y, pa_2)
    interval1 = quantile(unlist(pa1), c(0.025, 0.975))
    interval2 = quantile(unlist(pa2), c(0.025, 0.975))
    list(interval1, interval2)
  }
  i_temp = pa_interval(classify1, classify2)
  interval_second = append(interval_second, list(i_temp))
}
#####clean#####
rm(bootc1, bootc2)
rm(c1boombar, c2boombar, c1boovbar, c2boovbar, c1boorho, c2boorho)
rm(classify1, classify2)

#test 95% confidence interval
pa1interval = function(X){
  if(papop[1] >= X[[1]][1] && papop[1] <= X[[1]][2]){
    return(1)
  }else{
    return(0)
  }
}
pa2interval = function(X){
  if(papop[2] >= X[[2]][1] && papop[2] <= X[[2]][2]){
    return(1)
  }else{
    return(0)
  }
}
pa195per = (sum(sapply(interval_first, pa1interval))+sum(sapply(interval_second, pa1interval))) / (n*2)
pa295per = (sum(sapply(interval_first, pa2interval))+sum(sapply(interval_second, pa1interval))) / (n*2)
pa195per;pa295per
t4 = Sys.time()

par(mfrow = c(1, 2))
plot(c(papop[1], papop[1]), c(1,500),type = "l", col = "red", lwd = 3, xlim = c(0.8, 1.08), ylim = c(1,500), xlab = "Pruducer accuracy 1", ylab = "", main = "95% bootstrap interval(PA1)")
for(i in 1:500){
  lines(c(interval_first[[i]][[1]][1], interval_first[[i]][[1]][2]), c(i, i), col = rgb(190, 190 ,190, 200, maxColorValue = 255))
}
plot(c(papop[2], papop[2]), c(1,500),type = "l", col = "red", lwd = 3, xlim = c(0.8, 1.024), ylim = c(1,500), xlab = "Pruducer accuracy 2", ylab = "", main = "95% bootstrap interval(PA2)")
for(i in 1:500){
  lines(c(interval_first[[i]][[2]][1], interval_first[[i]][[2]][2]), c(i, i), col = rgb(190, 190 ,190, 200, maxColorValue = 255))
}
