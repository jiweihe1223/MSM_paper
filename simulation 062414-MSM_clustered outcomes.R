#simulation for the case with treatment interferece among cluster
library(Rlab)
library(nnet)
library(MASS)
library(geepack)
rm(list=ls(all=TRUE))

##################################
## function for simulating data ##
##################################

simulation<-function(parameters) { ## parameters contain K, param##
  
  param<-parameters[[2]] ## causal parameters
  K<-parameters[[1]] ## number of time points
  id<-1: N
  
  ## simulate time 1 (baseline) data ##
  
  ## 1. simulate baseline covariates L0
  sigma<-diag(3)
  sigma[1, 2]<-sigma[2, 1]<-0.3
  sigma[1, 3]<-sigma[3, 1]<-sigma[2, 3]<-sigma[3, 2]<-0.2
  
  l0<-matrix(0, N, 3)
  for (i in 1: N) {
    l0[i, ]<-mvrnorm(1, mu=numeric(3), Sigma=sigma)
  }
  
  ## 2. simulate 1st treatment A0 from multinomial logistic
  p01<-exp(0+0.5*l0[, 1]-0.5*l0[, 3])
  p10<-exp(0+0.5*l0[, 2]-0.5*l0[, 3])
  p11<-exp(0+0.5*l0[, 1]+0.5*l0[, 2]-l0[, 3])
  total<-1+p01+p10+p11
  
  probs<-cbind(1, p01, p10, p11)/total
  
  cat<-numeric(N) ##multinomial variable for joint treatment
  for (i in 1: N) {
    cat[i]<-sample(1: 4, 1, probs[i, ], replace=F)
  }
  
  ## treatment a01 for eye1 and a02 for eye2
  a01<-numeric(N)
  a01[cat==2 | cat==4]<-1
  a02<-numeric(N)
  a02[cat==3 | cat==4]<-1
  
  ## 3. observed outcome at time 1 ##
  
  q0<-1*(l0[,1:2]) - 0.5*cbind(l0[,3],l0[,3])
  
  if(parameters[[3]]=="no interference"){
    treat1<-a01
    treat2<-a02
    treat<-c(treat1, treat2)
  }else if(parameters[[3]]=="interference"){
    treat1<-cbind(a01, a02, a01*a02)
    treat2<-cbind(a02, a01, a01*a02)
    treat<-rbind(treat1, treat2)
  }
  
  sigma<-diag(2)
  sigma[1, 2]<-sigma[2, 1]<-0.4
  epsilon<-mvrnorm(N, mu=numeric(2), Sigma=sigma)
  
  Y1<-cbind(1, treat1)%*%param + q0[,1]+epsilon[,1]
  Y2<-cbind(1, treat2)%*%param + q0[,2]+epsilon[,2]
  Y<-c(Y1, Y2)
  
  unit<-c(rep(1, N), rep(2, N))
  id2<-c(id, id)
  data<-cbind(Y, treat, id2, unit, 0)
  
  data_total<-data
  
  ## t accumulates data over time; t0 refers to data from previous time ##
  t<-t0<-cbind(l0, a01, a02, Y1, Y2, id, 0)
  
  ## simulate data from time 2 to K ##
  for (j in 2: K) {
    
    ## 1. simulate L1 #
    sigma<-diag(3)
    sigma[1, 2]<-sigma[2, 1]<-0.3
    sigma[1, 3]<-sigma[3, 1]<-sigma[2, 3]<-sigma[3, 2]<-0.2
    
    l1<-matrix(0, N, 3)
    for (i in 1: N) {
      mu12<-0.3*t0[i, 1: 2]+0.2*t0[i, 4:5]+0.2*t0[i, 6:7]
      mu3<-0.3*t0[i, 3] + 0.1*t0[i, 4]+ 0.1*t0[i, 5]+0.1*t0[i, 6]+0.1*t0[i, 7]
      mu<-c(mu12, mu3)
      l1[i, ]<-mvrnorm(1, mu=mu, Sigma=sigma)
    }
    
    ## 2. effect of l1 on outcome minus its expected value; q0 sum over time ##
    q11<-1*(l1[, 1] - 0.3*t0[, 1]-0.2*t0[, 4]-0.2*t0[, 6]) - 0.5*(l1[, 3] - 0.3*t0[, 3] - 0.1*t0[, 4] - 0.1*t0[, 5]-0.1*t0[,6]-0.1*t0[,7])
    q12<-1*(l1[, 2] - 0.3*t0[, 2]-0.2*t0[, 5]-0.2*t0[, 7]) - 0.5*(l1[, 3] - 0.3*t0[, 3] - 0.1*t0[, 4] - 0.1*t0[, 5]-0.1*t0[,6]-0.1*t0[,7])
    q1<-cbind(q11,q12)
    q0<-q0+q1
    
    ## 3. simulate A1 from multinomial logistic ##
    p01<-exp(0+0.5*l1[, 1] -0.5*l1[, 3] + 0.2*t0[, 4]+0.2*t0[,6])
    p10<-exp(0+0.5*l1[, 2] -0.5*l1[, 3] + 0.2*t0[, 5]+0.2*t0[,7])
    p11<-exp(0+0.5*l1[, 1]+0.5*l1[, 2]-l1[, 3] + 0.2*t0[, 4]+ 0.2*t0[, 5]+0.1*t0[,6]+0.1*t0[,7])
    total<-1+p01+p10+p11
    probs<-cbind(1, p01, p10, p11)/total
    
    cat<-numeric(N)
    
    for (i in 1: N) {
      cat[i]<-sample(1: 4, 1, probs[i, ], replace=F)
    }
    
    a11<-numeric(N)
    a11[cat==2 | cat==4]<-1
    
    a12<-numeric(N)
    a12[cat==3 | cat==4]<-1
    
    # 4. observed outcome
    if(parameters[[3]]=="no interference"){
      treat1<-treat1+a11
      treat2<-treat2+a12
      treat<-c(treat1, treat2)
    }else if(parameters[[3]]=="interference"){
      treat1<-treat1+cbind(a11, a12, a11*a12)
      treat2<-treat2+cbind(a12, a11, a11*a12)
      treat<-rbind(treat1, treat2)
    }
    
    sigma<-diag(2)
    sigma[1, 2]<-sigma[2, 1]<-0.4
    epsilon<-mvrnorm(N, mu=numeric(2), Sigma=sigma)
    
    Y1<-cbind(1, treat1)%*%param + q0[,1]+epsilon[,1]
    Y2<-cbind(1, treat2)%*%param + q0[,2]+epsilon[,2]
    Y<-c(Y1, Y2)
    
    t0<-cbind(l1, a11, a12, Y1, Y2, id, j-1)
    t<-rbind(t, t0)
    
    unit<-c(rep(1, N), rep(2, N))
    id2<-c(id, id)
    data<-cbind(Y, treat, id2, unit, j-1)
    data_total<-rbind(data_total, data)  
  }
  
  ## output t (covariates and treatments) and observed outcomes
  t<-data.frame(t)
  data_total<-data.frame(data_total)
  return(list(t, data_total))
}

###############################################################
## function for calculating inverse probability weight (IPW) ##
###############################################################
ipw_with_unit<-function(history){
  
  ## cluster-specific weight ##
  history$A01<-c(numeric(N), history$A1[1: (N*(K-1))])
  history$A02<-c(numeric(N), history$A2[1: (N*(K-1))])
  history$Y01<-c(numeric(N), history$Y1[1: (N*(K-1))])
  history$Y02<-c(numeric(N), history$Y2[1: (N*(K-1))])
  
  history$treatment.cat<-NULL
  history$treatment.cat[history$A1==0 & history$A2==0]<-1
  history$treatment.cat[history$A1==1 & history$A2==0]<-2
  history$treatment.cat[history$A1==0 & history$A2==1]<-3
  history$treatment.cat[history$A1==1 & history$A2==1]<-4
  
  predictor<-history[, c(1:3,10:13)]
  m_treat<-multinom(factor(history$treatment.cat)~L1+L2+L+A01+A02+Y01+Y02, data=predictor)
  history$pred<-predict(m_treat,predictor,"probs")
  history$predicted<-NULL
  history$predicted[history$A1==0 & history$A2==0]<-history$pred[,1][history$A1==0 & history$A2==0]
  history$predicted[history$A1==1 & history$A2==0]<-history$pred[,2][history$A1==1 & history$A2==0]
  history$predicted[history$A1==0 & history$A2==1]<-history$pred[,3][history$A1==0 & history$A2==1]
  history$predicted[history$A1==1 & history$A2==1]<-history$pred[,4][history$A1==1 & history$A2==1]
  
  predictor<-history[, c(10,11)]
  m_treat2<-multinom(factor(history$treatment.cat)~A01+A02, data=predictor)
  history$pred2<-predict(m_treat2,predictor,"probs")
  history$predicted2<-NULL
  history$predicted2[history$A1==0 & history$A2==0]<-history$pred2[,1][history$A1==0 & history$A2==0]
  history$predicted2[history$A1==1 & history$A2==0]<-history$pred2[,2][history$A1==1 & history$A2==0]
  history$predicted2[history$A1==0 & history$A2==1]<-history$pred2[,3][history$A1==0 & history$A2==1]
  history$predicted2[history$A1==1 & history$A2==1]<-history$pred2[,4][history$A1==1 & history$A2==1]
  
  pred<-history$predicted[1:N]
  pred2<-history$predicted2[1:N]
  history$ipw[1:N]<-pred2/pred
  history$ipw_unstable[1:N]<-1/pred
  
  for (j in 2: K) {
    pred<-pred*history$predicted[(N*(j-1)+1): (N*j)]
    pred2<-pred2*history$predicted2[(N*(j-1)+1): (N*j)]
    
    history$ipw[(N*(j-1)+1): (N*j)]<-pred2/pred
    history$ipw_unstable[(N*(j-1)+1): (N*j)]<-1/pred
  }
  
  cluster_specific<-data.frame(history$id, history$ipw_unstable, history$ipw, history$Time)
  colnames(cluster_specific)<-c("id", "ipw_unstable", "ipw", "Time")
  
  ## unit-specific weight ##
  history$p1<-history$pred[,2]+history$pred[,4]
  history$p2<-history$pred[,3]+history$pred[,4]
  history$predicted.A1<-history$p1^history$A1*(1-history$p1)^(1-history$A1)
  history$predicted.A2<-history$p2^history$A2*(1-history$p2)^(1-history$A2)
  
  history$p1.2<-history$pred2[,2]+history$pred2[,4]
  history$p2.2<-history$pred2[,3]+history$pred2[,4]
  history$predicted2.A1<-history$p1.2^history$A1*(1-history$p1.2)^(1-history$A1)
  history$predicted2.A2<-history$p2.2^history$A2*(1-history$p2.2)^(1-history$A2)
  
  pred_unit1<-history$predicted.A1[1:N]
  pred_unit2<-history$predicted.A2[1:N]
  pred2_unit1<-history$predicted2.A1[1:N]
  pred2_unit2<-history$predicted2.A2[1:N]
  
  history$ipw.unstable_unit1[1:N]<-1/pred_unit1
  history$ipw.unstable_unit2[1:N]<-1/pred_unit2
  
  history$ipw_unit1[1:N]<-pred2_unit1/pred_unit1
  history$ipw_unit2[1:N]<-pred2_unit2/pred_unit2
  
  for (j in 2: K) {
    pred_unit1<-pred_unit1*history$predicted.A1[(N*(j-1)+1): (N*j)]
    pred_unit2<-pred_unit2*history$predicted.A2[(N*(j-1)+1): (N*j)]
    
    pred2_unit1<-pred2_unit1*history$predicted2.A1[(N*(j-1)+1): (N*j)]
    pred2_unit2<-pred2_unit2*history$predicted2.A2[(N*(j-1)+1): (N*j)]
    
    history$ipw.unstable_unit1[(N*(j-1)+1): (N*j)]<-1/pred_unit1
    history$ipw.unstable_unit2[(N*(j-1)+1): (N*j)]<-1/pred_unit2
    
    history$ipw_unit1[(N*(j-1)+1): (N*j)]<-pred2_unit1/pred_unit1
    history$ipw_unit2[(N*(j-1)+1): (N*j)]<-pred2_unit2/pred_unit2
  }
  
  sub1<-data.frame(history$id, 1, history$ipw.unstable_unit1, history$ipw_unit1, history$Time)
  sub2<-data.frame(history$id, 2, history$ipw.unstable_unit2, history$ipw_unit2, history$Time)
  colnames(sub1)<-colnames(sub2)<-c("id", "unit", "ipw.unstable_unit", "ipw_unit", "Time")
  unit_specific<-rbind(sub1, sub2)
  
  ipw<-merge(unit_specific, cluster_specific, by=c("Time", "id"))
  return(ipw)
}

###############################
# start of simulation (M=1000)#
###############################
K<-3
param<-c(0, 2, 1, 0.5)
parameters<-list(K, param)
N<-1000

set.seed(122384)
M<-1000

estimates.vec<-NULL
se.vec<-NULL
data.vec<-list(NULL)

for (k in 1: M){
  output<-simulation(parameters)
  data<-output[[2]]
  colnames(data)<-c("Y", "A",  "A_prime", "AA", "id", "unit", "Time")
  history<-output[[1]]
  colnames(history)<-c("L1", "L2", "L", "A1", "A2", "Y1", "Y2", "id", "Time")
  
  weight<-ipw_with_unit(history)
  weight_cluster.specific<-weight[weight$Time==K-1, ]
  data2<-merge(data, weight, by=c("id", "Time", "unit"))
  data2<-merge(data2, weight_cluster.specific, by=c("id", "unit"))
  data2<-data2[,-12]
  colnames(data2)<-c("id", "unit", "Time",  "Y", "A", "A_prime", "AA", "ipw.unstable_unit", "ipw_unit", "ipw_unstable", "ipw", "ipw.unstable_unit2", "ipw_unit2", "ipw_unstable2", "ipw2")
  data2<-data2[order(data2$id, data2$Time, data$unit), ] 
  
  ## no weight ##
  m1<-geeglm(Y~A+A_prime+AA, id=id, corstr="independence", data=data2)
  summary(m1)
  
  estimates<-m1$coef
  se<-summary(m1)$coefficients[,2]
  
  ## with unstablized IPW (cluster-specific) ##
  
  m2.1<-geeglm(Y~A+A_prime+AA, id=id, weights=ipw_unstable, corstr="independence", data=data2)
  summary(m2.1)
  
  ## with stablized IPW (cluster-specific) ##
  
  m2.2<-geeglm(Y~A+A_prime+AA, id=id, weights=ipw, corstr="independence", data=data2)
  summary(m2.2)
  
  estimates<-c(estimates, m2.1$coef, m2.2$coef)
  se<-c(se, summary(m2.1)$coefficients[,2], summary(m2.2)$coefficients[,2])
  
  ## repeat the above with unit-specific IPW ##
  
  m6.1<-geeglm(Y~A+A_prime+AA, id=id, weights=ipw.unstable_unit, corstr="independence", data=data2)
  summary(m6.1)
  
  m6.2<-geeglm(Y~A+A_prime+AA, id=id, weights=ipw_unit, corstr="independence", data=data2)
  summary(m6.2)
  
  estimates<-c(estimates, m6.1$coef, m6.2$coef)
  se<-c(se, summary(m6.1)$coefficients[,2], summary(m6.2)$coefficients[,2])
  
  estimates.vec<-rbind(estimates.vec, estimates)
  se.vec<-rbind(se.vec, se)
  
  #print(k)
  #print(estimates)
}

#######################
## summarize results ##
#######################
mean<-apply(estimates.vec, 2, mean)-param
mean<-matrix(mean, nrow=4)
sd<-apply(estimates.vec, 2, sd)
sd<-matrix(sd, nrow=4)
sdmean<-apply(se.vec, 2, mean)
sdmean<-matrix(sdmean, nrow=4)

summ<-function(k){
  result<-NULL
  for(i in k){
    result_i<-cbind(mean[,i], sd[,i], sdmean[,i])
    colnames(result_i)<-c("mean", "sd", "sdmean")
    result<-cbind(result, result_i)
  }
  return(result)
}

summ<-function(k){
  result<-NULL
  for(i in k){
    result_i<-cbind(mean[,i], sd[,i])
    colnames(result_i)<-c("mean", "sd")
    result<-cbind(result, result_i)
  }
  return(result)
}

CI_lower=estimates.vec - 1.96*se.vec
CI_upper=estimates.vec + 1.96*se.vec
beta_true<-rep(c(0,2,1,0.5), 5)

CI_ind<-matrix(NA, M, 20)
for (i in 1: M) {
  CI_ind[i, ]=ifelse( (CI_lower[i, ]<beta_true) & (CI_upper[i, ]>beta_true), 1, 0)
}

coverage=round(apply(CI_ind, 2, mean), 3)
coverage<-matrix(coverage, nrow=4)

summ<-function(k){
  result<-NULL
  for(i in k){
    result_i<-cbind(mean[,i], sd[,i], sdmean[,i], coverage[,i])
    colnames(result_i)<-c("mean", "sd", "sdmean", "95% coverage")
    result<-cbind(result, result_i)
  }
  return(result)
}

GEE<-summ(1)
ipw_independence<-summ(2:3)
ipw_independence_unit<-summ(4:5)

GEE
ipw_independence
ipw_independence_unit

data.vec[[1]]<-NULL
data_combined<-matrix(NA, 2*N*M, dim(data.vec[[1]])[2])

for (i in 1: M){
  data_combined[((i-1)*(2*N*K)+1) : (i*(2*N*K)), ]<-data.vec[[i]]
}


