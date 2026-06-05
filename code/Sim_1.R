####The simulation setting is the same in
####"Differentially Private False Discovery Rate Control"
####by Cynthia Dwork, Weijie J. Su and Li Zhang
library(VGAM) #generate Laplace distribution
library(parallel)
library(doParallel)
library(foreach)

one_iteration<-function(m = 1e5,
                        m1 = 500,
                        tr = 100,
                        q = 0.1,
                        eta = 1e-4,
                        delta = 0.001,
                        mu = 4,
                        epsilon = 0.5,
                        conserv = FALSE){
  nu <- 0.5*q/m
  
  #generate p values
  xi <- rnorm(n = tr, mean = 0, sd = 1)
  p <- unlist(lapply(xi, function(x) pnorm(q = x-mu,mean = 0,sd = 1)))
  #if (conserv) { p <- c(p,rbeta(n = m-tr, shape1 = 2, shape2 = 2))} else {
  if (conserv) { p <- c(p,runif(n = m-tr, min = 0, max = 1)^{1/3})} else {
  p <- c(p,runif(n = m-tr, min = 0, max = 1))}

  ########DP-BH
  lambda <- eta*sqrt(10*m1*log(1/delta))/(2*epsilon)
  theta <- log(pmax(p,nu))
  
  ##pealing arg min
  index<-c()
  new_theta<-c()
  for (i in 1:m1) {
    #add noisy
    noisy_theta <- theta + rlaplace(n = m, location = 0, scale = lambda)
    
    #select
    index[i] <- which.min(noisy_theta)
    new_theta[i] <- theta[index[i]]+rlaplace(n = 1, location = 0, scale = lambda)
    theta[index[i]] <- Inf
  }
  
  #BHq
  for (j in m1:1) {
    if (new_theta[j]<=(log(q*j/m)-2*lambda*log(6*m1/q))) {break}
  }
  
  if ((j==1)&(new_theta[1]>(log(q*1/m)-2*lambda*log(6*m1/q)))){
    dp_bh_fdp<-0
    dp_bh_power<-0
  } else{
  rej<-index[1:j]
  dp_bh_fdp <- length(which(rej>tr))/max(length(rej),1)
  dp_bh_power <- length(which(rej<=tr))/tr
  }
  
  ########DP-Bonf
  lambda_tilde <- eta*sqrt(10*m*log(1/delta))/(2*epsilon)
  theta_tilde <- log(pmax(p,nu)) + rlaplace(n = m, location = 0, scale = lambda_tilde)
  DP_Bonf_threshold <- log(q/m) - lambda_tilde*log(5*m/q)
  Bonf_rej <- which(theta_tilde <= DP_Bonf_threshold)
  dp_bonf_fdp <- length(which(Bonf_rej>tr))/max(length(Bonf_rej),1)
  dp_bonf_power <- length(which(Bonf_rej<=tr))/tr
  
  ####AdaPT
  control<-TRUE
  alpha <- q
  nullhp<-p[-c(1:tr)]
  #ini<-min(max(1-nullhp[order(nullhp,decreasing = TRUE)[floor(0.2*tr)]],p[1:tr]),0.4)
  ini<-max(1-nullhp[order(nullhp,decreasing = TRUE)[floor(0.2*tr)]])
  while (control) {
    At <- length(which(p>=(1-ini)))
    Rt <- length(which(p<=ini))
    fdpt <- (1+At)/max(Rt,1)
    
    if (fdpt<=alpha) {control <- FALSE} else {
      potent_down <- which(p<ini)
      potent_up <- which(p>(1-ini))
      seqq <- c(ini-p[potent_down],p[potent_up]+ini-1)
      ini <- ini - seqq[order(seqq)[2]]
    }
    Rt <- length(which(p<=ini))
    At <- length(which(p>=(1-ini)))
    if (Rt==0) {control<-FALSE}
    if (At==0) {control<-FALSE}
  }
  adapt_rej<-which(p<=ini)
  if (length(adapt_rej)==0) {
    adapt_True_fdr<-0
    adapt_power_est<-0} else {
      adapt_True_fdr<-length(which(adapt_rej>tr))/max(length(adapt_rej),1)
      adapt_power_est<-length(which(adapt_rej<=tr))/tr
    }
  if (fdpt>alpha) {adapt_True_fdr<-0; adapt_power_est<-0;}
  
  ####DP-AdaPT
  noisy_sd <- lambda*sqrt(2)
  
  ###peeling
  ind<-rep(TRUE,m)
  n_index<-c()
  newqp<-c()
  
  #transformation
  qp <- unlist(lapply(p, function(x) qnorm(p = x, mean = 0, sd = 1)))
  #select small (p,1-p)
  for (i in 1:m1) {
    #add noisy
    qptilde <- qp + rnorm(n = m, mean = 0, sd = noisy_sd)
    mirror_qptilde <- abs(qptilde)
    #select
    n_index[i] <- which(ind)[which.max(mirror_qptilde[ind])]
    newqp[i] <- qp[n_index[i]] + rnorm(n = 1, mean = 0, sd = noisy_sd)
    ind[n_index[i]] <- FALSE
  }
  
  #transform back
  newp <- unlist(lapply(newqp, function(x) pnorm(q = x, mean = 0, sd = 1)))
  
  control<-TRUE
  alpha <- q
  ini<-max(1-newp[order(newp,decreasing = TRUE)[floor(0.2*tr)]])
  while (control) {
    At <- length(which(newp>=(1-ini)))
    Rt <- length(which(newp<=ini))
    fdpt <- (1+At)/max(Rt,1)
    
    if (fdpt<=alpha) {control <- FALSE} else {
      potent_down <- which(newp<ini)
      potent_up <- which(newp>(1-ini))
      seqq <- c(ini-newp[potent_down],newp[potent_up]+ini-1)
      ini <- ini - seqq[order(seqq)[2]]
    }
    Rt <- length(which(newp<=ini))
    At <- length(which(newp>=(1-ini)))
    if (Rt==0) {control<-FALSE}
    if (At==0) {control<-FALSE}
  }
  new_rej<-n_index[which(newp<=ini)]
  if (length(new_rej)==0) {
    True_fdr<-0
    power_est<-0} else {
    True_fdr<-length(which(new_rej>tr))/max(length(new_rej),1)
    power_est<-length(which(new_rej<=tr))/tr
    }
  if (fdpt>alpha) {True_fdr<-0;power_est<-0;}
  return(list(dp_bh = c(dp_bh_fdp,dp_bh_power),dp_bonf = c(dp_bonf_fdp,dp_bonf_power),adapt = c(adapt_True_fdr,adapt_power_est),dp_adapt = c( True_fdr,  power_est)))
}

###Varying epsilon result
v_epsilon <- seq(0.1,1,by = 0.1)
v_epsilon_alpha<-matrix(data = NA,nrow = 4,ncol = length(v_epsilon))
v_epsilon_power<-matrix(data = NA,nrow = 4,ncol = length(v_epsilon))

cl_size = detectCores()
cl = makeCluster(cl_size)
registerDoParallel(cl)
for (vi in 1:length(v_epsilon)){
  set.seed(19951004)
  print(vi/length(v_epsilon))
  result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("VGAM")) %dopar% one_iteration(epsilon = v_epsilon[vi],m1 = 500)
  result1 <- cbind(Reduce("+",result[1,])/100,Reduce("+",result[2,])/100,Reduce("+",result[3,])/100,Reduce("+",result[4,])/100)
  v_epsilon_alpha[,vi] <- result1[1,]
  v_epsilon_power[,vi] <- result1[2,]
}
stopCluster(cl)

write.csv(v_epsilon_alpha,"./result/v_epsilon_alpha.csv")
write.csv(v_epsilon_power,"./result/v_epsilon_power.csv")

###Varying eta result
v_eta <- 10^(seq(-5,-3,by = 0.25))
v_eta_alpha<-matrix(data = NA,nrow = 4,ncol = length(v_eta))
v_eta_power<-matrix(data = NA,nrow = 4,ncol = length(v_eta))

cl_size = detectCores()
cl = makeCluster(cl_size)
registerDoParallel(cl)
for (vi in 1:length(v_eta)){
  set.seed(951004)
  print(vi/length(v_eta))
  result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("VGAM")) %dopar% one_iteration(eta = v_eta[vi],m1 = 500)
  result1 <- cbind(Reduce("+",result[1,])/100,Reduce("+",result[2,])/100,Reduce("+",result[3,])/100,Reduce("+",result[4,])/100)
  v_eta_alpha[,vi] <- result1[1,]
  v_eta_power[,vi] <- result1[2,]
}
stopCluster(cl)

write.csv(v_eta_alpha,"./result/v_eta_alpha.csv")
write.csv(v_eta_power,"./result/v_eta_power.csv")
sim1_plot(v_eta_alpha[1,],v_eta_alpha[2,],v_eta_alpha[3,],v_eta_alpha[4,],log10(v_eta),c("log$\\eta$","FDR"))
sim1_plot(v_eta_power[1,],v_eta_power[2,],v_eta_power[3,],v_eta_power[4,],log10(v_eta),c("log$\\eta$","Power"))

###Varying mu result
v_beta <- seq(1,6,by = 1)
v_beta_alpha<-matrix(data = NA,nrow = 4,ncol = length(v_beta))
v_beta_power<-matrix(data = NA,nrow = 4,ncol = length(v_beta))
cl_size = detectCores()
cl = makeCluster(cl_size)
registerDoParallel(cl)
for (vi in 1:length(v_beta)){
  set.seed(951004)
  print(vi/length(v_beta))
  result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("VGAM")) %dopar% one_iteration(mu = v_beta[vi],m1 = 500)
  result1 <- cbind(Reduce("+",result[1,])/100,Reduce("+",result[2,])/100,Reduce("+",result[3,])/100,Reduce("+",result[4,])/100)
  v_beta_alpha[,vi] <- result1[1,]
  v_beta_power[,vi] <- result1[2,]
}
stopCluster(cl)

write.csv(v_beta_alpha,"./result/v_beta_alpha.csv")
write.csv(v_beta_power,"./result/v_beta_power.csv")
sim1_plot(v_beta_alpha[1,],v_beta_alpha[2,],v_beta_alpha[3,],v_beta_alpha[4,],v_beta,c("$\\beta$","FDR"))
sim1_plot(v_beta_power[1,],v_beta_power[2,],v_beta_power[3,],v_beta_power[4,],v_beta,c("$\\beta$","Power"))

###Varying tr result
v_tr <- seq(50,300,by = 50)
v_tr_alpha<-matrix(data = NA,nrow = 4,ncol = length(v_tr))
v_tr_power<-matrix(data = NA,nrow = 4,ncol = length(v_tr))
cl_size = detectCores()
cl = makeCluster(cl_size)
registerDoParallel(cl)
for (vi in 1:length(v_tr)){
  set.seed(951004)
  print(vi/length(v_tr))
  result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("VGAM")) %dopar% one_iteration(tr = v_tr[vi],m1 = 500)
  result1 <- cbind(Reduce("+",result[1,])/100,Reduce("+",result[2,])/100,Reduce("+",result[3,])/100,Reduce("+",result[4,])/100)
  v_tr_alpha[,vi] <- result1[1,]
  v_tr_power[,vi] <- result1[2,]
}
stopCluster(cl)

write.csv(v_tr_alpha,"./result/v_tr_alpha.csv")
write.csv(v_tr_power,"./result/v_tr_power.csv")
sim1_plot(v_tr_alpha[1,],v_tr_alpha[2,],v_tr_alpha[3,],v_tr_alpha[4,],v_tr,c("number of true effects","FDR"))
sim1_plot(v_tr_power[1,],v_tr_power[2,],v_tr_power[3,],v_tr_power[4,],v_tr,c("number of true effects","Power"))


##################conserve p
###Varying epsilon result

###Varying epsilon result
v_epsilon <- seq(0.1,1,by = 0.1)
v_epsilon_alpha<-matrix(data = NA,nrow = 4,ncol = length(v_epsilon))
v_epsilon_power<-matrix(data = NA,nrow = 4,ncol = length(v_epsilon))

cl_size = detectCores()
cl = makeCluster(cl_size)
registerDoParallel(cl)
for (vi in 1:length(v_epsilon)){
  print(vi/length(v_epsilon))
  set.seed(951004)
  result<-NULL
  result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("VGAM")) %dopar% one_iteration(epsilon = v_epsilon[vi],conserv = TRUE,m1 = 500)
  result1 <- cbind(Reduce("+",result[1,])/100,Reduce("+",result[2,])/100,Reduce("+",result[3,])/100,Reduce("+",result[4,])/100)
  v_epsilon_alpha[,vi] <- result1[1,]
  v_epsilon_power[,vi] <- result1[2,]
}
stopCluster(cl)

write.csv(v_epsilon_alpha,"./result/v_epsilon_alpha_con.csv")
write.csv(v_epsilon_power,"./result/v_epsilon_power_con.csv")
sim1_plot(v_epsilon_alpha[1,],v_epsilon_alpha[2,],v_epsilon_alpha[3,],v_epsilon_alpha[4,],v_epsilon,c("$\\epsilon$","FDR"))
sim1_plot(v_epsilon_power[1,],v_epsilon_power[2,],v_epsilon_power[3,],v_epsilon_power[4,],v_epsilon,c("$\\epsilon$","Power"))

###Varying eta result
v_eta <- 10^(seq(-5,-3,by = 0.25))
v_eta_alpha<-matrix(data = NA,nrow = 4,ncol = length(v_eta))
v_eta_power<-matrix(data = NA,nrow = 4,ncol = length(v_eta))

cl_size = detectCores()
cl = makeCluster(cl_size)
registerDoParallel(cl)
for (vi in 1:length(v_eta)){
  set.seed(951004)
  print(vi/length(v_eta))
  result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("VGAM")) %dopar% one_iteration(eta = v_eta[vi],m1 = 500,conserv = TRUE)
  result1 <- cbind(Reduce("+",result[1,])/100,Reduce("+",result[2,])/100,Reduce("+",result[3,])/100,Reduce("+",result[4,])/100)
  v_eta_alpha[,vi] <- result1[1,]
  v_eta_power[,vi] <- result1[2,]
}
stopCluster(cl)

write.csv(v_eta_alpha,"./result/v_eta_alpha_con.csv")
write.csv(v_eta_power,"./result/v_eta_power_con.csv")
sim1_plot(v_eta_alpha[1,1:8],v_eta_alpha[2,1:8],v_eta_alpha[3,1:8],v_eta_alpha[4,1:8],log10(v_eta)[1:8],c("log$\\eta$","FDR"))
sim1_plot(v_eta_power[1,1:8],v_eta_power[2,1:8],v_eta_power[3,1:8],v_eta_power[4,1:8],log10(v_eta)[1:8],c("log$\\eta$","Power"))

###Varying mu result
v_beta <- seq(2,7,by = 1)
v_beta_alpha<-matrix(data = NA,nrow = 4,ncol = length(v_beta))
v_beta_power<-matrix(data = NA,nrow = 4,ncol = length(v_beta))
cl_size = detectCores()
cl = makeCluster(cl_size)
registerDoParallel(cl)
for (vi in 1:length(v_beta)){
  set.seed(951004)
  print(vi/length(v_beta))
  result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("VGAM")) %dopar% one_iteration(mu = v_beta[vi],conserv = TRUE,m1 = 500)
  result1 <- cbind(Reduce("+",result[1,])/100,Reduce("+",result[2,])/100,Reduce("+",result[3,])/100,Reduce("+",result[4,])/100)
  v_beta_alpha[,vi] <- result1[1,]
  v_beta_power[,vi] <- result1[2,]
}
stopCluster(cl)

write.csv(v_beta_alpha,"./result/v_beta_alpha_con.csv")
write.csv(v_beta_power,"./result/v_beta_power_con.csv")
sim1_plot(v_beta_alpha[1,],v_beta_alpha[2,],v_beta_alpha[3,],v_beta_alpha[4,],v_beta,c("$\\beta$","FDR"))
sim1_plot(v_beta_power[1,],v_beta_power[2,],v_beta_power[3,],v_beta_power[4,],v_beta,c("$\\beta$","Power"))

###Varying tr result
v_tr <- seq(50,300,by = 50)
v_tr_alpha<-matrix(data = NA,nrow = 4,ncol = length(v_tr))
v_tr_power<-matrix(data = NA,nrow = 4,ncol = length(v_tr))
cl_size = detectCores()
cl = makeCluster(cl_size)
registerDoParallel(cl)
for (vi in 1:length(v_tr)){
  set.seed(951004)
  result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("VGAM")) %dopar% one_iteration(tr = v_tr[vi],m1 = 500, conserv = TRUE)
  result1 <- cbind(Reduce("+",result[1,])/100,Reduce("+",result[2,])/100,Reduce("+",result[3,])/100,Reduce("+",result[4,])/100)
  v_tr_alpha[,vi] <- result1[1,]
  v_tr_power[,vi] <- result1[2,]
}
stopCluster(cl)

write.csv(v_tr_alpha,"./result/v_tr_alpha_con.csv")
write.csv(v_tr_power,"./result/v_tr_power_con.csv")
sim1_plot(v_tr_alpha[1,],v_tr_alpha[2,],v_tr_alpha[3,],v_tr_alpha[4,],v_tr,c("number of true effects","FDR"))
sim1_plot(v_tr_power[1,],v_tr_power[2,],v_tr_power[3,],v_tr_power[4,],v_tr,c("number of true effects","Power"))


