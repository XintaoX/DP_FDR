####The simulation setting is the same in
####"Differentially Private False Discovery Rate Control"
####This is the adaptive version to control alpha=0.05,0.10,0.15,0.20,0.25,0.30
####by Cynthia Dwork, Weijie J. Su and Li Zhang
library(VGAM) #generate Laplace distribution
set.seed(951004)

one_iteration<-function(m = 1e5,
                        m1 = 100,
                        tr = 100,
                        alphaa = seq(0.1,0.3,by = 0.01),
                        eta = 1e-4,
                        delta = 0.001,
                        mu = 4,
                        epsilon = 0.5,
                        conserv = FALSE){
  nu <- 0.5*0.05/m
  
  #generate p values
  xi <- rnorm(n = tr, mean = 0, sd = 1)
  p <- unlist(lapply(xi, function(x) pnorm(q = x-mu,mean = 0,sd = 1)))
  if (conserv) { p <- c(p,rbeta(n = m -tr,shape1 = 2,shape2 = 2))} else {
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
  t<-0
  dp_bh_fdp<-c()
  dp_bh_power<-c()
  for (q in alphaa) {
    t<-t+1
    for (j in m1:1) {
      if (new_theta[j]<=(log(q*j/m)-2*lambda*log(6*m1/q))) {break}
    }
    if ((j==1)&(new_theta[1]>(log(q*1/m)-2*lambda*log(6*m1/q)))){
      tdp_bh_fdp<-0
      tdp_bh_power<-0
    } else{
      rej<-index[1:j]
      tdp_bh_fdp <- length(which(rej>tr))/max(length(rej),1)
      tdp_bh_power <- length(which(rej<=tr))/tr
    }
    dp_bh_fdp[t] <- tdp_bh_fdp
    dp_bh_power[t] <- tdp_bh_power
  }

  
  ########DP-Bonf
  t<-0
  dp_bonf_fdp<-c()
  dp_bonf_power<-c()
  lambda_tilde <- eta*sqrt(10*m*log(1/delta))/(2*epsilon)
  theta_tilde <- log(pmax(p,nu)) + rlaplace(n = m, location = 0, scale = lambda_tilde)
  for (q in alphaa) {
    t <- t+1
    DP_Bonf_threshold <- log(q/m) - lambda_tilde*log(5*m/q)
    Bonf_rej <- which(theta_tilde <= DP_Bonf_threshold)
    dp_bonf_fdp[t] <- length(which(Bonf_rej>tr))/max(length(Bonf_rej),1)
    dp_bonf_power[t] <- length(which(Bonf_rej<=tr))/tr
  }
  
  ####DP-AdaPT
  ratio <- 0.5
  noisy_sd <- lambda*sqrt(2)*sqrt(1+ratio)
  
  ###peeling
  ind<-rep(TRUE,m)
  n_index<-c()
  newqp<-c()
  
  #transformation
  qp <- unlist(lapply(p, function(x) qnorm(p = x, mean = 0, sd = 1)))
  #select small p
  for (i in 1:m1) {
    #add noisy
    qptilde <- qp + rnorm(n = m, mean = 0, sd = noisy_sd)
    #select
    n_index[i] <- which(ind)[which.min(qptilde[ind])]
    newqp[i] <- qp[n_index[i]] + rnorm(n = 1, mean = 0, sd = noisy_sd)
    ind[n_index[i]] <- FALSE
  }
  
  #select large p
  for (i in (m1+1):(m1+floor(m1*ratio))) {
    #add noisy
    qptilde <- qp + rnorm(n = m, mean = 0, sd = noisy_sd)
    #select
    n_index[i] <- which(ind)[which.max(qptilde[ind])]
    newqp[i] <- qp[n_index[i]] + rnorm(n = 1, mean = 0, sd = noisy_sd)
    ind[n_index[i]] <- FALSE
  }
  
  #transform back
  newp <- unlist(lapply(newqp, function(x) pnorm(q = x, mean = 0, sd = 1)))
  control<-TRUE
  ini <- 0.45
  True_fdr <- rep(0,length(alphaa))
  power_est <- rep(0,length(alphaa))
  t<-1
  while (control) {
    At <- length(which(newp>=(1-ini)))
    Rt <- length(which(newp<=ini))
    fdpt <- (1+At)/max(Rt,1)
    
    if (fdpt<=alphaa[length(alphaa)+1-t]) {
      new_rej <- n_index[which(newp<=ini)]
      if (length(new_rej)>0) {
          True_fdr[length(alphaa)+1-t]<-length(which(new_rej>tr))/max(length(new_rej),1)
          power_est[length(alphaa)+1-t]<-length(which(new_rej<=tr))/tr
      }
      t<-t+1
      } else {
      potent_down <- which(newp<ini)
      potent_up <- which(newp>(1-ini))
      seqq <- c(ini-newp[potent_down],newp[potent_up]+ini-1)
      ini <- ini - (seqq[order(seqq)[2]]+seqq[order(seqq)[1]])/2
    }
    Rt <- length(which(newp<=ini))
    if (t==(length(alphaa)+1)) {control <- FALSE}
    if (Rt==0) {control<-FALSE}
  }

  return(list(dp_bh_1 = dp_bh_fdp,dp_bh_2 = dp_bh_power,dp_bonf_1 = dp_bonf_fdp,dp_bonf_2 = dp_bonf_power,dp_adapt_1 =  True_fdr,dp_adapt_2 = power_est))
}

###alpha
set.seed(951004)

result1 <- matrix(data = 0,nrow = 6,ncol = 21)
for (mc in 1:100) {
    print(mc)
    temp <- one_iteration()
    result1 <- result1 + rbind(temp$dp_bh_1,temp$dp_bh_2,temp$dp_bonf_1,temp$dp_bonf_2,temp$dp_adapt_1,temp$dp_adapt_2)
  }
result1 <- result1/100

plot(seq(0.1,0.3,by = 0.01),result1[1,])
plot(seq(0.1,0.3,by = 0.01),result1[3,])
plot(seq(0.1,0.3,by = 0.01),result1[5,])

plot(seq(0.1,0.3,by = 0.01),result1[2,])
plot(seq(0.1,0.3,by = 0.01),result1[4,])
plot(seq(0.1,0.3,by = 0.01),result1[6,])

write.csv(v_epsilon_alpha,"./result/v_epsilon_alpha.csv")
write.csv(v_epsilon_power,"./result/v_epsilon_power.csv")
sim1_plot(v_epsilon_alpha[1,],v_epsilon_alpha[2,],v_epsilon_alpha[3,],v_epsilon,c("$\\epsilon$","FDR"))
sim1_plot(v_epsilon_power[1,],v_epsilon_power[2,],v_epsilon_power[3,],v_epsilon,c("$\\epsilon$","Power"))
