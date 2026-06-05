nCores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
print(paste("nCores",nCores))
myCluster <- parallel::makeCluster(nCores)
doParallel::registerDoParallel(myCluster)

####The simulation setting is the same in
####AdaPT: An interactive procedure for multiple testing with side information
####by Lihua Lei and William Fithian
####This is the adaptive version to control alpha=0.05,0.10,0.15,0.20,0.25,0.30
source("AdaPT.R")
source('summarize_methods.R')
source("useful_functions.R")
source("All_q_est_functions.R")

library(VGAM) #generate Laplace distribution
library(foreach)
library(parallel)
library(doParallel)
library(mgcv)
##Example 1
one_iteration_example1<-function(location = 1,
                                 mu00 = 4,
                                 conserv = FALSE){
  #generate p values
  n <- 100^2
  x1 <- seq(-100, 100, length.out = 100)
  x2 <- seq(-100, 100, length.out = 100)
  x <- expand.grid(x1, x2)
  colnames(x) <- c("x1", "x2")
  pi.formula <- mu.formula <- "s(x1, x2)"
  alpha.list <- seq(0.01, 0.3, 0.01)
  
  if (location==1) {
    H0 <- apply(x, 1, function(coord){sum(coord^2) < 150})
    #120
    mu <- ifelse(H0, mu00, 0)
  }
  if (location==2) {
    H0 <- apply(x, 1, function(coord){sum((coord - 65)^2) < 150})
    #116
    mu <- ifelse(H0, mu00, 0)
  }
  if (location==3) {
    shape.fun <- function(coord){
      transform.coord <- c(coord[1] + coord[2], coord[2] - coord[1])/sqrt(2)
      transform.coord[1]^2 / 100^2 + transform.coord[2]^2 / 15^2 < 0.6/6
    }
    H0 <- apply(x, 1, shape.fun)
    #118
    mu <- ifelse(H0, mu00, 0)
  }
  
  z <- rnorm(n) + mu
  m1 <- length(which(mu == mu00))
  tr_index <- which(mu == mu00)
  peelm <- 500
  pvals <- 1 - pnorm(z)
  if (conserv == TRUE) {pvals[which(mu != mu00)] <- (runif(n = (n-m1),min = 0,max = 1))^{1/4}}
  t1<-Sys.time()
  res.AdaPT <- try(AdaPT.gam(x, pvals,
                             dist = beta.family(),
                             pi.formula ="s(x1, x2)",
                             mu.formula = "s(x1, x2)",
                             alphas = alpha.list))
  t2<-Sys.time()
  adapt.threshold <- res.AdaPT$s
  results <- apply(adapt.threshold, 2, function(xx){
    tmp <- which(pvals < xx)
    nfrej <- length(which(!tmp%in%tr_index))
    ntrej <- length(which(tmp%in%tr_index))
    return(c(nfrej, ntrej))
  })
  nfrej <- as.numeric(results[1, ])
  ntrej <- as.numeric(results[2, ])
  nrej <- nfrej + ntrej
  adapt.FDP <- nfrej / pmax(nrej, 1)
  adapt.power <- ntrej / max(m1,1)
  
  ########DP-BH
  DPBH_eta = 1e-4;DPBH_delta = 0.001;DPBH_epsilon = 0.5;DPBH_nu <- 0.5*0.01/n;
  lambda <- DPBH_eta*sqrt(10*peelm*log(1/DPBH_delta))/(2*DPBH_epsilon)
  theta <- log(pmax(pvals,DPBH_nu))
  
  ##pealing arg min
  index<-c()
  new_theta<-c()
  for (i in 1:peelm) {
    #add noisy
    noisy_theta <- theta + rlaplace(n = n, location = 0, scale = lambda)
    #select
    index[i] <- which.min(noisy_theta)
    new_theta[i] <- theta[index[i]]+rlaplace(n = 1, location = 0, scale = lambda)
    theta[index[i]] <- Inf
  }
  
  #BHq
  t<-0
  dp_bh_fdp<-c()
  dp_bh_power<-c()
  for (q in alpha.list) {
    t<-t+1
    for (j in peelm:1) {
      if (new_theta[j]<=(log(q*j/n)-2*lambda*log(6*peelm/q))) {break}
    }
    if ((j==1)&(new_theta[1]>(log(q*1/n)-2*lambda*log(6*peelm/q)))){
      tdp_bh_fdp<-0
      tdp_bh_power<-0
    } else{
      rej<-index[1:j]
      tdp_bh_fdp <- length(which(!rej%in%tr_index))/max(length(rej),1)
      tdp_bh_power <- length(which(rej%in%tr_index))/m1
    }
    dp_bh_fdp[t] <- tdp_bh_fdp
    dp_bh_power[t] <- tdp_bh_power
  }
  
  ####DP-AdaPT

  ###peeling
  ind<-rep(TRUE,n)
  n_index<-c()
  newqp<-c()
  
  #transformation
  qp <- unlist(lapply(pvals, function(x) qnorm(p = x, mean = 0, sd = 1)))
  #select small p
  for (i in 1:peelm) {
    #add noisy
    qptilde <- qp + rlaplace(n = n, location = 0, scale = lambda)
    #select p,1-p
    mirror_qptilde <- abs(qptilde)
    n_index[i] <- which(ind)[which.max(mirror_qptilde[ind])]
    newqp[i] <- qp[n_index[i]] + rnorm(n = 1, mean = 0, sd = noisy_sd)
    ind[n_index[i]] <- FALSE
  }
  
  
  #transform back
  newp <- unlist(lapply(newqp, function(x) pnorm(q = x, mean = 0, sd = 1)))
  t3<-Sys.time()
  #Adapt
  dp.res.AdaPT <- try(AdaPT.gam(x[n_index,], newp,
                                dist = beta.family(),
                                pi.formula ="s(x1, x2)",
                                mu.formula = "s(x1, x2)",
                                alphas = alpha.list))
  t4<-Sys.time()
  dp.adapt.threshold <- dp.res.AdaPT$s
  results <- apply(dp.adapt.threshold, 2, function(s){
    tmp <- (newp <= s)
    nfrej <- length(which(!n_index[tmp]%in%tr_index))
    ntrej <- length(which(n_index[tmp]%in%tr_index))
    return(c(nfrej, ntrej))
  })
  nfrej <- as.numeric(results[1, ])
  ntrej <- as.numeric(results[2, ])
  nrej <- nfrej + ntrej
  dp.adapt.FDP <- nfrej / pmax(nrej, 1)
  dp.adapt.power <- ntrej / max(m1,1)
  return(list(adapt_alpha = adapt.FDP,adapt_power = adapt.power,dp_bh_fdp = dp_bh_fdp,dp_bh_power =  dp_bh_power,dp_adapt_alpha = dp.adapt.FDP,dp_adapt_power = dp.adapt.power,t11 = difftime(t2, t1, units = "secs")[[1]],t22 = difftime(t4, t3, units = "secs")[[1]]))
}

result1<-NULL
result2<-NULL
result3<-NULL
###alpha case mu00=2.5,location=1
set.seed(951004)
l1m3.5result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("VGAM")) %dopar% one_iteration_example1(conserv = TRUE,mu00 = 2.5)
result1<-cbind(result1,Reduce("+",l1m3.5result[1,])/100,Reduce("+",l1m3.5result[3,])/100,Reduce("+",l1m3.5result[5,])/100)
result2<-cbind(result2,Reduce("+",l1m3.5result[2,])/100,Reduce("+",l1m3.5result[4,])/100,Reduce("+",l1m3.5result[6,])/100)
result3<-cbind(result3,Reduce("+",l1m3.5result[7,])/100,Reduce("+",l1m3.5result[8,])/100)

###alpha case mu00=3.5,location=1
set.seed(951004)
l1m4.5result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("VGAM")) %dopar% one_iteration_example1(conserv = TRUE,mu00 = 3.5)
result1<-cbind(result1,Reduce("+",l1m4.5result[1,])/100,Reduce("+",l1m4.5result[3,])/100,Reduce("+",l1m4.5result[5,])/100)
result2<-cbind(result2,Reduce("+",l1m4.5result[2,])/100,Reduce("+",l1m4.5result[4,])/100,Reduce("+",l1m4.5result[6,])/100)
result3<-cbind(result3,Reduce("+",l1m4.5result[7,])/100,Reduce("+",l1m4.5result[8,])/100)

###alpha case mu00=4.5,location=1
set.seed(951004)
l1m5.5result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("VGAM")) %dopar% one_iteration_example1(conserv = TRUE,mu00 = 4.5)
result1<-cbind(result1,Reduce("+",l1m5.5result[1,])/100,Reduce("+",l1m5.5result[3,])/100,Reduce("+",l1m5.5result[5,])/100)
result2<-cbind(result2,Reduce("+",l1m5.5result[2,])/100,Reduce("+",l1m5.5result[4,])/100,Reduce("+",l1m5.5result[6,])/100)
result3<-cbind(result3,Reduce("+",l1m5.5result[7,])/100,Reduce("+",l1m5.5result[8,])/100)

write.csv(result1,"/work/LAS/zhuz-lab/xintaox/result/let_seting_example1_fdr_con.csv")
write.csv(result2,"/work/LAS/zhuz-lab/xintaox/result/let_seting_example1_power_con.csv")

result1<-NULL
result2<-NULL
###alpha case mu00=2.5,location=2
set.seed(951004)
l1m3.5result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("VGAM")) %dopar% one_iteration_example1(conserv = TRUE,location = 2,mu00 = 2.5)
result1<-cbind(result1,Reduce("+",l1m3.5result[1,])/100,Reduce("+",l1m3.5result[3,])/100,Reduce("+",l1m3.5result[5,])/100)
result2<-cbind(result2,Reduce("+",l1m3.5result[2,])/100,Reduce("+",l1m3.5result[4,])/100,Reduce("+",l1m3.5result[6,])/100)
result3<-cbind(result3,Reduce("+",l1m3.5result[7,])/100,Reduce("+",l1m3.5result[8,])/100)

###alpha case mu00=3.5,location=2
set.seed(951004)
l1m4.5result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("VGAM")) %dopar% one_iteration_example1(conserv = TRUE,location = 2,mu00 = 3.5)
result1<-cbind(result1,Reduce("+",l1m4.5result[1,])/100,Reduce("+",l1m4.5result[3,])/100,Reduce("+",l1m4.5result[5,])/100)
result2<-cbind(result2,Reduce("+",l1m4.5result[2,])/100,Reduce("+",l1m4.5result[4,])/100,Reduce("+",l1m4.5result[6,])/100)
result3<-cbind(result3,Reduce("+",l1m4.5result[7,])/100,Reduce("+",l1m4.5result[8,])/100)

###alpha case mu00=4.5,location=2
set.seed(951004)
l1m5.5result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("VGAM")) %dopar% one_iteration_example1(conserv = TRUE,location = 2,mu00 = 4.5)
result1<-cbind(result1,Reduce("+",l1m5.5result[1,])/100,Reduce("+",l1m5.5result[3,])/100,Reduce("+",l1m5.5result[5,])/100)
result2<-cbind(result2,Reduce("+",l1m5.5result[2,])/100,Reduce("+",l1m5.5result[4,])/100,Reduce("+",l1m5.5result[6,])/100)
result3<-cbind(result3,Reduce("+",l1m5.5result[7,])/100,Reduce("+",l1m5.5result[8,])/100)

write.csv(result1,"/work/LAS/zhuz-lab/xintaox/result/let_seting_example2_fdr_con.csv")
write.csv(result2,"/work/LAS/zhuz-lab/xintaox/result/let_seting_example2_power_con.csv")

result1<-NULL
result2<-NULL
###alpha case mu00=2.5,location=3
set.seed(951004)
l1m3.5result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("VGAM")) %dopar% one_iteration_example1(conserv = TRUE,location = 3,mu00 = 2.5)
result1<-cbind(result1,Reduce("+",l1m3.5result[1,])/100,Reduce("+",l1m3.5result[3,])/100,Reduce("+",l1m3.5result[5,])/100)
result2<-cbind(result2,Reduce("+",l1m3.5result[2,])/100,Reduce("+",l1m3.5result[4,])/100,Reduce("+",l1m3.5result[6,])/100)
result3<-cbind(result3,Reduce("+",l1m3.5result[7,])/100,Reduce("+",l1m3.5result[8,])/100)

###alpha case mu00=3.5,location=3
set.seed(951004)
l1m4.5result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("VGAM")) %dopar% one_iteration_example1(conserv = TRUE,location = 3,mu00 = 3.5)
result1<-cbind(result1,Reduce("+",l1m4.5result[1,])/100,Reduce("+",l1m4.5result[3,])/100,Reduce("+",l1m4.5result[5,])/100)
result2<-cbind(result2,Reduce("+",l1m4.5result[2,])/100,Reduce("+",l1m4.5result[4,])/100,Reduce("+",l1m4.5result[6,])/100)
result3<-cbind(result3,Reduce("+",l1m4.5result[7,])/100,Reduce("+",l1m4.5result[8,])/100)

###alpha case mu00=4.5,location=3
set.seed(951004)
l1m5.5result <-  foreach(mc = 1:100, .combine = cbind,.packages = c("VGAM")) %dopar% one_iteration_example1(conserv = TRUE,location = 3,mu00 = 4.5)
result1<-cbind(result1,Reduce("+",l1m5.5result[1,])/100,Reduce("+",l1m5.5result[3,])/100,Reduce("+",l1m5.5result[5,])/100)
result2<-cbind(result2,Reduce("+",l1m5.5result[2,])/100,Reduce("+",l1m5.5result[4,])/100,Reduce("+",l1m5.5result[6,])/100)
result3<-cbind(result3,Reduce("+",l1m5.5result[7,])/100,Reduce("+",l1m5.5result[8,])/100)

write.csv(result1,"/work/LAS/zhuz-lab/xintaox/result/let_seting_example3_fdr_con.csv")
write.csv(result2,"/work/LAS/zhuz-lab/xintaox/result/let_seting_example3_power_con.csv")
write.csv(result3,"/work/LAS/zhuz-lab/xintaox/result/computation_time_con.csv")

stopCluster(myCluster)