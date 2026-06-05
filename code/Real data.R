source("AdaPT.R")
source('summarize_methods.R')
source("useful_functions.R")
source("All_q_est_functions.R")
library(VGAM) #generate Laplace distribution
library(mgcv)
#Bottomly data
#BiocManager::install()
library("IHWpaper")

set.seed(19951004)

bottomly <- analyze_dataset("bottomly")
dim(bottomly) #13932
n <- 13932
p <- bottomly$pvalue
p <- pnorm(bottomly$stat)
summary(bottomly$max)
summary(bottomly$min)
which(is.na(p))
x <- log(bottomly$baseMean)

peelm <- 2500
alphalist <- seq(0.01, 0.10, by=0.01) # target FDR level
pi.formulas <- paste0("ns(x, df = ", 6:10, ")")
mu.formulas <- paste0("ns(x, df = ", 6:10, ")")
formulas <- expand.grid(pi.formulas, mu.formulas)

x<-data.frame(x=x)
res.AdaPT <- try(AdaPT.glm(x, p,
                           dist = beta.family(),
                           pi.formula = "ns(x, df = 6)",
                           mu.formula = "ns(x, df = 6)",
                           alpha=alphalist))
res.AdaPT$num_rej

########DP-BH
DPBH_nu <- 0.5*0.01/13932/21
theta <- log(pmax(p,DPBH_nu))
lambda <- sqrt(16*8*peelm)/sqrt(2)/2e4
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
for (q in alphalist) {
  t<-t+1
  for (j in peelm:1) {
    if (new_theta[j]<=(log(q*j/n)-2*lambda*log(6*peelm/q))) {break}
  }
  if ((j==1)&(new_theta[1]>(log(q*1/n)-2*lambda*log(6*peelm/q)))){
    num_rej<-0
  } else{
    rej<-index[1:j]
    num_rej <- length(rej)
  }
  dp_bh_fdp[t] <- num_rej
}

####DP-AdaPT
noisy_sd <- lambda*sqrt(2)

###peeling
ind<-rep(TRUE,n)
n_index<-c()
newqp<-c()

#transformation
qp <- unlist(lapply(p, function(x) qnorm(p = x, mean = 0, sd = 1)))
#select small p
for (i in 1:peelm) {
  #add noisy
  qptilde <- qp + rnorm(n = n, mean = 0, sd = noisy_sd)
  #select p,1-p
  mirror_qptilde <- abs(qptilde)
  n_index[i] <- which(ind)[which.max(mirror_qptilde[ind])]
  newqp[i] <- qp[n_index[i]] + rnorm(n = 1, mean = 0, sd = noisy_sd)
  ind[n_index[i]] <- FALSE
}


#transform back
newp <- unlist(lapply(newqp, function(x) pnorm(q = x, mean = 0, sd = 1)))

#Adapt
dp.res.AdaPT <- try(AdaPT.gam(x[n_index,], newp,
                              dist = beta.family(),
                              pi.formula =formulas[, 1],
                              mu.formula = formulas[, 2],
                              alpha=alphalist))

dp.res.AdaPT$num_rej

#plot
big_table<-data.frame("alpha"=alphalist,"Adapt"=res.AdaPT$num_rej,"DP-BH"=dp_bh_fdp,"DP-Adapt"=dp.res.AdaPT$num_rej)
df <- big_table %>%
  select(alpha, Adapt, DP.BH,DP.Adapt) %>%
  gather(key = "variable", value = "value", -alpha)
colnames(df)[2]<-"method"
df$method[which(df$method=="DP.BH")]<-"DP-BH"
df$method[which(df$method=="DP.Adapt")]<-"DP-Adapt"
ggplot() + geom_line(aes(x=alpha,y=value,group=method,color=method,linetype=method),data=df) + 
  geom_point(aes(x=alpha,y=value,group=method,color=method,shape = method),data=df) + scale_x_continuous(breaks = seq(0, 0.1, by = 0.01))+
  labs(x="Target FDR level",y="Number of Rejections") +
  scale_shape_manual(values=c(15,19,0))+scale_linetype_manual(values=c("longdash","solid","dashed"))+
  scale_color_manual(values=c("blue","black","red"))


