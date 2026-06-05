library(ggplot2)
library(tidyverse)
library(latex2exp)
v_epsilon <- seq(0.1,1,by = 0.1)
v_eta <- 10^(seq(-5,-3,by = 0.25))
v_beta <- seq(2,7,by = 1)
v_tr <- seq(50,300,by = 50)
#sim1_fdr
epsilon <- read.csv("./result/v_epsilon_alpha.csv")
eta <- read.csv("./result/v_eta_alpha.csv")
beta <- read.csv("./result/v_beta_alpha.csv")
tr <- read.csv("./result/v_tr_alpha.csv")
table1 <- as.data.frame(t(epsilon))[-1,]
table2 <- as.data.frame(t(eta))[-1,]
table3 <- as.data.frame(t(beta))[-1,]
table4 <- as.data.frame(t(tr))[-1,]
table1 <- table1 %>% mutate(x=v_epsilon, class='epsilon') 
table2 <- table2 %>% mutate(x=log10(v_eta), class='log[10]~eta')
table3 <- table3 %>% mutate(x=v_beta, class='beta')
table4 <- table4 %>% mutate(x=v_tr, class='number~of~non~null~effects')
big_table <- bind_rows(table1,table2,table3,table4)
big_table <- big_table %>% reshape2::melt(id=c('class','x'))
colnames(big_table)[3] <- "method"
levels(big_table$method) <- c("DP-BH", "DP-Bonf", "Adapt", "DP-Adapt")
ggplot() + geom_line(aes(x=x,y=value,group=method,color=method,linetype=method),data=big_table) + facet_wrap(~class, scales = 'free_x',labeller = label_parsed) +
  ylim(0,0.12) + geom_point(aes(x=x,y=value,group=method,color=method,shape = method),data=big_table) + labs(x=NULL,y="FDR") +
  scale_shape_manual(values=c(0,1,15,19))+scale_linetype_manual(values=c("dashed","dotted","longdash","solid"))+
  scale_color_manual(values=c("red","green","blue","black"))

#sim1_power
epsilon <- read.csv("./result/v_epsilon_power.csv")
eta <- read.csv("./result/v_eta_power.csv")
beta <- read.csv("./result/v_beta_power.csv")
tr <- read.csv("./result/v_tr_power.csv")
table1 <- as.data.frame(t(epsilon))[-1,]
table2 <- as.data.frame(t(eta))[-1,]
table3 <- as.data.frame(t(beta))[-1,]
table4 <- as.data.frame(t(tr))[-1,]
table1 <- table1 %>% mutate(x=v_epsilon, class='epsilon') 
table2 <- table2 %>% mutate(x=log10(v_eta), class='log[10]~eta')
table3 <- table3 %>% mutate(x=v_beta, class='beta')
table4 <- table4 %>% mutate(x=v_tr, class='number~of~non~null~effects')
big_table <- bind_rows(table1,table2,table3,table4)
big_table <- big_table %>% reshape2::melt(id=c('class','x'))
colnames(big_table)[3] <- "method"
levels(big_table$method) <- c("DP-BH", "DP-Bonf", "Adapt", "DP-Adapt")
ggplot() + geom_line(aes(x=x,y=value,group=method,color=method,linetype=method),data=big_table) + facet_wrap(~class, scales = 'free_x',labeller = label_parsed) +
  ylim(0,1) + geom_point(aes(x=x,y=value,group=method,color=method,shape = method),data=big_table) + labs(x=NULL,y="Power") +
  scale_shape_manual(values=c(0,1,15,19))+scale_linetype_manual(values=c("dashed","dotted","longdash","solid"))+
  scale_color_manual(values=c("red","green","blue","black"))

#sim1_fdr_con
epsilon <- read.csv("./result/v_epsilon_alpha_con.csv")
eta <- read.csv("./result/v_eta_alpha_con.csv")
beta <- read.csv("./result/v_beta_alpha_con.csv")
tr <- read.csv("./result/v_tr_alpha_con.csv")
table1 <- as.data.frame(t(epsilon))[-1,]
table2 <- as.data.frame(t(eta))[-1,]
table3 <- as.data.frame(t(beta))[-1,]
table4 <- as.data.frame(t(tr))[-1,]
table1 <- table1 %>% mutate(x=v_epsilon, class='epsilon') 
table2 <- table2 %>% mutate(x=log10(v_eta), class='log[10]~eta')
table3 <- table3 %>% mutate(x=v_beta, class='beta')
table4 <- table4 %>% mutate(x=v_tr, class='number~of~non~null~effects')
big_table <- bind_rows(table1,table2,table3,table4)
big_table <- big_table %>% reshape2::melt(id=c('class','x'))
colnames(big_table)[3] <- "method"
levels(big_table$method) <- c("DP-BH", "DP-Bonf", "Adapt", "DP-Adapt")
ggplot() + geom_line(aes(x=x,y=value,group=method,color=method,linetype=method),data=big_table) + facet_wrap(~class, scales = 'free_x',labeller = label_parsed) +
  ylim(0,0.12) + geom_point(aes(x=x,y=value,group=method,color=method,shape = method),data=big_table) + labs(x=NULL,y="FDR") +
  scale_shape_manual(values=c(0,1,15,19))+scale_linetype_manual(values=c("dashed","dotted","longdash","solid"))+
  scale_color_manual(values=c("red","green","blue","black"))

#sim1_power_con
epsilon <- read.csv("./result/v_epsilon_power_con.csv")
eta <- read.csv("./result/v_eta_power_con.csv")
beta <- read.csv("./result/v_beta_power_con.csv")
tr <- read.csv("./result/v_tr_power_con.csv")
table1 <- as.data.frame(t(epsilon))[-1,]
table2 <- as.data.frame(t(eta))[-1,]
table3 <- as.data.frame(t(beta))[-1,]
table4 <- as.data.frame(t(tr))[-1,]
table1 <- table1 %>% mutate(x=v_epsilon, class='epsilon') 
table2 <- table2 %>% mutate(x=log10(v_eta), class='log[10]~eta')
table3 <- table3 %>% mutate(x=v_beta, class='beta')
table4 <- table4 %>% mutate(x=v_tr, class='number~of~non~null~effects')
big_table <- bind_rows(table1,table2,table3,table4)
big_table <- big_table %>% reshape2::melt(id=c('class','x'))
colnames(big_table)[3] <- "method"
levels(big_table$method) <- c("DP-BH", "DP-Bonf", "Adapt", "DP-Adapt")
ggplot() + geom_line(aes(x=x,y=value,group=method,color=method,linetype=method),data=big_table) + facet_wrap(~class, scales = 'free_x',labeller = label_parsed) +
  ylim(0,1) + geom_point(aes(x=x,y=value,group=method,color=method,shape = method),data=big_table) + labs(x=NULL,y="Power") +
  scale_shape_manual(values=c(0,1,15,19))+scale_linetype_manual(values=c("dashed","dotted","longdash","solid"))+
  scale_color_manual(values=c("red","green","blue","black"))


#sim2_fdr
ex1 <- read.csv("./result/let_seting_example1_fdr.csv")
ex2 <- read.csv("./result/let_seting_example2_fdr.csv")
ex3 <- read.csv("./result/let_seting_example3_fdr.csv")

alpha.list <- seq(0.01, 0.3, 0.01)

table1 <- (as.data.frame((ex1))[,-1])[,1:3]
table2 <- (as.data.frame((ex1))[,-1])[,4:6]
table3 <- (as.data.frame((ex1))[,-1])[,7:9]
table4 <- (as.data.frame((ex2))[,-1])[,1:3]
table5 <- (as.data.frame((ex2))[,-1])[,4:6]
table6 <- (as.data.frame((ex2))[,-1])[,7:9]
table7 <- (as.data.frame((ex3))[,-1])[,1:3]
table8 <- (as.data.frame((ex3))[,-1])[,4:6]
table9 <- (as.data.frame((ex3))[,-1])[,7:9]

table1 <- table1 %>% mutate(x=alpha.list, class=paste(as.roman(1),'~":"~beta==2.5')) 
table2 <- table2 %>% mutate(x=alpha.list, class=paste(as.roman(1),'~":"~beta==3.5'))
table3 <- table3 %>% mutate(x=alpha.list, class=paste(as.roman(1),'~":"~beta==4.5'))
table4 <- table4 %>% mutate(x=alpha.list, class=paste(as.roman(2),'~":"~beta==2.5'))
table5 <- table5 %>% mutate(x=alpha.list, class=paste(as.roman(2),'~":"~beta==3.5')) 
table6 <- table6 %>% mutate(x=alpha.list, class=paste(as.roman(2),'~":"~beta==4.5'))
table7 <- table7 %>% mutate(x=alpha.list, class=paste(as.roman(3),'~":"~beta==2.5'))
table8 <- table8 %>% mutate(x=alpha.list, class=paste(as.roman(3),'~":"~beta==3.5'))
table9 <- table9 %>% mutate(x=alpha.list, class=paste(as.roman(3),'~":"~beta==4.5')) 

colnames(table2)[1:3] <- colnames(table1)[1:3]
colnames(table3)[1:3] <- colnames(table1)[1:3]
colnames(table5)[1:3] <- colnames(table1)[1:3]
colnames(table6)[1:3] <- colnames(table1)[1:3]
colnames(table8)[1:3] <- colnames(table1)[1:3]
colnames(table9)[1:3] <- colnames(table1)[1:3]

big_table <- bind_rows(table1,table2,table3,table4,table5,table6,table7,table8,table9)
big_table <- big_table %>% reshape2::melt(id=c('class','x'))
colnames(big_table)[3] <- "method"
levels(big_table$method) <- c("Adapt","DP-BH","DP-Adapt")
ggplot() + geom_line(aes(x=x,y=value,group=method,color=method,linetype=method),data=big_table) + facet_wrap(~class,labeller = label_parsed) +
  ylim(0,0.31) + geom_point(aes(x=x,y=value,group=method,color=method,shape = method),data=subset(big_table,x %in% seq(0.01,0.30,by = 0.03))) + labs(x="Target FDR level",y="FDR") +
  scale_shape_manual(values=c(15,0,19))+scale_linetype_manual(values=c("longdash","dashed","solid"))+
  scale_color_manual(values=c("blue","red","black"))

#sim2_power
ex1 <- read.csv("./result/let_seting_example1_power.csv")
ex2 <- read.csv("./result/let_seting_example2_power.csv")
ex3 <- read.csv("./result/let_seting_example3_power.csv")

alpha.list <- seq(0.01, 0.3, 0.01)

table1 <- (as.data.frame((ex1))[,-1])[,1:3]
table2 <- (as.data.frame((ex1))[,-1])[,4:6]
table3 <- (as.data.frame((ex1))[,-1])[,7:9]
table4 <- (as.data.frame((ex2))[,-1])[,1:3]
table5 <- (as.data.frame((ex2))[,-1])[,4:6]
table6 <- (as.data.frame((ex2))[,-1])[,7:9]
table7 <- (as.data.frame((ex3))[,-1])[,1:3]
table8 <- (as.data.frame((ex3))[,-1])[,4:6]
table9 <- (as.data.frame((ex3))[,-1])[,7:9]

table1 <- table1 %>% mutate(x=alpha.list, class=paste(as.roman(1),'~":"~beta==2.5')) 
table2 <- table2 %>% mutate(x=alpha.list, class=paste(as.roman(1),'~":"~beta==3.5'))
table3 <- table3 %>% mutate(x=alpha.list, class=paste(as.roman(1),'~":"~beta==4.5'))
table4 <- table4 %>% mutate(x=alpha.list, class=paste(as.roman(2),'~":"~beta==2.5'))
table5 <- table5 %>% mutate(x=alpha.list, class=paste(as.roman(2),'~":"~beta==3.5')) 
table6 <- table6 %>% mutate(x=alpha.list, class=paste(as.roman(2),'~":"~beta==4.5'))
table7 <- table7 %>% mutate(x=alpha.list, class=paste(as.roman(3),'~":"~beta==2.5'))
table8 <- table8 %>% mutate(x=alpha.list, class=paste(as.roman(3),'~":"~beta==3.5'))
table9 <- table9 %>% mutate(x=alpha.list, class=paste(as.roman(3),'~":"~beta==4.5')) 

colnames(table2)[1:3] <- colnames(table1)[1:3]
colnames(table3)[1:3] <- colnames(table1)[1:3]
colnames(table5)[1:3] <- colnames(table1)[1:3]
colnames(table6)[1:3] <- colnames(table1)[1:3]
colnames(table8)[1:3] <- colnames(table1)[1:3]
colnames(table9)[1:3] <- colnames(table1)[1:3]

big_table <- bind_rows(table1,table2,table3,table4,table5,table6,table7,table8,table9)
big_table <- big_table %>% reshape2::melt(id=c('class','x'))
colnames(big_table)[3] <- "method"
levels(big_table$method) <- c("Adapt","DP-BH","DP-Adapt")
ggplot() + geom_line(aes(x=x,y=value,group=method,color=method,linetype=method),data=big_table) + facet_wrap(~class,labeller = label_parsed) +
  ylim(0,1) + geom_point(aes(x=x,y=value,group=method,color=method,shape = method),data=subset(big_table,x %in% seq(0.01,0.30,by = 0.03))) + labs(x="Target FDR level",y="Power") +
  scale_shape_manual(values=c(15,0,19))+scale_linetype_manual(values=c("longdash","dashed","solid"))+
  scale_color_manual(values=c("blue","red","black"))

#sim2_fdr_con
ex1 <- read.csv("./result/let_seting_example1_fdr_con.csv")
ex2 <- read.csv("./result/let_seting_example2_fdr_con.csv")
ex3 <- read.csv("./result/let_seting_example3_fdr_con.csv")

alpha.list <- seq(0.01, 0.3, 0.01)

table1 <- (as.data.frame((ex1))[,-1])[,1:3]
table2 <- (as.data.frame((ex1))[,-1])[,4:6]
table3 <- (as.data.frame((ex1))[,-1])[,7:9]
table4 <- (as.data.frame((ex2))[,-1])[,1:3]
table5 <- (as.data.frame((ex2))[,-1])[,4:6]
table6 <- (as.data.frame((ex2))[,-1])[,7:9]
table7 <- (as.data.frame((ex3))[,-1])[,1:3]
table8 <- (as.data.frame((ex3))[,-1])[,4:6]
table9 <- (as.data.frame((ex3))[,-1])[,7:9]

table1 <- table1 %>% mutate(x=alpha.list, class=paste(as.roman(1),'~":"~beta==2.5')) 
table2 <- table2 %>% mutate(x=alpha.list, class=paste(as.roman(1),'~":"~beta==3.5'))
table3 <- table3 %>% mutate(x=alpha.list, class=paste(as.roman(1),'~":"~beta==4.5'))
table4 <- table4 %>% mutate(x=alpha.list, class=paste(as.roman(2),'~":"~beta==2.5'))
table5 <- table5 %>% mutate(x=alpha.list, class=paste(as.roman(2),'~":"~beta==3.5')) 
table6 <- table6 %>% mutate(x=alpha.list, class=paste(as.roman(2),'~":"~beta==4.5'))
table7 <- table7 %>% mutate(x=alpha.list, class=paste(as.roman(3),'~":"~beta==2.5'))
table8 <- table8 %>% mutate(x=alpha.list, class=paste(as.roman(3),'~":"~beta==3.5'))
table9 <- table9 %>% mutate(x=alpha.list, class=paste(as.roman(3),'~":"~beta==4.5')) 

colnames(table2)[1:3] <- colnames(table1)[1:3]
colnames(table3)[1:3] <- colnames(table1)[1:3]
colnames(table5)[1:3] <- colnames(table1)[1:3]
colnames(table6)[1:3] <- colnames(table1)[1:3]
colnames(table8)[1:3] <- colnames(table1)[1:3]
colnames(table9)[1:3] <- colnames(table1)[1:3]

big_table <- bind_rows(table1,table2,table3,table4,table5,table6,table7,table8,table9)
big_table <- big_table %>% reshape2::melt(id=c('class','x'))
colnames(big_table)[3] <- "method"
levels(big_table$method) <- c("Adapt","DP-BH","DP-Adapt")
ggplot() + geom_line(aes(x=x,y=value,group=method,color=method,linetype=method),data=big_table) + facet_wrap(~class,labeller = label_parsed) +
  ylim(0,0.31) + geom_point(aes(x=x,y=value,group=method,color=method,shape = method),data=subset(big_table,x %in% seq(0.01,0.30,by = 0.03))) + labs(x="Target FDR level",y="FDR") +
  scale_shape_manual(values=c(15,0,19))+scale_linetype_manual(values=c("longdash","dashed","solid"))+
  scale_color_manual(values=c("blue","red","black"))

#sim2_power_con
ex1 <- read.csv("./result/let_seting_example1_power_con.csv")
ex2 <- read.csv("./result/let_seting_example2_power_con.csv")
ex3 <- read.csv("./result/let_seting_example3_power_con.csv")

alpha.list <- seq(0.01, 0.3, 0.01)

table1 <- (as.data.frame((ex1))[,-1])[,1:3]
table2 <- (as.data.frame((ex1))[,-1])[,4:6]
table3 <- (as.data.frame((ex1))[,-1])[,7:9]
table4 <- (as.data.frame((ex2))[,-1])[,1:3]
table5 <- (as.data.frame((ex2))[,-1])[,4:6]
table6 <- (as.data.frame((ex2))[,-1])[,7:9]
table7 <- (as.data.frame((ex3))[,-1])[,1:3]
table8 <- (as.data.frame((ex3))[,-1])[,4:6]
table9 <- (as.data.frame((ex3))[,-1])[,7:9]

table1 <- table1 %>% mutate(x=alpha.list, class=paste(as.roman(1),'~":"~beta==2.5')) 
table2 <- table2 %>% mutate(x=alpha.list, class=paste(as.roman(1),'~":"~beta==3.5'))
table3 <- table3 %>% mutate(x=alpha.list, class=paste(as.roman(1),'~":"~beta==4.5'))
table4 <- table4 %>% mutate(x=alpha.list, class=paste(as.roman(2),'~":"~beta==2.5'))
table5 <- table5 %>% mutate(x=alpha.list, class=paste(as.roman(2),'~":"~beta==3.5')) 
table6 <- table6 %>% mutate(x=alpha.list, class=paste(as.roman(2),'~":"~beta==4.5'))
table7 <- table7 %>% mutate(x=alpha.list, class=paste(as.roman(3),'~":"~beta==2.5'))
table8 <- table8 %>% mutate(x=alpha.list, class=paste(as.roman(3),'~":"~beta==3.5'))
table9 <- table9 %>% mutate(x=alpha.list, class=paste(as.roman(3),'~":"~beta==4.5')) 

colnames(table2)[1:3] <- colnames(table1)[1:3]
colnames(table3)[1:3] <- colnames(table1)[1:3]
colnames(table5)[1:3] <- colnames(table1)[1:3]
colnames(table6)[1:3] <- colnames(table1)[1:3]
colnames(table8)[1:3] <- colnames(table1)[1:3]
colnames(table9)[1:3] <- colnames(table1)[1:3]

big_table <- bind_rows(table1,table2,table3,table4,table5,table6,table7,table8,table9)
big_table <- big_table %>% reshape2::melt(id=c('class','x'))
colnames(big_table)[3] <- "method"
levels(big_table$method) <- c("Adapt","DP-BH","DP-Adapt")
ggplot() + geom_line(aes(x=x,y=value,group=method,color=method,linetype=method),data=big_table) + facet_wrap(~class,labeller = label_parsed) +
  ylim(0,1) + geom_point(aes(x=x,y=value,group=method,color=method,shape = method),data=subset(big_table,x %in% seq(0.01,0.30,by = 0.03))) + labs(x="Target FDR level",y="Power") +
  scale_shape_manual(values=c(15,0,19))+scale_linetype_manual(values=c("longdash","dashed","solid"))+
  scale_color_manual(values=c("blue","red","black"))


####computation time
computation_time <- read.csv("./result/computation_time.csv")
computation_time_con <- read.csv("./result/computation_time_con.csv")
computation_time <- as.numeric(computation_time)
computation_time_con <- as.numeric(computation_time_con)
xtable(rbind(computation_time[2:7],computation_time[8:13],computation_time[14:19]))
xtable(rbind(computation_time_con[2:7],computation_time_con[8:13],computation_time_con[14:19]))


######mirror density
np1 <- cbind(p1,rep("Uniform null before transformation",length(p1)))
np2 <- cbind(p2,rep("Uniform null after transformation",length(p2)))
np3 <- cbind(p3,rep("Quadratic null before transformation",length(p3)))
np4 <- cbind(p4,rep("Quadratic null after transformation",length(p4)))

den <- data.frame(rbind(np1,np2,np3,np4))
colnames(den) <- c("p","type")

ggplot(den, aes(x=p))+geom_density()+facet_grid(type ~ .)



