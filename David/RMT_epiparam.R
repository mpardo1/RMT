####### RANDOM MATRICES FOR METAPOPULATION MODELS #######
### 
### parent script
###
### generate, plot and integrate metapopulation
### epidemiological models
### 
rm(list = ls())
source("~/RMT/David/RMT_genrandom.R")
source("~/RMT/David/RMT_plotmobility.R")
source("~/RMT/David/d_functions_eigen_int.R")

####### GENERATE JACOBIAN ###############################
# number of patches
N <- 50

# epidemiological
#all rates must lie in (0,1) except for betas

Deltas <- rep(0.3, N) # birth rate
mub <- 0.1
sb <- 0.001
betas <- rep(mub, N) # transmission rates
# betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.3, N) # loss of immunity rates
mud <- 0.3
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.2
alphas <- rep(mua, N) # recovery rates
mudel <- 0.4
deltas <- rep(mudel, N) # disease-related death rates
gammas = deaths + alphas + deltas

# mobility
#commuting and migration networks
muw <- 0.2 
sw <- 0.05
rhow <- 0 #original rho (Gamma of baron et al)
Gammaw <- 0 #gamma of baron et al
rw <- 0
cw <- 0

muc <- 0.001
sc <- 0.0001
rhoc <- 0
Gammac <- 0
rc <- 0
cc <- 0

COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
diag(COMMUTING) <- 0
# COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
# COMMUTING[sample.int(N^2, round(p*N^2))] <- 0

MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
diag(MIGRATION) <- 0

sus_init <- rep(10000, N) # initial susceptibles
inf_init <- rep(100, N)    # initial infecteds

end_time <- 50
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)
plot_int(N,sol, "INF")

min_int <- 0
max_int <- 1
step <- 0.05

## Tethas Recovery rate ####
thetas_vec <- seq(min_int,max_int,step)
end_time <- 50
thetas <- rep(thetas_vec[1],N)
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)
df_sum_theta <- data.frame(time = sol[,1])
for(i in c(1:length(thetas_vec))){
  print(paste0("i: ", i))
  thetas <- rep(thetas_vec[i],N)
  sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
             COMMUTING,MIGRATION,
             sus_init,inf_init,end_time)
  sol_inf <- sol[,c(1,(N+2):(2*N+1))]
  sum_inf <- rowSums(sol_inf[,2:(N+1)])
  df_sum_theta[,ncol(df_sum_theta)+1] <- sum_inf
}

colnames(df_sum_theta) <- c("time",as.character(thetas_vec))
df_plot_theta <- reshape2::melt(df_sum_theta, id.vars = c("time"))

size_let <- 13

sum_tetha <- ggplot(data = df_plot_theta,
                    aes(x = time, y = value, 
                        color = as.numeric(as.character(variable)),
                        group = variable)) +
  geom_line() +
  scale_colour_gradient(name = ""*theta~","*delta~", d, "*alpha~"", 
                        low = "blue", high = "red",
                        breaks=c(0,0.5,1),
                        labels=c(0,0.5,1),
                        limits=c(0,1)) +
  ylab("Sum of infected individuals") +
  scale_y_continuous(labels = function(x)
    format(x, scientific = TRUE)) +
  theme_bw() +
  theme(text = element_text(size = size_let)) 

thetas <- rep(0.3, N) # loss of immunity rates

##### Disease mortality #####
deltas_vec <-  seq(min_int,max_int,step)
end_time <- 50
deltas <- rep(deltas_vec[1],N)
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)
df_sum_delt <- data.frame(time = sol[,1])
for(i in c(1:length(deltas_vec))){
  print(paste0("i: ", i))
  deltas <- rep(deltas_vec[i],N)
  sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
             COMMUTING,MIGRATION,
             sus_init,inf_init,end_time)
  sol_inf <- sol[,c(1,(N+2):(2*N+1))]
  sum_inf <- rowSums(sol_inf[,2:(N+1)])
  df_sum_delt[,ncol(df_sum_delt)+1] <- sum_inf
}

colnames(df_sum_delt) <- c("time",as.character(deltas_vec))
df_plot_delt <- reshape2::melt(df_sum_delt, id.vars = c("time"))

size_let <- 13

sum_deltas <- ggplot(data = df_plot_delt, aes(x = time, y = value,
                                      color = as.numeric(as.character(variable)),
                                      group = variable)) +
  geom_line() +
  scale_colour_gradient(name = ""*theta~","*delta~", d, "*alpha~"", 
                        low = "blue", high = "red",
                        breaks=c(0,0.5,1),
                        labels=c(0,0.5,1),
                        limits=c(0,1)) +
  ylab("Sum of infected individuals") +
  scale_y_continuous(labels = function(x)
    format(x, scientific = TRUE)) +
  theme_bw() +
  theme(text = element_text(size = size_let)) 

mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates

##### Natural mortility #####
death_vec <-  seq(min_int,max_int,step)
end_time <- 50
deaths <- rep(death_vec[1],N)
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)
df_sum_deat <- data.frame(time = sol[,1])
for(i in c(1:length(death_vec))){
  print(paste0("i: ", i))
  deaths <- rep(death_vec[i],N)
  sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
             COMMUTING,MIGRATION,
             sus_init,inf_init,end_time)
  sol_inf <- sol[,c(1,(N+2):(2*N+1))]
  sum_inf <- rowSums(sol_inf[,2:(N+1)])
  df_sum_deat[,ncol(df_sum_deat)+1] <- sum_inf
}

colnames(df_sum_deat) <- c("time",as.character(death_vec))
df_plot_deat <- reshape2::melt(df_sum_deat, id.vars = c("time"))

size_let <- 13

sum_deaths <- ggplot(data = df_plot_deat, aes(x = time, y = value,
                                      color = as.numeric(as.character(variable)),
                                      group = variable)) +
  geom_line() +
  scale_colour_gradient(name = ""*theta~","*delta~", d, "*alpha~"", 
                        low = "blue", high = "red",
                        breaks=c(0,0.5,1),
                        labels=c(0,0.5,1),
                        limits=c(0,1)) +
  ylab("Sum of infected individuals") +
  scale_y_continuous(labels = function(x)
    format(x, scientific = TRUE), limits = c(0, 300000)) +
  theme_bw() +
  theme(text = element_text(size = size_let)) 

mud <- 0.3
deaths <- rep(mud, N) # not disease-related death rates

##### Natural mortility #####
alphas_vec <-  seq(min_int,max_int,step) # recovery rates
end_time <- 50
alphas <- rep(alphas_vec[1],N)
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)
df_sum_alp <- data.frame(time = sol[,1])
for(i in c(1:length(alphas_vec))){
  print(paste0("i: ", i))
  alphas <- rep(alphas_vec[i],N)
  sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
             COMMUTING,MIGRATION,
             sus_init,inf_init,end_time)
  sol_inf <- sol[,c(1,(N+2):(2*N+1))]
  sum_inf <- rowSums(sol_inf[,2:(N+1)])
  df_sum_alp[,ncol(df_sum_alp)+1] <- sum_inf
}

colnames(df_sum_alp) <- c("time",as.character(alphas_vec))
df_plot_alp <- reshape2::melt(df_sum_alp, id.vars = c("time"))

size_let <- 13

sum_alphas <- ggplot(data = df_plot_alp, aes(x = time, y = value,
                                      color = as.numeric(as.character(variable)),
                                      group = variable)) +
  geom_line() +
  scale_colour_gradient(name = ""*theta~","*delta~", d, "*alpha~"", 
                        low = "blue", high = "red",
                        breaks=c(0,0.5,1),
                        labels=c(0,0.5,1),
                        limits=c(0,1)) +
  ylab("Sum of infected individuals") +
  scale_y_continuous(labels = function(x)
    format(x, scientific = TRUE)) +
  theme_bw() +
  theme(text = element_text(size = size_let), legend.position = "bottom") 



library("ggpubr")
library("latex2exp")
gg_param <- ggarrange(sum_tetha + ggtitle(TeX("$\\theta$")) + 
                        ylab("") + xlab(""),
                      sum_deltas + ggtitle(TeX("$\\delta$")) + 
                        ylab("") + xlab(""),
                      sum_deaths + ggtitle(TeX("$d$")) + 
                        ylab("") + xlab(""),
                      sum_alphas + ggtitle(TeX("$\\alpha$")) + 
                        ylab("") + xlab(""),
                      common.legend = TRUE,
                      ncol = 2, nrow = 2)

gg_param <- annotate_figure(gg_param,
                bottom = text_grob("Time", color = "black",
                                   size = 15),
                left = text_grob("Sum of infected individuals",
                                 color = "black", rot = 90,  size = 15))
gg_param

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/"
path <- paste0(Path,"Plot_epiparam1_b0,1_g0,5_muc_0,001_sc0,00001_muw0,2_sw0,05.png")
ggsave(path,
       plot = gg_param, device = "png")


###### RAMDOM PARAM #########
# Thetas ####
Deltas <- rep(0.3, N) # birth rate
mub <- 0.1
sb <- 0.001
betas <- rep(mub, N) # transmission rates
# betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.3, N) # loss of immunity rates
mud <- 0.3
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.3
alphas <- rep(mua, N) # recovery rates
mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates
gammas = deaths + alphas + deltas

mut <- 1
st <- 0.1
thetas <- rgamma(N, shape = (mut/st)^2, rate = mut/(st^2))
end_time <- 50
sol.rand <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

for(i in c(1:N)){
  colnames(sol.rand)[i+1] <-  paste0("S",i)
  colnames(sol.rand)[N+i+1] <-  paste0("I",i)
  colnames(sol.rand)[2*N+i+1] <-  paste0("R",i)
}

sol.rand <- as.data.frame(sol.rand)
df_plot <- reshape2::melt(sol.rand, id.vars = c("time"))
df_plot$type <- substr(df_plot$variable,1,1)
df_inf <- df_plot  %>% 
  filter( substr(df_plot$variable,1,1) == "I")
df_inf$thet <- 0
for(i in c(1:nrow(df_inf))){
  df_inf$thet[i] <- thetas[as.numeric(substr(df_plot$variable[i],2,3))]
}

max_int <- ceiling(max(thetas))
min_int <- floor(min(thetas))

library("latex2exp")
plot_stab.rand.hm.lv  <- ggplot(df_inf) + 
  geom_line(aes(time, value, colour = thet, group = variable))  +
  # geom_line(aes( colour =variable),size=0.5)  +
  ylab("Number of infected individuals")+
  scale_colour_gradient(name = ""*theta~" ", 
                        low = "blue", high = "red",
                        breaks=c(min_int,1,max_int),
                        labels=c(min_int,1,max_int),
                        limits=c(min_int,max_int)) + 
  theme_bw() + 
  ggtitle(TeX("$\\mu_{\\theta}: 1, \\sigma_{\\theta}:0.1$"))

plot_stab.rand.hm.lv

gg_thetas <- ggarrange(plot_stab.rand.hm + xlab(""),
          plot_stab.rand.lm + ylab(""),
          plot_stab.rand.hm.lv ,
          plot_stab.rand.lm.lv + ylab(""),
          common.legend = TRUE,
          nrow = 2, ncol = 2)

# Deaths ####
Deltas <- rep(0.3, N) # birth rate
mub <- 0.1
sb <- 0.001
betas <- rep(mub, N) # transmission rates
# betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.3, N) # loss of immunity rates
mud <- 0.3
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.3
alphas <- rep(mua, N) # recovery rates
mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates
gammas = deaths + alphas + deltas

mut <- 1
st <- 0.5
deaths <- rgamma(N, shape = (mut/st)^2, rate = mut/(st^2))
end_time <-200
sol.rand <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
                COMMUTING,MIGRATION,
                sus_init,inf_init,end_time)

plot_int(sol.rand, state = "INF")
for(i in c(1:N)){
  colnames(sol.rand)[i+1] <-  paste0("S",i)
  colnames(sol.rand)[N+i+1] <-  paste0("I",i)
  colnames(sol.rand)[2*N+i+1] <-  paste0("R",i)
}

sol.rand <- as.data.frame(sol.rand)
df_plot <- reshape2::melt(sol.rand, id.vars = c("time"))
df_plot$type <- substr(df_plot$variable,1,1)
df_inf <- df_plot  %>% 
  filter( substr(df_plot$variable,1,1) == "I")
df_inf$deaths <- 0
for(i in c(1:nrow(df_inf))){
  df_inf$deaths[i] <- deaths[as.numeric(substr(df_plot$variable[i],2,3))]
}

max_int <- ceiling(max(deaths))
min_int <- floor(min(deaths))

library("latex2exp")
plot_stab.rand.hm.hv_d  <- ggplot(df_inf) + 
  geom_line(aes(time, value, colour = deaths, group = variable))  +
  # geom_line(aes( colour =variable),size=0.5)  +
  ylab("Number of infected individuals")+
  scale_colour_gradient(name = ""*d~" ", 
                        low = "blue", high = "red",
                        breaks=c(min_int,1,max_int),
                        labels=c(min_int,1,max_int),
                        limits=c(min_int,max_int)) + 
  theme_bw() + 
  ggtitle(TeX("$\\mu_{d}: 1, \\sigma_{d}:0.5$"))

plot_stab.rand.hm.lv_d

gg_deaths <- ggarrange(plot_stab.rand.lm.lv_d +
                         xlab("") + ylim(c(0,20000)) + xlim(c(0,10)),
                       plot_stab.rand.lm.hv_d + 
                         ylab("") + ylim(c(0,20000)) + xlim(c(0,10)),
                       plot_stab.rand.hm.lv_d +
                         ylim(c(0,200)),
                       plot_stab.rand.hm.hv_d +
                         ylab("") + ylim(c(0,200)),
                       common.legend = TRUE,
                       nrow = 2, ncol = 2)
gg_deaths

# Induced deaths ####
Deltas <- rep(0.3, N) # birth rate
mub <- 0.1
sb <- 0.001
betas <- rep(mub, N) # transmission rates
# betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.3, N) # loss of immunity rates
mud <- 0.3
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.3
alphas <- rep(mua, N) # recovery rates
mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates
gammas = deaths + alphas + deltas

mut <- 0.1
st <- 0.5
deltas <- rgamma(N, shape = (mut/st)^2, rate = mut/(st^2))
end_time <- 50
sol.rand <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
                COMMUTING,MIGRATION,
                sus_init,inf_init,end_time)

for(i in c(1:N)){
  colnames(sol.rand)[i+1] <- paste0("S",i)
  colnames(sol.rand)[N+i+1] <- paste0("I",i)
  colnames(sol.rand)[2*N+i+1] <- paste0("R",i)
}

sol.rand <- as.data.frame(sol.rand)
df_plot <- reshape2::melt(sol.rand, id.vars = c("time"))
df_plot$type <- substr(df_plot$variable,1,1)
df_inf <- df_plot  %>% 
  filter( substr(df_plot$variable,1,1) == "I")
df_inf$delta <- 0
for(i in c(1:nrow(df_inf))){
  df_inf$delta[i] <- deltas[as.numeric(substr(df_plot$variable[i],2,3))]
}

max_int <- ceiling(max(deltas))
min_int <- floor(min(deltas))

library("latex2exp")
plot_stab.rand.lm.hv_delt  <- ggplot(df_inf) + 
  geom_line(aes(time, value, colour = delta, group = variable))  +
  # geom_line(aes( colour =variable),size=0.5)  +
  ylab("Number of infected individuals")+
  scale_colour_gradient(name = ""*deltas~" ", 
                        low = "blue", high = "red",
                        breaks=c(min_int,1,max_int),
                        labels=c(min_int,1,max_int),
                        limits=c(min_int,max_int)) + 
  theme_bw() + 
  ggtitle(TeX("$\\mu_{\\delta}: 0.1, \\sigma_{\\delta}:0.5$"))

plot_stab.rand.hm.lv_d

gg_deaths <- ggarrange(plot_stab.rand.lm.lv_d +
                         xlab("") + ylim(c(0,20000)) + xlim(c(0,10)),
                       plot_stab.rand.lm.hv_d + 
                         ylab("") + ylim(c(0,20000)) + xlim(c(0,10)),
                       plot_stab.rand.hm.lv_d +
                         ylim(c(0,200)),
                       plot_stab.rand.hm.hv_d +
                         ylab("") + ylim(c(0,200)),
                       common.legend = TRUE,
                       nrow = 2, ncol = 2)
gg_deaths