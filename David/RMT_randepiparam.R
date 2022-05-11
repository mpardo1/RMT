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

vec_mu <- c(0.1,1)
vec_s <- c(0.1,0.5)
list_plot <- list()
count <- 1
for(k in 1:2){
  for(j in 1:2){
    mut <- vec_mu[k]
    st <- vec_s[j]
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

    print(paste0("count:",count))
    list_plot[[count]]  <- ggplot(df_inf) +
      geom_line(aes(time, value, colour = thet, group = variable))  +
      # geom_line(aes( colour =variable),size=0.5)  +
      ylab("Number of infected individuals")+
      scale_colour_gradient(name = ""*theta~" ",
                            low = "blue", high = "red",
                            breaks=c(min_int,1,max_int),
                            labels=c(min_int,1,max_int),
                            limits=c(min_int,max_int)) +
      theme_bw() +
      ggtitle(TeX(sprintf('$\\mu_{\\theta }= %g, \\sigma_{\\theta }= %g$', mut,st)))
    count = count + 1
  }
}

library("ggpubr")
gg_thetas <- ggarrange(list_plot[[1]] + xlab("") + ylab(""),
                       list_plot[[2]] + xlab("") + ylab(""),
                       list_plot[[3]] + xlab("") + ylab(""),
                       list_plot[[4]] + xlab("") + ylab(""),
                       common.legend = TRUE,
                       nrow = 2, ncol = 2)

gg_thetas <- annotate_figure(gg_thetas,
                             bottom = text_grob("Time", color = "black",
                                                size = 15),
                             left = text_grob("Sum of infected individuals",
                                              color = "black", rot = 90,  size = 15))

gg_thetas

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/"
path <- paste0(Path,"rand_param_theta_muw",muw,"sw",sw,
               "muc",muc,"b",betas[1],"d",deltas[1],
               "D",Deltas[1],"a",alphas[1],"t",thetas[1],".png")
ggsave(path,
       plot = gg_thetas, device = "png")

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

vec_mu <- c(0.1,1)
vec_s <- c(0.1,0.5)
list_plot_d <- list()
count <- 1
for(k in 1:2){
  for(j in 1:2){
    mut <- vec_mu[k]
    st <- vec_s[j]
    deaths <- rgamma(N, shape = (mut/st)^2, rate = mut/(st^2))
    end_time <- 200
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
    df_inf$deaths <- 0
    for(i in c(1:nrow(df_inf))){
      df_inf$deaths[i] <- deaths[as.numeric(substr(df_plot$variable[i],2,3))]
    }
    
    max_int <- ceiling(max(deaths))
    min_int <- floor(min(deaths))
    
    library("latex2exp")
    list_plot_d[[count]]   <- ggplot(df_inf) + 
      geom_line(aes(time, value, colour = deaths, group = variable))  +
      # geom_line(aes( colour =variable),size=0.5)  +
      ylab("Number of infected individuals")+
      scale_colour_gradient(name = ""*d~" ", 
                            low = "blue", high = "red",
                            breaks=c(min_int,1,max_int),
                            labels=c(min_int,1,max_int),
                            limits=c(min_int,max_int)) + 
      theme_bw() + 
      ggtitle(TeX(sprintf('$\\mu_{d }= %g, \\sigma_{d }= %g$', mut,st)))
      count = count + 1
  }
}

gg_deaths <- ggarrange(list_plot_d[[1]] +
                         xlab("") +  ylab("") + 
                         ylim(c(0,20000)) + xlim(c(0,10)),
                       list_plot_d[[2]] + 
                         ylab("") + xlab("") + 
                         ylim(c(0,20000)) + xlim(c(0,10)),
                       list_plot_d[[3]] + xlab("") + ylab("") + 
                         xlim(c(0,10)),
                       list_plot_d[[4]] +
                         ylab("") + xlab("") + xlim(c(0,10)) ,
                       common.legend = TRUE,
                       nrow = 2, ncol = 2)

gg_deaths <- annotate_figure(gg_deaths,
                             bottom = text_grob("Time", color = "black",
                                                size = 15),
                             left = text_grob("Sum of infected individuals",
                                              color = "black", rot = 90,  size = 15))
gg_deaths

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/"
path <- paste0(Path,"rand_param_death_muw",muw,"sw",sw,
               "muc",muc,"b",betas[1],"d",deltas[1],
               "D",Deltas[1],"a",alphas[1],"t",thetas[1],".png")
ggsave(path,
       plot = gg_deaths, device = "png")

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


vec_mu <- c(0.1,1)
vec_s <- c(0.1,0.5)
list_plot_del <- list()
count <- 1
for(k in 1:2){
  for(j in 1:2){
    mut <- vec_mu[k]
    st <- vec_s[j]
    deltas <- rgamma(N, shape = (mut/st)^2, rate = mut/(st^2))
    end_time <- 500
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
    
    print(paste0("count:",count))
    list_plot_del[[count]]  <- ggplot(df_inf) + 
      geom_line(aes(time, value, colour = delta, group = variable))  +
      # geom_line(aes( colour =variable),size=0.5)  +
      ylab("Number of infected individuals")+
      scale_colour_gradient(name = ""*delta~" ", 
                            low = "blue", high = "red",
                            breaks=c(min_int,1,max_int),
                            labels=c(min_int,1,max_int),
                            limits=c(min_int,max_int)) + 
      theme_bw() + 
      ggtitle(TeX(sprintf('$\\mu_{\\delta}= %g, \\sigma_{\\delta}= %g$', mut,st)))
      count = count +1
  }
}


gg_deltas <- ggarrange(list_plot_del[[1]] +
                         xlab("") + ylab("") + xlim(c(0,100)),
                       list_plot_del[[2]] + 
                         xlab("") + ylab("") + xlim(c(0,100)),
                       list_plot_del[[3]] +
                         xlab("") + ylab("") + xlim(c(0,100)),
                       list_plot_del[[4]] +
                         xlab("") + ylab("") + xlim(c(0,100)) ,
                       common.legend = TRUE,
                       nrow = 2, ncol = 2)

gg_deltas <- annotate_figure(gg_deltas,
                             bottom = text_grob("Time", color = "black",
                                                size = 15),
                             left = text_grob("Sum of infected individuals",
                                              color = "black", rot = 90,  size = 15))
gg_deltas

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/"
path <- paste0(Path,"rand_param_delt_muw",muw,"sw",sw,
               "muc",muc,"b",betas[1],"d",deltas[1],
               "D",Deltas[1],"a",alphas[1],"t",thetas[1],".png")
ggsave(path,
       plot = gg_deltas, device = "png")

# alphas ####
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

vec_mu <- c(0.1,1)
vec_s <- c(0.1,0.5)
list_plot_al <- list()
count <- 1
for(k in 1:2){
  for(j in 1:2){
    mut <- vec_mu[k]
    st <- vec_s[j]
  alphas <- rgamma(N, shape = (mut/st)^2, rate = mut/(st^2))
  end_time <- 500
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
  df_inf$alphas <- 0
  for(i in c(1:nrow(df_inf))){
    df_inf$alphas[i] <- alphas[as.numeric(substr(df_plot$variable[i],2,3))]
  }
  
  max_int <- ceiling(max(alphas))
  min_int <- floor(min(alphas))
  
  library("latex2exp")
  
  print(paste0("count:",count))
  list_plot_al [[count]] <- ggplot(df_inf) + 
    geom_line(aes(time, value, colour = alphas, group = variable))  +
    # geom_line(aes( colour =variable),size=0.5)  +
    ylab("Number of infected individuals")+
    scale_colour_gradient(name = ""*alpha~" ", 
                          low = "blue", high = "red",
                          breaks=c(min_int,1,max_int),
                          labels=c(min_int,1,max_int),
                          limits=c(min_int,max_int)) + 
    theme_bw() + 
    ggtitle(TeX(sprintf('$\\mu_{\\alpha}= %g, \\sigma_{\\alpha}= %g$', mut,st)))
  count = count +1
  }
}

gg_alphas <- ggarrange(list_plot_al[[1]] +
                         xlab("") + ylab("") + xlim(c(0,100)),
                       list_plot_al[[2]] + 
                         xlab("") + ylab("") + xlim(c(0,100)),
                       list_plot_al[[3]] +
                         xlab("") + ylab("") + xlim(c(0,100)),
                       list_plot_al[[4]] +
                         xlab("") + ylab("") + xlim(c(0,100)) ,
                       common.legend = TRUE,
                       nrow = 2, ncol = 2)

gg_alphas <- annotate_figure(gg_alphas,
                             bottom = text_grob("Time", color = "black",
                                                size = 15),
                             left = text_grob("Sum of infected individuals",
                                              color = "black", rot = 90,  size = 15))

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/"
path <- paste0(Path,"rand_param_alp_muw",muw,"sw",sw,
               "muc",muc,"b",betas[1],"d",deltas[1],
               "D",Deltas[1],"a",alphas[1],"t",thetas[1],".png")
ggsave(path,
       plot = gg_alphas, device = "png")


ggarrange(gg_alphas,gg_deltas,gg_deaths,gg_thetas)

# ERROR OUTLIER PRED ####
# Alpha
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

mum = 1
sm = 0.5
mu_vec <- rgamma(100, shape = (mum/sm)^2, rate = mum/(sm^2))
s_vec <- seq(0,1,0.01)
df_err <- data.frame(muw = muw, sw = sw, muc = muc, sc = sc, 
                     mual = 0, sal = 0, max_eig = 0, out_mean = 0 ,
                     count_re = 0)

count = 1
while(count < 1000){
  # for(i in c(1:length(mu_vec))){
    for(j in c(1:length(s_vec))){
      mut <- 0.5
      st <- s_vec[j]
      alphas <- rgamma(N, shape = (mut/st)^2, rate = mut/(st^2))
      gammas = deaths + alphas + deltas
      # Compute Jacobian
      jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
        diag(deaths + alphas + deltas + colSums(MIGRATION))
      # Compute eigenvalues
      eig <- eigen_mat(jacobian)
      max_eig <- max(eig$re)
      mub <- betas[1]
      mug <- mean(gammas)
      out_mean <- mub - mug + mub*muw*(N-1)
      df_err[nrow(df_err)+1,] <- c(muw, sw, muc, sc, mut,
                                   st, max_eig, out_mean , count) 
    # }
  }
  print(paste0("count:", count))
  count = count + 1 
}

df_err$err_mean <- (df_err$out_mean - df_err$max_eig)^2/df_err$out_mean
df_err_g <- df_err %>% group_by(sal) %>%
  summarise(mean_err_mean = mean(err_mean))
df_err_g <- df_err_g[-1,]

library("latex2exp")
ggplot(df_err_g) + 
  geom_line(aes(sal, mean_err_mean), color ="#A40E4C", size = 0.4) + 
  geom_point(aes(sal, mean_err_mean), color ="#2C2C54", size = 0.9 ) + 
  xlim(c(0,0.8)) + ylim(c(-1.5,1)) +
  theme_bw() + xlab(TeX("$\\sigma_{\\alpha}$")) + ylab("Mean squared error")

# Deaths
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

mum = 1
sm = 0.5
mu_vec <- rgamma(100, shape = (mum/sm)^2, rate = mum/(sm^2))
s_vec <- seq(0,1,0.01)
df_err <- data.frame(muw = muw, sw = sw, muc = muc, sc = sc, 
                     mual = 0, sal = 0, max_eig = 0, out_mean = 0 ,
                     count_re = 0)

count = 1
while(count < 1000){
  # for(i in c(1:length(mu_vec))){
  for(j in c(1:length(s_vec))){
    mut <- 0.5
    st <- s_vec[j]
    deaths <- rgamma(N, shape = (mut/st)^2, rate = mut/(st^2))
    gammas = deaths + alphas + deltas
    # Compute Jacobian
    jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      diag(deaths + alphas + deltas + colSums(MIGRATION))
    # Compute eigenvalues
    eig <- eigen_mat(jacobian)
    max_eig <- max(eig$re)
    mub <- betas[1]
    mug <- mean(gammas)
    out_mean <- mub - mug + mub*muw*(N-1)
    df_err[nrow(df_err)+1,] <- c(muw, sw, muc, sc, mut,
                                 st, max_eig, out_mean , count) 
    # }
  }
  print(paste0("count:", count))
  count = count + 1 
}

df_err$err_mean <- (df_err$out_mean - df_err$max_eig)^2/df_err$out_mean
df_err_g <- df_err %>% group_by(sal) %>%
  summarise(mean_err_mean = mean(err_mean))
df_err_g <- df_err_g[-1,]

library("latex2exp")
gg_deaths <- ggplot(df_err_g) + 
  geom_line(aes(sal, mean_err_mean), color ="#A40E4C", size = 0.4) + 
  geom_point(aes(sal, mean_err_mean), color ="#2C2C54", size = 0.9 ) + 
  xlim(c(0,0.8)) + ylim(c(-0.2,0.5)) +
  theme_bw() + xlab(TeX("$\\sigma_{d}$")) + ylab("Mean squared error")

Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/Gen/"
path <- paste0(Path,"rand_deaths",format(muw,decimal.mark=","),
               "sw",format(sw,decimal.mark=","),
               "muc",format(muc,decimal.mark=","),
               "b",format(betas[1],decimal.mark=","),
               "d",format(deltas[1],decimal.mark=","), 
               "D",format(Deltas[1],decimal.mark=","),
               "a",format(alphas[1],decimal.mark=","),
               "t",format(thetas[1],decimal.mark=","),".png")
ggsave(path,
       plot = gg_deaths, device = "png")

# Deltas
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

mum = 1
sm = 0.5
mu_vec <- rgamma(100, shape = (mum/sm)^2, rate = mum/(sm^2))
s_vec <- seq(0,1,0.01)
df_err <- data.frame(muw = muw, sw = sw, muc = muc, sc = sc, 
                     mual = 0, sal = 0, max_eig = 0, out_mean = 0 ,
                     count_re = 0)

count = 1
while(count < 1000){
  # for(i in c(1:length(mu_vec))){
  for(j in c(1:length(s_vec))){
    mut <- 0.5
    st <- s_vec[j]
    deltas <- rgamma(N, shape = (mut/st)^2, rate = mut/(st^2))
    gammas = deaths + alphas + deltas
    # Compute Jacobian
    jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      diag(deaths + alphas + deltas + colSums(MIGRATION))
    # Compute eigenvalues
    eig <- eigen_mat(jacobian)
    max_eig <- max(eig$re)
    mub <- betas[1]
    mug <- mean(gammas)
    out_mean <- mub - mug + mub*muw*(N-1)
    df_err[nrow(df_err)+1,] <- c(muw, sw, muc, sc, mut,
                                 st, max_eig, out_mean , count) 
    # }
  }
  print(paste0("count:", count))
  count = count + 1 
}

df_err$err_mean <- (df_err$out_mean - df_err$max_eig)^2/df_err$out_mean
df_err_g <- df_err %>% group_by(sal) %>%
  summarise(mean_err_mean = mean(err_mean))
df_err_g <- df_err_g[-1,]

library("latex2exp")
gg_deltas <- ggplot(df_err_g) + 
  geom_line(aes(sal, mean_err_mean), color ="#A40E4C", size = 0.4) + 
  geom_point(aes(sal, mean_err_mean), color ="#2C2C54", size = 0.9 ) + 
  xlim(c(0,0.75)) + 
  ylim(c(-5,5)) +
  theme_bw() + xlab(TeX("$\\sigma_{\\delta}$")) + ylab("Mean squared error")

Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/Gen/"
path <- paste0(Path,"rand_deltas",format(muw,decimal.mark=","),
               "sw",format(sw,decimal.mark=","),
               "muc",format(muc,decimal.mark=","),
               "b",format(betas[1],decimal.mark=","),
               "d",format(deltas[1],decimal.mark=","), 
               "D",format(Deltas[1],decimal.mark=","),
               "a",format(alphas[1],decimal.mark=","),
               "t",format(thetas[1],decimal.mark=","),".png")
ggsave(path,
       plot = gg_deltas, device = "png")


