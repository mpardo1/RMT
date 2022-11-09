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
library("ggpubr")
library("ggforce")
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

# Gammas
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
df_err_gam <- data.frame(muw = muw, sw = sw, muc = muc, sc = sc, 
                         mual = 0, sal = 0, max_eig = 0, out_mean = 0 ,
                         count_re = 0)

count = 1
while(count < 100){
  # for(i in c(1:length(mu_vec))){
  for(j in c(1:length(s_vec))){
    mut <- 0.5
    st <- s_vec[j]
    gammas = rgamma(N, shape = (mut/st)^2, rate = mut/(st^2))
    # Compute Jacobian
    jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      diag(gammas + colSums(MIGRATION))
    # Compute eigenvalues
    eig <- eigen_mat(jacobian)
    max_eig <- max(eig$re)
    gammas = mean(gammas)
    # Compute Jacobian
    jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      diag(gammas + colSums(MIGRATION))
    # Compute eigenvalues
    eig <- eigen_mat(jacobian)
    max_eig_mean <- max(eig$re)
    df_err_gam[nrow(df_err_gam)+1,] <- c(muw, sw, muc, sc, mut,
                                     st, max_eig, max_eig_mean , count) 
    # }
  }
  print(paste0("count:", count))
  count = count + 1 
}

df_err_gam$err_mean <- (df_err_gam$out_mean - df_err_gam$max_eig)^2/df_err_gam$out_mean
df_err_g_gam <- df_err_gam %>% group_by(sal) %>%
  summarise(mean_err_mean = mean(err_mean))
df_err_g_gam <- df_err_g_gam[-1,]

library("latex2exp")
gg_gammas <- ggplot(df_err_g_gam) + 
  geom_line(aes(sal, mean_err_mean), color ="#3293CE", size = 0.4) + 
  geom_point(aes(sal, mean_err_mean), color ="#2C2C54", size = 0.9 ) + 
  theme_bw() + xlab(TeX("$\\sigma_{\\gamma}$")) + 
  ylab("Relative error") +
  theme(text = element_text(size = 15), legend.position = "bottom") 
gg_gammas

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/"
path <- paste0(Path,"rand_gammas","N", N,
               "muw",format(muw,decimal.mark=","),
               "sw",format(sw,decimal.mark=","),
               "muc",format(muc,decimal.mark=","),
               "b",format(betas[1],decimal.mark=","),
               "d",format(deltas[1],decimal.mark=","), 
               "D",format(Deltas[1],decimal.mark=","),
               "a",format(alphas[1],decimal.mark=","),
               "t",format(thetas[1],decimal.mark=","),".pdf")
ggsave(path,
       plot = gg_gammas, device = "pdf")



# Betas
N = 200
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

MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
diag(MIGRATION) <- 0

mum = 1
sm = 0.5
mu_vec <- rgamma(101, shape = (mum/sm)^2, rate = mum/(sm^2))
s_vec <- seq(0,1,0.01)
df_err_bet <- data.frame(muw = muw, sw = sw, muc = muc, sc = sc, 
                         mual = 0, sal = 0, max_eig = 0, out_mean = 0 ,
                         count_re = 0)

count = 1
while(count < 50){
  for(j in c(1:length(s_vec))){
    mut <- mu_vec[j]
    st <- s_vec[j]
    if(mut <= st){
      st <- 0.2
      mut <- 0.5
    }
    print(paste0("j:", j))
    betas = rgamma(N, shape = (mut/st)^2, rate = mut/(st^2))
    # Compute Jacobian
    jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      diag(gammas + colSums(MIGRATION))
    # Compute eigenvalues
    eig <- eigen_mat(jacobian)
    max_eig <- max(eig$re)
    betas = rep(mean(betas),N)
    # Compute Jacobian
    jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      diag(gammas + colSums(MIGRATION))
    # Compute eigenvalues
    eig <- eigen_mat(jacobian)
    max_eig_mean <- max(eig$re)
    df_err_bet[nrow(df_err_bet)+1,] <- c(muw, sw, muc, sc, mut,
                                     st, max_eig, max_eig_mean , count) 
  }
  print(paste0("count:", count))
  count = count + 1 
}

df_err_bet$err_mean <- (df_err_bet$out_mean - df_err_bet$max_eig)^2/df_err_bet$out_mean
df_err_bet$var_bet <- (df_err_bet$sal)/(df_err_bet$mual)
df_err_bet_g <- df_err_bet %>% group_by(var_bet) %>%
  summarise(mean_err_mean = mean(err_mean))
df_err_bet_g <- df_err_bet_g[-1,]

library("latex2exp")
x_inter <- 0.7650498
gg_betas <- ggplot(df_err_bet_g) + 
  geom_line(aes(var_bet, mean_err_mean), color ="#3293CE", size = 0.4) + 
  geom_point(aes(var_bet, mean_err_mean), color ="#2C2C54", size = 0.9 ) + 
  theme_bw() + xlab(TeX("CV")) + 
  ylab("Relative error")  +
  # geom_vline(xintercept = x_inter, color = "blue", linetype = "longdash") +
  theme(text = element_text(size = 15), legend.position = "bottom",
        plot.margin = margin(1, 1, 1, 1, "cm")) 
gg_betas

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/"
path <- paste0(Path,"rand_betas_no_dash","N", N,
               "muw",format(muw,decimal.mark=","),
               "sw",format(sw,decimal.mark=","),
               "muc",format(muc,decimal.mark=","),
               "b",format(betas[1],decimal.mark=","),
               "d",format(deltas[1],decimal.mark=","), 
               "D",format(Deltas[1],decimal.mark=","),
               "a",format(alphas[1],decimal.mark=","),
               "t",format(thetas[1],decimal.mark=","),".pdf")
ggsave(path, plot = gg_betas, device = "pdf")


####CONVERGENCE TEST###

## BETA #
# Test whether the mean square error for the predicted outlier
# decrease when increasing the size of the network
# The j position for the sigma vector.
x_inter <- 0.75
j <- which(df_err_bet$var_bet <= (x_inter + 0.001) & 
             df_err_bet$var_bet >= (x_inter - 0.001))[1]

mu <- df_err_bet$mual[j]
sig <- df_err_bet$sal[j]
count = 1
df_err_bet_N <- data.frame(N = 0, err = 0)
while(count < 25){
  for(i in c(50:300)){
    N = i
    COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
    diag(COMMUTING) <- 0
    
    MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
    diag(MIGRATION) <- 0
    
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
    
    mut <- mu
    st <- sig
    betas = rgamma(N, shape = (mut/st)^2, rate = mut/(st^2))
    # Compute Jacobian
    jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      diag(gammas + colSums(MIGRATION))
    # Compute eigenvalues
    eig <- eigen_mat(jacobian)
    max_eig <- max(eig$re)
    betas = rep(mean(betas),N)
    # Compute Jacobian
    jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      diag(gammas + colSums(MIGRATION))
    # Compute eigenvalues
    eig <- eigen_mat(jacobian)
    max_eig_mean <- max(eig$re)
    
    err <- (max_eig_mean - max_eig)^2/max_eig_mean
    df_err_bet_N[nrow(df_err_bet_N)+1,] <- c(i, err)
  }
  count = count + 1
  print(paste0("count:", count))
}

df_err_bet_N_g <- df_err_bet_N %>% group_by(N) %>%
  summarise(mean_err_mean = mean(err))
df_err_bet_N_g <- df_err_bet_N_g[-1,]

library("latex2exp")
gg_betas_err_N <- ggplot(df_err_bet_N_g) + 
  geom_line(aes(N, mean_err_mean), color ="#3293CE", size = 0.4) + 
  geom_point(aes(N,mean_err_mean), color ="#2C2C54", size = 0.9 ) + 
  theme_bw() +
  xlab(TeX("$N$")) + 
  ylab("Relative error") +
  geom_vline(xintercept = 200, color = "blue", linetype = "longdash") +
  theme(text = element_text(size = 15), legend.position = "bottom") 

gg_betas_err_N

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/"
path <- paste0(Path,"rand_bet_N_",
               "muw",format(muw,decimal.mark=","),
               "sw",format(sw,decimal.mark=","),
               "muc",format(muc,decimal.mark=","),
               "b",format(betas[1],decimal.mark=","),
               "d",format(deltas[1],decimal.mark=","), 
               "D",format(Deltas[1],decimal.mark=","),
               "a",format(alphas[1],decimal.mark=","),
               "t",format(thetas[1],decimal.mark=","),".pdf")
ggsave(path, plot = gg_betas_err_N, device = "pdf")

gg_arr <- ggarrange(gg_betas +
          geom_vline(xintercept = 0.75, color = "blue", linetype = "longdash")  +
            ggtitle(TeX("$N = 200$")),
          gg_betas_err_N + ylab("") + 
            ggtitle(TeX("$CV = 0.75$")),
          labels = c("a", "b"))

gg_arr

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/"
path <- paste0(Path,"arr_bet_rand",
               "muw",format(muw,decimal.mark=","),
               "sw",format(sw,decimal.mark=","),
               "muc",format(muc,decimal.mark=","),
               "b",format(betas[1],decimal.mark=","),
               "d",format(deltas[1],decimal.mark=","), 
               "D",format(Deltas[1],decimal.mark=","),
               "a",format(alphas[1],decimal.mark=","),
               "t",format(thetas[1],decimal.mark=","),".pdf")
ggsave(path, plot = gg_arr, device = "pdf")

## GAMMA #
# Test whether the mean square error for the predicted outlier
# decrease when increasing the size of the network
# The j position for the sigma vector.
j <- 100
s_vec[j]
count = 1
df_err_gam_N <- data.frame(N = 0, err = 0)
while(count < 50){
  for(i in c(50:250)){
    N = i
    COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
    diag(COMMUTING) <- 0
    
    MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
    diag(MIGRATION) <- 0
    
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
    
    mut <- 0.5
    st <- s_vec[j]
    gammas = rgamma(N, shape = (mut/st)^2, rate = mut/(st^2))
    # Compute Jacobian
    jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      diag(gammas + colSums(MIGRATION))
    # Compute eigenvalues
    eig <- eigen_mat(jacobian)
    max_eig <- max(eig$re)
    gammas = rep(mean(gammas),N)
    # Compute Jacobian
    jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      diag(gammas + colSums(MIGRATION))
    # Compute eigenvalues
    eig <- eigen_mat(jacobian)
    max_eig_mean <- max(eig$re)
    
    err <- (max_eig_mean - max_eig)^2/max_eig_mean
    df_err_gam_N[nrow(df_err)+1,] <- c(i, err)
  }
  count = count + 1
  print(paste0("count:", count))
}

df_err_g_N <- df_err_gam_N %>% group_by(N) %>%
  summarise(mean_err_mean = mean(err))
df_err_g_N <- df_err_g_N[-1,]

library("latex2exp")
gg_gammas_err_N <- ggplot(df_err_g_N) + 
  geom_line(aes(N, mean_err_mean), color ="#3293CE", size = 0.4) + 
  geom_point(aes(N,mean_err_mean), color ="#2C2C54", size = 0.9 ) + 
  theme_bw() +
  xlab(TeX("$N$")) + 
  ylab("Relative error")

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/"
path <- paste0(Path,"rand_gamm_N_",
               "muw",format(muw,decimal.mark=","),
               "sw",format(sw,decimal.mark=","),
               "muc",format(muc,decimal.mark=","),
               "b",format(betas[1],decimal.mark=","),
               "d",format(deltas[1],decimal.mark=","), 
               "D",format(Deltas[1],decimal.mark=","),
               "a",format(alphas[1],decimal.mark=","),
               "t",format(thetas[1],decimal.mark=","),".pdf")
ggsave(path, plot = gg_gammas_err_N, device = "pdf")

####TEST Sigma or mu importance####
## Sigma
mut = 1
sm = 0.5
# mu_vec <- rgamma(101, shape = (mum/sm)^2, rate = mum/(sm^2))
s_vec <- seq(0,0.9,0.01)
df_err_bet_fixed_fm <- data.frame(muw = muw, sw = sw, muc = muc, sc = sc, 
                         mual = 0, sal = 0, max_eig = 0, out_mean = 0 ,
                         count_re = 0)

count = 1
while(count < 50){
  for(j in c(1:length(s_vec))){
    st <- s_vec[j]
    if(mut <= st){
      st <- 0.2
      mut <- 0.5
    }
    print(paste0("j:", j))
    betas = rgamma(N, shape = (mut/st)^2, rate = mut/(st^2))
    # Compute Jacobian
    jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      diag(gammas + colSums(MIGRATION))
    # Compute eigenvalues
    eig <- eigen_mat(jacobian)
    max_eig <- max(eig$re)
    betas = rep(mean(betas),N)
    # Compute Jacobian
    jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      diag(gammas + colSums(MIGRATION))
    # Compute eigenvalues
    eig <- eigen_mat(jacobian)
    max_eig_mean <- max(eig$re)
    df_err_bet_fixed_fm[nrow(df_err_bet_fixed_fm)+1,] <- c(muw, sw, muc, sc, mut,
                                         st, max_eig, max_eig_mean , count) 
  }
  print(paste0("count:", count))
  count = count + 1 
}

df_err_bet_fixed_fm$err_mean <- (df_err_bet_fixed_fm$out_mean - df_err_bet_fixed_fm$max_eig)^2/df_err_bet_fixed_fm$out_mean
df_err_bet_fixed_fm$var_bet <- (df_err_bet_fixed_fm$sal)/(df_err_bet_fixed_fm$mual)
df_err_bet_g_fm <- df_err_bet_fixed_fm %>% group_by(sal) %>%
  summarise(mean_err_mean = mean(err_mean))
df_err_bet_g_fm <- df_err_bet_g_fm[-1,]

library("latex2exp")
x_inter <- 0.7650498
gg_betas_fm <- ggplot(df_err_bet_g_fm) + 
  geom_line(aes(sal, mean_err_mean), color ="#3293CE", size = 0.4) + 
  geom_point(aes(sal, mean_err_mean), color ="#2C2C54", size = 0.9 ) + 
  theme_bw() + xlab(TeX("$\\sigma_{\\beta}$")) + 
  ylab("Relative error")  +
  # geom_vline(xintercept = x_inter, color = "blue", linetype = "longdash") +
  theme(text = element_text(size = 15), legend.position = "bottom",
        plot.margin = margin(1, 1, 1, 1, "cm")) 
gg_betas_fm

## mu
mut = 1
sm = 0.1
mu_vec <- seq(0.2,1,0.01)
s_vec <- seq(0,0.9,0.01)
df_err_bet_fixed_mu <- data.frame(muw = muw, sw = sw, muc = muc, sc = sc, 
                                   mual = 0, sal = 0, max_eig = 0, out_mean = 0 ,
                                   count_re = 0)

count = 1
while(count < 50){
  for(j in c(1:length(mu_vec))){
    mut <- mu_vec[j]
    if(mut <= st){
      st <- 0.2
      mut <- 0.5
    }
    print(paste0("j:", j))
    betas = rgamma(N, shape = (mut/st)^2, rate = mut/(st^2))
    # Compute Jacobian
    jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      diag(gammas + colSums(MIGRATION))
    # Compute eigenvalues
    eig <- eigen_mat(jacobian)
    max_eig <- max(eig$re)
    betas = rep(mean(betas),N)
    # Compute Jacobian
    jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      diag(gammas + colSums(MIGRATION))
    # Compute eigenvalues
    eig <- eigen_mat(jacobian)
    max_eig_mean <- max(eig$re)
    df_err_bet_fixed_mu[nrow(df_err_bet_fixed_mu)+1,] <- c(muw, sw, muc, sc, mut,
                                                             st, max_eig, max_eig_mean , count) 
  }
  print(paste0("count:", count))
  count = count + 1 
}

df_err_bet_fixed_mu$err_mean <- (df_err_bet_fixed_mu$out_mean - df_err_bet_fixed_mu$max_eig)^2/df_err_bet_fixed_mu$out_mean
df_err_bet_fixed_mu$var_bet <- (df_err_bet_fixed_mu$sal)/(df_err_bet_fixed_mu$mual)
df_err_bet_g <- df_err_bet_fixed_mu %>% group_by(mual) %>%
  summarise(mean_err_mean = mean(err_mean))
df_err_bet_g <- df_err_bet_g[-1,]

library("latex2exp")
x_inter <- 0.7650498
gg_betas_fs <- ggplot(df_err_bet_g) + 
  geom_line(aes(mual, mean_err_mean), color ="#3293CE", size = 0.4) + 
  geom_point(aes(mual, mean_err_mean), color ="#2C2C54", size = 0.9 ) + 
  theme_bw() + xlab(TeX("$\\mu_{\\beta}$")) + 
  ylab("Relative error")  +
  # geom_vline(xintercept = x_inter, color = "blue", linetype = "longdash") +
  theme(text = element_text(size = 15), legend.position = "bottom",
        plot.margin = margin(1, 1, 1, 1, "cm")) 
gg_betas_fs


ggarrange(gg_betas_fm + ggtitle(TeX("$\\mu_{\\beta} = 1$")),
          gg_betas_fs + rremove("ylab") + ggtitle(TeX("$\\sigma_{\\beta} = 0.1$")))

## fixed CV vary sigma
CV = 0.5
s_vec <- seq(0.1,0.9,0.01)
# mu = sigma/CV
df_err_bet_fixed_CV <- data.frame(muw = muw, sw = sw, muc = muc, sc = sc, 
                                  mual = 0, sal = 0, max_eig = 0, out_mean = 0 ,
                                  CV = CV)

count = 1
while(count < 50){
  for(j in c(1:length(s_vec))){
    st <- s_vec[j]
    mut <-  s_vec[j]/CV
    if(mut <= st){
      st <- 0.2
      mut <- 0.5
    }
    print(paste0("j:", j))
    betas = rgamma(N, shape = (mut/st)^2, rate = mut/(st^2))
    # Compute Jacobian
    jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      diag(gammas + colSums(MIGRATION))
    # Compute eigenvalues
    eig <- eigen_mat(jacobian)
    max_eig <- max(eig$re)
    betas = rep(mean(betas),N)
    # Compute Jacobian
    jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      diag(gammas + colSums(MIGRATION))
    # Compute eigenvalues
    eig <- eigen_mat(jacobian)
    max_eig_mean <- max(eig$re)
    df_err_bet_fixed_CV[nrow(df_err_bet_fixed_CV)+1,] <- c(muw, sw, muc, sc, mut,
                                                           st, max_eig, max_eig_mean , CV) 
  }
  print(paste0("count:", count))
  count = count + 1 
}

df_err_bet_fixed_CV$err_mean <- (df_err_bet_fixed_CV$out_mean -
                                   df_err_bet_fixed_CV$max_eig)^2/df_err_bet_fixed_CV$out_mean
df_err_bet_fixed_CV$var_bet <- (df_err_bet_fixed_CV$sal)/(df_err_bet_fixed_CV$mual)
df_err_bet_g <- df_err_bet_fixed_CV %>% group_by(sal) %>%
  summarise(mean_err_mean = mean(err_mean))
df_err_bet_g <- df_err_bet_g[-1,]

library("latex2exp")
x_inter <- 0.7650498
gg_betas_fCV_s <- ggplot(df_err_bet_g) + 
  geom_line(aes(sal, mean_err_mean), color ="#3293CE", size = 0.4) + 
  geom_point(aes(sal, mean_err_mean), color ="#2C2C54", size = 0.9 ) + 
  theme_bw() + xlab(TeX("$\\sigma_{\\beta}$")) + 
  ylab("Relative error")  +
  # geom_vline(xintercept = x_inter, color = "blue", linetype = "longdash") +
  theme(text = element_text(size = 15), legend.position = "bottom",
        plot.margin = margin(1, 1, 1, 1, "cm")) 
gg_betas_fCV_s

## fixed CV vary mu
CV = 0.5
mu_vec <- seq(0.1,0.9,0.01)
# sigma = mu*CV
df_err_bet_fixed_CV <- data.frame(muw = muw, sw = sw, muc = muc, sc = sc, 
                                  mual = 0, sal = 0, max_eig = 0, out_mean = 0 ,
                                  CV = CV)

count = 1
while(count < 50){
  for(j in c(1:length(mu_vec))){
    mut <- mu_vec[j]
    st <-  mu_vec[j]*CV
    if(mut <= st){
      st <- 0.2
      mut <- 0.5
    }
    print(paste0("j:", j))
    betas = rgamma(N, shape = (mut/st)^2, rate = mut/(st^2))
    # Compute Jacobian
    jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      diag(gammas + colSums(MIGRATION))
    # Compute eigenvalues
    eig <- eigen_mat(jacobian)
    max_eig <- max(eig$re)
    betas = rep(mean(betas),N)
    # Compute Jacobian
    jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      diag(gammas + colSums(MIGRATION))
    # Compute eigenvalues
    eig <- eigen_mat(jacobian)
    max_eig_mean <- max(eig$re)
    df_err_bet_fixed_CV[nrow(df_err_bet_fixed_CV)+1,] <- c(muw, sw, muc, sc, mut,
                                                           st, max_eig, max_eig_mean , CV) 
  }
  print(paste0("count:", count))
  count = count + 1 
}

df_err_bet_fixed_CV$err_mean <- (df_err_bet_fixed_CV$out_mean -
                                   df_err_bet_fixed_CV$max_eig)^2/df_err_bet_fixed_CV$out_mean
df_err_bet_fixed_CV$var_bet <- (df_err_bet_fixed_CV$sal)/(df_err_bet_fixed_CV$mual)
df_err_bet_g <- df_err_bet_fixed_CV %>% group_by(mual) %>%
  summarise(mean_err_mean = mean(err_mean))
df_err_bet_g <- df_err_bet_g[-1,]

library("latex2exp")
x_inter <- 0.7650498
gg_betas_fCV <- ggplot(df_err_bet_g) + 
  geom_line(aes(mual, mean_err_mean), color ="#3293CE", size = 0.4) + 
  geom_point(aes(mual, mean_err_mean), color ="#2C2C54", size = 0.9 ) + 
  theme_bw() + xlab(TeX("$\\mu_{\\beta}$")) + 
  ylab("Relative error")  +
  # geom_vline(xintercept = x_inter, color = "blue", linetype = "longdash") +
  theme(text = element_text(size = 15), legend.position = "bottom",
        plot.margin = margin(1, 1, 1, 1, "cm")) 
gg_betas_fCV


gg1 <- ggarrange(gg_betas_fCV_s,
                 gg_betas_fCV + rremove("ylab"),
                 widths = c(1.1,1)) 

gg1 <- annotate_figure(gg1, top = text_grob("CV = 0.5", 
                                      color = "black", face = "bold", size = 14))

gg2 <- ggarrange(gg_betas_fm ,
          gg1, 
          nrow = 2,
          heights = c(1.3,2),
          labels = c("a","b"))

gg2 <- annotate_figure(gg2, top = text_grob(TeX("$\\mu_\\beta = 1$"), 
                                            color = "black", face = "bold", size = 14))

###### PLOT GRID #####
library("cowplot")
gg1 <- plot_grid(gg_betas_fCV_s + ggtitle("CV = 0.5"),
                 gg_betas_fCV + ggtitle("CV = 0.5") + rremove("ylab"),
                 rel_widths = c(1.1,1)) 

gg2 <- plot_grid(gg_betas_fm  + ggtitle(TeX("$\\mu_\\beta = 1$")),
                 gg1, 
                 nrow = 2,
                 rel_heights = c(1.6,2),
                 labels = c("a","b"))

gg2
