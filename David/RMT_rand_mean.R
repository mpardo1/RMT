####### RANDOM MATRICES FOR METAPOPULATION MODELS #######
### 
### parent script
###
### generate, plot and integrate metapopulation
### epidemiological models
### 
# rm(list = ls())
source("~/RMT/David/RMT_genrandom.R")
source("~/RMT/David/RMT_plotmobility.R")
source("~/RMT/David/d_functions_eigen_int.R")

####### GENERATE JACOBIAN ###############################

# number of patches
N <- 100

# epidemiological
#all rates must lie in (0,1) except for betas

Deltas <- rep(0.3, N) # birth rate
mub <- 0.02
sb <- 0.001
betas <- rep(mub, N) # transmission rates
# betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.3, N) # loss of immunity rates
mud <- 0.3
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.2
alphas <- rep(mua, N) # recovery rates
mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates
gammas = deaths + alphas + deltas

# mobility
#commuting and migration networks
muw <- 0.1 
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

mub <- 0.5
df_epi <- data.frame(mean = mub, mean_arit = 0, sigma = 0, max_rand = 0, max_mean = 0)
sig_vec <- seq(0,1,0.01)
len <- length(sig_vec)
count = 0
while( count < 1000){
  for(i in c(1:len)){
    ##### Random beta ####
    sb <- sig_vec[i]
    betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2)) 
    
    jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      diag(deaths + alphas + deltas + colSums(MIGRATION))
    eig <- eigen_mat(jacobian)
    max_rand <- max(eig$re)
    
    # Mean(rand(betas)):
    betas <- rep(mean(betas),N)
    jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      diag(deaths + alphas + deltas + colSums(MIGRATION))
    eig <- eigen_mat(jacobian)
    max_mean <- max(eig$re)
    
    df_epi[nrow(df_epi)+1,] <- c(mub, mean(betas), sb, max_rand, max_mean)
  }
  count = count + 1
}

df_epi$sq_err <- ((df_epi$max_rand - df_epi$max_mean)/df_epi$max_rand)^2



