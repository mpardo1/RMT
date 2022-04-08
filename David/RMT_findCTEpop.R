####### RANDOM MATRICES FOR METAPOPULATION MODELS #######
### 
### parent script
###
### generate, plot and integrate metapopulation
### epidemiological models
### 
#########################################################
rm(list=ls()) 
source("~/RMT/David/RMT_genrandom.R")
source("~/RMT/David/RMT_plotmobility.R")
source("~/RMT/David/d_functions_eigen_int.R")
library(matlib)
####### GENERATE JACOBIAN ###############################
N <- 5
# epidemiological
#all rates must lie in (0,1) except for betas

Deltas <- rep(0.3, N) # birth rate
mub <- 0.05
sb <- 0.001
betas <- rep(mub, N) # transmission rates
# betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.1, N) # loss of immunity rates
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
sw <- 0.02
rhow <- 0 #original rho (Gamma of baron et al)
Gammaw <- 0 #gamma of baron et al
rw <- 0
cw <- 0

muc <- 0.1
sc <- 0.001
rhoc <- 0
Gammac <- 0
rc <- 0
cc <- 0

COMMUTING <- rand_mat(N, muw, sw, distrib = "gamma")
COMMUTING <- matrix(0, N,N)
diag(COMMUTING) <- 0

MIGRATION <- rand_mat(N, muc, sc, distrib = "gamma")
diag(MIGRATION) <- 0

mat <- MIGRATION
diag(mat) <- - colSums(MIGRATION)

init_pop <- c()
N0 <- 100
# Cramer's Rule:
for(i in c(1:(N))){
  mat_1 <- matrix(1,1,(N-1))
  if(i > 1 & i < N){
    print("i entre 2 y N-1")
    sub_mat <- mat[c((1:(i-1)),((i+1):N)),c((1:(i-1)),((i+1):N))]
    mat_i <- -mat[c(1:(i-1),(i+1):N),i]
    mat_i <- mat_i%*%mat_1
  }else if(i == 1){
    print("i = 1")
    sub_mat <- mat[(2:N),(2:N)]
    mat_i <- -mat[(2:N),i]
    mat_i <- mat_i%*%mat_1
   }
  det_1 <- det(sub_mat)
  det_2 <- det(sub_mat + mat_i)
  init_pop[i] <- det_1*N0/(det_2)
}

init_pop_1 <- init_pop

# Integrate the system:
sus_init <- snit_pop
inf_init <- rep(0,N)
end_time <- 50
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
                   COMMUTING,MIGRATION,
                   sus_init,inf_init,end_time)
# library("ggpubr")
plot_no_com_tot <- plot_int(N, sol, state = "TOT")

