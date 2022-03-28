####### RANDOM MATRICES FOR METAPOPULATION MODELS #######
### 
### parent script
###
### generate, plot and integrate metapopulation
### epidemiological models
### 
#########################################################
source("~/RMT/David/RMT_genrandom.R")
source("~/RMT/David/RMT_plotmobility.R")
source("~/RMT/David/d_functions_eigen_int.R")
library(matlib)
####### GENERATE JACOBIAN ###############################

# number of patches
N <- 3

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
diag(COMMUTING) <- 0

MIGRATION <- rand_mat(N, muc, sc, distrib = "gamma")
diag(MIGRATION) <- 0


# list <- births_func2(MIGRATION)
# deaths <- list[[1]]
# Deltas <- list[[2]]

# ### TRY with N = inv(D)*Deltas
init_pop <- DFE_func(MIGRATION, Deltas, deaths)
sus_init <- init_pop # initial susceptibles
inf_init <- rep(0, N)    # initial infecteds

sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
                   COMMUTING,MIGRATION,
                   sus_init,inf_init,end_time)
plot_int(N, sol, state = "TOT")

### TRY with Gamma = D*N
sus_init <- rep(100, N) # initial susceptibles
inf_init <- rep(0, N)    # initial infecteds
init_pop <- sus_init+inf_init
# Deltas <- births_func(MIGRATION, init_pop, deaths)

if(all(Deltas > 0)){
  print("All good bro")
}else{
  print("Problem with the births rates, NEGATIVE!!")
}

end_time <- 10
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

sol_1 <- rowSums(sol[, c(2,4,6)])
sol_2 <- rowSums(sol[, c(3,5,7)])
plot_int(N, sol, state = "TOT")

# Check if the pop is cte along time at each patch:
pop_patch <- matrix(0, nrow = nrow(sol), ncol = N+1)
pop_patch[,1] <- sol[,1]
for(i in c(1: nrow(sol))){
  for(j in c(1:N)){
    pop_patch[i,j+1] <- sol[i,(j+1)] + sol[i,(1+j+N)] + sol[i,(1+j+2*N)]
  }
}

########################################################################
####### EIGENVALUES #######
sus_init <- init_pop_func(MIGRATION, Deltas, deaths)
# jacobian
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

# plot the eigenvalues of the system
library("ggforce")
plot_eigen_rmt(jacobian,N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,muc,sc,rhoc,Gammac,
               tau = 0, alp = 0, K = 0) +
               scale_y_continuous( breaks=c(0)) 

print(plot_eigen(jacobian))
plot_int(N, sol, state = "INF")

