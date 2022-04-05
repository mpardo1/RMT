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

####### GENERATE JACOBIAN ###############################

# number of patches
 N <- 100

# epidemiological
#all rates must lie in (0,1) except for betas

Deltas <- rep(0.3, N) # birth rate
mub <- 0.1
sb <- 0.001
betas <- rep(mub, N) # transmission rates
# betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.1, N) # loss of immunity rates
mud <- 0.3
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.65
alphas <- rep(mua, N) # recovery rates
mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates
gammas = deaths + alphas + deltas

# mobility
#commuting and migration networks

muw <- 0.6
sw <- 0.05
rhow <- 0 #original rho (Gamma of baron et al)
Gammaw <- 0 #gamma of baron et al
rw <- 0
cw <- 0

muc <- 0.01
sc <- 0.002
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

# jacobian

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

# plot the eigenvalues of the system
library("ggforce")
eigen_stab <- plot_eigen_rmt(jacobian,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc,sc,rhoc,Gammac,
               tau = 0, alp = 0, K = 0) +
               scale_y_continuous( breaks=c(0)) +
              xlim(c(-102,-99)) + ylim(c(-1,1))

eig <- eigen_mat(jacobian)
eig <- eig[-nrow(eig),] 
ggplot(eig, aes(re, im)) + 
  geom_point(size = 0.5) +
  coord_fixed() + 
  xlim(c(-2,-1.8))

print(plot_eigen(jacobian))
####### INTEGRATE SYSTEM ################################

# initial populations
# for constant populations, set deltas = 0, Deltas = deaths

sus_init <- rep(100000, N) # initial susceptibles
inf_init <- rep(100, N)    # initial infecteds

end_time <- 20

# integro el sistema con condiciones iniciales 
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

# plot SUS, INF, REC or TOT population
plot_inf_stab <- plot_int(N, sol, state = "INF")  + theme_bw()+ theme(legend.position = "none")
plot_inf_stab
plot_int(N, sol, state = "TOT")


# sol_df <-  as.data.frame(sol)
# for(i in c(1:N)){
#   colnames(sol_df)[i+1] <-  paste0("S",i)
#   colnames(sol_df)[N+i+1] <-  paste0("I",i)
#   colnames(sol_df)[2*N+i+1] <-  paste0("R",i)
# }

####### MOBILITY MATRIX ###########################################

# plot the mobility network
#legend for plotmobility2 can be found in RMT_plotmobility
plotmobility(COMMUTING)
plotmobility2(MIGRATION, COMMUTING)



eigen_df <- eigen_mat(jacobian)
####### PERTURBATIONS ###################################
alp <- 1.4
K = 1
ind <- sample(1:N,K)
betas[ind] <- alp + betas[ind]

# jacobian

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

# plot the eigenvalues of the system

plot_eigen_rmt(jacobian,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc,sc,rhoc,Gammac,
               tau = 0, alp, K)

