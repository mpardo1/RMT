####### RANDOM MATRICES FOR METAPOPULATION MODELS #######
### 
### parent script
###
### generate, plot and integrate metapopulation
### epidemiological models
### 
#########################################################

setwd("~/D/bmb/r0/miR_R0")

source("./RMT_genrandom.R")
source("./RMT_plotmobility.R")
source("./d_functions_eigen_int.R")

####### GENERATE JACOBIAN ###############################

# number of patches
N <- 100

# epidemiological
#all rates must lie in (0,1) except for betas

Deltas <- rep(0.6, N) # birth rate
mub <- 0.9
sb <- 0.001
betas <- rep(mub, N) # transmission rates
betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.3, N) # loss of immunity rates
mud <- 0.3
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.5
alphas <- rep(mua, N) # recovery rates
mudel <- 0.15
deltas <- rep(mudel, N) # disease-related death rates
#gammas = deaths + alphas + deltas

# mobility
#commuting and migration networks

muw <- 0.2 
sw <- 0.05
rhow <- 0 #original rho (Gamma of baron et al)
Gammaw <- 0 #gamma of baron et al
rw <- 0
cw <- 0

muc <- 0.001
sc <- 0.0005
rhoc <- 0
Gammac <- 0
rc <- 0
cc <- 0

COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
COMMUTING[sample.int(N^2, round(p*N^2))] <- 0

MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")

# jacobian

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(deaths + alphas + deltas)

####### INTEGRATE SYSTEM ################################

# initial populations
# for constant populations, set deltas = 0, Deltas = deaths

sus_init <- rep(100000, N) # initial susceptibles
inf_init <- rep(100, N)    # initial infecteds

end_time <- 10

# integro el sistema con condiciones iniciales 
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

####### PLOTS ###########################################

# plot the mobility network
#legend for plotmobility2 can be found in RMT_plotmobility
plotmobility(COMMUTING)
plotmobility2(MIGRATION, COMMUTING)

# plot the eigenvalues of the system

plot_eigen_rmt(jacobian,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc,sc,rhoc,Gammac,
               tau = 0)

eigen_df <- eigen_mat(jacobian)

# plot SUS, INF, REC or TOT population
plot_int(N, sol, state = "INF")

sol_df <-  as.data.frame(sol)
for(i in c(1:N)){
  colnames(sol_df)[i+1] <-  paste0("S",i)
  colnames(sol_df)[N+i+1] <-  paste0("I",i)
  colnames(sol_df)[2*N+i+1] <-  paste0("R",i)
}

####### PERTURBATIONS ###################################

# improve

errores <- vector()
for (j in c(1:1000)) {
  
  N <- 300
  beta <- 2
  alpha <- 0
  gamma <- 4
  muw <- 3
  muws <- 20
  sw <- 3
  
  betas <- matrix(rep(0,N^2), nrow = N)
  diag(betas) <- rep(beta,N)
  betas[1,1] <- beta+alpha
  
  BIGT <- rand_mat(N, muw, sw, distrib = "gamma")
  diag(BIGT) <- rep(1,N)
  #BIGT[1,2] <- muws # single edge, single direction
  #BIGT[1,2] <- BIGT[2,1] <- muws # single edge, both directions
  BIGT[1,c(2:N)] <- BIGT[2,1] <- BIGT[2,c(3:N)] <- BIGT[2,1]+muws-muw # two rows
  #BIGT[1,c(2:N)] <- BIGT[2,1] <- BIGT[2,c(3:N)] <- BIGT[3,c(1:2)] <- BIGT[3,c(4:N)] <- BIGT[2,1]+muws-muw # three edges, all directions
  BIGT <- BIGT%*%betas
  
  BIGS <- matrix(rep(0,N^2), nrow = N)
  diag(BIGS) <- rep(-gamma, N)
  
  jacobian <- BIGT+BIGS
  ngm <- -BIGT%*%solve(BIGS)
  
  outlier <- beta*(1-muw +(muw*(N-1))/2 +muws/2 + sqrt(muw^2*((N-3)^2)/4 - (muw*muws*(3*N-5))/2 + (muws^2)/4) ) - gamma
  #outlier <- beta*(1-muw +muw*N/2 + (1/2)*(alpha/beta) + sqrt((muw^2)*(N^2)/4+(N-1)*muw^2*alpha/beta-(N/2-1)*muw*alpha/beta+((1/2)*(alpha/beta))^2))-gamma # single edge, single direction
  
  plot_eigen(jacobian) + coord_fixed() + theme_bw() +
    geom_vline(xintercept = 0, color = "blue") +
    #geom_circle(aes(x0 = beta-gamma, y0 = 0, r = beta*(N-1)*muw), color = "green") +
    #geom_circle(aes(x0 = beta-gamma, y0 = 0, r = beta*((N-2)*muw+muws)), color = "green") +
    geom_point(aes(x = trueout, y = 0), color = "red") +
    geom_ellipse(aes(x0 = beta*(1-muw)-gamma, y0 = 0, a = sqrt(N)*beta*sw, b = sqrt(N)*beta*sw, angle = 0), color = "red")
  
  #outlier
  #eigen_mat(jacobian)$re[1]
  
  errores[j] <- outlier-eigen_mat(jacobian)$re[1]
  
}
#########################################################