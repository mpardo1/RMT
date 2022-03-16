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

Deltas <- rep(0.6, N) # birth rate
mub <- 0.6
sb <- 0.001
betas <- rep(mub, N) # transmission rates
# betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.3, N) # loss of immunity rates
mud <- 0.6
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.5
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

muc <- 0.1
sc <- 0.05
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

# plot_eigen_rmt(jacobian,
               # N,mub,mug = mud + mua + mudel,
               # muw,sw,rhow,Gammaw,
               # muc,sc,rhoc,Gammac,
               # tau = 0, alp = 0, K = 0)

####### INTEGRATE SYSTEM ################################

# initial populations
# for constant populations, set deltas = 0, Deltas = deaths
mu_sus <- 10000
s_sus <- 100
# alphag <- alphagamma(mu_sus,s_sus)
# betag <- betagamma(mu_sus,s_sus)
# sus_init <- rgamma(N,alphag,betag) 
sus_init <- rep(10000, N) # initial susceptibles

mu_inf <- 200
s_inf <- 150
# alphag <- alphagamma(mu_inf,s_inf)
# betag <- betagamma(mu_inf,s_inf)
# inf_init <- rgamma(N,alphag,betag) 
df_inf <- data.frame(inf_init=0,max_inf=0,max_inf1=0,max_inf2=0)
list_sol <- list()
vec_init_inf <- seq(1,3000,20)
for(i in c(1:length(vec_init_inf))){
  inf_init <- rep(vec_init_inf[i], N)    # initial infecteds
  end_time <- 100
  Deltas <- rep(0.3, N) # birth rate
  sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
             COMMUTING,MIGRATION,
             sus_init,inf_init,end_time)
  max_inf <- max(sol[nrow(sol),])
  
  Deltas <- rep(0.8, N) # birth rate
  sol1 <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
             COMMUTING,MIGRATION,
             sus_init,inf_init,end_time)
  max_inf1 <- max(sol1[nrow(sol1),])
  
  Deltas <- rep(0.1, N) # birth rate
  sol2 <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
              COMMUTING,MIGRATION,
              sus_init,inf_init,end_time)
  max_inf2 <- max(sol1[nrow(sol2),])
  df_inf[nrow(df_inf) +1,] <- c(vec_init_inf[i], max_inf, max_inf1, max_inf2)
 
  
}

path <- paste0("~/RMT/Integration/ci_inf",Sys.Date(),".csv")
write.csv(df_inf, path, row.names = TRUE)

# ggplot(df_inf) + 
#   geom_line(aes(inf_init, max_inf)) +
#   geom_line(aes(inf_init, max_inf1)) +
#   xlab("Number of initial infected individuals")+
#   ylab("Number of infected individuals at equilibrium")+
#   theme_bw()
