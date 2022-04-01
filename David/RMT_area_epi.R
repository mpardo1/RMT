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
muw <- 0.6 
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

end_time <- 100


step <- 0.001
beta_vec <- seq(0,0.1,step)
alp_vec <- seq(0,2.5,step)
df_sol <- data.frame(beta = 0, gamma = 0, N = 0, alp = 0, state = FALSE)
N = 100
mug <- gammas[1]
for(i in c(1:length(beta_vec))){
  print(paste0(" i : ", i))
  for(j in c(1:length(alp_vec))){
    a <- beta_vec[i]*muw + muc
    b <- alp_vec[j]
    c <- alp_vec[j]*muw
    outl <- (1/2)*(N*a + b + sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
    outl <- outl + (mub*(1-muw) - N*muc - mug)
    state <- ifelse(outl < 0, TRUE, FALSE)
    df_sol[nrow(df_sol) + 1,1:4] <- list(beta_vec[i], gammas[1], N, alp_vec[j])
    df_sol[nrow(df_sol) ,5] <- state
  }
}

df_sol <- df_sol[-1,]
Path <- "~/RMT/David/OUTPUT/"
path <- paste0(Path,"Areaepi_g0,5_muc_0,001_sc0,0001_muw0,6_sw0,05_",Sys.Date(), ".csv")
write.csv(df_sol, path,row.names = TRUE)