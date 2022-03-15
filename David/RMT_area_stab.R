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
N <- 10

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

muc <- 0.01
sc <- 0.00001
rhoc <- 0
Gammac <- 0
rc <- 0
cc <- 0

N_vec <- seq(5,500,1)
muw_vec <- seq(0.01,1,0.001)
muw_vec <- muw_vec[1:length(N_vec)]
df_sol <- data.frame(beta = 0, gamma = 0, N = 0, muw = 0, state = FALSE)
for(i in c(1:length(N_vec))){
  for(j in c(1:length(N_vec))){
    COMMUTING <- rand_mat(N_vec[i], muw_vec[j], sw, distrib = "beta")
    diag(COMMUTING) <- 0
    MIGRATION <- rand_mat(N_vec[i], muc, sc, distrib = "beta")
    diag(MIGRATION) <- 0
    
    # EPI param:
    betas <- rep(mub, N_vec[i])
    deltas <- rep(mudel, N_vec[i])
    deaths <- rep(mud, N_vec[i]) 
    alphas <- rep(mua, N_vec[i])
    
    jacobian <- (COMMUTING + diag(N_vec[i])) %*% diag(betas) + MIGRATION -
      diag(deaths + alphas + deltas + colSums(MIGRATION))
    eigen <- eigen_mat(jacobian)
    max_eig <- max(eigen$re)
    state <- ifelse(max_eig >0, FALSE, TRUE)
    df_sol[nrow(df_sol) + 1,1:4] <- list(betas[1], gammas[1], N_vec[i], muw_vec[j])
    df_sol[nrow(df_sol) ,5] <- state
  }
}

df_sol <- df_sol[-1,]

path <- paste0("~/RMT/Integration/area_gen_",Sys.Date(),".csv")
write.csv(df_sol, row.names = TRUE)

# library(latex2exp)
# ggplot(df_sol) +
#   geom_point(aes(N,muw, colour = state)) + theme_bw()  +
#   scale_color_manual(values=c("#6622CC", "#A755C2")) + 
#   ylab(""*mu[c]~"") + 
#   # ggtitle(""*gamma/beta~": 4")
#   ggtitle(TeX("\\frac{\\gamma}{\\beta}: 0.4"))
