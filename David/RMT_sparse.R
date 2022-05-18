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
mug <- gammas[1]
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
MIGRATION <- matrix(0, N, N)
# jacobian

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

### SPARSE MATRIX ####
p <- 0.5
vec_p <- which(rbinom(N^2, 1, p)==0)
COM_SPA <- COMMUTING
COM_SPA[vec_p] <- 0
MIGRATION[vec_p] <- 0

jacobian_sparse <- (COM_SPA + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

# Plot the eigenvalues of the system
sigma <- sqrt(mub^2*sw^2  + sc^2)
radius <- sigma*sqrt(N)

library("ggpubr")
outl <- mub*(muw*(N-1)+1) - mug
center <- mub*(1-muw) - mug - N*muc

plot_eig <- plot_eigen(jacobian) + 
  geom_point(aes(outl,0)) + 
  geom_circle(aes(x0 = center, y0 = 0, r = radius), color = "red")

outl_spa <- mub*(p*muw*(N-1)+1) - mug
center_spa <- mub*(1-p*muw) - mug - N*p*muc
sigma_spa <- sqrt((p*(1-p) - p^2)*(sigma - (mub*muw + muc)^2) - (p*(mub*muw+muc)))
radius_spa <- sqrt(N)*sigma_spa

plot_eig_spa <- plot_eigen(jacobian_sparse) + 
  geom_point(aes(outl_spa,0)) + 
  geom_circle(aes(x0 = center_spa, y0 = 0, r = radius_spa), color = "red")

ggarrange(plot_eig, plot_eig_spa)

#### MEAN = p*muw ###  
muw <- p*muw
muc <- p*muc
center_mean <- mub*(1-muw) - mug - N*muc
outl_mean <- mub*(muw*(N-1)+1) - mug

COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
diag(COMMUTING) <- 0
MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
diag(MIGRATION) <- 0

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

plot_eig_mean <- plot_eigen(jacobian) + 
  geom_point(aes(outl_mean,0)) + 
  geom_circle(aes(x0 = center_mean, y0 = 0, r = radius), color = "red")

ggarrange(plot_eig_spa,plot_eig_mean, nrow = 2)
