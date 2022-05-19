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
N <- 50

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
# MIGRATION <- matrix(0, N, N)
# jacobian

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

### SPARSE MATRIX ####
p <- 0.4
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
library("ggforce")
outl <- mub*(muw*(N-1)+1) - mug
center <- mub*(1-muw) - mug - N*muc

plot_eig <- plot_eigen(jacobian) + 
  geom_point(aes(outl,0), color = "#5448C8") + 
  geom_circle(aes(x0 = center, y0 = 0, r = radius), color = "#5448C8")
plot_eig

outl_spa <- mub*(p*muw*(N-1)+1) - mug
center_spa <- mub*(1-p*muw) - mug - N*p*muc
# Variance of the product of Bernouilli and the initial jacobian:
sigma_spa <- sqrt((p*(1-p)*sigma^2)+(p*(1-p)*(mub*muw + muc)^2) + sigma^2*p^2)
radius_spa <- sqrt(N)*sigma_spa

plot_eig_spa <- plot_eigen(jacobian_sparse) + 
  geom_point(aes(outl_spa,0), color = "#5448C8") + 
  geom_circle(aes(x0 = center_spa, y0 = 0, r = radius_spa), color = "#5448C8") + 
  theme_bw()
plot_eig_spa

ggarrange(plot_eig, plot_eig_spa)

### MOBILITY MATRIX ##

# plot the mobility network
#legend for plotmobility2 can be found in RMT_plotmobility
plotmobility(COM_SPA)
plot_mob <- plotmobility2(MIGRATION, COM_SPA)
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
  geom_point(aes(outl_mean,0), color = "#5448C8") + 
  geom_circle(aes(x0 = center_mean, y0 = 0, r = radius), color = "#5448C8") +
  theme_bw()

### MOBILITY MATRIX ##

# plot the mobility network
#legend for plotmobility2 can be found in RMT_plotmobility
plotmobility(COMMUTING)
plot_mob_mean <- plotmobility2(MIGRATION, COMMUTING)

library("latex2exp")
gg_arr <- ggarrange(plot_eig_spa + ggtitle("Sparse matrix with p probability"),
          plot_eig_mean + ggtitle(TeX("$Mean:p(\\beta\\mu_w + \\mu_c)$")) +
            ylab(""),
          plot_mob,
          plot_mob_mean,
          nrow = 2, ncol = 2)
gg_arr

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/"
path <- paste0(Path,"Plot_spa_mean1_b0,1_g0,95_muc_0,004_sc0,002_muw0,24_sw0,05.png")
ggsave(path,
       plot = gg_arr, device = "png")
# Compute the difference between the right most eigenvalue with sparse
# and with the matix with mean p*muc, p*muw
p_vec <- seq(0.1,1,0.01)
dim <- length(p_vec)
df_spa <- data_frame(p = 0, max_eig_m = 0, max_eig_spa = 0)
for(i in c(1:dim)){
  COMMUTING <- rand_mat(N, p_vec[i]*muw, sw, distrib = "beta")
  diag(COMMUTING) <- 0
  MIGRATION <- rand_mat(N, p_vec[i]*muc, sc, distrib = "beta")
  diag(MIGRATION) <- 0
  
  jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
    diag(deaths + alphas + deltas + colSums(MIGRATION))
  
  eig_mean <- eigen_mat(jacobian)
  max_eig_mean <- max(eig_mean$re)
  
  #### Sparse matrix ###
  COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
  diag(COMMUTING) <- 0
  MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
  diag(MIGRATION) <- 0
  vec_p <- which(rbinom(N^2, 1, p_vec[i])==0)
  COMMUTING[vec_p] <- 0
  MIGRATION[vec_p] <- 0
  
  ### MOBILITY MATRIX #
  
  # plot the mobility network
  #legend for plotmobility2 can be found in RMT_plotmobility
  plotmobility(COMMUTING)
  plotmobility2(MIGRATION, COMMUTING)
  jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
    diag(deaths + alphas + deltas + colSums(MIGRATION))
  
  eig_spa <- eigen_mat(jacobian)
  max_eig_spa <- max(eig_spa$re)
  
  df_spa[nrow(df_spa)+1,] <- list(p_vec[i], max_eig_mean, max_eig_spa)
}

df_spa <- df_spa[-1,]
df_plot <- reshape2::melt(df_spa, id.vars="p")
gg_max_eig <- ggplot(df_plot) +
  geom_line(aes(p, value, colour=variable)) +
  xlab("Right most eigenvalue real part") +
  theme_bw()

Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/"
path <- paste0(Path,"Plot_spa_max_eig1_b0,1_g0,95_muc_0,004_sc0,002_muw0,24_sw0,05.png")
ggsave(path,
       plot = gg_max_eig, device = "png")
