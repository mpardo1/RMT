####### RANDOM MATRICES FOR METAPOPULATION MODELS #######
### 
### parent script
###
### generate, plot and integrate metapopulation
### epidemiological models
### 
#########################################################

setwd("~/RMT/David/")

source("~/RMT/David/RMT_genrandom_1.R")
source("~/RMT/David/RMT_plotmobility.R")
source("~/RMT/David/d_functions_eigen_int.R")

####### GENERATE JACOBIAN ###############################

# number of patches
N <- 300

muw <- 0.2
sw <- 0.05
rhow <- 0 #original rho (Gamma of baron et al)
Gammaw <- .15 #gamma of baron et al
rw <- .1
cw <- .3

(Gammaw/sqrt(rw*cw) < 1) & ((N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)
# Outliers:
outl_BG <- -1 + muw + (muw/2)*(1+rhow/Gammaw)*(sqrt(1+(4*Gammaw*sw/muw^2))-1)
outl <- muw

# Matrix MPA form:
com_mat <- rand_mat_cor_norm_MPA(N,muw,sw,rhow,Gammaw,rw,cw)
plot_eig <- plot_eigen(com_mat[[1]])
com_mat[[2]]

plot_eig + 
  geom_point(aes(outl_BG,0), colour = "green") +
  geom_point(aes(outl,0), colour = "blue") +
  theme_bw()

# Matrix DGG form:
com_mat_D <- rand_mat_cor_norm_DGG(N,muw,sw,rhow,Gammaw,rw,cw)
plot_eig <- plot_eigen(com_mat_D[[1]])
com_mat_D[[2]]

plot_eig + 
  geom_point(aes(outl_BG,0), colour = "green") +
  geom_point(aes(outl,0), colour = "blue") +
  theme_bw()
