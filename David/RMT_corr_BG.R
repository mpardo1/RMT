####### RANDOM MATRICES FOR METAPOPULATION MODELS #######
### 
### parent script
###
### generate, plot and integrate metapopulation
### epidemiological models
### 
#########################################################
rm(list=ls()) 

source("~/RMT/David/RMT_genrandom_1.R")
source("~/RMT/David/RMT_plotmobility.R")
source("~/RMT/David/d_functions_eigen_int.R")

####### GENERATE JACOBIAN ###############################

# number of patches
N <- 400

muw <- 1
sw <- 0.7
rhow <- 0.4 #original rho (Gamma of baron et al)
Gammaw <- 0.1 #gamma of baron et al
rw <- 0.4
cw <- 0.3

(Gammaw/sqrt(rw*cw) < 1) & ((N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)
# Outliers:
outl_BG <- -1 + muw + (muw/2)*(1+rhow/Gammaw)*(sqrt(1+(4*Gammaw*sw^2/muw^2))-1)
outl <- -1 + muw + rhow*sw^2/muw

# Matrix MPA form:
com_mat <- rand_mat_cor_norm_MPA(N,muw,sw,rhow,Gammaw,rw,cw)
plot_eig <- plot_eigen(com_mat[[1]])
com_mat[[2]]

plot_eig + 
  geom_point(aes(outl_BG,0), colour = "green") 

eig_MPA <- eigen_mat(com_mat[[1]])
max_eig <- max(eig_MPA$re)
eig_MPA_circ <- eig_MPA[-which(eig_MPA$re == max(eig_MPA$re)),]

library(latex2exp)
plot_eig + 
  geom_point(aes(outl_BG,0), colour = "green") +
  geom_point(aes(outl,0), colour = "blue") +
  ggtitle(TeX("$outBG:$")) +
  theme_bw()

eig_MPA <- eigen_mat(com_mat[[1]])
eig_MPA_circ <- eig_MPA[-which(eig_MPA$re == max(eig_MPA$re)),]

center <-  -1
radius <- sw
ggplot(eig_MPA) + 
  geom_point(aes(re, im), size = 0.1) + 
  geom_circle(aes(outl_BG,0), colour = "red") +
  geom_point(aes(outl,0), colour = "green") +
  geom_ellipse(aes(x0 = center, y0 = 0, a = (1+rhow)*radius,
                   b = (1-rhow)*radius, angle = 0), color = "red") 

# Matrix DGG form:
com_mat_D <- rand_mat_cor_norm_DGG(N,muw,sw,rhow,Gammaw,rw,cw)
plot_eig <- plot_eigen(com_mat_D[[1]])
com_mat_D[[2]]

plot_eig + 
  geom_point(aes(outl_BG,0), colour = "green")

plot_eig + 
  geom_point(aes(outl_BG,0), colour = "green") +
  geom_point(aes(outl,0), colour = "blue", size= 0.2) +
  theme_bw() 

eig_DD <- eigen_mat(com_mat_D[[1]])
eig_DD_circ <- eig_DD[-which(eig_DD$re == max(eig_DD$re)),]

center <-  -1
radius <- sw/sqrt(N)
ggplot(eig_DD_circ) + 
  geom_point(aes(re, im), size = 0.1) + 
  geom_ellipse(aes(x0 = center, y0 = 0, a = (1+rhow)*radius,
                   b = (1-rhow)*radius, angle = 0), color = "red") 
  
# Reescale variables:
com_mat_resc <- rand_mat_cor_norm_MPA_resc(N,muw,sw,rhow,Gammaw,rw,cw)
plot_eig <- plot_eigen(com_mat_resc[[1]])
com_mat_resc[[2]]

outl_BG <- -1 + N*muw + ((N*muw/2)*(1+(rhow/(N*Gammaw)))*(sqrt(1+(4*Gammaw*sw^2/muw^2))-1))
outl <- -1 + N*muw + rhow*sw^2/muw

plot_eig + 
  geom_point(aes(outl_BG,0), colour = "green") +
  geom_point(aes(outl,0), colour = "blue") +
  theme_bw()

plot_eig + 
  geom_point(aes(outl,0), colour = "blue") +
  theme_bw()

eig_MPA_resc <- eigen_mat(com_mat_resc[[1]])
eig_MPA_circ <- eig_MPA[-which(eig_MPA_resc$re == max(eig_MPA_resc$re)),]

center <-  -1
radius <- sqrt(N)*sw
ggplot(eig_MPA_resc) + 
  geom_point(aes(re, im), size = 0.1) + 
  geom_point(aes(outl_BG,0), colour = "green") +
  geom_point(aes(outl,0), colour = "blue") +
  geom_ellipse(aes(x0 = center, y0 = 0, a = (1+ rhow)*radius,
                   b = (1-rhow)*radius, angle = 0), color = "red") 
