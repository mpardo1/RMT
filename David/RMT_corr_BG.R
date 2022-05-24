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
N <- 4000

muw <- 1
sw <- 0.7
rhow <- -0.4 #original rho (Gamma of baron et al)
Gammaw <- 0.7 #gamma of baron et al
rw <- 1.4
cw <- 1.6

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
  geom_point(aes(outl_BG,0), colour = "red") +
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
N = 1000
rw <- 0.34
cw <- 0.6
Gammaw <- 0.1

df_err <- data.frame(it = 0, out_re = 0)
count=0
max_it <- 2000
while(count < max_it){
  count = count +1
  com_mat_resc <- rand_mat_cor_norm_MPA_resc1(N,muw,sw,rhow,Gammaw,rw,cw)
  eig_resc <- eigen_mat(com_mat_resc[[1]])
  max_eig <- max(eig_resc$re)
  df_err[nrow(df_err) + 1, ] <-  c(count, max_eig)
}

outl_BG <- -1 + N*muw + ((N*muw/2)*(1+(rhow/(N*Gammaw)))*(sqrt(1+(4*Gammaw*sw^2/muw^2))-1))
outl <- -1 + N*muw + rhow*sw^2/muw

df_err$out_RMT <- outl
df_err$out_BG <- outl_BG

mean(df_err$out_BG - df_err$out_re)
mean(df_err$out_RMT - df_err$out_re)

N = 300
muw <- 0.1
sw <- 0.01
rhow <- 0.5
Gammaw <- 0.01
rw <- 0.75
cw <- 0.6

# outl_BG <- -1 + N*muw + ((N*muw/2)*(1+(rhow/(N*Gammaw)))*(sqrt(1+(4*Gammaw*sw^2/muw^2))-1))
outl_BG <- -1 + N*muw + ((N*muw/2)*(1+(rhow/Gammaw))*(sqrt(1+(4*Gammaw*sw^2/(N*muw^2)))-1))
outl <- -1 + N*muw 

com_mat_resc <- rand_mat_cor_norm_MPA_resc1(N,muw,sw,rhow,Gammaw,rw,cw)
plot_eig <- plot_eigen(com_mat_resc[[1]])
com_mat_resc[[2]]

plot_eig + 
  geom_point(aes(outl_BG,0), colour = "green") +
  geom_point(aes(outl,0), colour = "blue") +
  theme_bw()

plot_eig + 
  geom_point(aes(outl,0), colour = "blue") +
  theme_bw()

com_mat_resc[[2]]
eig_MPA_resc <- eigen_mat(com_mat_resc[[1]])
eig_MPA_circ <- eig_MPA[-which(eig_MPA_resc$re == max(eig_MPA_resc$re)),]

center <-  -1
radius <- sqrt(N)*sw
plot_eig + 
  geom_point(aes(outl_BG,0), colour = "green") +
  geom_point(aes(outl,0), colour = "blue") +
  geom_ellipse(aes(x0 = center, y0 = 0, a = (1+rhow)*radius,
                   b = (1-rhow)*radius, angle = 0), color = "red")
