rm(list = ls())
library("tidyverse")
library("deSolve")
library("ggplot2")
library("gifski")

#----------------------------------------------------------------------------#
source("~/RMT/Integration/functions_eigen_int.R")
#-------------------EPIDEMIOLOGICAL----------------------
N = 100 # Number of patches
# CTE parameters:
del_N <- rep(0.6, N) # Birth rate
# bet_cte <- 6
bet_cte <-  0.03
bet <- rep(bet_cte, N)  # Transmission rate
# bet <- abs(rnorm(N,1,1))  # Transmission rate
d_vec <- rep(0.2, N) # Natural mortality rate
thet <- rep(0.1, N) # Rate of loss of immunity
alp <- rep(0.02, N) # Rate of disease overcome
delt <- rep(0, N) # Diseases related mortality rate
beta_ct = bet_cte  # beta
gamma_ct = alp[1] + delt[1] + d_vec[1]     # gamma  
# -------------------CHECK VISUALLY ----------------------#

# # Random matrix general, not bounded by [0,1]:
commut_mat <- rand_mat(N,60 ,12,"gamma")
migrate_mat <- rand_mat(N,60 ,0.01,"gamma")
jac <- jacobian(N,bet,gamma_ct, commut_mat, migrate_mat,mu_m, 3)
eig <- eigen_mat(jac)

ind <-  which(eig$re == max(eig$re))
eig <-  eig[-ind,]
plot_eig <- ggplot(eig) + geom_point(aes(re,im), size = 0.05) 
plot_eig + coord_fixed() 
# plot_eig + xlim(c(-5,5))+ ylim(c(-5,5))

#-----------------CHECK OUTLIER---------------#
# Compute mean and sd for migration coefficients:
s_w <- seq(0,(1/4),0.001)
l <- length(s_w)
df = data.frame(mean = 0, sigma = 0, cond = FALSE)
for(i in c(1:l)){
  mu_w <- seq(0,1,0.01)
  list <- sapply(mu_w, validate_mu_s, sigma = s_w[i])
  df1 <- data.frame(mean = mu_w, sigma = s_w[i], cond = list)
  df <- rbind(df,df1)
}

# Color map TRUE posible values for mu and sigma in a beta distribution.
ggplot(df) + 
  geom_point(aes(mean,sigma, colour = cond))

# Compute sigma of the migration as lamba*s_w:
df$s_c <- scale_sigma(N,beta_ct,df$sigma)
