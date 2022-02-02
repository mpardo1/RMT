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

# Filter the values where sigma = 0 and the beta dist properties holds.
df_filt <- df %>% filter(sigma != 0 & cond == TRUE)

beta_a_b_mat <- function(mat){
  beta_a_b(mat[1], mat[2])
}

mat_ab <- as.matrix(df_filt[,c(1,2)])
vec_ab <- apply(mat_ab,1, beta_a_b_mat)
# df_filt a df with the possible values for the mean and variance of the beta distribution
# with the respectives values for a and b (parameters of the beta distribution). 
df_filt$a <- vec_ab[1,]
df_filt$b <- vec_ab[2,]

N = 100
beta_ct = 0.3
gamma_ct = 0.4

mu_com <- 0.45
sigma_com <- 0.24
com <- df_filt[(df_filt$mean == mu_com & df_filt$sigma == sigma_com),]
alp_c <- com$a
bet_c <- com$b
mu_w <-  com$mean

check_mat <- function(mat){
  # print(paste0("a:", mat[1]))
  # print(paste0("b:", mat[2]))
  diff <- check_outl(N,beta_ct,gamma_ct,alp_c,bet_c, mu_w,mat[1],mat[2])
  # print(paste0("diff:", diff))
  return(diff)
}
mig <- as.matrix(df_filt[df_filt$mean == mu_w,])
vec <- apply(mig, 1, check_mat)
df <- as.data.frame(mig)
df$diff <- vec*100
high_diff <- df[df$diff > 0.5,]
low_diff <-  df[df$diff < 0.5,]
mean(high_diff$sigma)
mean(low_diff$sigma)
var(high_diff$sigma)
var(low_diff$sigma)
mean(high_diff$mean)
mean(low_diff$mean)
var(high_diff$mean)
var(low_diff$mean)
