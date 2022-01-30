rm(list = ls())
library("tidyverse")
library("deSolve")
library("ggplot2")
library("gifski")

#----------------------------------------------------------------------------#
source("~/RMT/Integration/functions_eigen_int.R")

# -------------------CHECK ISSUE SUM CIJ ----------------------#

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
beta_ct = bet_cte  # beta
gamma_ct = alp[1] + delt[1] + d_vec[1]     # gamma   

a <-  2
b <-  3
diff <-  check_outl(N,a,b,beta_ct,gamma_ct)

dim <- 1000
a <- abs(rnorm(dim, 5.0, 10))
b <- abs(rnorm(dim, 5.0, 10))

# Compute mean and sd for migration coefficients:
mu_w <- seq(0,1,0.01)

appfunc <- function(mu){
  s_w <- seq(0,0.25,0.01)
  sapply(mu, beta_a_b, sigma = s_w)
}

sapply(mu_w, appfunc)

alp_m <- 1.03
bet_m <- 0.01
mu_m <- alp_m/(alp_m + bet_m)
s_m <-  sqrt((alp_m*bet_m)/(((alp_m + bet_m)^2)*(1+alp_m+bet_m)))
print(paste0("mu :", mu_m))
print(paste0("sigma :", s_m))

vec <-  sapply(a,check_outl,N=N, b = b[1],beta_ct =beta_ct, gamma_ct= gamma_ct, a_m=a_m, b_m=bm)
for(i in c(2:dim)){
  vec1 <-  sapply(a,check_outl,N=N, b = b[i],beta_ct =beta_ct, gamma_ct= gamma_ct)
  vec <-  c(vec,vec1)
}
sort(vec, decreasing = TRUE)
mean(vec)
