rm(list = ls())
library("tidyverse")
library("deSolve")
library("ggplot2")
library("gifski")
library("ggpubr")
library("ggforce")
library("stats")

#----------------------------------------------------------------------------#
source("~/RMT/Integration/functions_eigen_int.R")
#----------------PARAMETERS-----------------
#-------------------LOGIC----------------------
# Mobility parameter, 0: Just commuting, 1: just migration 2: migration & commuting.
MOB <- 2
# Integration parameter, 0: No integration, 1: integration.
INT <- 1
# Parameter for initial population. 0: No cte, 1: cte.
CTE_POP <- 1
# Parameter for constant population at each patch. 0: No cte, 1: cte.
CTE_POP_Patch <- 0
# Parameter for transmission rate. 0: No cte, 1: cte.
BETA_CTE <- 0
# Parameter for initial infected ind. 0: No cte , 1: cte.
CTE_INF <- 1
# Parameter for distribution, "normal", "beta", "gamma":
DIST <-  "beta"

#-------------------EPIDEMIOLOGICAL----------------------
N = 100 # Number of patches
# CTE parameters:
del_N <- rep(0.6, N) # Birth rate
# bet_cte <- 6
bet_cte <- 9.12014
bet <- rep(bet_cte, N)  # Transmission rate
# bet <- abs(rnorm(N,1,1))  # Transmission rate
d_vec <- rep(1, N) # Natural mortality rate
thet <- rep(0.6, N) # Rate of loss of immunity
alp <- rep(1, N) # Rate of disease overcome
delt <- rep(0, N) # Diseases related mortality rate
gamma_ct <-  alp[1] + delt[1] + d_vec[1]
print(paste0("gamma:", alp[1] + delt[1] + d_vec[1]))
print(paste0("beta - gamma:", bet[1] - (alp[1] + delt[1] + d_vec[1])))

#-------------------- MOBILITY ---------------------------
### Migration:
mu_c <- 0.01
s_c <- 0.00001
alp_c <- beta_a_b(mu_c, s_c)[1]
bet_c <- beta_a_b(mu_c, s_c)[2]
# Compute mean and sd:
print(paste0("mu :", mu_c))
print(paste0("sigma :", s_c))

migrate_mat <- mat_conect(N,alp_c,bet_c,DIST)
### Commuting
mu_w <- 0.4
s_w <- 0.2
alp_w <- beta_a_b(mu_w, s_w)[1]
bet_w <- beta_a_b(mu_w, s_w)[2]

# Compute mean and sd:
print(paste0("mu :", mu_w))
print(paste0("sigma :",s_w))

commut_mat <- mat_conect(N,alp_w,bet_w,DIST)

tau_ct <- 0
# Initial populations:
init_pop <- matrix(100, nrow = N)
print(paste0("beta - gamma:", bet[1] - (alp[1] + delt[1] + d_vec[1])))

#-----------------POPULATION INIT----------------------#
# Number of initial individuals by compartments:
SUS_INIT <- 100000 #Susceptible
INF_INIT <- 100    #Infected

# End time integration:
end_time <- 10

#-----------------------RAND(BETA) VS MEAN(BETA)---------------------------
#-------------------- MOBILITY ------------------#
d <- 50000
alphag <- alphagamma(1,4)
betag <- betagamma(1,4)
mu_bet <- rgamma(d,alphag,betag) 
sig_bet <- rgamma(d,alphag,betag) 
mat_rand_mean <-  matrix(0, ncol = 8, nrow = d)
for(i in c(1:d)){
  alphag <- alphagamma(mu_bet[i], sig_bet[i])
  betag <- betagamma(mu_bet[i], sig_bet[i])
  bet <- rgamma(N,alphag,betag)  # Transmission rate
  ### Migration:
  mu_c <- 0.01
  s_c <- 0.00001
  alp_c <- beta_a_b(mu_c, s_c)[1]
  bet_c <- beta_a_b(mu_c, s_c)[2]
  migrate_mat <- mat_conect(N,alp_c,bet_c,DIST)
  ### Commuting
  mu_w <- 0.1
  s_w <- 0.01
  alp_w <- beta_a_b(mu_w, s_w)[1]
  bet_w <- beta_a_b(mu_w, s_w)[2]
  commut_mat <- mat_conect(N,alp_w,bet_w,DIST)
  
  #-------------------------------RANDOM------------------------------------#
  # Compute the Jacobian matrix:
  
  jac_rand <- jacobian(N, bet, gamma_ct, commut_mat, migrate_mat, mu_c, MOB)
  eig_rand <- eigen_mat(jac_rand)
  out_rand <-  max(eig_rand$re)
  
  sol_rand <- int(N, del_N,bet,d_vec,thet,alp,delt,
                  commut_mat,migrate_mat,50,MOB, CTE_POP,
                  CTE_INF,SUS_INIT, INF_INIT,init_pop)
  
  print(paste0("i:",i))
  max_inf_rand <- max(sol_rand[,c((N+2):(2*N+1))])
  time_max_rand <- sol_rand[which(max(sol_rand[,c((N+2):(2*N+1))]) == max_inf_rand),1]
  
  #------------------------------MEAN(BETA)---------------------------------#
  bet <- rep(mean(bet),N)
  # Compute the Jacobian matrix:
  jac_mean <- jacobian(N, bet, gamma_ct, commut_mat, migrate_mat, mu_c, MOB)
  eig_mean <- eigen_mat(jac_mean)
  out_mean <-  max(eig_mean$re)
  
  sol_mean <- int(N, del_N,bet,d_vec,thet,alp,delt,
                  commut_mat,migrate_mat,50,MOB, CTE_POP,
                  CTE_INF,SUS_INIT, INF_INIT,init_pop)
  
  print(paste0("i:",i))
  max_inf_mean <- max(sol_mean[,c((N+2):(2*N+1))])
  time_max_mean <- sol_mean[which(max(sol_mean[,c((N+2):(2*N+1))]) == max_inf_mean),1]
  
  mat_rand_mean[i,] = c(mu_bet[i], sig_bet[i], 
                        out_rand, max_inf_rand, time_max_rand,
                        out_mean, max_inf_mean, time_max_mean)
}

df.bet <- as.data.frame(mat_rand_mean)
path <- paste0("~/RMT/Integration/random_mean_bet_",Sys.Date(),".csv")
write.csv(df.bet,path, row.names = TRUE)