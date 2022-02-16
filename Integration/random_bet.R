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
  end_time <- 100
#--------------------------Integration-----------------------------------
  mu <-  9 
  sigma <-  30
  alp <-  (mu^2/sigma)
  bet <- mu/sigma
  bet <-  rgamma(N,alp, bet)
  sol <- int(N, del_N,bet,d_vec,thet,alp,delt,
             commut_mat,migrate_mat,end_time,
             MOB, CTE_POP, CTE_INF,SUS_INIT, INF_INIT,init_pop)
  
  sol.df <-  as.data.frame(sol)
  
  # Plot the infected:
  state <- "INF"
  plot.inf.k <- plot_int(N, sol, state) +
    theme_bw() + xlim(c(0,20))
  plot.inf.k
  
  # Compute the jacobian matrix:
  jac <- jacobian(N,bet,gamma_ct, commut_mat, migrate_mat,mu_c, MOB)
  eig <- eigen_mat(jac)
  
  plot.eig.rand <- ggplot(eig) + geom_point(aes(re,im), size = 0.05)
  plot.eig.rand + coord_fixed()
  
  rad <- pred_radius(N, mean(bet), gamma_ct, tau_ct,
                     mu_c, sqrt(s_c), mu_w, sqrt(s_w), MOB)
  cent <- pred_center(N, mean(bet), gamma_ct, tau_ct,
                      mu_c, sqrt(s_c), mu_w,sqrt(s_w), MOB)
  outl <- pred_outlier(N, mean(bet), gamma_ct, mu_w, MOB)
  
  plot.eigen.rand <- plot_eigen(eig, cent, rad, outl, MOB)
  
  plot.eigen.rand
  
  sol <- int(N, del_N,bet,d_vec,thet,alp,delt,
             commut_mat,migrate_mat,end_time,
             MOB, CTE_POP, CTE_INF,SUS_INIT, 
             INF_INIT,init_pop)
  
  sol.df <-  as.data.frame(sol)
  # Plot the susceptible, infected and recovered:
  state <- "INF"
  plot.inf <- plot_int(N, sol, state) + 
    theme_bw() 
  
  plot.inf
  
  