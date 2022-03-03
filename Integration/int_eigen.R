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
  
  # bet <- abs(rnorm(N,1,1))  # Transmission rate
  d_vec <- rep(1, N) # Natural mortality rate
  thet <- rep(2.6, N) # Rate of loss of immunity
  alp <- rep(1, N) # Rate of disease overcome
  delt <- rep(0, N) # Diseases related mortality rate
  gamma_ct <-  alp[1] + delt[1] + d_vec[1]
  print(paste0("gamma:", alp[1] + delt[1] + d_vec[1]))
  print(paste0("beta - gamma:", bet[1] - (alp[1] + delt[1] + d_vec[1])))

#-------------------- MOBILITY ------------------
### Migration:
  mu_c <- 0.0001
  s_c <- 0.00005
  alp_c <- beta_a_b(mu_c, s_c)[1]
  bet_c <- beta_a_b(mu_c, s_c)[2]
  # Compute mean and sd:
  print(paste0("mu :", mu_c))
  print(paste0("sigma :", s_c))

  migrate_mat <- mat_conect(N,alp_c,bet_c,DIST)
### Commuting
  mu_w <- 0.2
  s_w <- 0.05
  alp_w <- beta_a_b(mu_w, s_w)[1]
  bet_w <- beta_a_b(mu_w, s_w)[2]
  
  # Compute mean and sd:
  print(paste0("mu :", mu_w))
  print(paste0("sigma :",s_w))

  commut_mat <- mat_conect(N,alp_w,bet_w,DIST)

tau_ct <- 0
  # Initial populations:
  init_pop <- matrix(100, nrow = N)
  init_pop <- abs(round(rnorm(N, 1000000, 800000),0))
  if(CTE_POP_Patch == 1){
    list_param <- diff_f(migrate_mat,d_vec[1,1],init_pop)
    d_vec[1,1] <- list_param[1]
    del_N <- list_param[2:(N+1)]
    d_vec[1,1] = 1.7
  }
  # init_pop <- matrix(rgamma(N,shape = (mu_p/s_p)^2,rate = mu_p/(s_p^2)), nrow = N)
  
  
  print(paste0("beta - gamma:", bet[1] - (alp[1] + delt[1] + d_vec[1])))
  
#-----------------POPULATION INIT----------------------#
  # Number of initial individuals by compartments:
  SUS_INIT <- 100000 #Susceptible
  INF_INIT <- 500    #Infected

# End time integration:
end_time <- 50

#----------------------------CTE BETA------------------------------------
bet_cte <- 0.1
bet <- rep(bet_cte, N)  # Transmission rate
# alphag <- alphagamma(1,4)
# betag <- betagamma(1,4)
# bet <- rgamma(N,alphag,betag)  # Transmission rate
d_vec <- rep(0.3, N) # Natural mortality rate
del_N <- rep(0.6, N) # Birth rate
thet <-  rep(0.3, N)# Rate of loss of immunity
alp <-  rep(0.5, N)# Rate of disease overcome
delt <- rep(0.15, N) # Diseases related mortality rate
gamma_ct <-  alp[1] + delt[1] + d_vec[1]
print(paste0("gamma:", alp[1] + delt[1] + d_vec[1]))
print(paste0("beta - gamma:", bet[1] - (alp[1] + delt[1] + d_vec[1])))
beta_ct = bet_cte  # beta
gamma_ct = alp[1] + delt[1] + d_vec[1]     # gamma   

sol <- int(N, del_N,bet,d_vec,thet,alp,delt,
           commut_mat,migrate_mat,end_time,
           MOB, CTE_POP, CTE_INF,SUS_INIT, INF_INIT,init_pop)

sol_df <-  as.data.frame(sol)
# Plot the susceptible, infected and recovered:
state <- "INF"

# Parameters:
print(paste0("beta - gamma", beta_ct - gamma_ct))
# 
cond_gen(N, mu_c, s_c, mu_w, s_w, gamma_ct, beta_ct, tau_ct)
# Make distribution:
jac <- jacobian(N,bet,gamma_ct, commut_mat, migrate_mat,mu_c, MOB)
eig <- eigen_mat(jac)

plot_eig_cte <- ggplot(eig) + geom_point(aes(re,im), size = 0.05) 
plot_eig_cte + coord_fixed() 

# Predicted distribution for constant beta: 
rad <- pred_radius(N, beta_ct, gamma_ct, tau_ct, mu_c, sqrt(s_c), mu_w, sqrt(s_w), MOB)
cent <- pred_center(N, beta_ct, gamma_ct, tau_ct, mu_c, sqrt(s_c), mu_w,sqrt( s_w), MOB)
outl <- pred_outlier(N, beta_ct, gamma_ct, mu_w, MOB)

plot_eigen_cte_pred <- plot_eigen(eig, cent, rad, outl, 1) +
  geom_point(aes(outl,0), colour =  "blue",
             show.legend = NA)
plot_eigen_cte_pred + ggtitle(paste0("N: ",N,"\n", "gamma: ", gamma_ct, ", beta:", bet_cte,"\n",
                              "mu_w: ", mu_w, ", s_w: ", s_w,"\n",
                              "mu_c: ", mu_c, ", s_c: ", s_c))

# Plot integration:
state <- "INF"
cte_bet <- plot_int(N, sol, state) + 
  theme_bw() + xlim(c(0,10))

low_var_eigen <- low_var_eigen + labs(title="a" )
high_var_eigen <- high_var_eigen + labs(title="b" )
ggeigen <- ggarrange(low_var_eigen, high_var_eigen, ncol = 1, nrow = 2)
low_var_inf <- low_var_inf + labs(title="c" ) + xlim(c(0,3))
high_var_inf <- high_var_inf + labs(title="d" ) + ylab("")  + xlim(c(0,3))
gginf <- ggarrange( low_var_inf, high_var_inf, ncol = 2, nrow = 1)
ggfull <- ggarrange(ggeigen, gginf, ncol = 1, nrow=2)
low_var_eigen <- low_var_eigen + scale_y_continuous(breaks = c(0))
ggarrange(cte_bet + ylim(c(0, 100000)), rand_bet + ylim(c(0, 100000)), align = "h")
#-----------------------RANDOM PARAMETERS-------------------------#
plot_inf_d <- plot_int(N, sol, state) + 
  theme_bw() +
  xlim(c(0,4)) + ggtitle(""*mu[d]~": , "*sigma[d]~": 4 ") 

plot_inf_theta <- plot_int(N, sol, state) + 
  theme_bw() +
  xlim(c(0,4)) + ggtitle(""*mu[theta]~": , "*sigma[theta]~": 4 ") 

plot_inf_alp <- plot_int(N, sol, state) + 
  theme_bw() +
  xlim(c(0,6)) + ggtitle(""*mu[alpha]~": , "*sigma[alpha]~": 4 ") 

plot_inf_delt <- plot_int(N, sol, state) + 
  theme_bw() +
  xlim(c(0,6)) + ggtitle(""*mu[delta]~": , "*sigma[delta]~": 4 ")

plot_inf_beta <- plot_int(N, sol, state) + 
  theme_bw() +
  xlim(c(0,10)) + ggtitle(""*mu[beta]~": , "*sigma[beta]~": 4 ")

gg_random <- ggarrange(plot_inf_d + xlab(""),plot_inf_theta  + xlab("") +ylab(""),
          plot_inf_alp   ,plot_inf_delt +ylab("")) +
  ggtitle("Random parameters")
#-----------------------PLOT EXIT RATES----------------------------##
plot_inf_cte_delt_4.6 <-  plot_inf_cte +
  ggtitle(paste0("delta",":4.6")) +
  xlim(c(0,20))  
# + xlim(c(0,10))

# plot_inf_cte_del_4.4 <-  plot_inf_cte
#  
# gg.alp1 <- ggarrange(plot_inf_cte_alp_0.6,plot_inf_cte_alp_1.6)
# gg.alp2 <- ggarrange(plot_inf_cte_alp_10.6,plot_inf_cte_alp_11.6)
# gg_alp_inf <-  ggarrange(gg.alp1,gg.alp2, ncol = 1)

# gg.delt1 <- ggarrange(plot_inf_cte_delt_0,plot_inf_cte_delt_0.6)
# gg.delt2 <- ggarrange(plot_inf_cte_delt_6.6,plot_inf_cte_delt_11.6)
# gg_delt_inf <-  ggarrange(gg.delt1,gg.delt2, ncol = 1)

# gg.d1 <- ggarrange(plot_inf_cte_d_0.6,plot_inf_cte_d_1.6)
# gg.d2 <- ggarrange(plot_inf_cte_d_4.6,plot_inf_cte_d_11.6)
# gg_d_inf <-  ggarrange(gg.d1,gg.d2, ncol = 1)

gg.theta1 <- ggarrange(plot_inf_cte_theta_0.006,
                       plot_inf_cte_theta_0.6)
gg.theta2 <- ggarrange(plot_inf_cte_theta_1.6,
                       plot_inf_cte_theta_4.6)
gg_theta_inf <-  ggarrange(gg.theta1,gg.theta2, ncol = 1)

#STOP#
plot_inf_cte_d_1 <- plot_inf_cte_d_0.6 + 
  xlim(c(0,5)) + ggtitle("d: 0.6") + xlab("")
plot_inf_cte_d_2 <- plot_inf_cte_d_4.6 + 
  xlim(c(0,5)) + ggtitle("d: 4.6") + xlab("") + ylab("")

plot_inf_cte_alp_1 <- plot_inf_cte_alp_0.6 + 
  xlim(c(0,5)) + ggtitle(""*alpha~": 0.6") + xlab("")
plot_inf_cte_alp_2 <- plot_inf_cte_alp_4.6 + 
  xlim(c(0,5)) + ggtitle(""*alpha~": 4.6") + xlab("") + ylab("")

plot_inf_cte_delt_1 <- plot_inf_cte_delt_0.6 + 
  xlim(c(0,5)) + ggtitle(""*delta~": 0.6")
plot_inf_cte_delt_2 <- plot_inf_cte_delt_4.6  + ylab("") +
  xlim(c(0,5)) + ggtitle(""*delta~": 4.6")


gg_exit <- ggarrange(plot_inf_cte_d_1,plot_inf_cte_d_2,
          plot_inf_cte_alp_1,plot_inf_cte_alp_2,
          plot_inf_cte_delt_1,plot_inf_cte_delt_2,
        ncol = 2, nrow = 3)

#---------------------THETA-------------------------------#
plot_inf_cte_thet_1 <- plot_inf_cte_theta_0.006 + 
  xlim(c(0,5)) + ggtitle(""*theta~":0.006") 
plot_inf_cte_thet_2 <- plot_inf_cte_theta_1.6 + 
  xlim(c(0,5)) + ggtitle(""*theta~":1.6") + ylab("")

ggtheta <- ggarrange(plot_inf_cte_thet_1,plot_inf_cte_thet_2 )

#---------------------------------------------------------#
gg.d <- ggarrange(plot_inf_cte_del_0.6,plot_inf_cte_del_4.4)
gg_del_inf <-  ggarrange(gg.d,plot_inf_cte_del_8.4, ncol = 1)


plot_eigen_cte_pred_sw_0.05 <- plot_eigen_cte_pred_sw_0.05 + ggtitle("") +
  labs(title="a" )+ theme(legend.position = "none") + xlab("real")
plot_inf_cte_sw_0.05 <- plot_inf_cte_sw_0.05 + ggtitle("") +
  labs(title="c" )+ theme(legend.position = "none") + xlab("time")


plot_eigen_cte_pred <- plot_eigen_cte_pred + ggtitle("") +
  labs(title="b" )+ theme(legend.position = "none") + xlab("real") + ylab("")
plot_inf_cte <- plot_inf_cte + ggtitle("") +
  labs(title="d" )+ theme(legend.position = "none")+ ylab("")

ggarr1 <-  ggarrange(plot_inf_cte_sw_0.05,plot_inf_cte, ncol = 2, nrow = 1)
ggarr2 <-  ggarrange(plot_eigen_cte_pred_sw_0.05,plot_eigen_cte_pred, ncol = 2, nrow = 1)
gg_full <-  ggarrange(ggarr2,ggarr1, ncol = , nrow = 2)

#--------------------------SAVE FILE--------------------------------#
gamma_ct_w <- format(round(gamma_ct,2), decimal.mark = ',')
beta_ct_w <- format(round(beta_ct,2), decimal.mark = ',')
mu_w_w <- format(round(mu_w,2), decimal.mark = ',')
s_w_w <- format(round(s_w,2), decimal.mark = ',')
mu_c_w <- format(round(mu_c,2), decimal.mark = ',')
s_c_w <- format(round(s_c,2), decimal.mark = ',')

Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/"
Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/"

path <- paste0(Path,"rand_param","N",N,"b",beta_ct_w,"mw",
               mu_w_w,"sw1","0.01","sw2","0.2","mm",mu_c_w,"sm",s_c_w,".png")
ggsave(path,
       plot =ggfull, device = "png")

#-----------------------SAVE FILE---------------------------#
Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/Gen/"
Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/"
gamma_ct_w <- format(round(gamma_ct,2), decimal.mark = ',')
beta_ct_w <- format(round(beta_ct,2), decimal.mark = ',')
mu_w_w <- format(round(mu_w,2), decimal.mark = ',')
s_w_w <- format(round(s_w,2), decimal.mark = ',')
mu_c_w <- format(round(mu_c,2), decimal.mark = ',')
s_c_w <- format(round(s_c,2), decimal.mark = ',')
path <- paste0(Path,"gen","N",N_w,"g",gamma_ct_w,"b",beta_ct_w,"mw",
       mu_w_w,"b_new_1p", i,"sw",s_w_w,"mm",mu_c_w,"sm",s_c_w,".png")
# path <- paste0(Path,"mod1patchtrans","N",N,"g",gamma_ct,"b",
               # beta_ct,"mw","bnew","14",
               # mu_w,"sw",s_w,"mm",mu_c,"sm",s_c,".png")
png(file = path, width = 8000, height = 6000, res = 1100)
plot_outl_inf
dev.off()
