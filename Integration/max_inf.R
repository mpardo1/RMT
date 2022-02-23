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
N = 50 # Number of patches
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

##--------------------------------------------------------------------------------#

#
#-----------------------TIME FOR MAX INFECTED---------------------------
alpha_r0_1 <-  function(alp_p){
  a <- bet_cte*mu_w + mu_c
  b <- alp_p
  c <- alp_p*mu_w
  outl <- (1/2)*(N*a + b + sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
  outl <- outl + (bet_cte*(1-mu_w) - N*mu_c - gamma_ct)
  return(outl[1])
}

# Compute the value of alpha for the bifurcationn point,
# The numbers inside de c() are the interval to look for roots:
uni <- uniroot(alpha_r0_1, c(0, 20))$root

bet_cte = 0.1
bet <- rep(bet_cte, N)  # Transmission rate
beta_ct = bet_cte  # beta
d_vec <- rep(1, N) # Natural mortality rate
alp <- rep(1, N) # Rate of disease overcome
delt <- rep(0, N) # Diseases related mortality rate

gamma_ct = d_vec[1] + alp[1] +  delt[1]   # gamma   
count = 1
alp_vec <- seq(0,10,0.1)
l <-  length(alp_vec) + 1
mat_max_inf <-  matrix(0, ncol = 3, nrow = l)
ind <- sample(1:N,1)
for(i in alp_vec){
  print(paste0("New beta : ",i))
  bet_new <- i
  bet[ind] <- bet_cte + bet_new
  
  sol <- int(N, del_N,bet,d_vec,thet,alp,delt,
             commut_mat,migrate_mat,end_time,
             MOB, CTE_POP, CTE_INF,SUS_INIT, INF_INIT,init_pop)
  
  sol_df <-  as.data.frame(sol)
  # Compute the time where the sum of infected individuals is maximum,
  inf_df <- sol_df[, c(1,(N+2):(2*N+1))]
  inf_df$sum <- rowSums(inf_df[,2:(N+1)])
  t_max_inf <- inf_df$time[inf_df$sum==max(inf_df$sum)]
  if(length(t_max_inf) > 1){
    t_max_inf <- t_max_inf[1]
  }
  mat_max_inf[count,] <- c(i, t_max_inf,max(inf_df$sum))
  print(paste0("alp , time:",mat_max_inf[count,]))
  count = count + 1
}

len <- length(alp_vec)
mat_max_inf_df <-  data.frame(alp = alp_vec ,time = mat_max_inf[1:len,2]
                              ,max_inf = mat_max_inf[1:len,3])

write.csv(mat_max_inf_df,"~/Documentos/PHD/2022/RMT_SIR/Plots/epi_param/Alpha.csv", row.names = TRUE)

plot_time <- ggplot() + 
  geom_line(data = mat_max_inf_df,aes(alp,time)) + theme_bw() +
  ylab("Time of maximum infected individuals") + xlab(expression(alpha)) +
  ggtitle(paste0("N: ",N,"\n", "gamma: ", gamma_ct, ", beta:", bet_cte,"\n", 
                 "mu_w: ",mu_w, ", s_w: ", s_w,"\n",
                 "mu_c: ", mu_c, ", s_c: ", s_c))

plot_time

plot_max_inf <- ggplot(mat_max_inf_df) + 
  geom_line(aes(alp,max_inf)) + theme_bw() +
  ylab("Number of maximum infected individuals") + xlab(expression(alpha))

plot_max_inf

#---------------------S_w vs MAX INF ----------------------------------
count = 1
s_w_vec <- seq(0.01,0.24,0.01)
l <-  length(s_w_vec) + 1
mat_max_inf <-  matrix(0, ncol = 3, nrow = l)
plot_list <- list()
d_vec <- rep(1, N) # Natural mortality rate
alp <- rep(1, N) # Rate of disease overcome
delt <- rep(0, N) # Diseases related mortality rate

gamma_ct = d_vec[1] + alp[1] +  delt[1]   # gamma   
for(i in c(1:l)){
  alp_w <- beta_a_b(mu_w, s_w_vec[i])[1]
  bet_w <- beta_a_b(mu_w, s_w_vec[i])[2]
  
  # Compute mean and sd:
  print(paste0("mu :", mu_w))
  print(paste0("sigma :", s_w_vec[i]))
  
  commut_mat <- mat_conect(N,alp_w,bet_w,DIST)
  sol <- int(N, del_N,bet,d_vec,thet,alp,delt,
             commut_mat,migrate_mat,end_time,
             MOB, CTE_POP, CTE_INF,SUS_INIT, INF_INIT,init_pop)
  
  sol_df <-  as.data.frame(sol)
  # Compute the time where the sum of infected individuals is maximum,
  inf_df <- sol_df[, c(1,(N+2):(2*N+1))]
  inf_df$sum <- rowSums(inf_df[,2:(N+1)])
  t_max_inf <- inf_df$time[inf_df$sum==max(inf_df$sum)]
  if(length(t_max_inf) > 1){
    t_max_inf <- t_max_inf[1]
  }
  mat_max_inf[i,] <- c(s_w_vec[i], t_max_inf,max(inf_df$sum))
  print(paste0("sigma:",mat_max_inf[i,1], "  "," time:",mat_max_inf[i,2], 
               " "," max_inf:",mat_max_inf[i,3]))
  
  jac <- jacobian(N,bet,gamma_ct, commut_mat, migrate_mat,mu_c, MOB)
  eig <- eigen_mat(jac)
  
  plot_eig_cte <- ggplot(eig) + geom_point(aes(re,im), size = 0.05) 
  plot_eig_cte + coord_fixed() 
  
  # Predicted distribution for constant beta: 
  rad <- pred_radius(N, beta_ct, gamma_ct, tau_ct, mu_c, sqrt(s_c), mu_w, sqrt(s_w_vec[i]), MOB)
  cent <- pred_center(N, beta_ct, gamma_ct, tau_ct, mu_c, sqrt(s_c), mu_w,sqrt( s_w_vec[i]), MOB)
  outl <- pred_outlier(N, beta_ct, gamma_ct, mu_w, MOB)
  
  plot_eigen_pred <- plot_eigen(eig, cent, rad, outl, 1) +
    geom_point(aes(outl,0), colour =  "blue",
               show.legend = NA) +
    ggtitle((paste("mu[w]", ": ",s_w_vec[i])))
  plot_eigen_pred + ggtitle(paste0("N: ",N,"\n", "gamma: ", gamma_ct, ", beta:", bet_cte,"\n",
                                   "mu_w: ", s_w_vec[i], ", s_w: ", s_w,"\n",
                                   "mu_c: ", mu_c, ", s_c: ", s_c))
  
  state <- "INF"
  plot_inf <- plot_int(N, sol, state) + 
    theme_bw() + xlim(c(0,20)) + ylim(c(0,90000))
  
  plot_list[[i]] <- ggarrange(plot_eigen_pred, plot_inf, nrow = 2, ncol = 1) 
  
  count = count + 1
}

ggarrange(plot_list[[1]],  plot_list[[24]])
len <- length(s_w_vec)
mat_max_inf_df <-  data.frame(sig_w = s_w_vec ,time = mat_max_inf[1:len,2]
                              ,max_inf = mat_max_inf[1:len,3])
plot_time <- ggplot() + 
  geom_line(data = mat_max_inf_df,aes(sig_w,time)) +
  # # geom_line(data = mat_max_inf_df,aes(alp,max_inf)) +
  # geom_line(data = alp_df, aes( sig_w, outl), colour = "blue")  +
  # geom_segment(aes(x = 0, y = 0, xend = max_alp, yend =0,
  #                  colour = "segment"), linetype=2,
  #              show.legend = FALSE) +
  ylab("Time") + xlab(expression(alpha)) + theme_bw() 
plot_time

plot_max_inf <- ggplot(mat_max_inf_df) + 
  geom_line(aes(sig_w,max_inf)) + theme_bw() 

plot_max_inf

#-------------------EIGEN vs MAGNITUDE EPI----------------------
N = 100 # Number of patches
# CTE parameters:
del_N <- rep(0.6, N) # Birth rate
# bet_cte <- 6
d_vec <- rep(0.6, N) # Natural mortality rate
thet <- rep(2.6, N) # Rate of loss of immunity
alp <- rep(1, N) # Rate of disease overcome
delt <- rep(0, N) # Diseases related mortality rate
gamma_ct <-  alp[1] + delt[1] + d_vec[1]
print(paste0("gamma:", alp[1] + delt[1] + d_vec[1]))
print(paste0("beta - gamma:", bet[1] - (alp[1] + delt[1] + d_vec[1])))

#-------------------- MOBILITY ------------------#
### Migration:
mu_c <- 0.85
s_c <- 0.1
alp_c <- beta_a_b(mu_c, s_c)[1]
bet_c <- beta_a_b(mu_c, s_c)[2]
migrate_mat <- mat_conect(N,alp_c,bet_c,DIST)
### Commuting
mu_w <- 0.1
s_w <- 0.01
alp_w <- beta_a_b(mu_w, s_w)[1]
bet_w <- beta_a_b(mu_w, s_w)[2]
commut_mat <- mat_conect(N,alp_w,bet_w,DIST)

max_bet <- 10
bet_vec <- seq(0,max_bet,0.1)
leng <-  length(bet_vec)
mat.bet <-  matrix(0, ncol = 4, nrow = leng)
for(i in c(1:N)){
  # Compute  right most eigenvalue for 1 patch:
  bet_cte <-  bet_vec[i]
  bet <-  rep(bet_cte,N)
  # Compute the jacobian matrix:
  jac <- jacobian(N, bet, gamma_ct, commut_mat, migrate_mat, mu_c, MOB)
  eig <- eigen_mat(jac)
  out1 <-  max(eig$re)
  
  sol <- int(N, del_N,bet,d_vec,thet,alp,delt,
             commut_mat,migrate_mat,50,MOB, CTE_POP,
             CTE_INF,SUS_INIT, INF_INIT,init_pop)
  
  print(paste0("i:",i))
  max_inf <- max(sol[,c((N+2):(2*N+1))])
  time_max <- sol[which(max(sol[,c((N+2):(2*N+1))]) == max_inf),1]
  mat.bet[i,] = c(bet_cte, out1, max_inf, time_max)
}

df.bet <-  as.data.frame(mat.bet)
colnames(df.bet) <- c("beta", "outlier","max_inf", "time_max")

POP_INIT <- SUS_INIT + INF_INIT
plot_eigen <- ggplot(df.bet) + geom_line(aes(beta,outlier))
plot_inf <- ggplot(df.bet) + geom_line(aes(beta,max_inf)) +
  geom_segment(aes(x = 0, y = POP_INIT, xend = max_bet, yend = POP_INIT,
                   colour = "segment"), linetype = 2, show.legend = FALSE)
ggarrange(plot_eigen,plot_inf)

gg_outl_inf <- ggplot(df.bet) + geom_line(aes(outlier,max_inf))

path <- paste0("~/Documents/PHD/2022/RMT_SIR/Plots/","N: ",N,"\n", "gamma: ", gamma_ct, ", beta:", bet_cte,"\n", 
       "mu_w: ",mu_w, ", s_w: ", s_w,"\n",
       "mu_c: ", mu_c, ", s_c: ", s_c)
#------------------CI vs MAGNITUDE EPI----------------------
N = 100 # Number of patches
# CTE parameters:
del_N <- rep(0.6, N) # Birth rate
# bet_cte <- 6
d_vec <- rep(0.6, N) # Natural mortality rate
thet <- rep(2.6, N) # Rate of loss of immunity
alp <- rep(1, N) # Rate of disease overcome
delt <- rep(0, N) # Diseases related mortality rate
gamma_ct <-  alp[1] + delt[1] + d_vec[1]
print(paste0("gamma:", alp[1] + delt[1] + d_vec[1]))
print(paste0("beta - gamma:", bet[1] - (alp[1] + delt[1] + d_vec[1])))
CTE_POP <- 1

#-------------------- MOBILITY ------------------#
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

max_inf <- 500000
inf_vec <- seq(10,max_inf,1000)
leng <-  length(inf_vec)
mat.inf <-  matrix(0, ncol = 3, nrow = leng)
bet_cte <-  1
bet <-  rep(bet_cte,N)
plot_list <- list()
for(i in c(1:leng)){
  INF_INIT <- inf_vec[i]
  # Compute  right most eigenvalue for 1 patch:
  sol <- int(N, del_N,bet,d_vec,thet,alp,delt,
             commut_mat,migrate_mat,50,MOB, CTE_POP,
             CTE_INF,SUS_INIT, INF_INIT,init_pop)
  
  print(paste0("i:",i))
  max_inf <- max(sol[,c((N+2):(2*N+1))])
  time_max <- sol[which(max(sol[,c((N+2):(2*N+1))]) == max_inf)]
  mat.inf[i,] = c(INF_INIT, max_inf, time_max)
  plot_list[[i]] <-  plot_int(N, sol, state) + 
    theme_bw() 
}

df.inf <-  as.data.frame(mat.inf)
colnames(df.inf) <- c("inf_init","max_inf", "time_max")

POP_INIT <- rowSums(sol[,2:(ncol(sol))])[1]
plot_eigen <- ggplot(df.inf) + geom_line(aes(inf_init,time_max))
plot_inf <- ggplot(df.inf) + geom_line(aes(inf_init,max_inf)) +
  theme_bw() + xlab("Initial number of infected individuals") +
  ylab("Maximum number of infected individuals")
ggarrange(plot_eigen,plot_inf)

