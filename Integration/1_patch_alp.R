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

#--------------------------CHANGING BETA-------------------------------
# Change parameters to characters with decimal coma:
beta_ct = bet_cte  # beta
gamma_ct = alp[1] + delt[1] + d_vec[1]     # gamma   
gamma_ct_w <- format(round(gamma_ct,2), decimal.mark = ',')
beta_ct_w <- format(round(beta_ct,2), decimal.mark = ',')
mu_w_w <- format(round(mu_w,2), decimal.mark = ',')
s_w_w <- format(round(s_w,2), decimal.mark = ',')
mu_c_w <- format(round(mu_c,2), decimal.mark = ',')
s_c_w <- format(round(s_c,2), decimal.mark = ',')

# Reinitialize the beta to cte vector:
d_vec <- rep(1.2, N) # Natural mortality rate
alp <- rep(0.5, N) # Rate of disease overcome
delt <- rep(0, N) # Diseases related mortality rate
gamma_ct = d_vec[1] + alp[1] +  delt[1]   # gamma   

bet <- rep(bet_cte, N)  # Transmission rate

print(paste0("beta-gamma:", bet_cte - gamma_ct))
count = 1
alp_vec <- seq(0,10,0.1)
l <-  length(alp_vec) + 1
plot_list <- list()
ind <-  sample(1:N,1)
for(i in alp_vec){
  print(paste0("New beta : ",i))
  bet_new <- i
  bet[ind] <- bet_cte + bet_new
  
  sol <- int(N, del_N,bet,d_vec,thet,alp,delt,
             commut_mat,migrate_mat,end_time,
             MOB, CTE_POP, CTE_INF,SUS_INIT, INF_INIT,init_pop)
  
  sol_df <-  as.data.frame(sol)
  
  
  # Plot the susceptible, infected and recovered:
  state <- "INF"
  plot_inf_1 <- plot_int(N, sol, state) +
    theme_bw() + xlim(c(0,20))
  plot_inf_1
  
  # Parameters:
  print(paste0("beta - gamma", beta_ct - gamma_ct))
  
  cond_gen(N, mu_c, s_c, mu_w, s_w, gamma_ct, beta_ct, tau_ct)
  # Make distribution:
  jac <- jacobian(N,bet,gamma_ct, commut_mat, migrate_mat,mu_c, MOB)
  eig <- eigen_mat(jac)
  
  plot_eig <- ggplot(eig) + geom_point(aes(re,im), size = 0.05)
  plot_eig + coord_fixed()
  
  # Compute the outliers, center and radius:
  a <- bet_cte*mu_w + mu_c
  b <- bet_new
  c <- bet_new*mu_w
  rad <- pred_radius(N, beta_ct, gamma_ct, tau_ct, mu_c, sqrt(s_c), mu_w, sqrt(s_w), MOB)
  cent <- pred_center(N, beta_ct, gamma_ct, tau_ct, mu_c, sqrt(s_c), mu_w,sqrt( s_w), MOB)
  
  outl <- (1/2)*(N*a + b + sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
  outl2 <- (1/2)*(N*a + b - sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
  outl <- outl + (bet_cte*(1-mu_w) - N*mu_c - gamma_ct)
  outl2 <- outl2 + (bet_cte*(1-mu_w) - N*mu_c - gamma_ct)
  
  
  # If last parameter is 1 he outlier is not computed:
  plot_eigen_1 <- plot_eigen(eig, cent, rad, outl, 1)+
    geom_point(aes(outl,0), colour =  "blue",
               show.legend = NA) +
    geom_point(aes(outl2,0), colour =  "blue",
               show.legend = NA)
  
  plot_eigen_1
  
  # Create a vector to colour differently the time series related to the patch with
  # different beta:
  vec_col <-  vector(mode="character", length=N)
  vec_col[1:N] <- "red4"
  # vec_col[vec_rand] <- "royalblue3"
  vec_col[ind] <- "royalblue3"
  
  plot_inf_1 <-  plot_inf_1 +
    scale_colour_manual(values = vec_col)
  
  plot_inf_1
  plot_inf_1_lim <- plot_inf_1 +
    xlim(c(0,30)) +
    # scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    ggtitle(paste(expression(alpha),":",i))
  
  plot_list[[count]] <-  plot_inf_1_lim
  count = count + 1
}

png_files <-  c()
for(i in c(1:(count-1))){
  # png_files[i] <- paste0("/home/marta/Documentos/PHD/2022/RMT_SIR/Plots/1patch/High_both///genN100g0,82b0mw0,5sw0,46mm0,5sm0,46_",i+1,".png")
  plot_list[[i]] +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE), limits = c(0,20000))
  ggsave(paste0("~/Documents/PHD/2022/RMT_SIR/Plots/kpatch/High_both/plot_",i,".png" ),
         plot = last_plot(), device = "png")
  png_files[i] <-  paste0("~/Documents/PHD/2022/RMT_SIR/Plots/kpatch/High_both/plot_",i,".png" )
}

# png_files <- list.files("~/Documents/PHD/2022/RMT_SIR/Plots/1patch/test/", pattern = ".*png$", full.names = TRUE)
# gifski(png_files, gif_file = "~/Documentos/PHD/2022/RMT_SIR/Plots/1patch/High_both/animation.gif", width = 800, height = 600, delay = 0.3)
gifski(png_files, gif_file = "~/Documents/PHD/2022/RMT_SIR/Plots/kpatch/High_both/animation.gif", width = 800, height = 600, delay = 0.3)

