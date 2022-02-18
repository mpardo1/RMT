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
  INF_INIT <- 100    #Infected

# End time integration:
end_time <- 50
#----------------------------CTE BETA------------------------------------
bet_cte <-  0.3
bet <- rep(bet_cte, N)  # Transmission rate
d_vec <- rep(0.8, N) # Natural mortality rate
del_N <- rep(0.6, N) # Birth rate
thet <- rep(0.6, N) # Rate of loss of immunity
alp <- rep(0.5, N) # Rate of disease overcome
delt <- rep(0, N) # Diseases related mortality rate
gamma_ct <-  alp[1] + delt[1] + d_vec[1]
print(paste0("gamma:", alp[1] + delt[1] + d_vec[1]))
print(paste0("beta - gamma:", bet[1] - (alp[1] + delt[1] + d_vec[1])))
beta_ct = bet_cte  # beta
gamma_ct = alp[1] + delt[1] + d_vec[1]     # gamma   
bet <- rep(bet_cte, N)  # Transmission rate
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
plot_int(N, sol, state) + 
  theme_bw() 
plot_inf_cte <- plot_int(N, sol, state) + 
  theme_bw()
# + xlim(c(0,10))

# plot_inf_cte_del_4.4 <-  plot_inf_cte
#  
# gg.d <- ggarrange(plot_inf_d_0.4,plot_inf_d_1.8)
# gg_d_inf <-  ggarrange(gg.d,plot_inf_d_8.8, ncol = 1)

#STOP#

plot_inf_cte_alp_0.5 <-  plot_inf_cte +
  ggtitle("alpha:0.5") +
  xlim(c(0,5))

plot_inf_ct_alp_2.4 <-  plot_inf_cte_del_2.4 +
  ggtitle("alpha:2.4") +
  xlim(c(0,5))

plot_inf_cte_alp_4.4 <-  plot_inf_cte_del_4.4 +
  ggtitle("alpha:4.4") +
  xlim(c(0,5))

plot_inf_cte_alp_8.4 <-  plot_inf_cte_del_8.4 +
  ggtitle("alpha:8.4") +
  xlim(c(0,5))


plot_inf_cte_del_0.6 <-  plot_inf_cte_del_0.6 +
  ggtitle("delta:0.6") +
  xlim(c(0,5))

plot_inf_cte_del_4.4 <-  plot_inf_cte_del_4.4 +
  ggtitle("delta:4.4") +
  xlim(c(0,5))

plot_inf_cte_del_8.4 <-  plot_inf_cte_del_8.4 +
  ggtitle("delta:8.4") +
  xlim(c(0,5))

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
Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/"

path <- paste0(Path,"sigma_radius","N",N,"g",gamma_ct_w,"b",beta_ct_w,"mw",
               mu_w_w,"sw1",s_w_w,"sw2","0.05","mm",mu_c_w,"sm",s_c_w,".png")
ggsave(path,
       plot =gg_full, device = "png")

#---------------------------------------------------------------------#
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


#---------------------RANDOM BETAS------------------------
bet <- abs(rnorm(N,1,1))  # Transmission rate
sol <- int(N, del_N,bet,d_vec,thet,alp,delt,
           commut_mat,migrate_mat,end_time,
           MOB, CTE_POP, CTE_INF,SUS_INIT, INF_INIT,init_pop)
sol_df <- as.data.frame(sol)
sol_df <-  sol_df[,c(1,(N+2):(2*N+1))]
vec <- as.character(c(1:N))
colnames(sol_df) <-  c("time",vec)
head(sol_df)
df_plot <- reshape2::melt(sol_df, id.vars = c("time"))
bet_df <-  data.frame(variable = c(1:length(bet)), bet = bet)
df_plot <- merge(df_plot,bet_df, by = "variable", all.x = TRUE )
colnames(df_plot)[4] <- expression(beta)
df_plot$norm <-  df_plot$bet/max(df_plot$bet)
mid<-mean(df_plot$beta)
# Plot gradient lines infected:
plot  <- ggplot(df_plot,aes(time, value)) + 
  geom_line(aes( colour = beta, group = beta), size = 1)  +
  ylab("Number of infected individuals")+
  theme(legend.position="right") + 
  scale_color_viridis(discrete = FALSE)

plot
# Plot matrix commuting:
bet_n <-  as.character(bet[1:5]/max(bet))
heatmap(commut_mat,Rowv=NA, Colv =NA, RowSideColors = bet_n)
mat_max_inf <-  as.data.frame(mat_max_inf)
colnames(mat_max_inf) <-  c("Beta","time")
ggplot(mat_max_inf) + 
  geom_line(aes(Beta, time))

ggarrange(plot_inf_1_low, plot_inf_1_HALF, plot_inf_1_high)
ggarrange(plot_inf_1_high, plot_inf_1_lim)
#---------------------------------------------------------------------

plot_inf_2_lim <- plot_inf_2_lim + labs(title="a")
plot_inf_1 <- plot_inf_1 + labs(title="a")
plot_inf_2 <- plot_inf_2 + labs(title="c")
plot_eigen_1 <- plot_eigen_1 + labs(title="e")
plot_eigen_2 <- plot_eigen_2 + labs(title="e")


#-------------------SAVE FILE-------------------
Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/"
Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/"
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
plot_1
dev.off()

