rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")
library("ggpubr")
library("ggsci")
library("ggforce")
library("ggplot2")
library("gifski")

# Set up plots theme:
ggplot2::theme_set(ggplot2::theme_bw() %+replace% 
                     ggplot2::theme(axis.ticks =
                                      element_line(color = 'black'),
                                    # axis.title = element_text(color = 'black', size = 15),
                                    # axis.text = element_text(color = 'black', size = 15),
                                    # legend.text = element_text(color = 'black', size = 15),
                                    # legend.title = element_text(color = 'black', size = 16),
                                    # plot.title = element_text(color = 'black', size = 18, face = 'bold'),
                                    # strip.text = element_text(color = 'black', size = 15),
                                    legend.position = "none"
                     ))


#----------------------------------------------------------------------------#
source("~/RMT/Integration/functions_eigen_int.R")
#----------------PARAMETERS-----------------
#-------------------LOGIC----------------------
# Mobility parameter, 0: Just commuting, 1: just migration 2: migration & commuting.
MOB <- 2
# Integration parameter, 0: No integration, 1: integration.
INT <- 1
# Parameter for initial population. 0: No cte, 1: cte.
CTE_POP <- 0
# Parameter for transmission rate. 0: No cte, 1: cte.
BETA_CTE <- 0
# Parameter for initial infected ind. 0: No cte , 1: cte.
CTE_INF <- 1


#-------------------EPIDEMIOLOGICAL----------------------
N = 100 # Number of patches
# CTE parameters:
del_N <- rep(0.6, N) # Birth rate
bet_cte <- 6
bet_cte <-  0.001
bet <- rep(bet_cte, N)  # Transmission rate
# bet <- abs(rnorm(N,2,5))  # Transmission rate
d_vec <- rep(0.8, N) # Natural mortality rate
thet <- rep(0.1, N) # Rate of loss of immunity
alp <- rep(0.02, N) # Rate of disease overcome
delt <- rep(0, N) # Diseases related mortality rate

# Changing transmission rate:
mu = 0.5
sig = 4
ind <- sample(1:N,1)
bet_new <- 1
# bet_new <-  0.001
bet[ind] <- bet_new
# size <- sample(1:N,1)
# vec_rand <- sample(1:N,size)
# vec_rand <- seq(1,round(N/4),1)
# bet[vec_rand] <- bet_new
# 
print(paste0("gamma:", alp[1] + delt[1] + d_vec[1]))
print(paste0("beta - gamma:", bet[1] - (alp[1] + delt[1] + d_vec[1])))

#-------------------- MOBILITY ------------------#
### Migration:
alp_m <- 0.1
bet_m <- 0.1

# Compute mean and sd:
mu_m <- alp_m/(alp_m + bet_m)
s_m <-  sqrt((alp_m*bet_m)/(((alp_m + bet_m)^2)*(1+alp_m+bet_m)))
print(paste0("mu :", mu_m))
print(paste0("sigma :", s_m))

migrate_mat <- mat_conect(N,alp_m,bet_m,MOB)
### Commuting
alp_c <- 0.1
bet_c <- 0.1

# Compute mean and sd:
mu_w <- alp_c/(alp_c + bet_c)
s_w <- sqrt((alp_c*bet_c)/((alp_c + bet_c)^2*(1+alp_c+bet_c)))
print(paste0("mu :", mu_w))
print(paste0("sigma :",s_w))

commut_mat <- mat_conect(N,alp_c,bet_c,MOB)

tau_ct <- 0
# Initial populations:
init_pop <- matrix(10000, nrow = N)
if(CTE_POP == 1){
  list_param <- diff_f(migrate_mat,d_vec[1,1],init_pop)
  d_vec[1,1] <- list_param[1]
  del_N <- list_param[2:(N+1)]
  d_vec[1,1] = 1.7
}
# init_pop <- matrix(rgamma(N,shape = (mu_p/s_p)^2,rate = mu_p/(s_p^2)), nrow = N)


print(paste0("beta - gamma:", bet[1] - (alp[1] + delt[1] + d_vec[1])))
# End time integration:
end_time <- 50

#-------------------------------------------------------------------------------#
beta_ct = bet_cte  # gamma
gamma_ct = alp[1] + delt[1] + d_vec[1]        # beta
gamma_ct_w <- format(round(gamma_ct,2), decimal.mark = ',')
beta_ct_w <- format(round(beta_ct,2), decimal.mark = ',')
mu_w_w <- format(round(mu_w,2), decimal.mark = ',')
s_w_w <- format(round(s_w,2), decimal.mark = ',')
mu_m_w <- format(round(mu_m,2), decimal.mark = ',')
s_m_w <- format(round(s_m,2), decimal.mark = ',')
# Integrate the system:
count = 1
for(i in seq(0,5,0.05)){
  
  print(paste0("New beta : ",i))
  bet_new <- i
  bet[ind] <- bet_cte + bet_new
  sol <- int(N, del_N,bet,d_vec,thet,alp,delt,
             commut_mat,migrate_mat,end_time,
             MOB, CTE_POP, CTE_INF)
  
  sol_df <- as.data.frame(sol)
  
  # Plot the susceptible, infected and recovered:
  state <- "INF"
  # Parameters:
  print(paste0("beta - gamma", beta_ct - gamma_ct))
  
  cond_gen(N, mu_m, s_m, mu_w, s_w, gamma_ct, beta_ct, tau_ct)
  # Make distribution: 
  jac <- jacobian(N,bet,gamma_ct, commut_mat, migrate_mat,mu_m, MOB)
  eig <- eigen_mat(jac)
  
  # Predicted distribution:
  rad <- pred_radius(N, beta_ct, gamma_ct, tau_ct, mu_m, s_m, mu_w, s_w, MOB)
  cent <- pred_center(N, beta_ct, gamma_ct, tau_ct, mu_m, s_m, mu_w, s_w, MOB)
  outl <- pred_outlier(N, beta_ct, gamma_ct, tau_ct, mu_m, s_m, mu_w, s_w, MOB)
  
  a <- bet_cte*mu_w +mu_m
  b <- bet_new
  c <- bet_new*mu_w
  
  outl <- (1/2)*(N*a + b + sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
  outl2 <- (1/2)*(N*a + b - sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
  outl <- outl + (bet_cte*(1-mu_w) - N*mu_m - gamma_ct)
  outl2 <- outl2 + (bet_cte*(1-mu_w) - N*mu_m - gamma_ct)
  
  plot_inf_1 <- plot_int(N, sol, state)
  # If last parameter is 1 he outlier is jnoatj acjoajmputed:
  plot_eigen_1 <- plot_eigen(eig, cent, rad, outl, 1) +
    geom_point(aes(outl2,0), colour = "blue",
               show.legend = NA) + 
    geom_point(aes(outl,0), colour = "blue",
               show.legend = NA)
  
  
  plot_eigen_1
  
  vec_col <-  vector(mode="character", length=N)
  vec_col[1:N] <- "red4"
  # vec_col[vec_rand] <- "royalblue3"
  vec_col[ind] <- "royalblue3"
  
  plot_inf_1 <-  plot_inf_1 +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + 
    scale_colour_manual(values = vec_col)
  
  plot_inf_1_lim <- plot_inf_1 + xlim(c(0,30)) + ggtitle(paste(expression(alpha),":",i))
  
  count = count + 1
  
  print(paste0("N*(bet_cte*mu_w +mu_m): ", N*(bet_cte*mu_w +mu_m)))
  Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/1patch/High_both/"
  # Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/1patch/High_com/"
  path <- paste0(Path,"gen","N",N,"g",gamma_ct_w,"b",beta_ct_w,"mw",
                 mu_w_w,"sw",s_w_w,"mm",mu_m_w,"sm",s_m_w,"_",count,".png")
  # path <- paste0(Path,"mod1patchtrans","N",N,"g",gamma_ct,"b",
  # beta_ct,"mw","bnew","14",
  # mu_w,"sw",s_w,"mm",mu_m,"sm",s_m,".png")
  ggsave(file = path, plot = plot_inf_1_lim,width = 5, height = 5)
  # plot_inf_1_lim
  # dev.off()
}

# png_files <- list.files("~/Documents/PHD/2022/RMT_SIR/Plots/1patch/High_com/", pattern = ".*png$", full.names = TRUE)
# gifski(png_files, gif_file = "~/Documents/PHD/2022/RMT_SIR/Plots/1patch/High_com/animation.gif", width = 800, height = 600, delay = 0.3)

png_files <- list.files("~/Documentos/PHD/2022/RMT_SIR/Plots/1patch/High_both/", pattern = ".*png$", full.names = TRUE)

len <- length(png_files)
for(i in c(1:len)){
  png_files[i] <- paste0("/home/marta/Documentos/PHD/2022/RMT_SIR/Plots/1patch/High_both///genN100g0,82b0mw0,5sw0,46mm0,5sm0,46_",i+1,".png")
}
gifski(png_files, gif_file = "~/Documentos/PHD/2022/RMT_SIR/Plots/1patch/High_both/animation.gif", width = 800, height = 600, delay = 0.3)


ggarrange(plot_inf_1_low, plot_inf_1_HALF, plot_inf_1_high)
ggarrange(plot_inf_1_high, plot_inf_1_lim)
#--------------------------------------------------------------------------#
plot_inf_1_lim <- plot_inf_1_lim + labs(title="a")
plot_1 <- ggarrange(plot_inf_1_lim, plot_inf_2_lim)
# plot_1 <- ggarrange(plot_inf_1,plot_inf_2, labels = c("a","b","c","d"))
# plot_2 <- ggarrange(plot_eigen_1,plot_eigen_2, ncol = 2, labels = c("a","b","c","d"))
plot_1 <- ggarrange(plot_inf_1_lim,
                    plot_inf_2,
                    plot_eigen_1,
                    plot_eigen_2,
                    nrow = 2,ncol = 2,
                    label.x = 0,
                    label.y = 1)

plot_2 <- ggarrange(plot_inf_2,
                    plot_inf_2_lim,
                    nrow = 1,ncol = 2)

plot_bet <- ggarrange(plot_2,plot_eigen_1,
                      nrow = 2,ncol = 1)
plot_bet

plot_full <- ggarrange(plot,
                       plot_bet,
                       nrow = 1,
                       ncol = 2)
plot_inf_2_lim <- plot_inf_2_lim + labs(title="a")
plot_inf_1 <- plot_inf_1 + labs(title="a")
plot_inf_2 <- plot_inf_2 + labs(title="c")
plot_eigen_1 <- plot_eigen_1 + labs(title="e")
plot_eigen_2 <- plot_eigen_2 + labs(title="e")

plot <- ggarrange(plot_inf_1, plot_inf_2,
                  plot_eigen_1, plot_eigen_2,
                  nrow = 2,ncol = 2,
                  label.x = 0,
                  label.y = 1)
plot

plot <- ggarrange(plot_inf_2_lim,
                  plot_eigen_2,
                  nrow = 2,ncol = 1,
                  label.x = 0,
                  label.y = 1)
plot

#-------------------SAVE FILE-------------------
Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/"
Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/"
gamma_ct_w <- format(round(gamma_ct,2), decimal.mark = ',')
beta_ct_w <- format(round(beta_ct,2), decimal.mark = ',')
mu_w_w <- format(round(mu_w,2), decimal.mark = ',')
s_w_w <- format(round(s_w,2), decimal.mark = ',')
mu_m_w <- format(round(mu_m,2), decimal.mark = ',')
s_m_w <- format(round(s_m,2), decimal.mark = ',')
path <- paste0(Path,"gen","N",N_w,"g",gamma_ct_w,"b",beta_ct_w,"mw",
               mu_w_w,"b_new_1p", i,"sw",s_w_w,"mm",mu_m_w,"sm",s_m_w,".png")
# path <- paste0(Path,"mod1patchtrans","N",N,"g",gamma_ct,"b",
# beta_ct,"mw","bnew","14",
# mu_w,"sw",s_w,"mm",mu_m,"sm",s_m,".png")
png(file = path, width = 8000, height = 6000, res = 1100)
plot_1
dev.off()


