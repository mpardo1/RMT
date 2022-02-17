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
bet_cte <-  1
bet <- rep(bet_cte, N)  # Transmission rate
# bet <- abs(rnorm(N,1,1))  # Transmission rate
d_vec <- rep(1, N) # Natural mortality rate
thet <- rep(0.6, N) # Rate of loss of immunity
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
end_time <- 100
#----------------------AREA of STAB muw vs alp-------------------------
  # Influence of alpha in outlier:
  # alp_vet <- 3
  gamma_ct <- 5
  bet_vec <- seq(0.01,5,0.01)
  len_vec <- length(bet_vec)
  plot_list_area <- list()
  ind <- sample(1:N,1)
  # gamma_ct = alp[1] + delt[1] + d_vec[1] # gamma   
  for(i in c(1:len_vec)){
    # beta_ct = bet_cte  # beta
    # bet[ind] <- bet_cte + alp_vet
    max_alp <- 5
    vec_alp <- seq(0,max_alp,0.01)
    mu_w_vec <- seq(0.01,1.01,0.01)
    len <- length(mu_w_vec)
    df_mu_w_slp <- data.frame(alp = 0, mu_w = 0, outl1 = 0, outl2 = 0)
    for(j in c(1:len)){
      vec <- sapply(vec_alp, outl_1patch, bet_vec[i] , N, mu_w_vec[j], mu_c, gamma_ct)
      df <- data.frame(alp =vec_alp, mu_w =mu_w_vec[j], outl1 = vec[1,], outl2 = vec[2,])
      df_mu_w_slp <- rbind(df_mu_w_slp,df)
    }
    
    df_mu_w_slp <- df_mu_w_slp[-1,]
    df_mu_w_slp$stability <- ifelse((df_mu_w_slp$outl1 > 0 | df_mu_w_slp$outl2 > 0) , FALSE, TRUE)
    plot_list_area[[i]] <- ggplot(df_mu_w_slp) + 
      geom_point(aes(alp, mu_w, colour = stability)) +
      xlab(expression(alpha)) + 
      ylab(~ paste( mu [w])) +
      theme_bw() +
      ggtitle(paste0("beta_cte*mu_w - gamma :",  bet_vec[i] - gamma_ct ))
  }
  
  st <- 1
  end <- 18
  plot_bet_0.01 <- plot_list_area[[st]] + ggtitle("") +
    labs(title="a")+ theme(legend.position = "none")
  plot_bet_0.18 <- plot_list_area[[end]]+ labs(title="b")+ theme(legend.position = "none")
  
  ggarra <- ggarrange(plot_bet_0.01, plot_bet_0.18, common.legend = FALSE)
  
  gg_full <-  ggarrange(ggarra, ggarra1, common.legend = TRUE, ncol = 1, nrow = 2)
  
  #--------------------------SAVE FILE--------------------------------#
  gamma_ct_w <- format(round(gamma_ct,2), decimal.mark = ',')
  beta_ct_1 <- format(round(bet_vec[st],2), decimal.mark = ',')
  beta_ct_2 <- format(round(bet_vec[end],2), decimal.mark = ',')
  mu_w_w <- format(round(mu_w,2), decimal.mark = ',')
  s_w_w <- format(round(s_w,2), decimal.mark = ',')
  mu_c_w <- format(round(mu_c,2), decimal.mark = ',')
  s_c_w <- format(round(s_c,2), decimal.mark = ',')
  
  
  Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/epi_param/Area/"
  Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/epi_param/Area/"
  
  path <- paste0(Path,"Area","N",N,"g1_2","g2_",gamma_ct_w,"b1",beta_ct_1,"b2",beta_ct_2,"mw",
                 mu_w_w,"sw",s_w_w,"mm",mu_c_w,"sm",s_c_w,".png")
  ggsave(path,
         plot =gg_full, device = "png")
  
  #---------------------------------------------------------------------#
  png_files <-  c()
  for(i in c(1:len_vec)){
    # png_files[i] <- paste0("/home/marta/Documentos/PHD/2022/RMT_SIR/Plots/1patch/High_both///genN100g0,82b0mw0,5sw0,46mm0,5sm0,46_",i+1,".png")
    
    ggsave(paste0("~/Documentos/PHD/2022/RMT_SIR/Plots/epi_param/Area/Mid_gam/plot_gam_",gamma_ct,"_",i,".png" ),
           plot = plot_list_area[[i]], device = "png")
    png_files[i] <-  paste0("~/Documentos/PHD/2022/RMT_SIR/Plots/epi_param/Area/Mid_gam/plot_gam_",gamma_ct,"_",i,".png" )
  }
  
  # png_files <- list.files("~/Documents/PHD/2022/RMT_SIR/Plots/1patch/test/", pattern = ".*png$", full.names = TRUE)
  # gifski(png_files, gif_file = "~/Documentos/PHD/2022/RMT_SIR/Plots/1patch/High_both/animation.gif", width = 800, height = 600, delay = 0.3)
  gifski(png_files, gif_file = "~/Documentos/PHD/2022/RMT_SIR/Plots/epi_param/Area/Mid_gam/animation_gam_",gamma_ct,".gif", width = 800, height = 600, delay = 0.3)
  
  # Compute the right most eigenvalue vs alpha:
  vec <- sapply(vec_alp, outl_1patch, bet_cte , N, mu_w, mu_c, gamma_ct)
  alp_df <- data.frame(alp_bet = vec_alp, outl <- vec[1,])
  
  plot_alp_max_eigen  <- ggplot(alp_df) + 
    geom_line(aes( alp_bet, outl))  +
    ylab("Rigth most eigenvalue") + xlab(expression(alpha)) +
    geom_segment(aes(x = 0, y = 0, xend = max_alp, yend =0,
                     colour = "segment"), linetype=2, 
                 show.legend = FALSE) +
    theme_bw() 
  
  plot_alp_max_eigen 
  #----------------------AREA of STAB muw vs N-----------------------------------
  # Influence of alpha in outlier:
  # alp_vet <- 3
  gamma_ct <- 10.4
  N_vec = seq(50,250,1)
  len_vec <- length(N_vec)
  plot_list_area <- list()
  ind <- sample(1:N,1)
  # gamma_ct = alp[1] + delt[1] + d_vec[1] # gamma   
  for(i in c(1:len_vec)){
    max_bet <- 5
    vec_bet <- seq(0,max_bet,0.01)
    mu_w_vec <- seq(0.01,1.01,0.01)
    len <- length(mu_w_vec)
    df_mu_w_slp <- data.frame(N = 0, mu_w = 0, outl = 0)
    for(j in c(1:len)){
      vec <- sapply(N_vec, out_gen, vec_bet[i], mu_w_vec[j], gamma_ct)
      df <- data.frame(N =N_vec, mu_w =mu_w_vec[j], outl= vec)
      df_mu_w_slp <- rbind(df_mu_w_slp,df)
    }
    
    df_mu_w_slp <- df_mu_w_slp[-1,]
    df_mu_w_slp$stability <- ifelse((df_mu_w_slp$outl > 0) , FALSE, TRUE)
    plot_list_area[[i]] <- ggplot(df_mu_w_slp) + 
      geom_point(aes(N, mu_w, colour = stability)) +
      xlab("Number of patches") + 
      ylab(~ paste( mu [w])) +
      theme_bw() +
      ggtitle(paste("gamma: ",gamma_ct))
  }
  
  st <- 1
  end <- 18
  plot_bet_0.01 <- plot_list_area[[st]] + ggtitle("") +
    labs(title="a")+ theme(legend.position = "none")
  plot_bet_0.18 <- plot_list_area[[end]]+ labs(title="b")+ theme(legend.position = "none")
  
  ggarra <- ggarrange(plot_bet_0.01, plot_bet_0.18, common.legend = FALSE)
  
  gg_full <-  ggarrange(ggarra, ggarra1, common.legend = TRUE, ncol = 1, nrow = 2)
  
  gg_arr <-  ggarrange(plot_gen_1,plot_gen_2,
                       plot_gen_3,plot_gen_4,
                       nrow = 2, ncol = 2, common.legend = TRUE)
  #--------------------------SAVE FILE--------------------------------#
  gamma_ct_w <- format(round(gamma_ct,2), decimal.mark = ',')
  beta_ct_1 <- format(round(bet_vec[st],2), decimal.mark = ',')
  beta_ct_2 <- format(round(bet_vec[end],2), decimal.mark = ',')
  mu_w_w <- format(round(mu_w,2), decimal.mark = ',')
  s_w_w <- format(round(s_w,2), decimal.mark = ',')
  mu_c_w <- format(round(mu_c,2), decimal.mark = ',')
  s_c_w <- format(round(s_c,2), decimal.mark = ',')
  
  Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/Gen/Area/"
  Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/Area/"
  
  path <- paste0(Path,"Area","N",N,"g1_2","g2_",gamma_ct_w,"b1",beta_ct_1,"b2",beta_ct_2,"mw",
                 mu_w_w,"sw",s_w_w,"mm",mu_c_w,"sm",s_c_w,".png")
  ggsave(path,
         plot =gg_full, device = "png")
  
  #---------------------------------------------------------------------#
  png_files <-  c()
  for(i in c(1:len_vec)){
    # png_files[i] <- paste0("/home/marta/Documentos/PHD/2022/RMT_SIR/Plots/1patch/High_both///genN100g0,82b0mw0,5sw0,46mm0,5sm0,46_",i+1,".png")
    
    ggsave(paste0("~/Documents/PHD/2022/RMT_SIR/Plots/Gen/Area/High_gam/plot_gam_",gamma_ct,"_",i,".png" ),
           plot = plot_list_area[[i]], device = "png")
    png_files[i] <-  paste0("~/Documents/PHD/2022/RMT_SIR/Plots/Gen/Area/High_gam/plot_gam_",gamma_ct,"_",i,".png" )
  }
  
  # png_files <- list.files("~/Documents/PHD/2022/RMT_SIR/Plots/1patch/test/", pattern = ".*png$", full.names = TRUE)
  # gifski(png_files, gif_file = "~/Documentos/PHD/2022/RMT_SIR/Plots/1patch/High_both/animation.gif", width = 800, height = 600, delay = 0.3)
  gifski(png_files, gif_file = "~/Documents/PHD/2022/RMT_SIR/Plots/Gen/Area/High_gam/animation_gam_",gamma_ct,".gif", width = 800, height = 600, delay = 0.3)
  
  # Compute the right most eigenvalue vs alpha:
  vec <- sapply(vec_alp, outl_1patch, bet_cte , N, mu_w, mu_c, gamma_ct)
  alp_df <- data.frame(alp_bet = vec_alp, outl <- vec[1,])
  
  plot_alp_max_eigen  <- ggplot(alp_df) + 
    geom_line(aes( alp_bet, outl))  +
    ylab("Rigth most eigenvalue") + xlab(expression(alpha)) +
    geom_segment(aes(x = 0, y = 0, xend = max_alp, yend =0,
                     colour = "segment"), linetype=2, 
                 show.legend = FALSE) +
    theme_bw() 
  
  plot_alp_max_eigen 
  