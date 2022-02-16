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
CTE_POP <- 1
# Parameter for transmission rate. 0: No cte, 1: cte.
BETA_CTE <- 0
# Parameter for initial infected ind. 0: No cte , 1: cte.
CTE_INF <- 1
# Parameter for distribution, "normal", "beta", "gamma":
DIST <-  "beta"

# Number of initial individuals by compartments:
SUS_INIT <- 100000 #Susceptible
INF_INIT <- 100    #Infected
#-------------------EPIDEMIOLOGICAL----------------------
N = 50 # Number of patches
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
# ind <- sample(1:N,1)
bet_new <- 6
# bet_new <-  0.001
# bet[ind] <- bet_new
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

migrate_mat <- mat_conect(N,alp_m,bet_m,DIST)
### Commuting
alp_c <- 0.1
bet_c <- 0.1

# Compute mean and sd:
mu_w <- alp_c/(alp_c + bet_c)
s_w <- sqrt((alp_c*bet_c)/((alp_c + bet_c)^2*(1+alp_c+bet_c)))
print(paste0("mu :", mu_w))
print(paste0("sigma :",s_w))

commut_mat <- mat_conect(N,alp_c,bet_c,DIST)

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
end_time <- 30

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
max_n_patch <-  6
k_vec <- seq(0,max_n_patch,1)
l <-  length(k_vec) + 1
mat_max_inf <-  matrix(0, ncol = 2, nrow = l)
plot_list <- list()
bet_new <- 1
for(i in k_vec){
  print(paste0("New beta : ",i))
  bet[1:i] <- bet_cte + bet_new
  sol <- int(N, del_N,bet,d_vec,thet,alp,delt,
             commut_mat,migrate_mat,end_time,
             MOB, CTE_POP, CTE_INF,SUS_INIT, INF_INIT,init_pop)
  
  sol_df <- as.data.frame(sol)
  
  # inf_df <- sol_df[, c(1,(N+2):(2*N+1))]
  # inf_df$sum <- rowSums(inf_df[,2:(N+1)])
  # t_max_inf <- inf_df$time[inf_df$sum==max(inf_df$sum)]
  # mat_max_inf[count,] <- c(i, t_max_inf)
  # Plot the susceptible, infected and recovered:
  state <- "INF"
  # Predicted distribution:
  plot_inf_1 <- plot_int(N, sol, state)
  
  vec_col <-  vector(mode="character", length=N)
  vec_col[1:N] <- "red4"
  vec_col[1:i] <- "royalblue3"
  
  plot_inf_1 <-  plot_inf_1 +
    scale_colour_manual(values = vec_col)
  
  plot_inf_1_lim <- plot_inf_1 +
    xlim(c(0,30)) +
    # scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
    ggtitle(paste("Numer of modified patches:",i)) 

  plot_list[[count]] <-  plot_inf_1_lim
  count = count + 1
  
  print(paste0("N*(bet_cte*mu_w +mu_m): ", N*(bet_cte*mu_w +mu_m)))
 
}

png_files <-  c()
for(i in c(1:N)){
  # png_files[i] <- paste0("/home/marta/Documentos/PHD/2022/RMT_SIR/Plots/1patch/High_both///genN100g0,82b0mw0,5sw0,46mm0,5sm0,46_",i+1,".png")
  plot_list[[i]] +
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE), limits = c(0,450000))
  ggsave(paste0("~/Documents/PHD/2022/RMT_SIR/Plots/Kpaches/low_both/plot_",i,".png" ),
         plot = last_plot(), device = "png")
  png_files[i] <-  paste0("~/Documents/PHD/2022/RMT_SIR/Plots/Kpaches/High_both/plot_",i,".png" )
}

# gifski(png_files, gif_file = "~/Documentos/PHD/2022/RMT_SIR/Plots/1patch/High_both/animation.gif", width = 800, height = 600, delay = 0.3)
gifski(png_files, gif_file = "~/Documents/PHD/2022/RMT_SIR/Plots/Kpaches/High_both/animation.gif", width = 800, height = 600, delay = 0.3)

mat_max_inf <-  as.data.frame(mat_max_inf)
colnames(mat_max_inf) <-  c("Beta","time")
ggplot(mat_max_inf) + 
  geom_line(aes(Beta, time))

ggarrange(plot_inf_1_low, plot_inf_1_HALF, plot_inf_1_high)
ggarrange(plot_inf_1_high, plot_inf_1_lim)

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

#--------------------------------------------------------------------#
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


