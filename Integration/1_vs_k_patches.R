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
  del_N <- rep(0.7, N) # Birth rate

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
  SUS_INIT <- 10000 #Susceptible
  INF_INIT <- 100    #Infected
  
  # End time integration:
  end_time <- 100

#----------Change 1 patch or change K patches with alpha/k----------
#--------------------------1 PATCH----------------------------------
  thet <- rep(0.6, N) # Rate of loss of immunity
  # Exit Rates:
  d_vec <- rep(0.7, N) # Natural mortality rate
  alp.vec <- rep(0.72, N) # Rate of disease overcome
  delt <- rep(0, N) # Diseases related mortality rate
  gamma_ct <-  alp.vec[1] + delt[1] + d_vec[1]
  print(paste0("gamma_ct:", gamma_ct))
  
  #Transmission rate:
  bet_cte <-  1
  bet <- rep(bet_cte, N)  # Transmission rate
  bet_new <- 5
  ind <-  sample(1:N,1)
  K <- 48 # Number of patches to change transmission rate.
  bet[ind] <-( bet_cte + K*bet_new)
  
  sol <- int(N, del_N,bet,d_vec,thet,alp.vec,delt,
             commut_mat,migrate_mat,end_time,
             MOB, CTE_POP, CTE_INF,SUS_INIT, INF_INIT,init_pop)

  sol_df <-  as.data.frame(sol)
  
  
  # Plot the susceptible, infected and recovered:
  state <- "INF"
  plot_inf_1 <- plot_int(N, sol, state) +
    theme_bw() + xlim(c(0,20))
  plot_inf_1
  
  # Make distribution:
  jac <- jacobian(N,bet,gamma_ct, commut_mat, migrate_mat,mu_c, MOB)
  eig <- eigen_mat(jac)
  
  plot_eig <- ggplot(eig) + geom_point(aes(re,im), size = 0.05)
  plot_eig + coord_fixed()
  
  # Compute the outliers, center and radius:
  a <- bet_cte*mu_w + mu_c
  b <- K*bet_new
  c <- K*bet_new*mu_w
  beta_ct <-  bet_cte
  rad <- pred_radius(N, beta_ct, gamma_ct, tau_ct, mu_c, sqrt(s_c), mu_w, sqrt(s_w), MOB)
  cent <- pred_center(N, beta_ct, gamma_ct, tau_ct, mu_c, sqrt(s_c), mu_w, sqrt( s_w), MOB)
  
  outl <- (1/2)*(N*a + b + sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
  outl2 <- (1/2)*(N*a + b - sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
  outl <- outl + (bet_cte*(1-mu_w) - N*mu_c - gamma_ct)
  outl2 <- outl2 + (bet_cte*(1-mu_w) - N*mu_c - gamma_ct)
  
  
  # If last parameter is 1 he outlier is not computed:
  plot_eigen_1 <- plot_eigen(eig, cent, rad, outl, 1) +
    geom_point(aes(outl,0), colour =  "blue",
               show.legend = NA) +
    geom_point(aes(outl2,0), colour =  "blue",
               show.legend = NA)
  
  plot_eigen_1
  
  #--------------------------K PATCH-----------------------------------
  ind <- sample(1:N,K)
  bet.k <- rep(bet_cte, N)  # Transmission rate
  bet.k[ind] <- bet_cte + bet_new
  
  sol.k <- int(N, del_N,bet.k,d_vec,thet,alp.vec,delt,
             commut_mat,migrate_mat,end_time,
             MOB, CTE_POP, CTE_INF,SUS_INIT, INF_INIT,init_pop)

  sol.k.df <-  as.data.frame(sol.k)

  
  # Plot the infected :
  state <- "INF"
  plot_inf_k <- plot_int(N, sol.k, state) +
    theme_bw() + xlim(c(0,20))
  plot_inf_k
  
  # Compute the jacobian matrix:
  jac.k <- jacobian(N,bet.k,gamma_ct, commut_mat, migrate_mat,mu_c, MOB)
  eig.k <- eigen_mat(jac.k)
  
  plot_eig_k <- ggplot(eig.k) + geom_point(aes(re,im), size = 0.05)
  plot_eig_k + coord_fixed()
  
  # Compute outliers:
  a <- bet_cte*mu_w + mu_c
  b <- bet_new
  c <- bet_new*mu_w
  
  outl <- (1/2)*(N*a + b + (K-1)*c + sqrt((N*a)^2 - (2*N-4*K)*a*b +
                                            (2*(K+1)*N-4*K)*a*c + b^2 + (2*K-2)*b*c + (K-1)^2*c^2))
  outl2 <- (1/2)*(N*a + b + (K-1)*c - sqrt((N*a)^2 - (2*N-4*K)*a*b +
                                             (2*(K+1)*N-4*K)*a*c + b^2 + (2*K-2)*b*c + (K-1)^2*c^2))
  outl3 <- b - c
  
  outl <- outl + (bet_cte*(1-mu_w) - N*mu_c - gamma_ct)
  outl2 <- outl2 + (bet_cte*(1-mu_w) - N*mu_c - gamma_ct)
  outl3 <- outl3 + (bet_cte*(1-mu_w) - N*mu_c - gamma_ct)
  
  rad.k <- pred_radius(N, bet_cte, gamma_ct[1], tau_ct, mu_c, sqrt(s_c), mu_w, sqrt(s_w), MOB)
  cent.k <- pred_center(N, bet_cte, gamma_ct[1], tau_ct, mu_c, sqrt(s_c), mu_w,sqrt(s_w), MOB)
  
  rad.k <- pred_radius(N, bet_new, gamma_ct[1], tau_ct, mu_c, sqrt(s_c), mu_w, sqrt(s_w), MOB)
  cent.k <- pred_center(N, bet_new, gamma_ct[1], tau_ct, mu_c, sqrt(s_c), mu_w,sqrt(s_w), MOB)
  
  
  plot_eigen_k <- plot_eigen(eig.k, cent.k, rad.k, outl, 1) +
    geom_point(aes(outl,0), colour =  "blue",
               show.legend = NA) +
    geom_point(aes(outl2,0), colour =  "blue",
               show.legend = NA) +
    geom_point(aes(outl3,0), colour =  "blue",
               show.legend = NA)
  
  plot_eigen_k
#--------------------------------------------------------------------#
  ggarrange(plot_eigen_1, plot_eigen_k,
            plot_inf_1, plot_inf_k,
            ncol = 2, nrow= 2)
#--------------------------------------------------------------------#
  Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/"
  Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/"
  gamma_ct_w <- format(round(gamma_ct,2), decimal.mark = ',')
  beta_ct_w <- format(round(beta_ct,2), decimal.mark = ',')
  mu_w_w <- format(round(mu_w,2), decimal.mark = ',')
  s_w_w <- format(round(s_w,2), decimal.mark = ',')
  mu_c_w <- format(round(mu_c,2), decimal.mark = ',')
  s_c_w <- format(round(s_c,2), decimal.mark = ',')
  path <- paste0(Path,"gen","N",N,"g",gamma_ct_w,"b",beta_ct_w,"mw",
                 mu_w_w,"b_new_1p",K*bet_new,"b_new_kp", "n_patch",K,
                 bet_new,"sw",s_w_w,"mm",mu_c_w,"sm",s_c_w,".png")
  # path <- paste0(Path,"mod1patchtrans","N",N,"g",gamma_ct,"b",
  # beta_ct,"mw","bnew","14",
  # mu_w,"sw",s_w,"mm",mu_c,"sm",s_c,".png")
  png(file = path, width = 8000, height = 6000, res = 1100)
  ggarrange(plot_eigen_1, plot_eigen_k,
            plot_inf_1, plot_inf_k,
            ncol = 2, nrow= 2)
  dev.off()
  
#--------------------------------------------------------------------#
  # Check that the outlier for 1 patch with K*beta^* is the same as 
  # K patches with transmission rate beta^*
  # Varying K:
  # Compute outliers:
  N = 100
  bet_cte <-  3
  dim = 1000
  mat.pred <- matrix(0, ncol = 3, nrow = dim) 
  for(i in c(1:dim)){
    print(paste0("i: ",i))
    print(paste0("alp: ",alp))
    alp <- abs(rnorm(1,10,9))
    out.pred.1 <- max(outl_1patch(alp, bet_cte, N,
                                  mu_w, mu_c, gamma_ct))
    out.pred.k <- outl_Kpatch(i, alp/i, bet_cte,
                              N, mu_w, mu_c, gamma_ct)
    mat.pred[i,] <-  c(i,out.pred.1, out.pred.k)
  }
  
  df.Kpatch <-  data.frame(num_patches = mat.pred[,1],
                           patch1 = mat.pred[,2],  patchk = mat.pred[,3])
  head(df.Kpatch)
  
  df.plot <- reshape2::melt(df.Kpatch, id.vars = c("num_patches"))
  
  plot  <- ggplot(df.plot,aes(num_patches, value)) + 
    geom_line(aes( colour = variable))  +
    ylab("Real part") + xlab("Number of patches") 
  
  plot_k.alp <- plot + theme_bw() +
    theme(legend.position="none")
  plot_k.alp 
  
  plot_k_alp <- plot + theme_bw() +
    theme(legend.position="none")
  plot_k_alp
  # ------------------COMPARISON RIGHT MOST----------------------------
  d_vec <- rep(0.7, N) # Natural mortality rate
  alp.vec <- rep(0.72, N) # Rate of disease overcome
  delt <- rep(0, N) # Diseases related mortality rate
  # Check if the right most eigenvalue is the same for many iterations
  thet <- rep(0.6, N) # Rate of loss of immunity
  gamma_ct <-  alp.vec[1] + delt[1] + d_vec[1]
  
d <- 100
alphag <- alphagamma(2,6)
betag <- betagamma(2,6)
bet_vec <- rgamma(d,alphag,betag) 
alp_vec <- rgamma(d,alphag,betag) 
for(i in c(1:d)){
  ind <-  sample(1:N,1)
  bet_cte <- bet_vec[i]
  bet_new <- alp_vec[i]
  alp <- bet_new
  mat.comp <-  matrix(0, ncol = 8, nrow = N)
  df.comp <- data.frame(N.patches=0,out1=0,outk=0,diff=0,
                          out.pred.1=0, out.pred.k=0, max.inf.1=0,max.inf.k=0)
  for(i in c(1:N)){
    # Compute  right most eigenvalue for 1 patch:
    bet <-  rep(bet_cte,N)
    bet[ind] <- bet[ind] + alp
    # Compute the jacobian matrix:
    jac <- jacobian(N, bet, gamma_ct, commut_mat, migrate_mat, mu_c, MOB)
    eig <- eigen_mat(jac)
    out1 <-  max(eig$re)
    out.pred.1 <- max(outl_1patch(alp, bet_cte, N, mu_w, mu_c, gamma_ct))
    
    sol <- int(N, del_N,bet,d_vec,thet,alp.vec,delt,
               commut_mat,migrate_mat,5,
               MOB, CTE_POP, CTE_INF,SUS_INIT, INF_INIT,init_pop)
    
    # Compute  right most eigenvalue for K patches:
    bet <-  rep(bet_cte,N)
    bet[1:i] <- bet[ind] + alp/i
    # Compute the jacobian matrix:
    jac <- jacobian(N,bet,gamma_ct, commut_mat, migrate_mat,mu_c, MOB)
    eig <- eigen_mat(jac)
    outk <-  max(eig$re)
    out.pred.k <- outl_Kpatch(i, alp/i, bet_cte, N, mu_w, mu_c, gamma_ct)
    
    sol.k <- int(N, del_N,bet,d_vec,thet,alp.vec,delt,
               commut_mat,migrate_mat,5,
               MOB, CTE_POP, CTE_INF,SUS_INIT, INF_INIT,init_pop)
    
    print(paste0("i:",i))
    mat.comp[i,] = c(i, out1, outk,
                     abs((outk -out1)/outk),out.pred.1,
                     out.pred.k, max(sol[,c((N+2):(2*N+1))]),
                     max(sol.k[,c((N+2):(2*N+1))]))
  }
  
df.comp1 <-  as.data.frame(mat.comp)
colnames(df.comp1) <-  c("N.patches","out1","outk","diff",
                        "out.pred.1", "out.pred.k", "max.inf.1","max.inf.k")
df.comp <- rbind(df.comp, df.comp1)
}

write.csv(df.comp,"~/RMT/Integration/1_vs_k_patch.csv", row.names = TRUE)
df.comp$diff.1 <- abs(df.comp$out1 - df.comp$out.pred.1)/df.comp$out.pred.1
df.comp$diff.k <- abs(df.comp$outk - df.comp$out.pred.k)/df.comp$out.pred.k
df.comp$diff.pred <- abs(df.comp$out.pred.1 - df.comp$out.pred.k)/df.comp$out.pred.k
df.comp$K.alp <- df.comp$N.patches * alp
df.comp$diff.inf.abs <- df.comp$max.inf.1 - df.comp$max.inf.k
df.comp$diff.inf <- abs(df.comp$diff.inf.abs)/max(df.comp$max.inf.1,df.comp$max.inf.k)

df.comp.round <- round(df.comp,2)
max(df.comp$diff)

df.comp.filt <- df.comp[,c(1,2,3,5,6)]
colnames(df.comp.filt) <- c("N.patches","Real outlier 1 patch", "Real outlier k patches",
                            "Predicted outlier 1 patch",
                            "Predicted outlier k patches")
df.plot <- reshape2::melt(df.comp.filt, id.vars = c("N.patches"))
plot.pred.real  <- ggplot(df.plot,aes(N.patches, value)) + 
  geom_line(aes( colour = variable))  +
  ylab("Real part") + xlab("Number of patches") + 
  scale_fill_continuous(name = FALSE,
                      labels = ) + theme_bw() 
plot.pred.real

ggarr1 <- ggarrange(plot_k.alp, plot_k_alp)
ggarrange(ggarr1,plot, ncol=1, nrow =2)

gamma_ct_w <- format(round(gamma_ct,2), decimal.mark = ',')
beta_ct_w <- format(round(beta_ct,2), decimal.mark = ',')
mu_w_w <- format(round(mu_w,2), decimal.mark = ',')
s_w_w <- format(round(s_w,2), decimal.mark = ',')
mu_c_w <- format(round(mu_c,2), decimal.mark = ',')
s_c_w <- format(round(s_c,2), decimal.mark = ',')


Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/epi_param/"
Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/epi_param/"

path <- paste0(Path,"pred_real_k_1","N",N,"g",gamma_ct_w,"b",beta_ct_w,
               "bet_mew", bet_new,
               "mw",mu_w_w,"sw",s_w_w,"mm",mu_c_w,"sm",s_c_w,".png")
gg_full <- ggarrange(ggarr1,plot, ncol=1, nrow =2)
ggsave(path,
       plot =gg_full, device = "png")
