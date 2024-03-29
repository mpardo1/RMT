rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")
library("ggpubr")
library("ggsci")
library("ggforce")
library("dplyr")

# -------------------------PARAMETERS ---------------------------------
N = 200 # Number of patches
mu = 1 
sig = 4
# Random generation of parameters:
# del_N <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2)) # Birth rate
# bet <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2))   # Transmission rate
# d_vec <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2)) # Natural mortality rate
# thet <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2))  # Rate of loss of immunity 
# alp <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2)) # Rate of disease overcome
# delt <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2)) # Diseases related mortality rate

# CTE parameters:
del_N <- matrix(0.6, ncol = N, nrow = 1) # Birth rate
bet <- matrix(0.3, ncol = N, nrow = 1)  # Transmission rate
d_vec <- matrix(1.3, ncol = N, nrow = 1) # Natural mortality rate
thet <- matrix(0.1, ncol = N, nrow = 1) # Rate of loss of immunity
alp <- matrix(0.31, ncol = N, nrow = 1) # Rate of disease overcome
delt <- matrix(0, ncol = N, nrow = 1) # Diseases related mortality rate

print(paste0("gamma:", alp[1,1] + delt[1,1] + d_vec[1,1]))
print(paste0("beta - gamma:", bet[1,1] - (alp[1,1] + delt[1,1] + d_vec[1,1])))

# Integration parameter, 0: No integration, 1: integration.
INT <- 1
#--------------------MOBILITY PARAMETERS --------------------------
# Mobility parameter, 0: Just commuting, 1: just migration 2: migration & commuting.
MOB <- 2

#-----------MIGRATION--------------
# Beta distribution:
alp_m <- 2
bet_m <- 0.5
mu_m <- alp_m/(alp_m + bet_m)
s_m <-  sqrt((alp_m*bet_m)/(((alp_m + bet_m)^2)*(1+alp_m+bet_m)))
print(paste0("mu :", mu_m))
print(paste0("sigma :", s_m))
  
if( MOB == 1 | MOB ==2){
  print("Migration added")
  # connect_mat <- matrix(runif(N^2, min , max), nrow = N)
  connect_mat <- matrix(rbeta(N^2, alp_m, bet_m, ncp = 0), nrow = N)
}else{
  print("No migration")
  connect_mat <- matrix(0, nrow = N, ncol = N) # No migration
}

diag(connect_mat) <- 0
#-----------COMMUTING--------------
# Beta distribution:
alp_c <- 0.3
bet_c <- 0.9
mu_w <- alp_c/(alp_c + bet_c)
s_w <- sqrt((alp_c*bet_c)/((alp_c + bet_c)^2*(1+alp_c+bet_c)))
print(paste0("mu :", mu_w))
print(paste0("sigma :",s_w))

if( MOB == 0 | MOB ==2){
  print("Commuting added")
  # commut_mat <- matrix(runif(N^2, min, max, nrow = N)
  commut_mat <- matrix(rbeta(N^2, alp_c, bet_c, ncp = 0), nrow = N)
}else{
  print("No commuting")
  commut_mat <- matrix(0, nrow = N, ncol = N) # No commuting
}

diag(commut_mat) <- 0
# Parameter for initial population. 0: No cte pop, 1: cte pop.
CTE_POP <- 1
# Constant population at each patch (ie. no mortality induced by the disease):
if(CTE_POP == 0){
  print("No constant population")
}else{
  print("Constant population")
  # mu_p = 10000
  # s_p = 5000
  # Initial populations:
  init_pop <- matrix(10000, nrow = N)
  # init_pop <- matrix(rgamma(N,shape = (mu_p/s_p)^2,rate = mu_p/(s_p^2)), nrow = N)
  delt <- matrix(0, ncol = N, nrow = 1) # Diseases related mortality rate
  del_N <- c()
  for(i in c(1:N)){
    del_N[i] <- sum(connect_mat[,i]) - sum(connect_mat[i,]*init_pop)/init_pop[i] 
    if(del_N[i] < 0 && max(connect_mat[i,]) < 1){
      print("Error delta_N less than zero")
    }
  }
}


# Create vector of parameters for ode function:
parameters <- list(
  dim = N,
  delta_N = del_N,
  beta_r = bet,
  d = d_vec,
  theta_r = thet,
  C = connect_mat,
  W = commut_mat,
  alpha_r = alp,
  delta_r = delt )

# ----------------------------MODELS ----------------------------------
#--------------------ND: SIR metapopulation model-----------------------
SIR <- function(t, y, parameters) {
  with(as.list(c(y, parameters)),{
    dy <- c()
    dim1 <- dim + 1
    dim2 <- dim*2
    dim3 <- (dim*2)+1
    dim4 <- dim*3
    for(i in c(1:dim)){
      # Total population (N = S+I+R)
      N <- y[i] + y[i+dim] + y[i+dim2]
      
      # Susceptible individuals:
      q1 <- delta_N[i]*N
      q2 <- beta_r[i]*(y[i]/N)*y[i+dim]
      q3 <- d[i]*y[i]
      q4 <- theta_r[i]*y[i+dim2]
      q5 <- y[i]*sum(C[,i])
      q6 <- sum(y[1:dim]*C[i,]) 
      q7 <- (y[i]/N)*sum(beta_r*W[,i]*y[dim1:dim2])
      # dS/dt
      dy[i] <- q1 - q2 - q3 + q4 - q5 + q6 - q7
      
      # Infected individuals:
      q1 <- beta_r[i]*(y[i]/N)*y[i+dim]
      q2 <- (alpha_r[i] + delta_r[i] + d[i])*y[i+dim] 
      q3 <- y[dim+i]*sum(C[,i]) 
      q4 <- sum(y[dim1:dim2]*C[i,])
      q5 <- (y[i]/N)*sum(beta_r*W[,i]*y[dim1:dim2])
      # dI/dt
      dy[i+dim] <- q1 - q2 - q3 + q4 + q5
      
      # Recovered individuals:
      q1 <- alpha_r[i]*y[i+dim] 
      q2 <- (theta_r[i] + d[i])*y[i+2*dim] 
      q3 <- y[i+2*dim]*sum(C[,i])
      q4 <- sum(y[dim3:dim4]*C[i,])
      # dR/dt
      dy[i+dim2] <- q1 - q2 - q3 + q4
        
    }
    list(dy)
  }) 
}

end_time <- 30
times = seq(0,end_time, 0.1)

# Vector of initial values:
population <- c(matrix(0,ncol=3*N,nrow =1 ))
# Initial value for infected individuals:
population[(N+1):(2*N)] <- ceiling(abs(rnorm(N,100,60)))

if(CTE_POP == 0){
  print("No constant population")
  population[1:N] <- 100000
}else{
  print("Constant population")
  # Susceptible individuals:
  population[1:N] <- init_pop - population[(N+1):(2*N)]
}

if(INT == 1){
  print("Doing integration")
  # Run integration:
  z <- ode(population, times, SIR, parameters)
  
  head(z)
  
  # Change labels:
  for(i in c(1:N)){
    colnames(z)[i+1] <-  paste0("S",i)
    colnames(z)[N+i+1] <-  paste0("I",i)
    colnames(z)[2*N+i+1] <-  paste0("R",i)
  }
  
  z <- as.data.frame(z)
  # z  <- z   %>% filter( z$time < 1) 
  df_plot <- reshape2::melt(z, id.vars = c("time"))
  
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
  
  head(df_plot)
  # Plot number of individuals at each time.
  plot_tot  <- ggplot(df_plot,aes(time, value)) + 
    geom_line(aes( colour = variable))  +
    ylab("Number of individuals") 
  
  # plot_tot
  
  # Filter Susceptibles:
  df_sus <- df_plot  %>% filter( substr(df_plot$variable,1,1) == "S")
  plot_sus  <- ggplot(df_sus,aes(time, value)) + 
    geom_line(aes( colour = variable))  +
    ylab("Number of individuals")  + 
    ggtitle("Susceptible individuals")
  # plot_sus
  
  # Filter Infected:
  df_inf <- df_plot  %>% filter( substr(df_plot$variable,1,1) == "I")
  plot_inf  <- ggplot(df_inf,aes(time, value)) + 
    geom_line(aes( colour = variable))  +
    ylab("Number of infected individuals") 
  # +
  #   ggtitle("Infected individuals")
  
  # plot_inf
  
  # Filter Recovered:
  df_rec <- df_plot  %>% filter( substr(df_plot$variable,1,1) == "R")
  plot_rec  <- ggplot(df_rec,aes(time, value)) + 
    geom_line(aes( colour = variable))  +
    ylab("Number of individuals") +
    ggtitle("Recovered individuals")
  
  # plot_rec
  
  # ggarrange(plot_tot,plot_inf)
}

#---------------EIGENVALUE DISTRIBUTION------------------
#-------------------FUNCTIONS----------------------------
# Function which plot the eigenvalues:
eigen_mat <- function(mat){
  eigen_m <- as.complex(eigen(mat, only.values = TRUE)$values)
  df <- data.frame(re = Re(eigen_m), im = Im(eigen_m))
  return(df)
}

# Create a random matrix
rand_mat <- function(N,mu,sig,distrib){
  muln <- log(mu^2/sqrt(mu^2 + sig^2))
  sdln <- sqrt(log(1+sig^2/mu^2))
  rmatrix <- dplyr::case_when(
    distrib == "gamma" ~ matrix(rgamma(N^2,shape = (mu/sig)^2,rate = mu/(sig^2)), nrow = N),
    distrib == "lognormal" ~ matrix(rlnorm(N^2,meanlog = muln,sdlog = sdln), nrow = N),
    TRUE ~ matrix(rep(0, N^2), nrow = N)
  )
  rmatrix <- matrix(rmatrix, nrow = N)
  #return(rmatrix)
}

# Parameters:
beta_ct = bet[1,1]  # gamma
gamma_ct = alp[1,1] + delt[1,1] + d_vec[1,1]        # beta

#------------------- COMMUTING ----------------------
if(MOB == 0){
  print("Just commuting")
  # Generate de Jacobian:
  betas <- matrix(rep(0,N^2), nrow = N)
  diag(betas) <- rep(beta_ct,N)
  BIGT <- commut_mat
  diag(BIGT) <- rep(1,N)
  BIGT <- BIGT%*%betas
  BIGS <- matrix(rep(0,N^2), nrow = N)
  diag(BIGS) <- rep(-gamma_ct, N)
  jacobian <- BIGT+BIGS
  
  # Compute the eigenvalues:
  eig <- eigen_mat(jacobian)
  
  # Compute the center and radius for the circular law:
  center = beta_ct*(1-mu_w)-gamma_ct
  radius = beta_ct*s_w*sqrt(N)
  outlier <- beta_ct*(mu_w*(N-1)+1)-gamma_ct
}

#------------------- MIGRATION ----------------------
if(MOB == 1){
  print("Just migration")
  # Generate de Jacobian:
  BIGT <- connect_mat
  diag(BIGT) <- rep(beta_ct - gamma_ct - mu_m*(N-1),N)
  # diag(BIGT) <- 0
  # diag(BIGT) <- rep(beta_ct - gamma_ct,N) - colSums(BIGT)
  jacobian <- BIGT
  
  # Compute the eigenvalues:
  eig <- eigen_mat(jacobian)
  
  # Compute the center and radius for the Gershgorin th:
  center = beta_ct-gamma_ct-mu_m*(N-1)
  radius = mu_m*(N-1)
}

#--------------- MIGRATION & COMMUTING -------------------------
if(MOB == 2){
  print("Migration and commuting")
  # Generate de Jacobian:
  betas <- matrix(rep(0,N^2), nrow = N)
  diag(betas) <- rep(beta_ct,N)
  BIGT <- commut_mat
  diag(BIGT) <- rep(1,N)
  BIGT <- BIGT%*%betas
  BIGS <- connect_mat
  diag(BIGS) <- rep(-gamma_ct - mu_m*(N-1), N)
  jacobian <- BIGT+BIGS
  
  # Compute the eigenvalues:
  eig <- eigen_mat(jacobian)
  
  # Compute the center and radius for the circular law:
  tau_ct <- 0
  center = beta_ct*(1-mu_w) - N*mu_m - gamma_ct
  radius = sqrt(N*(beta_ct^2*s_w^2 + 2*beta_ct*tau_ct + s_m^2))
  outlier <- beta_ct*(mu_w*(N-1) + 1) - gamma_ct
}

#-----------------------PLOTS------------------------------------
# Generate plots:
plot_eig <- ggplot(eig) + geom_point(aes(re,im), size = 0.05) 

if(MOB == 0 | MOB == 2){
  print("Commuting included")
  plot_eig <- plot_eig + 
    geom_circle(aes(x0 = center,
                    y0 = 0,
                    r = radius), colour = "blue",
                show.legend = NA,size = 0.2) +
    geom_point(aes(outlier,0), colour =  "blue",
               show.legend = NA) +
    coord_fixed() +
    theme_bw() 
}else{
  print("Just migration")
  plot_eig <- plot_eig + 
    geom_circle(aes(x0 = center,
                    y0 = 0,
                    r = radius), colour = "blue",
                show.legend = NA,size = 0.2) +
    coord_fixed() +
    theme_bw() 
}

# plot_eig

# Plot with the segment defining x = 0:
max_im <- max(eig$im) + max(eig$im)/4
df <- data.frame(x1 = 0, x2 = 0, y1 =-max_im, y2 = max_im)
plot_eig <- plot_eig +
  geom_segment(aes(x = 0, y = -max_im, xend = 0, yend = max_im,
                   colour = "segment"), data = df) +
  theme(legend.position = "none")

plot_eig

if(INT == 1){
  print("Doing plots")
  top_row = ggarrange(plot_tot, plot_inf, ncol = 2)
  bottom_row = ggarrange(NULL, plot_eig, NULL, ncol = 3, widths = c(1,4,1), heights = 2)
  final_plot = ggarrange(top_row, bottom_row, ncol = 1)
  # final_plot
  
  plot = ggarrange(plot_inf, plot_eig, nrow = 2)
   
  
  Path <- "/home/marta/Documentos/PHD/2021/R0_SIR_RMT/Plots/"
  gamma_ct <- format(round(gamma_ct,2), decimal.mark = ',')
  beta_ct <- format(round(beta_ct,2), decimal.mark = ',')
  mu_w <- format(round(mu_w,2), decimal.mark = ',')
  s_w <- format(round(s_w,2), decimal.mark = ',')
  mu_m <- format(round(mu_m,2), decimal.mark = ',')
  s_m <- format(round(s_m,2), decimal.mark = ',')
  
  Path <- paste0(Path,"gen","N",N,"g",gamma_ct,"b",beta_ct,"mw",
                 mu_w,"sw",s_w,"mm",mu_m,"sm",s_m,".png")
}

png(file = Path, width = 8000, height = 6000, res = 1100)
plot
dev.off()

plot_uns_int <- plot_inf
plot_uns_eig <- plot_eig

# plot_s_int <- plot_inf
# plot_s_eig <- plot_eig

plot_1 <- ggarrange(plot_s_int,plot_uns_int)
plot_2 <- ggarrange(plot_s_eig,plot_uns_eig, nrow = 2)
