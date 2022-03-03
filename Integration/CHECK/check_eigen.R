rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")
library("ggpubr")
library("ggsci")
library("ggplot2")


#------------------------DETERMINISTIC MATRIX-------------------------#
# 1 patch:
# Parameters: 
dim <- 400
bet <- 1
mu_w <- 0.6
mu_c <- 0.01
alp <- 700*45
mat <- matrix(0, nrow = dim, ncol = dim)

a <- bet*mu_w +mu_c
b <- alp
c <- alp*mu_w

mat[,] <- a
mat[1,1] <- a+b
mat[2:dim,1] <- a+c

eigen_m <- as.complex(eigen(mat, only.values = TRUE)$values)
df <- data.frame(re = Re(eigen_m), im = Im(eigen_m))
plot_eig <- ggplot(df) + geom_point(aes(re,im), size = 0.05) 

bet_cte <- bet
mu_m <- mu_c
N <- dim
bet_new <- bet_cte + alp


outl <- (1/2)*(N*a + b + sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
outl2 <- (1/2)*(N*a + b - sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))

plot_eig + 
  geom_point(aes(outl,0), colour =  "blue",
             show.legend = NA) +
  geom_point(aes(outl2,0), colour =  "blue",
           show.legend = NA) +
  geom_point(aes(0,0), colour =  "blue",
             show.legend = NA)

# K patches:
# Parameters: 
dim <- 100
bet <- 1
mu_w <- 0.6
mu_c <- 0.01
alp <- 700
mat <- matrix(0, nrow = dim, ncol = dim)
k <-  45

a <- bet*mu_w +mu_c
b <- alp
c <- alp*mu_w

mat[,] <- a
for(i in c(1:k)){
  mat[,i] <- a+c
  mat[i,i] <- a+b
}

eigen_m <- as.complex(eigen(mat, only.values = TRUE)$values)
df <- data.frame(re = Re(eigen_m), im = Im(eigen_m))
plot_eig <- ggplot(df) + geom_point(aes(re,im), size = 0.05) 

bet_cte <- bet
mu_m <- mu_c
N <- dim
bet_new <- bet_cte + alp


outl <- (1/2)*(N*a + b + (k-1)*c + sqrt((N*a)^2 - (2*N-4*k)*a*b +
                                          (2*(k+1)*N-4*k)*a*c + b^2 + (2*k-2)*b*c + (k-1)^2*c^2))
outl2 <- (1/2)*(N*a + b + (k-1)*c - sqrt((N*a)^2 - (2*N-4*k)*a*b +
                                           (2*(k+1)*N-4*k)*a*c + b^2 + (2*k-2)*b*c + (k-1)^2*c^2))
outl3 <- b - c

plot_eig + 
  geom_point(aes(outl,0), colour =  "blue",
             show.legend = NA) +
  geom_point(aes(outl2,0), colour =  "blue",
             show.legend = NA) +
  geom_point(aes(outl3,0), colour =  "blue",
             show.legend = NA) +
  geom_point(aes(0,0), colour =  "blue",
             show.legend = NA)

png_files <- list.files("~/Documents/PHD/2022/RMT_SIR/Plots/1patch/OLD", pattern = ".*png$", full.names = TRUE)
gifski(png_files, gif_file = "~/Documents/PHD/2022/RMT_SIR/Plots/1patch/animation.gif", width = 800, height = 600, delay = 0.25)

#------------------------DETERMINISTIC MATRIX-------------------------#
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
N = 200 # Number of patches
# CTE parameters:
del_N <- rep(0.6, N) # Birth rate
bet_cte <- 6
bet_cte <-  0.001
bet <- rep(bet_cte, N)  # Transmission rate
# bet <- abs(rnorm(N,2,5))  # Transmission rate
d_vec <- rep(0.8, N) # Natural mortality rate
thet <- rep(0.1, N) # Rate of loss of immunity
alp <- rep(0.6, N) # Rate of disease overcome
delt <- rep(0, N) # Diseases related mortality rate

# Changing transmission rate:
mu = 0.5
sig = 4
ind <- sample(1:N,1)
bet_new <- 13
# bet_new <-  0.001
bet[ind] <- bet_new + bet_cte
# size <- sample(1:N,1)
# vec_rand <- sample(1:N,size)
# vec_rand <- seq(1,round(N/4),1)
# bet[vec_rand] <- bet_new
# 
print(paste0("gamma:", alp[1] + delt[1] + d_vec[1]))
print(paste0("beta - gamma:", bet[1] - (alp[1] + delt[1] + d_vec[1])))

#-------------------- MOBILITY ------------------#
### Migration:
alp_m <- 0.01
bet_m <- 0.1

# Compute mean and sd:
mu_m <- alp_m/(alp_m + bet_m)
s_m <-  sqrt((alp_m*bet_m)/(((alp_m + bet_m)^2)*(1+alp_m+bet_m)))
print(paste0("mu :", mu_m))
print(paste0("sigma :", s_m))

migrate_mat <- mat_conect(N,alp_m,bet_m,MOB)
### Commuting
alp_c <- 0.01
bet_c <- 0.1

# Compute mean and sd:
mu_w <- alp_c/(alp_c + bet_c)
s_w <- sqrt((alp_c*bet_c)/((alp_c + bet_c)^2*(1+alp_c+bet_c)))
print(paste0("mu :", mu_w))
print(paste0("sigma :",s_w))

commut_mat <- mat_conect(N,alp_c,bet_c,MOB)

tau_ct <- 0

print(paste0("beta - gamma:", bet[1] - (alp[1] + delt[1] + d_vec[1])))

#-------------------------------------------------------------------------------#
beta_ct = bet_cte  # gamma
gamma_ct = alp[1] + delt[1] + d_vec[1]        # beta
gamma_ct_w <- format(round(gamma_ct,2), decimal.mark = ',')
beta_ct_w <- format(round(beta_ct,2), decimal.mark = ',')
mu_w_w <- format(round(mu_w,2), decimal.mark = ',')
s_w_w <- format(round(s_w,2), decimal.mark = ',')
mu_m_w <- format(round(mu_m,2), decimal.mark = ',')
s_m_w <- format(round(s_m,2), decimal.mark = ',')

vec_bet <- seq(0.001,40.001,0.1)
len <-  length(vec_bet)
mat_bet <-  matrix(0,nrow = len, ncol= 3 )
# Integrate the system:

count = 1
for(i in vec_bet){
  print(paste0("New beta : ",i))
  bet_new <- i
  bet[ind] <- bet_new
  
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
  
  mat_bet[count,] <-  c(bet_new, outl, outl2)
  count = count + 1
}

df_bet <-  as.data.frame(mat_bet)
colnames(df_bet) <-  c("alpha", "outlier1", "outlier2")
ggplot(df_bet) + 
  geom_line(aes(alpha,outlier1), color = "red")+ 
  geom_line(aes(alpha,outlier2), color = "blue")


