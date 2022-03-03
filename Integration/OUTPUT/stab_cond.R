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
#--------------------BETA DISTRIBUTION PARAM SPACE----------------------
mu <- seq(0, 1,0.01)
sig <- seq(0,1,0.001)
len <- length(sig)
df.beta <- data.frame(mu =0 ,sigma = 0, val = TRUE)
for(i in c(1:len)){
  func <- function(a){
    validate_mu_s(a, sig[i])
  }
  log.vec <- sapply(mu,func)
  df.prov <- data.frame(mu =mu ,sigma = sig[i], val = log.vec)
  df.beta <- rbind(df.beta,df.prov)
}

ggplot(df.beta) + geom_point(aes(x = mu, y = sigma, colour = val)) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1,
                   colour = "segment"))

#-------------------------STAB COND -------------------------------------
# Randomly generate parameters and check the stability conditions.
# If cond_1 or cond_2 = 0 then Unstable, stable otherwise
len <- 100000
df.param<- data.frame(N = 0, mu_w= 0, mu_c= 0, sig_w= 0,
                      sig_c= 0, tau= 0, gam= 0, bet= 0,
                      cond_1 = 0, cond_2 = 0)
while(nrow(df.param) < 10){
  N_rand <- sample(40:100,1)
  mu_w_rand <- runif(1,0,0.3)
  mu_c_rand <- runif(1,0,0.3)
  sig_w_rand <- runif(1,0,0.25)
  sig_c_rand <- runif(1,0,0.25)
  while(mu_w_rand < sig_w_rand | mu_c_rand < sig_c_rand){
    mu_w_rand <- runif(1,0,0.3)
    mu_c_rand <- runif(1,0,0.3)
    sig_w_rand <- runif(1,0,0.25)
    sig_c_rand <- runif(1,0,0.25)
  }
  gam_rand <- runif(1,0,10)
  bet_rand <- runif(1,0,1)
  cond.stab.1 <- ifelse(cond1(N_rand, mu_w_rand, gam_rand, bet_rand) == TRUE, 1, 0)
  cond.stab.2 <- ifelse(cond2(N_rand, mu_w_rand, mu_c_rand, sig_w_rand,
                              sig_c_rand, 0, gam_rand, bet_rand) == TRUE, 1, 0)
  if(cond.stab.1 == 1 & cond.stab.2 == 0){
    df.param[nrow(df.param) + 1,] <- c(N_rand, mu_w_rand, mu_c_rand,sig_w_rand,
                                       sig_c_rand, 0,  gam_rand, bet_rand, cond.stab.1,
                                       cond.stab.2)
  }
}

which(df.param$cond_1 == 0 && df.param$cond_2 == 1 )
df.param <- df.param[2:nrow(df.param),] 
#-----------------SOL WITH FIRST STAB SECOND NO--------------------
i <- 2
N = df.param$N[i] # Number of patches
# CTE parameters:
del_N <- rep(0.6, N) # Birth rate
bet <- rep(df.param$bet[i],N)
# bet <- abs(rnorm(N,1,1))  # Transmission rate
d_vec <- rep(1, N) # Natural mortality rate
thet <- rep(2.6, N) # Rate of loss of immunity
alp <- rep(df.param$gam[i] - d_vec[1], N) # Rate of disease overcome
delt <- rep(0, N) # Diseases related mortality rate

gamma_ct <-  alp[1] + delt[1] + d_vec[1]
print(paste0("gamma:", alp[1] + delt[1] + d_vec[1]))
print(paste0("beta - gamma:", bet[1] - (alp[1] + delt[1] + d_vec[1])))

#-------------------- MOBILITY ------------------
### Migration:
mu_c <- mu_c_rand
s_c <- sig_c_rand
alp_c <- beta_a_b(mu_c, s_c)[1]
bet_c <- beta_a_b(mu_c, s_c)[2]
migrate_mat <- mat_conect(N,alp_c,bet_c,DIST)
### Commuting
mu_w <- mu_w_rand
s_w <- sig_w_rand
alp_w <- beta_a_b(mu_w, s_w)[1]
bet_w <- beta_a_b(mu_w, s_w)[2]
commut_mat <- mat_conect(N,alp_w,bet_w,DIST)

# Create jacobian
jac <- jacobian(N,bet,gamma_ct, commut_mat, migrate_mat,mu_c, MOB)
eig <- eigen_mat(jac)

plot_eig_cte <- ggplot(eig) + geom_point(aes(re,im), size = 0.05) 
plot_eig_cte + coord_fixed() 


state <- "INF"
plot_int(N, sol, state) + 
  theme_bw() 