rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")
library("ggpubr")
library("ggsci")
library("ggforce")

# -------------------------PARAMETERS ---------------------------------
N = 100 # Number of patches
mu = 2 
sig = 3

# Random generation of parameters:
# del_N <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2)) # Birth rate
# bet <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2))   # Transmission rate
# d_vec <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2)) # Natural mortality rate
# thet <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2))  # Rate of loss of immunity 
# alp <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2)) # Rate of disease overcome
# delt <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2)) # Diseases related mortality rate

# CTE parameters:
del_N <- matrix(0.1, ncol = N, nrow = 1) # Birth rate
bet <- matrix(1, ncol = N, nrow = 1)  # Transmission rate
d_vec <- matrix(0.8, ncol = N, nrow = 1) # Natural mortality rate
thet <- matrix(0.001, ncol = N, nrow = 1) # Rate of loss of immunity
alp <- matrix(0.31, ncol = N, nrow = 1) # Rate of disease overcome
delt <- matrix(0.19, ncol = N, nrow = 1) # Diseases related mortality rate

print(paste0("gamma:", alp[1,1] + delt[1,1] + d_vec[1,1]))
#--------------------MOBILITY PARAMETERS --------------------------
# Migration matrix:
mu_m = 0.4
s_m = 0.25
connect_mat <- matrix(rgamma(N^2,shape = (mu_m/s_m)^2,rate = mu_m/(s_m^2)), nrow = N)
# connect_mat <- matrix(0, nrow = N, ncol = N) # No migration

#Commuting matrix:
mu_w = 5
s_w = 30
commut_mat <- matrix(rgamma(N^2,shape = (mu_w/s_w)^2,rate = mu_w/(s_w^2)), nrow = N)
# commut_mat <- matrix(0, nrow = N, ncol = N) # No commuting

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
#------------------------- 1D SIR model--------------------------------
SIR_1 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dS1 <- delta_N[1]*(S1+I1+R1) - beta_r[1]*(S1/(S1+I1+R1))*I1  -
      d[1]*S1 + theta_r[1]*R1
    dI1 <-beta_r[1]*(S1/(S1+I1+R1))*I1 -(alpha_r[1] + delta_r[1] + d[1])*I1
    dR1 <- alpha_r[1]*I1 - (theta_r[1] + d[1])*R1 
    list(c(dS1, dI1, dR1))
  }) 
}

times = seq(0, 1, 0.1)
# Toy model 2 patches:
population <- c(S1 = 100, I1 = 10, R1= 0)
z <- ode(population, times, SIR_1, parameters)

head(z)
#---------------------------------------------------------------------
# 2D: SIR model in two patches:
SIR_2 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    
    # First patch: 
    # Susceptible individuals:
    q1 <- delta_N[1]*(S1+I1+R1)
    q2 <-  beta_r[1]*(S1/(S1+I1+R1))*I1
    q3 <-  d[1]*S1 
    q4 <- theta_r[1]*R1 
    q5 <- S1*sum(C[,1])
    q6 <- (S1*C[1,1] + S2*C[1,2])
    q7 <- (S1/(S1+I1+R1))*(beta_r[1]*W[1,1]*I1 + beta_r[2]*W[2,1]*I2)
    
    dS1 <- q1 - q2  - q3 +  q4 - q5  +  q6 - q7
    
    # Infected individuals:
    q1 <- beta_r[1]*(S1/(S1+I1+R1))*I1
    q2 <- (alpha_r[1] + delta_r[1] + d[1])*I1
    q3 <-  I1*(C[1,1]+C[2,1]) 
    q4 <- (I1*C[1,1] + I2*C[1,2])
    q5 <- (S1/(S1+I1+R1))*(beta_r[1]*W[1,1]*I1 + beta_r[2]*W[2,1]*I2)
    
    dI1 <- q1 - q2  - q3 + q4 + q5
    
    # Recovered individuals:
    q1 <- alpha_r[1]*I1
    q2 <- (theta_r[1] + d[1])*R1
    q3 <- R1*(C[1,1]+C[2,1])
    q4 <- (R1*C[1,1] + R2*C[1,2])
    dR1 <-  q1 - q2  - q3 + q4
    
    #---------------------------------------------------------------------#
    # Second patch:
    
    q1 <- delta_N[2]*(S2+I2+R2)
    q2 <- beta_r[2]*(S2/(S2+I2+R2))*I2 
    q3 <- d[2]*S2 
    q4 <- theta_r[2]*R2
    q5 <- S2*sum(C[,2]) 
    q6 <- (S1*C[2,1] + S2*C[2,2])
    q7 <- (S2/(S2+I2+R2))*(beta_r[1]*W[1,2]*I1 + beta_r[2]*W[2,2]*I2)
    dS2 <- q1 - q2 - q3 + q4 - q5 + q6 - q7
    
    # Infected individuals:
    q1 <- beta_r[2]*(S2/(S2+I2+R2))*I2 
    q2 <- (alpha_r[2] + delta_r[2] + d[2])*I2
    q3 <-  I2*(C[1,2]+C[2,2]) 
    q4 <-  (I1*C[2,1] + I2*C[2,2]) #sum(I2*C[2,])
    q5 <- (S2/(S2+I2+R2))*(beta_r[1]*W[1,2]*I1 + beta_r[2]*W[2,2]*I2)
    dI2 <- q1 - q2 - q3 + q4 + q5
    
    # Recovered individuals:
    q1 <- alpha_r[2]*I2 
    q2 <-  (theta_r[2] + d[2])*R2 
    q3 <-  R2*(C[1,2]+C[2,2])
    q4 <-  (R1*C[2,1] + R2*C[2,2])
    dR2 <-  q1- q2  - q3 + q4
    
    list(c(dS1, dS2, dI1, dI2, dR1, dR2))
  }) 
}

times = seq(0, 0.5, 0.1)
# Toy model 2 patches:
population <- c(S1 = 100, S2 = 100, I1 = 10, I2 = 10, R1= 0, R2 = 0 )
z <- ode(population, times, SIR_2, parameters)

head(z)
rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")
library("ggpubr")
library("ggsci")
library("ggforce")

# -------------------------PARAMETERS ---------------------------------
N = 100 # Number of patches
mu = 2 
sig = 3

# Random generation of parameters:
# del_N <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2)) # Birth rate
# bet <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2))   # Transmission rate
# d_vec <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2)) # Natural mortality rate
# thet <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2))  # Rate of loss of immunity 
# alp <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2)) # Rate of disease overcome
# delt <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2)) # Diseases related mortality rate

# CTE parameters:
N = 2
del_N <- matrix(0.6, ncol = N, nrow = 1) # Birth rate
bet <- matrix( 0.00001, ncol = N, nrow = 1)  # Transmission rate
d_vec <- matrix(0.8, ncol = N, nrow = 1) # Natural mortality rate
thet <- matrix(0.1, ncol = N, nrow = 1) # Rate of loss of immunity
alp <- matrix(0.02, ncol = N, nrow = 1) # Rate of disease overcome
delt <- matrix(0, ncol = N, nrow = 1) # Diseases related mortality rate

print(paste0("gamma:", alp[1,1] + delt[1,1] + d_vec[1,1]))
#--------------------MOBILITY PARAMETERS --------------------------
# Migration matrix:
mu_m = 0.4
s_m = 0.25
connect_mat <- matrix(rgamma(N^2,shape = (mu_m/s_m)^2,rate = mu_m/(s_m^2)), nrow = N)
# connect_mat <- matrix(0, nrow = N, ncol = N) # No migration

#Commuting matrix:
mu_w = 5
s_w = 30
commut_mat <- matrix(rgamma(N^2,shape = (mu_w/s_w)^2,rate = mu_w/(s_w^2)), nrow = N)
# commut_mat <- matrix(0, nrow = N, ncol = N) # No commuting

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
#------------------------- 1D SIR model--------------------------------
SIR_1 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dS1 <- delta_N[1]*(S1+I1+R1) - beta_r[1]*(S1/(S1+I1+R1))*I1  -
      d[1]*S1 + theta_r[1]*R1
    dI1 <-beta_r[1]*(S1/(S1+I1+R1))*I1 -(alpha_r[1] + delta_r[1] + d[1])*I1
    dR1 <- alpha_r[1]*I1 - (theta_r[1] + d[1])*R1 
    list(c(dS1, dI1, dR1))
  }) 
}

times = seq(0, 1, 0.1)
# Toy model 2 patches:
population <- c(S1 = 100, I1 = 10, R1= 0)
z <- ode(population, times, SIR_1, parameters)

head(z)
#---------------------------------------------------------------------
# 2D: SIR model in two patches:
SIR_2 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    
    # First patch: 
    # Susceptible individuals:
    q1 <- delta_N[1]*(S1+I1+R1)
    q2 <-  beta_r[1]*(S1/(S1+I1+R1))*I1
    q3 <-  d[1]*S1 
    q4 <- theta_r[1]*R1 
    q5 <- S1*sum(C[,1])
    q6 <- (S1*C[1,1] + S2*C[1,2])
    q7 <- (S1/(S1+I1+R1))*(beta_r[1]*W[1,1]*I1 + beta_r[2]*W[2,1]*I2)
    
    dS1 <- q1 - q2  - q3 +  q4 - q5  +  q6 - q7
    
    # Infected individuals:
    q1 <- beta_r[1]*(S1/(S1+I1+R1))*I1
    q2 <- (alpha_r[1] + delta_r[1] + d[1])*I1
    q3 <-  I1*(C[1,1]+C[2,1]) 
    q4 <- (I1*C[1,1] + I2*C[1,2])
    q5 <- (S1/(S1+I1+R1))*(beta_r[1]*W[1,1]*I1 + beta_r[2]*W[2,1]*I2)
    
    dI1 <- q1 - q2  - q3 + q4 + q5
    
    # Recovered individuals:
    q1 <- alpha_r[1]*I1
    q2 <- (theta_r[1] + d[1])*R1
    q3 <- R1*(C[1,1]+C[2,1])
    q4 <- (R1*C[1,1] + R2*C[1,2])
    dR1 <-  q1 - q2  - q3 + q4
    
    #---------------------------------------------------------------------#
    # Second patch:
    
    q1 <- delta_N[2]*(S2+I2+R2)
    q2 <- beta_r[2]*(S2/(S2+I2+R2))*I2 
    q3 <- d[2]*S2 
    q4 <- theta_r[2]*R2
    q5 <- S2*sum(C[,2]) 
    q6 <- (S1*C[2,1] + S2*C[2,2])
    q7 <- (S2/(S2+I2+R2))*(beta_r[1]*W[1,2]*I1 + beta_r[2]*W[2,2]*I2)
    dS2 <- q1 - q2 - q3 + q4 - q5 + q6 - q7
    
    # Infected individuals:
    q1 <- beta_r[2]*(S2/(S2+I2+R2))*I2 
    q2 <- (alpha_r[2] + delta_r[2] + d[2])*I2
    q3 <-  I2*(C[1,2]+C[2,2]) 
    q4 <-  (I1*C[2,1] + I2*C[2,2]) #sum(I2*C[2,])
    q5 <- (S2/(S2+I2+R2))*(beta_r[1]*W[1,2]*I1 + beta_r[2]*W[2,2]*I2)
    dI2 <- q1 - q2 - q3 + q4 + q5
    
    # Recovered individuals:
    q1 <- alpha_r[2]*I2 
    q2 <-  (theta_r[2] + d[2])*R2 
    q3 <-  R2*(C[1,2]+C[2,2])
    q4 <-  (R1*C[2,1] + R2*C[2,2])
    dR2 <-  q1- q2  - q3 + q4
    
    list(c(dS1, dS2, dI1, dI2, dR1, dR2))
  }) 
}

times = seq(0, 25, 0.1)
# Toy model 2 patches:
population <- c(S1 = 100000, S2 = 100000, I1 = 100, I2 = 100, R1= 0, R2 = 0 )
z <- ode(population, times, SIR_2, parameters)
 sol <-  as.data.frame(z)
 plot_inf_1 <- plot_int(N, sol, state)
head(z)

