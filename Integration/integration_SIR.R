rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")
library("ggpubr")
library("ggsci")

# -------------------------PARAMETERS ---------------------------------
N = 2 # Number of patches
mu = 1 
sig = 0.5

# Random generation of parameters:
del_N <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2)) # Birth rate
bet <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2))   # Transmission rate
d_vec <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2)) # Natural mortality rate
thet <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2))  # Rate of loss of immunity 
connect_mat <- matrix(rgamma(N^2,shape = (mu/sig)^2,rate = mu/(sig^2)), nrow = N) # Connectivity matrix
commut_mat <- matrix(rgamma(N^2,shape = (mu/sig)^2,rate = mu/(sig^2)), nrow = N) # Commuting matrix
alp <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2)) # Rate of disease overcome
delt <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2)) # Diseases related mortality rate

# CTEf parameters:
del_N <- 1 # Birth rate
bet <- 1   # Transmission rate
d_vec <- 1 # Natural mortality rate
thet <- 1  # Rate of loss of immunity 
connect_mat <- 0 # Connectivity matrix
commut_mat <- 0 # Commuting matrix
alp <- 1 # Rate of disease overcome
delt <- 1 # Diseases related mortality rate

# Parameter vector
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

# ----------------------------MODEL ----------------------------------
# 1D SIR model:
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
      dS1 <- delta_N[1]*(S1+I1+R1) - 
        beta_r[1]*(S1/(S1+I1+R1))*I1  -
        d[1]*S1 + theta_r[1]*R1 - S1*sum(C[,1]) +
        # (S1*C[1,1] + S2*C[1,2]) 
      # - (S1/(S1+I1+R1))*(beta_r[1]*W[1,1]*I1 +
                                                     # beta_r[2]*W[2,1]*I2)
      dI1 <- beta_r[1]*(S1/(S1+I1+R1))*I1 - 
        (alpha_r[1] + delta_r[1] + d[1])*I1 - I1*(C[1,1]+C[2,1]) 
      # +
        # sum(I1*C[1,])
        # + (S1/(S1+I1+R1))*(beta_r[1]*W[1,1]*I1 +
      #                                                 beta_r[2]*W[2,1]*I2)
      # print("C[1,]:")
      # print(C[1,])
      # print(paste("(I1*C[1,1] + I2*C[1,2]):", (I1*C[1,1] + I2*C[1,2])))
      dR1 <- alpha_r[1]*I1 - (theta_r[1] + d[1])*R1 
      # - R1*(C[1,1]+C[2,1]) 
      # +  (R1*C[1,1] + R2*C[1,2])
      #---------------------------------------------------------------------#
      # Second patch:
      dS2 <- delta_N[2]*(S2+I2+R2) - 
        beta_r[2]*(S2/(S2+I2+R2))*I2  -
        d[2]*S2 + theta_r[2]*R2 - S2*sum(C[,2]) +
        # (S2*C[2,1] + S2*C[2,2]) 
      # - (S2/(S2+I2+R2))*(beta_r[2]*W[1,2]*I2 +
                                                     # beta_r[2]*W[2,2]*I2)
      dI2 <- beta_r[2]*(S2/(S2+I2+R2))*I2 - 
        (alpha_r[2] + delta_r[2] + d[2])*I2 - I2*(C[1,2]+C[2,2]) 
      # +
      #   sum(I2*C[2,]) 
      # + (S2/(S2+I2+R2))*(beta_r[2]*W[1,2]*I2 +
      #                                                 beta_r[2]*W[2,2]*I2)
      # print("C[2,]:")
      # print(C[2,])
      # print(paste(" (I2*C[2,1] + I2*C[2,2]):",  (I2*C[2,1] + I2*C[2,2])))
      dR2 <- alpha_r[2]*I2 - (theta_r[2] + d[2])*R2 
      # - R2*(C[1,2]+C[2,2]) 
      # + (R2*C[2,1] + R2*C[2,2])
    list(c(dS1, dS2, dI1, dI2, dR1, dR2))
  }) 
}

times = seq(0, 1, 0.1)
# Toy model 2 patches:
population <- c(S1 = 100, S2 = 100, I1 = 10, I2 = 10, R1= 0, R2 = 0 )
z <- ode(population, times, SIR_2, parameters)

head(z)
#---------------------------------------------------------------------
# ND: SIR metapopulation model:
SIR <- function(t, y, parameters) {
  with(as.list(c(y, parameters)),{
    dy <- c()
    for(i in c(1:dim)){
      # Total population (N = S+I+R)
      N <- y[i] + y[i+dim] + y[i+2*dim]
      
      # dS/dt
      dy[i] <- delta_N[i]*N - beta_r[i]*(y[i]/N)*y[i+dim] -
                d[i]*y[i] + theta_r[i]*y[i+2*dim]  
      #- y[i]*sum(C[,i])
      # +
                # sum(y[1:dim]*C[i,]) 
      # - (y[i]/N)*sum(beta_r*W[,i]*y[dim+1:2*dim])
      
      # dI/dt
      dy[i+dim] <- beta_r[i]*(y[i]/N)*y[i+dim] -
                    (alpha_r[i] + delta_r[i] + d[i])*y[i+dim] - y[i+dim]*sum(C[,i]) 
      # + sum(y[dim+1:2*dim]*C[i,])
      # + (y[i]/N)*sum(beta_r*W[,i]*y[dim+1:2*dim])
      # print("i:")
      # print(i)
      # print("C[i,]:")
      # print(C[i,])
      # print(paste( "sum(y[dim1:dim2]*C[i,])", sum(y[dim1:dim2]*C[i,])))
      # dR/dt
      dy[i+2*dim] <- alpha_r[i]*y[i+dim] - (theta_r[i] + d[i])*y[i+2*dim]  
      # - y[i+2*dim]*sum(C[,i])
      # + sum(y[2*dim+1:i+3*dim]*C[i,])
    }
    list(dy)
  }) 
}

population <- c(S <- matrix(100,ncol=N,nrow =1 ), I <- matrix(10,ncol=N,nrow =1 ),
                R <- matrix(0,ncol=N,nrow =1 ))
z <- ode(population, times, SIR, parameters)

head(z)

