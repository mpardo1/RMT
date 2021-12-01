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
# del_N <- 1 # Birth rate
# bet <- 1   # Transmission rate
# d_vec <- 1 # Natural mortality rate
# thet <- 1  # Rate of loss of immunity 
# connect_mat <- 0 # Connectivity matrix
# commut_mat <- 0 # Commuting matrix
# alp <- 1 # Rate of disease overcome
# delt <- 1 # Diseases related mortality rate

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
      
      # print("------------------------------------------------------------")
      # print(paste("t :",t))
      # print(paste("q1 :",q1))
      # print(paste("q2 :",q2))
      # print(paste("q3 :",q3))
      # print(paste("q4 :",q4))
      # print(paste("q5 :",q5))
      # print(paste("dI1:", dI1))
      
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
      
      # print("")
      # print(paste("q1 :",q1))
      # print(paste("q2 :",q2))
      # print(paste("q3 :",q3))
      # print(paste("q4 :",q4))
      # print(paste("q5 :",q5))
      # print(paste("dI2:", dI2))
      
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
#---------------------------------------------------------------------
# ND: SIR metapopulation model:
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
      
      # print("------------------------------------------------------------")
      # print(paste("t :",t))
      # print(paste("i :",i))
      # print(paste("q1 :",q1))
      # print(paste("q2 :",q2))
      # print(paste("q3 :",q3))
      # print(paste("q4 :",q4))
      # print(paste("q5 :",q5))
      # print(paste("dy[i+dim] :",dy[i+dim]))
      
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

population <- c(matrix(0,ncol=3*N,nrow =1 ))
population[1:N] <- 100
N1 <- N+1
N2 <- 2*N
population[N1:N2] <- 10
z <- ode(population, times, SIR, parameters)

head(z)
