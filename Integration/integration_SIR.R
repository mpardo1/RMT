rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")
library("ggpubr")
library("ggsci")

# -------------------------PARAMETERS ---------------------------------
N = 100 # Number of patches
mu = 1 
sig = 0.5

# Random generation of parameters:
# del_N <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2)) # Birth rate
bet <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2))   # Transmission rate
# d_vec <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2)) # Natural mortality rate
thet <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2))  # Rate of loss of immunity 
alp <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2)) # Rate of disease overcome
delt <- rgamma(N,shape = (mu/sig)^2,rate = mu/(sig^2)) # Diseases related mortality rate

# CTE parameters:
del_N <- matrix(1.2, ncol = N, nrow = 1) # Birth rate
# bet <- 1   # Transmission rate
d_vec <- matrix(0.5, ncol = N, nrow = 1) # Natural mortality rate
# thet <- 1  # Rate of loss of immunity 
# alp <- 1 # Rate of disease overcome
# delt <- 1 # Diseases related mortality rate

#--------------------MOBILITY PARAMETERS --------------------------
# Migration matrix:
connect_mat <- matrix(rgamma(N^2,shape = (mu/sig)^2,rate = mu/(sig^2)), nrow = N) 
# connect_mat <- matrix(0, nrow = N, ncol = N) # No migration

#Commuting matrix:
commut_mat <- matrix(rgamma(N^2,shape = (mu/sig)^2,rate = mu/(sig^2)), nrow = N)
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

times = seq(0, 100, 0.1)

# Vector of initial values:
population <- c(matrix(0,ncol=3*N,nrow =1 ))

# Susceptible initial values:
population[1:N] <- 100
population[(N+1):(2*N)] <- 10

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
df_plot <- reshape2::melt(z, id.vars = c("time"))

ggplot2::theme_set(ggplot2::theme_bw() %+replace% 
                    ggplot2::theme(axis.ticks = 
                    element_line(color = 'black'),axis.title = element_text(color = 'black', size = 15),
                    axis.text = element_text(color = 'black', size = 15),
                    legend.text = element_text(color = 'black', size = 15),
                    legend.title = element_text(color = 'black', size = 16),
                    plot.title = element_text(color = 'black', size = 18, face = 'bold'),
                    strip.text = element_text(color = 'black', size = 15),
                    legend.position = "none"
                    ))

head(df_plot)
# Plot number of individuals at each time.
plot_1  <- ggplot(df_plot,aes(time, value)) + 
  geom_line(aes( colour = variable))  +
  ylab("Number of individuals") 

plot_1

# Filter Susceptibles:
df_sus <- df_plot  %>% filter( substr(df_plot$variable,1,1) == "S")
plot_sus  <- ggplot(df_sus,aes(time, value)) + 
  geom_line(aes( colour = variable))  +
  ylab("Number of individuals")  + 
  ggtitle("Susceptible individuals")
plot_sus

# Filter Infected:
df_inf <- df_plot  %>% filter( substr(df_plot$variable,1,1) == "I")
plot_sus  <- ggplot(df_inf,aes(time, value)) + 
  geom_line(aes( colour = variable))  +
  ylab("Number of individuals") +
  ggtitle("Infected individuals")

plot_sus

# Filter Recovered:
df_rec <- df_plot  %>% filter( substr(df_plot$variable,1,1) == "R")
plot_sus  <- ggplot(df_rec,aes(time, value)) + 
  geom_line(aes( colour = variable))  +
  ylab("Number of individuals") +
  ggtitle("Recovered individuals")

plot_sus
