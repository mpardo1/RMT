##### LOAD PACKAGES #####

library("tidyverse")
library("parallel")
library("deSolve")
library("ggplot2")
library("ggpubr")
library("ggsci")
library("ggforce")
library("gifski")

##### FUNCTIONS #####

##### 1.SIR: system of equations #####

SIR <- function(t, y, parameters) {
  with(as.list(c(y, parameters)),{
    dy <- c()
    dim1 <- dim + 1
    dim2 <- dim*2
    dim3 <- (dim*2)+1
    dim4 <- dim*3
    for(i in c(1:dim)){
      # Total population (N = S+I+R)
      N <- y[i] + y[i+dim] + y[i+dim2]      # Susceptible individuals:
      q1 <- delta_N[i]*N
      q2 <- beta_r[i]*(y[i]/N)*y[i+dim]
      q3 <- d[i]*y[i]
      q4 <- theta_r[i]*y[i+dim2]
      q5 <- y[i]*sum(C[,i])
      q6 <- sum(y[1:dim]*C[i,])
      q7 <- (y[i]/N)*sum(beta_r*W[,i]*y[dim1:dim2])
      # dS/dt
      dy[i] <- q1 - q2 - q3 + q4 - q5 + q6 - q7      # Infected individuals:
      q1 <- beta_r[i]*(y[i]/N)*y[i+dim]
      q2 <- (alpha_r[i] + delta_r[i] + d[i])*y[i+dim]
      q3 <- y[dim+i]*sum(C[,i])
      q4 <- sum(y[dim1:dim2]*C[i,])
      q5 <- (y[i]/N)*sum(beta_r*W[,i]*y[dim1:dim2])
      # dI/dt
      dy[i+dim] <- q1 - q2 - q3 + q4 + q5      # Recovered individuals:
      q1 <- alpha_r[i]*y[i+dim]
      q2 <- (theta_r[i] + d[i])*y[i+2*dim]
      q3 <- y[i+2*dim]*sum(C[,i])
      q4 <- sum(y[dim3:dim4]*C[i,])
      # dR/dt
      dy[i+dim2] <- q1 - q2 - q3 + q4    }
    list(dy)
  })
}

##### 2.integrate the system #####

int <- function(N,Deltas,betas,deaths,thetas,alphas,deltas,COMMUTING,MIGRATION,sus_init,inf_init,end_time){
  
  # create vector of parameters for ode function:
  parameters <- list(
    dim = N,
    delta_N = Deltas,
    beta_r = betas,
    d = deaths,
    theta_r = thetas,
    alpha_r = alphas,
    delta_r = deltas,
    C = COMMUTING,
    W = MIGRATION)

  # time steps for integration:
  times = seq(0,end_time, 0.1)

  # initial values:
  pops <- c(sus_init, inf_init, rep(0,N))

  # run integration:
  z <- ode(pops, times, SIR, parameters)
  return(z)
}

##### 3.plot the result #####

plot_int <- function(N, z, state){
  # Change labels:
  for(i in c(1:N)){
    colnames(z)[i+1] <-  paste0("S",i)
    colnames(z)[N+i+1] <-  paste0("I",i)
    colnames(z)[2*N+i+1] <-  paste0("R",i)
  }
  
  z <- as.data.frame(z)
  # z  <- z   %>% filter( z$time < 1) 
  df_plot <- reshape2::melt(z, id.vars = c("time"))
  
  head(df_plot)
  # Plot number of individuals at each time.
  if( state == "TOT"){
    plot  <- ggplot(df_plot,aes(time, value)) + 
      geom_line(aes( colour = variable))  +
      ylab("Number of individuals") 
  }else if( state == "SUS"){
    # Filter Susceptibles:
    df_sus <- df_plot  %>% filter( substr(df_plot$variable,1,1) == "S")
    plot  <- ggplot(df_sus,aes(time, value)) + 
      geom_line(aes( colour = variable))  +
      ylab("Number of individuals")  + 
      ggtitle("Susceptible individuals")
  }else if( state == "INF"){
    # Filter Infected:
    df_inf <- df_plot  %>% filter( substr(df_plot$variable,1,1) == "I")
    plot  <- ggplot(df_inf,aes(time, value)) + 
      geom_line(aes( colour = variable),size=1)  +
      ylab("Number of infected individuals") 
  }else if( state == "REC"){
    # Filter Recovered:
    df_rec <- df_plot  %>% filter( substr(df_plot$variable,1,1) == "R")
    plot  <- ggplot(df_rec,aes(time, value)) + 
      geom_line(aes( colour = variable))  +
      ylab("Number of individuals") +
      ggtitle("Recovered individuals")
  }
  return(plot + theme(legend.position = "none"))
}