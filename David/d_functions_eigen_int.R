##### LOAD PACKAGES #####

library("tidyverse")
library("parallel")
library("deSolve")
library("ggplot2")

##### FUNCTIONS #####

##### 1.SIR: system of equations #####

SIR <- function(t, y, parameters) {
  with(as.list(c(y, parameters)),{
    dy <- c()
    dim1 <- dime + 1
    dim2 <- dime*2
    dim3 <- (dime*2)+1
    dim4 <- dime*3
    for(i in c(1:dime)){
      # Total population (N = S+I+R)
      N <- y[i] + y[i+dime] + y[i+dim2]      # Susceptible individuals:
      q1 <- delta_N[i]*N
      q2 <- beta_r[i]*(y[i]/N)*y[i+dime]
      q3 <- d[i]*y[i]
      q4 <- theta_r[i]*y[i+dim2]
      q5 <- y[i]*sum(C[,i])
      q6 <- sum(y[1:dime]*C[i,])
      q7 <- (y[i]/N)*sum(beta_r*W[,i]*y[dim1:dim2])
      # dS/dt
      dy[i] <- q1 - q2 - q3 + q4 - q5 + q6 - q7      # Infected individuals:
      q1 <- beta_r[i]*(y[i]/N)*y[i+dime]
      q2 <- (alpha_r[i] + delta_r[i] + d[i])*y[i+dime]
      q3 <- y[dime+i]*sum(C[,i])
      q4 <- sum(y[dim1:dim2]*C[i,])
      q5 <- (y[i]/N)*sum(beta_r*W[,i]*y[dim1:dim2])
      # dI/dt
      dy[i+dime] <- q1 - q2 - q3 + q4 + q5      # Recovered individuals:
      q1 <- alpha_r[i]*y[i+dime]
      q2 <- (theta_r[i] + d[i])*y[i+dim2]
      q3 <- y[i+dim2]*sum(C[,i])
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
    dime = N,
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
    df_plot$type <- substr(df_plot$variable,1,1)
    plot  <- ggplot(df_plot,aes(time, value)) + 
      geom_line(aes(group =variable, colour = type),size=0.5)  +
      ylab("Number of individuals") 
  }else if( state == "SUS"){
    # Filter Susceptibles:
    df_sus <- df_plot  %>% filter( substr(df_plot$variable,1,1) == "S")
    plot  <- ggplot(df_sus,aes(time, value)) + 
      geom_line(aes( colour = variable),size=0.5)  +
      ylab("Number of individuals")  + 
      ggtitle("Susceptible individuals")
  }else if( state == "INF"){
    # Filter Infected:
    df_plot$type <- substr(df_plot$variable,1,1)
    df_inf <- df_plot  %>% filter( substr(df_plot$variable,1,1) == "I")
    plot  <- ggplot(df_inf,aes(time, value)) + 
      geom_line(aes( group =variable, colour = type),size=0.5)  +
      ylab("Number of infected individuals") 
  }else if( state == "REC"){
    # Filter Recovered:
    df_rec <- df_plot  %>% filter( substr(df_plot$variable,1,1) == "R")
    plot  <- ggplot(df_rec,aes(time, value)) + 
      geom_line(aes( colour = variable),size=0.5)  +
      ylab("Number of individuals") +
      ggtitle("Recovered individuals")
  }else if( state == "REC_INF" ){
    # Filter Recovered:
    df_rec_inf <- df_plot  %>% filter( substr(df_plot$variable,1,1) == "R" |  substr(df_plot$variable,1,1) == "I")
    df_rec_inf$type <- substr(df_rec_inf$variable,1,1)
    plot  <- ggplot(df_rec_inf,aes(time, value), show.legend = NA) + 
      geom_line(aes( colour = type),size=0.5)  +
      ylab("Number of individuals") +
      ggtitle("Recovered individuals")
  }
  if(state != "TOT"){
    plot <- plot + theme(legend.position = "none")
  }
  return(plot)
}

#SIR N PATCH DISCONECTED#
SIR1 <- function(t, y, parameters) {
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
      # dS/dt
      dy[i] <- q1 - q2 - q3 + q4 
      
      # Infected individuals:
      q1 <- beta_r[i]*(y[i]/N)*y[i+dim]
      q2 <- (alpha_r[i] + delta_r[i] + d[i])*y[i+dim] 
      # dI/dt
      dy[i+dim] <- q1 - q2 
      
      # Recovered individuals:
      q1 <- alpha_r[i]*y[i+dim] 
      q2 <- (theta_r[i] + d[i])*y[i+2*dim] 
      # dR/dt
      dy[i+dim2] <- q1 - q2 
      
    }
    list(dy)
  }) 
}

int1 <- function(N,Deltas,betas,deaths,thetas,alphas,deltas,COMMUTING,MIGRATION,sus_init,inf_init,end_time){
  
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
  z <- ode(pops, times, SIR1, parameters)
  return(z)
}

births_func <- function(MIGRATION, init_pop, deaths){
  mat <- - MIGRATION
  diag(mat) <- deaths + colSums(MIGRATION)
  mat <- mat%*%init_pop
  return(mat)
}

init_pop_func <- function(MIGRATION, Deltas, deaths){
  mat <- - MIGRATION
  diag(mat) <- deaths + colSums(MIGRATION)
  mat1 <- inv(mat)%*%Deltas
  return(mat1)
}

DFE_func <- function(MIGRATION, deaths, Deltas){
  mat <- - MIGRATION
  diag(mat) <- deaths + colSums(MIGRATION)
  mat1 <- inv(mat)%*%Deltas
  return(mat1)
}
##### LOAD PACKAGES #####

##### 1.SIR: system of equations #####

SIR_cte_pop <- function(t, y, parameters) {
  with(as.list(c(y, parameters)),{
    dy <- c()
    dim1 <- dime + 1
    dim2 <- dime*2
    dim3 <- (dime*2)+1 
    dim4 <- dime*3
    for(i in c(1:dime)){
      # Total population (N = S+I+R)
      N <- tot_pop      # Susceptible individuals:
      # print(paste0("N: ",N))
      q1 <- delta_N[i]
      q2 <- beta_r[i]*(y[i]/N)*y[i+dime]
      q3 <- d[i]*y[i]
      q4 <- theta_r[i]*y[i+dim2]
      q5 <- y[i]*sum(C[,i])
      q6 <- sum(y[1:dime]*C[i,])
      q7 <- (y[i]/N)*sum(beta_r*W[,i]*y[dim1:dim2])
      # dS/dt
      dy[i] <- q1 - q2 - q3 + q4 - q5 + q6 - q7      # Infected individuals:
      q1 <- beta_r[i]*(y[i]/N)*y[i+dime]
      q2 <- (alpha_r[i] + delta_r[i] + d[i])*y[i+dime]
      q3 <- y[dime+i]*sum(C[,i])
      q4 <- sum(y[dim1:dim2]*C[i,])
      q5 <- (y[i]/N)*sum(beta_r*W[,i]*y[dim1:dim2])
      # dI/dt
      dy[i+dime] <- q1 - q2 - q3 + q4 + q5      # Recovered individuals:
      q1 <- alpha_r[i]*y[i+dime]
      q2 <- (theta_r[i] + d[i])*y[i+dim2]
      q3 <- y[i+dim2]*sum(C[,i])
      q4 <- sum(y[dim3:dim4]*C[i,])
      # dR/dt
      dy[i+dim2] <- q1 - q2 - q3 + q4    }
    list(dy)
  })
}

##### 2.integrate the system #####

int_cte_pop <- function(N_dim,Deltas,betas,deaths,thetas,alphas,deltas,COMMUTING,MIGRATION,sus_init,inf_init,end_time){
  
  # create vector of parameters for ode function:
  parameters <- list(
    dime = N_dim,
    delta_N = Deltas,
    beta_r = betas,
    d = deaths,
    theta_r = thetas,
    alpha_r = alphas,
    delta_r = deltas,
    C = COMMUTING,
    W = MIGRATION,
    tot_pop = sus_init[1]+inf_init[1])
  
  # time steps for integration:
  times = seq(0,end_time, 0.1)
  
  # initial values:
  pops <- c(sus_init, inf_init, rep(0,N))
  
  # run integration:
  z <- ode(pops, times, SIR_cte_pop, parameters)
  return(z)
}
