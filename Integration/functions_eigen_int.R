#----------------FUNCTIONS-----------------
# SIR multipatch model:
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

mat_conect <- function(N,alp,bet,DIST){
  rmatrix <- dplyr::case_when(
    DIST == "normal" ~ matrix(rnorm(N^2, alp, bet), nrow = N, ncol = N),
    DIST == "gamma" ~ matrix(rgamma(N^2, alp, bet), nrow = N),
    DIST == "beta" ~ matrix(rbeta(N^2, alp, bet, ncp = 0), nrow = N),
    TRUE ~ matrix(rep(0, N^2), nrow = N)
  )
  rmatrix <- matrix(rmatrix, nrow = N)
  return(rmatrix)
}

diff_f <- function(migrate_mat,d,init_pop){
  # Constant population at each patch (ie. no mortality induced by the disease):
  diff <- c()
  for(i in c(1:N)){
    diff[i] <- sum(migrate_mat[,i]) - sum(migrate_mat[i,]*init_pop)/init_pop[i] 
  }
  if(min(diff) < 0){
    # Change mortality to not have del_N negative
    d <- abs(min(diff)) + 0.01
    print(paste0("Change d:",d))
  }
  del_N <- diff + d
  return(c(d,del_N))
}

# function which integrate:
int <- function(N, del_N,bet,d_vec,thet,alp,delt, commut_mat,migrate_mat ,end_time, MOB, CTE_POP, CTE_INF,SUS_INIT, INF_INIT, init_pop){
  # Create vector of parameters for ode function:
  parameters <- list(
    dim = N,
    delta_N = del_N,
    beta_r = bet,
    d = d_vec,
    theta_r = thet,
    C = migrate_mat,
    W = commut_mat,
    alpha_r = alp,
    delta_r = delt )
  
  times = seq(0,end_time, 0.1)
  
  # Vector of initial values:
  pops <- c(matrix(0,ncol=3*N,nrow =1 ))
  # Initial value for infected individuals:
  if(CTE_INF == 1){
    pops[(N+1):(2*N)] <- INF_INIT
  }else{
    pops[(N+1):(2*N)] <- ceiling(abs(rnorm(N,100,60)))
  }
  
  
  if(CTE_POP == 1){
    print("No constant population")
    pops[1:N] <- SUS_INIT
  }else{
    print("Constant population")
    # Susceptible individuals:
    pops[1:N] <- init_pop - pops[(N+1):(2*N)]
  }
  # Run integration:
  z <- ode(pops, times, SIR, parameters)
  return(z)
}

# Functions plots integration:
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
      geom_line(aes( colour = variable), show.legend = FALSE)  +
      ylab("Number of individuals")  + 
      ggtitle("Susceptible individuals")
  }else if( state == "INF"){
    # Filter Infected:
    df_inf <- df_plot  %>% filter( substr(df_plot$variable,1,1) == "I")
    plot  <- ggplot(df_inf,aes(time, value)) + 
      geom_line(aes( colour = variable), show.legend = FALSE)  +
      ylab("Number of infected individuals") 
  }else if( state == "REC"){
    # Filter Recovered:
    df_rec <- df_plot  %>% filter( substr(df_plot$variable,1,1) == "R")
    plot  <- ggplot(df_rec,aes(time, value)) + 
      geom_line(aes( colour = variable), show.legend = FALSE)  +
      ylab("Number of individuals") +
      ggtitle("Recovered individuals")
  }
  return(plot)
}

# Function which plot the eigenvalues:
eigen_mat <- function(mat){
  eigen_m <- as.complex(eigen(mat, only.values = TRUE)$values)
  df <- data.frame(re = Re(eigen_m), im = Im(eigen_m))
  return(df)
}

pred_radius <- function(N, beta_ct, gamma_ct, tau_ct, mu_c, s_c, mu_w, s_w, MOB){
  radius <- dplyr::case_when(
    MOB == 0 ~ beta_ct*s_w*sqrt(N),
    MOB == 1 ~ mu_c*(N-1),
    MOB == 2 ~ sqrt(N*(beta_ct^2*s_w^2 + 2*beta_ct*tau_ct + s_c^2)),
    TRUE ~ 0
  )
  return(radius)
}

pred_center <- function(N, beta_ct, gamma_ct, tau_ct, mu_c, s_c, mu_w, s_w, MOB){
  radius <- dplyr::case_when(
    MOB == 0 ~ beta_ct*(1-mu_w)-gamma_ct,
    MOB == 1 ~ beta_ct-gamma_ct-mu_c*(N-1),
    MOB == 2 ~ beta_ct*(1-mu_w) - N*mu_c - gamma_ct,
    TRUE ~ 0
  )
  return(radius)
}

pred_outlier <- function(N, beta_ct, gamma_ct, mu_w, MOB){
  radius <- dplyr::case_when(
    MOB == 0 ~ beta_ct*(mu_w*(N-1)+1)-gamma_ct,
    MOB == 1 ~ 0,
    MOB == 2 ~ beta_ct*(mu_w*(N-1) + 1) - gamma_ct,
    TRUE ~ 0
  )
  return(radius)
}

# Function which compute eigen distribution predicted:
jacobian <- function(N,beta_ct,gamma_ct, commut_mat, migrate_mat,mu_c, MOB){
  if(MOB == 0){
    print("Just commuting")
    # Generate de Jacobian:
    betas <- matrix(rep(0,N^2), nrow = N)
    diag(betas) <- beta_ct
    BIGT <- commut_mat
    diag(BIGT) <- rep(1,N)
    BIGT <- BIGT%*%betas
    BIGS <- matrix(rep(0,N^2), nrow = N)
    diag(BIGS) <- -gamma_ct
    jacobian <- BIGT+BIGS
  }else if(MOB == 1){
    print("Just migration")
    # Generate de Jacobian:
    BIGT <- migrate_mat
    diag(BIGT) <- beta_ct + rep(- gamma_ct - mu_c*(N-1),N)
    jacobian <- BIGT
  }else if(MOB == 2){
    print("Migration and commuting")
    # Generate de Jacobian:
    betas <- matrix(rep(0,N^2), nrow = N)
    diag(betas) <- beta_ct
    BIGT <- commut_mat
    diag(BIGT) <- rep(1,N)
    BIGT <- BIGT%*%betas
    BIGS <- migrate_mat
    diag(BIGS) <- -gamma_ct - colSums(migrate_mat)
    jacobian <- BIGT+BIGS
  }
  else if(MOB == 3){
    print("Migration and commuting with aprox sum cij")
    # Generate de Jacobian:
    betas <- matrix(rep(0,N^2), nrow = N)
    diag(betas) <- beta_ct
    BIGT <- commut_mat
    diag(BIGT) <- rep(1,N)
    BIGT <- BIGT%*%betas
    BIGS <- migrate_mat
    diag(migrate_mat) <- 0
    vec <-  -gamma_ct - (N-1)*mu_c
    diag(BIGS) <- vec
    jacobian <- BIGT+BIGS
  }
  return(jacobian)
}

# Function which plots the eig and the predicted dist:
plot_eigen <- function(eig, center, radius, outlier, MOB){
  # Simplest plot:
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
  
  # Plot with the segment defining x = 0:
  max_im <- max(eig$im) + max(eig$im)/4
  df <- data.frame(x1 = 0, x2 = 0, y1 =-max_im, y2 = max_im)
  plot_eig <- plot_eig +
    geom_segment(aes(x = 0, y = -max_im, xend = 0, yend = max_im,
                     colour = "segment"), data = df) +
    theme(legend.position = "none")
  
  return(plot_eig)
}

cond_gen <- function(N, mu_c,s_c, mu_w,s_w, gam, bet, tau){
  cond1 <- (N-1)*mu_w - (gam/bet) + 1
  cond2 <- mu_w + sqrt(N*(s_w^2 + 2*(tau/bet)+(s_c/bet)^2))-(gam/bet) + 1 - N*(mu_c/bet)
  return(c(cond1,cond2))
}


# Function that validates the mu and sigma for a beta distribution:
# ¡¡ Sigma siempre está al cuadrado en mis cálculos !!
validate_mu_s <-  function(mu,sigma){
  l <-  TRUE
  c <- mu*(1-mu)
  if( mu < 0 | mu > 1 | sigma > 0.25 | c < sigma){
    l <-  FALSE
  }
  return(l)
}

# Compute a and b from mu and sigma of a beta distribution:
# ¡¡ Sigma siempre está al cuadrado en mis cálculos !!
beta_a_b <-  function(mu, sigma){
  l <-  validate_mu_s(mu, sigma)
  if( l == TRUE){
    a <-  mu*((mu/sigma)*(1-mu)-1)
    print(paste0("a: ",a))
    b <- ((1-mu)/mu)*a
  }else{
    print("Problem with mu and sigma")
  }
  return(c(a,b))
}

# Compute mean and variance from a and b beta distribution:
comp_mean_var <- function(vec){
  a <- vec[1]
  b <- vec[2]
  mean <-a/(a+b)
  var <- a*b/((a+b)^2*(a+b+1))
  return(c(mean, var))
}

# Functions which compute the outlier when varying beta in 1 patch:
outl_1patch <- function(alp, bet_cte, N, mu_w, mu_c, gamma_ct){
  a <- bet_cte*mu_w + mu_c
  b <- alp
  c <- alp*mu_w
  
  # print(paste("bet_cte: ",bet_cte))
  # print(paste("alp: ", alp))
  # print(paste("mu_w: ", mu_w))
  
  outl <- (1/2)*(N*a + b + sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
  outl2 <- (1/2)*(N*a + b - sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
  outl <- outl + (bet_cte*(1-mu_w) - N*mu_c - gamma_ct)
  outl2 <- outl2 + (bet_cte*(1-mu_w) - N*mu_c - gamma_ct)
  return(c(outl,outl2))
}

out_gen <-  function(N,beta_ct,  mu_w, gamma_ct){
  out <-  beta_ct*(mu_w*(N-1) + 1) - gamma_ct
  return(out)
}

outl_Kpatch <- function(K, alp, bet_cte, N, mu_w, mu_c, gamma_ct){
  a <- bet_cte*mu_w +mu_c
  b <- alp
  c <- alp*mu_w
  
  outl <- (1/2)*(N*a + b + (K-1)*c + sqrt((N*a)^2 - (2*N-4*K)*a*b +
                                            (2*(K+1)*N-4*K)*a*c + b^2 + 
                                            (2*K-2)*b*c + (K-1)^2*c^2))
  outl <- outl + (bet_cte*(1-mu_w) - N*mu_c - gamma_ct)
  
  return(outl)
}

# Functions which compute the outlier when varying beta in 1 patch:
outl_1patch <- function(alp, bet_cte, N, mu_w, mu_c, gamma_ct){
  a <- bet_cte*mu_w + mu_c
  b <- alp
  c <- alp*mu_w
  
  outl <- (1/2)*(N*a + b + sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
  outl2 <- (1/2)*(N*a + b - sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
  outl <- outl + (bet_cte*(1-mu_w) - N*mu_c - gamma_ct)
  outl2 <- outl2 + (bet_cte*(1-mu_w) - N*mu_c - gamma_ct)
  return(c(outl,outl2))
}

alphagamma <- function(mu,sig) {
  alphagamma <- (mu/sig)^2
}
betagamma <- function(mu,sig) {
  betagamma <- mu/(sig^2)
}

cond1 <- function(N, mu_w, gam, bet){
  cond <- (N-1)*mu_w - (gam/bet) + 1
  log <- TRUE
  if(cond > 0 ){
    log <- FALSE
  }
  return(log)
}

cond2 <- function(N, mu_w, mu_c, sig_w, sig_c, tau, gam, bet){
  sig_hat <- sqrt(sig_w^2 + 2*(tau/bet) + (sig_c/bet)^2 )
  cond <- mu_w + sqrt(N)*sig_hat - (gam/bet) + 1 - N*(mu_c/bet)
  log <- TRUE
  if(cond > 0 ){
    log <- FALSE
  }
  return(log)
}
