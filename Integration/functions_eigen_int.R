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

mat_conect <- function(N,alp,bet,MOB){
  rmatrix <- dplyr::case_when(
    MOB == 0 ~ matrix(0, nrow = N, ncol = N),
    MOB == 1 ~ matrix(rbeta(N^2, alp, bet, ncp = 0), nrow = N),
    MOB == 2 ~ matrix(rbeta(N^2, alp, bet, ncp = 0), nrow = N),
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
      geom_line(aes( colour = variable))  +
      ylab("Number of individuals")  + 
      ggtitle("Susceptible individuals")
  }else if( state == "INF"){
    # Filter Infected:
    df_inf <- df_plot  %>% filter( substr(df_plot$variable,1,1) == "I")
    plot  <- ggplot(df_inf,aes(time, value)) + 
      geom_line(aes( colour = variable))  +
      ylab("Number of infected individuals") 
  }else if( state == "REC"){
    # Filter Recovered:
    df_rec <- df_plot  %>% filter( substr(df_plot$variable,1,1) == "R")
    plot  <- ggplot(df_rec,aes(time, value)) + 
      geom_line(aes( colour = variable))  +
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
    diag(BIGS) <- rep(-gamma_ct, N)
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
    # else if(MOB == 3){
  #   # print("Migration and commuting with sum cij")
  #   # Generate de Jacobian:
  #   betas <- matrix(rep(0,N^2), nrow = N)
  #   diag(betas) <- beta_ct
  #   BIGT <- commut_mat
  #   diag(BIGT) <- rep(1,N)
  #   BIGT <- BIGT%*%betas
  #   BIGS <- migrate_mat
  #   diag(migrate_mat) <- 0
  #   vec <-  -gamma_ct - colSums(migrate_mat)
  #   diag(BIGS) <- vec
  #   jacobian <- BIGT+BIGS
  # }else{
  #   print("Migration and commuting with 1/N")
  #   # Generate de Jacobian:
  #   betas <- matrix(rep(0,N^2), nrow = N)
  #   diag(betas) <- beta_ct
  #   BIGT <- commut_mat
  #   diag(BIGT) <- rep(1,N)
  #   BIGT <- (1/N)*(BIGT%*%betas)
  #   BIGS <- (1/N)*migrate_mat
  #   diag(migrate_mat) <- 0
  #   vec <-  -(1/N)*gamma_ct - (1/N)*colSums(migrate_mat)
  #   diag(BIGS) <- vec
  #   jacobian <- BIGT+BIGS
  # }
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

# Check the difference between the outlier and the predicted one. 3 is because of
# the model with the sum cij.
# And 2: is to do the prediction with the commuting and migration.
check_outl <-  function(N,beta_ct,gamma_ct,alp_w,bet_w, mu_w,alp_c,bet_c, MOB, df_filt){
  # the 2 is to use the model with commuting and migration:
  mig_mat <- mat_conect(N,alp_c,bet_c,MOB)
  l <- length(which(is.na(mig_mat)))
  count = 1
  while(l > 0){
    print("Migration matrix with NAN")
    mig_mat <- mat_conect(N,alp_c,bet_c,MOB)
    l <- length(which(is.na(mig_mat)))
    if(count > 100 ){
      print("Set mig_mt to 10000")
      print(paste0("alp_m",alp_c))
      print(paste0("bet_m",bet_c))
      mig_mat <- matrix(10000, ncol = N, nrow= N)
      break
    }
    count = count + 1
  }
  com_mat <- mat_conect(N,alp_w,bet_w,MOB)
  
  # The 0 its because we dont use the mean of migration in the jacobian as before.
  # The 3 its to use the sum cij in the diagonal terms
  jac <- jacobian(N,beta_ct,gamma_ct, com_mat, mig_mat,0, MOB)
  eig <- eigen_mat(jac)
  
  # the 2 is to use the model with commuting and migration:
  outl <- pred_outlier(N, beta_ct, gamma_ct, mu_w, MOB)
  max_eig <-  eig$re[which(eig$re == max(eig$re))]   
  diff <- abs(outl - max_eig)/abs(max_eig)
  
  
  # Plot eigenvalues:
  ind <-  which(eig$re == max(eig$re))
  eig_filt <-  eig[-ind,]
  mig <- df_filt[df_filt$mean == mu_w,]
  
  # Compute the mean and sigma for the migration
  mu <- mig[which(mig$a == alp_c | mig$b == bet_c  ),1]
  sig <- mig[which(mig$a == alp_c | mig$b == bet_c  ),2]
  
  # Compute the predicted radius and center
  rad <- pred_radius(N, beta_ct, gamma_ct, 0, mu, sig, mu_w, sigma_w, MOB)
  cent <- pred_center(N, beta_ct, gamma_ct, 0, mu, sig, mu_w, sigma_w, MOB)
  
  plot_eig_noout <- ggplot(eig) +
    geom_point(aes(re,im), size = 0.05) + 
    geom_circle(aes(x0 = cent,
                    y0 = 0,
                    r = rad), colour = "blue",
                show.legend = NA,size = 0.2)+
    coord_fixed() + 
    ggtitle(paste0("Commuting param:", expression(mu),":", mu,", ", expression(sigma),":", sig))
  
  print("Param:")
  print(paste0("N",N))
  print(paste0("alp_c",alp_c))
  print(paste0("bet_c",bet_c))
  print(paste0("mu",mu))
  print(paste0("sig",sig))
  print(paste0("alp_w",alp_w))
  print(paste0("bet_w",bet_w))
  print(paste0("mu_w",mu_w))
  print(paste0("sigma_w",sigma_w))
  
  
  return(list(alp_c = alp_c, bet_c = bet_c, diff= diff, eigen = eig, plot = plot_eig_noout, com_mat = com_mat, mig_mat = mig_mat))
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
