rm(list = ls())
library(easypackages)
libraries("gdata", "ggExtra","ggplot2", "numbers",
          "tidyverse", "MASS", "bivariate", "barsurf",
          "ggforce", "MethylCapSig")

############FUNCTIONS#############
# Function which create a circle:
circle <- function(t, N, mu, sigma){
  x <- mu + sigma*sqrt(N)*cos(t)
  y <- sigma*sqrt(N)*sin(t)
  vec <- c(x,y)
  return(vec)
}

# Function which create a normal random matrix.
rand_norm_mat <- function(N, mu, sigma){
  rmatrix <- matrix(rnorm(N^2, mu, sigma),nrow = N)
  return(rmatrix)
}

# Function which create a elliptic (bivariate) random matrix.
rand_ellip_mat_norm <- function(N, mu1, s1,rho){
  mu2 = mu1
  s2 = s1
  mu <- c(mu1,mu2) # Mean
  sigma <- matrix(c(s1^2, s1*s2*rho, s1*s2*rho, s2^2),
                  2) # Covariance matrix
  biv_norm <- mvrnorm(N^2, mu = mu, Sigma = sigma )
  rmatrix <- matrix(0,N,N)
  for(i in c(1:N)){
    for(j in c(i:N)){
      ind = j + (i-1)*N
      rmatrix[i,j] = biv_norm[ind,1]
      rmatrix[j,i] = biv_norm[ind,2]
    }
  }
  
  return(rmatrix)
}

# Function which create a elliptic (bivariate) random matrix.
rand_ellip_mat_lognorm <- function(N, mu1, s1,rho){
  mu2 = mu1
  s2 = s1
  mu <- c(mu1,mu2) # Mean
  sigma <- c(s1*s2*rho, s1*s2*rho) # Covariance matrix
  rho_mat <- matrix(c(1,rho,rho,1),2)
  biv_norm <- mvlognormal(n = N^2, Mu = mu, 
                          Sigma = sigma, R = rho_mat)
  rmatrix <- matrix(0,N,N)
  for(i in c(1:N)){
    for(j in c(i:N)){
      ind = j + (i-1)*N
      rmatrix[i,j] = biv_norm[ind,1]
      rmatrix[j,i] = biv_norm[ind,2]
    }
  }
  
  return(rmatrix)
}
# Function which change the diagonal of a matrix:
change_diag_norm <- function(N, mat, mu, sigma){
  diag(mat) <- rnorm(N, mu, sigma)
    return(mat)
}

change_diag_lognorm <- function(N, mat, mu, sigma){
  diag(mat) <- rlnorm(N, mu, sigma)
  return(mat)
}

# Create the NGM:
NGM_matrix <- function(N, mu, sigma, mu_d,sigma_d, gam,rho){
  ellip <- rand_ellip_mat_norm(N, mu, sigma, rho)
  mu_diag <- mu_d - mu
  mat_1 <- change_diag_norm(N, ellip, mu_diag, sigma_d)
  mat_2 <- (1/gam)*mat_1
  return(mat_2)
}

NGM_matrix_lognorm <- function(N, mu, sigma, mu_d,sigma_d, gam,rho){
  ellip <- rand_ellip_mat_lognorm(N, mu, sigma, rho)
  mu_diag <- mu_d - mu
  mat_1 <- change_diag_lognorm(N, ellip, mu_diag, sigma_d)
  mat_2 <- (1/gam)*mat_1
  return(mat_2)
}


# Create the Jacobian:
J_matrix <- function(N, mu, sigma, mu_d,sigma_d, gam,rho){
  ellip <- rand_ellip_mat_norm(N, mu, sigma, rho)
  mu_diag <- mu_d - mu
  mat_1 <- change_diag_norm(N, ellip, mu_diag, sigma_d)
  mat_2 <- gam*diag(N) + mat_1
  return(mat_2)
}  

J_matrix_lognorm <- function(N, mu, sigma, mu_d,sigma_d, gam,rho){
  ellip <- rand_ellip_mat_lognorm(N, mu, sigma, rho)
  mu_diag <- mu_d - mu
  mat_1 <- change_diag_lognorm(N, ellip, mu_diag, sigma_d)
  mat_2 <- gam*diag(N) + mat_1
  return(mat_2)
}  

# Function which plot the eigenvalues:
eigen_mat <- function(mat){
  eigen_m <- as.complex(eigen(mat)$values)
  df <- data.frame(re = Re(eigen_m), im = Im(eigen_m))
  return(df)
}


##### Stability conditions ######
# NGM:
NGM_stability <- function(N, mu, sigma, mu_d, rho, gam, eps){
  cond1 = (sigma*sqrt(N)*(1+rho) + (mu_d - mu))*(1/gam)
  cond2 = (-sigma*sqrt(N)*(1+rho) + (mu_d - mu))*(1/gam)
  cond3 = (sigma*sqrt(N)*(1-rho) + (mu_d - mu))*(1/gam)
  cond4 = (-sigma*sqrt(N)*(1-rho) + (mu_d - mu))*(1/gam)
  cond5 = ((mu_d - mu) + mu*N)*(1/gamma)
  stab = "FALSE"
  tol <- 1 + eps
  if(cond1 < tol){
    if(abs(cond2) < tol){
      if(cond3 < tol){
        if(abs(cond4) < tol){
          if(cond5 < tol){
            stab = "TRUE"
          }else{
            print("Condition 5 (outlier) failed")
          }
        }else{
          print("Condition 4 (down y) failed")
        }
      }else{
        print("Condition 3 (up y) failed")
      }
    }else{
      print("Condition 2 (left x) failed")
    }
  }else{
    print("Condition 1 (right x) failed")
  }
  return(stab)
}

# Jacobian:
J_stability <- function(N, mu, sigma, mu_d, rho, gam, eps){
  cond = sigma*sqrt(N)*(1+rho) + (mu_d - mu) - gam
  stab = "FALSE"
  if( cond < eps){
    stab = "TRUE"
  }else{
    print("Condition of stability failed")
  }
  return(stab)
}


############COMPUTATIONS#############
# Normal distribution ######
  mu = 0
  sigma = 0.00001
  N = 100
  mat <- rand_norm_mat(N, mu, sigma)
  eig <- eigen_mat(mat)
  ggplot(eig) +
    geom_point(aes(re,im), size = 0.05)
  
# Circle:
  min_seq <- 0
  max_seq <- 2*pi
  x <- seq(min_seq, max_seq, 0.01)
  circ <- lapply(x, circle, N, mu, sigma)
  l = length(x)
  mat <- matrix(unlist(circ,recursive =TRUE), ncol = l, byrow = FALSE)
  df_mat <- data.frame(t(mat))
  out <- mu*N
  outlier <- c(out,0)
  df_mat <- rbind(df_mat, outlier)
  head(df_mat)
  colnames(df_mat) <- c("X","Y")
  ggplot(df_mat) +
    geom_point(aes(X,Y), size = 0.05, color = "red")+
    geom_point(data = eig, aes(re,im) , size = 0.05)

  
# Bivariate distribution #######
  ## NGM Matrix:
  mu_bi = 3        # Bivariate mean
  sigma_bi = 0.08  # Bivariate variance
  rho = 0.8        # Correlation
  N = 300          # Size matrix
  gam = 0.0003        # gamma
  mu_d = 0.0008 # Media diagonal
  sigma_d = 2   # Variance diagonal
  
  # with normal distribution:
  mat_2 <- NGM_matrix(N, mu_bi, sigma_bi, mu_d,sigma_d, gam,rho)
  eig <- eigen_mat(mat_2)
  outlier <- data.frame(x = ((mu_d - mu_bi) + mu_bi*N)*(1/gam), y = 0)
  plot_J <- ggplot(eig) +
    geom_point(aes(re,im), size = 0.05) +
    geom_ellipse(aes(x0 = (mu_d - mu_bi)*(1/gam), y0 = 0,
                     a = sigma_bi*sqrt(N)*(1+rho)*(1/gam),
                     b = sigma_bi*sqrt(N)*(1-rho)*(1/gam),
                     angle = 0, colour = "red")) +
    geom_point(data = outlier, aes(x,y), color = "red",size = 0.3) +
    ggtitle("NGM matrix normal distribution") +
    theme_bw() 
  
  plot_J
  xmin = (mu_d - mu_bi)*(1/gam) + sigma_bi*sqrt(N)*(1+rho)*(1/gam)
  xmax = (mu_d - mu_bi)*(1/gam) - sigma_bi*sqrt(N)*(1+rho)*(1/gam)
  plot_J + xlim(c(xmin,xmax))  
  
  
  # with lognormal distribution:
  mat_2 <- NGM_matrix_lognorm(N, mu_bi, sigma_bi, mu_d,sigma_d, gam,rho)
  eig <- eigen_mat(mat_2)
  outlier <- data.frame(x = ((mu_d - mu_bi) + mu_bi*N)*(1/gam), y = 0)
  plot_J <- ggplot(eig) +
    geom_point(aes(re,im), size = 0.05) +
    geom_ellipse(aes(x0 = (mu_d - mu_bi)*(1/gam), y0 = 0,
                     a = sigma_bi*sqrt(N)*(1+rho)*(1/gam),
                     b = sigma_bi*sqrt(N)*(1-rho)*(1/gam),
                     angle = 0, colour = "red")) +
    geom_point(data = outlier, aes(x,y), color = "red",size = 0.3) +
    ggtitle("NGM matrix lognormal distribution") +
    theme_bw()
  
  plot_J
  xmin = (mu_d - mu_bi)*(1/gam) + sigma_bi*sqrt(N)*(1+rho)*(1/gam)
  xmax = (mu_d - mu_bi)*(1/gam) - sigma_bi*sqrt(N)*(1+rho)*(1/gam)
  plot_J + xlim(c(xmin,xmax))
  
  ### Jacobian Matrix: 
  # With normal distribution:
  mat_2 <- J_matrix(N, mu_bi, sigma_bi, mu_d,sigma_d, gam,rho)
  eig <- eigen_mat(mat_2)
  outlier <- data.frame(x = ((mu_d - mu_bi) + mu_bi*N) + gam, y = 0)
  plot_J <- ggplot(eig) +
    geom_point(aes(re,im), size = 0.05)+
    geom_ellipse(aes(x0 = (mu_d - mu_bi) - gam, y0 = 0,
                     a = sigma_bi*sqrt(N)*(1+rho),
                     b = sigma_bi*sqrt(N)*(1-rho),
                     angle = 0, colour = "red")) +
    geom_point(data = outlier, aes(x,y), color = "red",size = 0.3) +
    ggtitle("Jacobian matrix normal distribution")+
    theme_bw()
  
  plot_J
  xmin = (mu_d - mu_bi) + sigma_bi*sqrt(N)*(1+rho) - gam
  xmax = (mu_d - mu_bi) - sigma_bi*sqrt(N)*(1+rho) - gam
  plot_J + xlim(c(xmin,xmax))
  
  xmin = -10
  xmax = 10
  plot_J + xlim(c(xmin,xmax))
  
  # with lognormal distribution:
  mat_2 <- J_matrix_lognorm(N, mu_bi, sigma_bi, mu_d,sigma_d, gam,rho)
  eig <- eigen_mat(mat_2)
  outlier <- data.frame(x = ((mu_d - mu_bi) + mu_bi*N) + gam, y = 0)
  plot_J <- ggplot(eig) +
    geom_point(aes(re,im), size = 0.05)+
    geom_ellipse(aes(x0 = (mu_d - mu_bi) - gam, y0 = 0,
                     a = sigma_bi*sqrt(N)*(1+rho),
                     b = sigma_bi*sqrt(N)*(1-rho),
                     angle = 0, colour = "red")) +
    geom_point(data = outlier, aes(x,y), color = "red",size = 0.3)+
    ggtitle("Jacobian matrix lognormal distribution") +
    theme_bw()
  
  plot_J
  xmin = (mu_d - mu_bi) + sigma_bi*sqrt(N)*(1+rho) - gam
  xmax = (mu_d - mu_bi) - sigma_bi*sqrt(N)*(1+rho) - gam
  plot_J + xlim(c(xmin,xmax))
  
  plot_J + xlim(c(-250,250))
# Stability computations:
  mu_vec <- runif(10,0,3)
  sigma_vec <- runif(10,0,3)  
  mud_vec <- runif(10,0,3)
  rho_vec <- runif(10,0,1)
  gam_vec <- runif(10,0,3)
    
  N = 500
  mu = 5
  sigma = 0.001
  mu_d = 1
  rho = 0.2
  gam = 5
  eps = 0.00001
  
  J_stability(N, mu, sigma, mu_d, rho, gam, eps)
  mu_bi = mu
  sigma_bi = sigma
  mat_2 <- J_matrix_lognorm(N, mu_bi, sigma_bi, mu_d,sigma_d, gam,rho)
  eig <- eigen_mat(mat_2)
  outlier <- data.frame(x = ((mu_d - mu_bi) + mu_bi*N) + gam, y = 0)
  plot_J <- ggplot(eig) +
    geom_point(aes(re,im), size = 0.05)+
    geom_ellipse(aes(x0 = (mu_d - mu_bi) +gam, y0 = 0,
                     a = sigma_bi*sqrt(N)*(1+rho),
                     b = sigma_bi*sqrt(N)*(1-rho),
                     angle = 0, colour = "red")) +
    geom_point(data = outlier, aes(x,y), color = "red",size = 0.3)+
    ggtitle("Jacobian matrix lognormal distribution") + 
    theme_bw()
  
    plot_J
    plot_J + xlim(c(-50,50))
    