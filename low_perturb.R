rm(list =ls())
library(easypackages)
libraries("gdata", "ggExtra","ggplot2", "numbers",
          "tidyverse", "MASS", "bivariate", "barsurf",
          "ggforce", "MethylCapSig", "simstudy")

############FUNCTIONS#############
# Function which create a normal random matrix.
rand_norm_mat <- function(N, mu, sigma){
  rmatrix <- matrix(rnorm(N^2, mu, sigma),nrow = N)
  return(rmatrix)
}

# Function which create a elliptic (bivariate) random matrix.
rand_ellip_mat_norm <- function(N, mu1, s1,c){
  mu <- c(mu1,mu1) # Mean
  sig <- matrix(c(s1^2, s1^2*c, s1^2*c,s1^2),
                2) # Covariance matrix
  print(paste("sigma:", sig))
  biv_norm <- mvrnorm(N^2, mu = mu, Sigma = sig )
  rmatrix <- matrix(0,N,N)
  for(i in c(1:N)){
    for(j in c(i:N)){
      ind = j + (i-1)*N
      rmatrix[i,j] = biv_norm[ind,1]
      rmatrix[j,i] = biv_norm[ind,2]
    }
  }
  rmatrix <- rmatrix * s1
  return(rmatrix)
}

# Function which create a elliptic (bivariate) random matrix.
rand_ellip_mat_lognorm <- function(N, mu1, s1,rho){
  mu2 = mu1
  s2 = s1
  mu <- c(mu1,mu2) # Mean
  sig <- c(s1*s2*rho, s1*s2*rho) # Covariance matrix
  print(paste("sigma:", sigma))
  rho_mat <- matrix(c(1,rho,rho,1),2)
  biv_norm <- mvlognormal(n = N^2, Mu = mu, 
                          Sigma = sig, R = rho_mat)
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
  mat_1 <- change_diag_norm(N, ellip, mu_d, sigma_d)
  mat_2 <- (1/(gam+(N-1)*mu))*mat_1
  return(mat_2)
}

NGM_matrix_lognorm <- function(N, mu, sigma, mu_d,sigma_d, gam,rho){
  ellip <- rand_ellip_mat_lognorm(N, mu, sigma, rho)
  mat_1 <- change_diag_lognorm(N, ellip, mu_d, sigma_d)
  mat_2 <- (1/(gam+(N-1)*mu))*mat_1
  return(mat_2)
}


# Create the Jacobian:
J_matrix <- function(N, mu, sigma, mu_d,sigma_d, gam,rho){
  ellip <- rand_ellip_mat_norm(N, mu, sigma, rho)
  mat_1 <- change_diag_norm(N, ellip, mu_d, sigma_d)
  mat_2 <- mat_1 - (gam+(N-1)*mu)*diag(N)
  return(mat_2)
}  

J_matrix_lognorm <- function(N, mu, sigma, mu_d,sigma_d, gam,rho){
  ellip <- rand_ellip_mat_lognorm(N, mu, sigma, rho)
  mat_1 <- change_diag_lognorm(N, ellip, mu_d, sigma_d)
  mat_2 <- (gam+(N-1)*mu)*diag(N) + mat_1
  return(mat_2)
}  

# Function which plot the eigenvalues:
eigen_mat <- function(mat){
  eigen_m <- as.complex(eigen(mat, only.values = TRUE)$values)
  df <- data.frame(re = Re(eigen_m), im = Im(eigen_m))
  return(df)
}

############COMPUTATIONS#############
# Normal distribution ######
mu = 0
sigma = 0.0000000000001
N = 100
mat <- rand_norm_mat(N, mu, sigma)
eig <- eigen_mat(mat)
ggplot(eig) +
  geom_point(aes(re,im), size = 0.05) +
  geom_ellipse(aes(x0 = mu_d, y0 = 0,
                   a = sigma_bi*sqrt(N),
                   b = sigma_bi*sqrt(N),
                   angle = 0, colour = "red"))
