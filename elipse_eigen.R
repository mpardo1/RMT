rm(list = ls())
library(easypackages)
libraries("gdata", "ggExtra","ggplot2", "numbers",
          "tidyverse", "MASS", "bivariate", "barsurf")

############FUNCTIONS#############
# Function which create a circle:
circle <- function(t, N, mu, sigma){
  x <- mu + sigma*sqrt(N)*cos(t)
  y <- sigma*sqrt(N)*sin(t)
  vec <- c(x,y)
  return(vec)
}
# Function which create an ellipse.
# For the NGM:
elipse_NGM <- function(t, N, mu, mu_d, sigma, rho, gamma){
  x <- (mu_d - mu)/gamma + (sigma*sqrt(N)*(1+rho)/gamma)*cos(t)
  y <- (sigma*sqrt(N)*(1-rho)/gamma)*sin(t)
  vec <- c(x,y)
  return(vec)
}
# For the Jacobian:
elipse_J <- function(t, N, mu, mu_d, sigma, rho, gamma){
  x <- (mu_d - mu + gamma) + (sigma*sqrt(N)*(1+rho))*cos(t)
  y <- (sigma*sqrt(N)*(1-rho))*sin(t)
  vec <- c(x,y)
  return(vec)
}

# Function which create a normal random matrix.
rand_norm_mat <- function(N, mu, sigma){
  rmatrix <- matrix(rnorm(N^2, mu, sigma),nrow = N)
  return(rmatrix)
}

# Bivariate distribution
rbvn<-function (n, m1, s1, m2, s2, rho)
{
  X1 <- rnorm(n, mu1, s1)
  X2 <- rnorm(n, mu2 + (s2/s1) * rho *
                (X1 - mu1), sqrt((1 - rho^2)*s2^2))
  cbind(X1, X2)
}

# Function which create a elliptic (bivariate) random matrix.
rand_ellip_mat <- function(N, mu1, s1,rho){
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

# Function which change the diagonal of a matrix:
change_diag <- function(N, mat,sigma, mu){
  diag(mat) <- rnorm(N, mu, sigma)
    return(mat)
}

NGM_matrix <- function(N, mu, sigma, mu_d,sigma_d, gamma){
  mat <- rand_ellip_mat(N, mu, sigma)
  mat_1 <- change_diag(N, mat, mu_d, sigma_d)
  id_gamma <- 1/gamma*diag(N)
  mat_1 <- mat_1 %*% id_gamma
  return(mat_1)
}
# Function which plot the eigenvalues:
eigen_mat <- function(mat){
  eigen_m <- as.complex(eigen(mat)$values)
  df <- data.frame(re = Re(eigen_m), im = Im(eigen_m))
  return(df)
}




############COMPUTATIONS#############
# Normal distribution:
  mu = 0.01
  sigma = 0.009
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

# Bivariate distribution:
  N = 20
  mu1 = 1
  mu2 = 1
  s1 = 0.4
  s2 = 0.4
  rho = 0.6
  bvn3 <- rbvn(N,mu1,s1,mu2,s2,rho)
  
  
  
  
  
  
# Variables for the bivariate distribution:
  # sigma_bi <-rbind(c(1.2,0.1), c(0.1,1.2)) # Covariance matrix
  # mu_bi <-c(3,3) 
  # rho <- cov2cor(sigma_bi)[1,2]
  # gamma = 0.4
# Diagonal distribution
  mu_d <- 2
  sigma_d <- 0.4
# Size of the NxN matrix
  N = 300

# Bivariate distribution:
  mu_bi = 1.2
  sigma_bi = 0.2
  rho = 0.8
  N = 800
  mat <- rand_ellip_mat(N,mu_bi,sigma_bi,rho)
  eig <- eigen_mat(mat)
  ggplot(eig) +
    geom_point(aes(re,im), size = 0.05)+
    xlim(c(-10,10))
  
  
  mat <- change_diag(N,mat, mu,sigma)
  eig <- eigen_mat(mat)
  ggplot(eig) +
    geom_point(aes(re,im), size = 0.05)+
    xlim(c(-20,20))
  
  mat <- NGM_matrix(N, mu, sigma, mu_d,sigma_d, gamma)
  eig <- eigen_mat(mat)
  ggplot(eig) +
    geom_point(aes(re,im), size = 0.05)

############# Draw the ellipse  ################:
  min_seq <- 0
  max_seq <- 2*pi
  x <- seq(min_seq, max_seq, 0.01)
  sigma <- sigma_bi[1,1]
  lap <- lapply(x, elipse_NGM, N, mu, mu_d, sigma, rho, gamma)
  l = length(x)
  mat <- matrix(unlist(lap,recursive =TRUE), ncol = l, byrow = FALSE)
  
  df_mat <- data.frame(t(mat))
  out <- (mu* sqrt(N) + (mu_d - mu))*gamma^(-1)
  outlier <- c(out,0)
  df_mat <- rbind(df_mat, outlier)
  head(df_mat)
  colnames(df_mat) <- c("X","Y")
# df_1 <- data.frame(x,vec = -vec)
# df_join <- rbind(df,df_1)
  df_center <- data.frame(x = ((mu_d - mu))*gamma^(-1), y=0)
  ggplot(df_mat) +
    geom_point(aes(X,Y), size = 0.05)+
    geom_point(data = df_center, aes(x,y), color = "red")

