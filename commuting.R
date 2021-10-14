rm(list =ls())
library(easypackages)
libraries("gdata", "ggExtra","ggplot2", "numbers",
          "tidyverse", "MASS", "bivariate", "barsurf",
          "ggforce", "MethylCapSig", "simstudy")

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

# Function which change the diagonal of a matrix:
change_diag_norm <- function(N, mat, mu, sigma){
  diag(mat) <- rnorm(N, mu, sigma)
  return(mat)
}


# Create the NGM:
NGM_matrix <- function(N, mu, sigma, mu_d,sigma_d, beta,gam){
  circ <- rand_norm_mat(N, mu, sigma)
  mat1 <- change_diag_norm(N, circ, mu_d, sigma_d)
  mat2 <- (1/gam)*beta*mat1
  return(mat2)
}

# Create the Jacobian:
J_matrix <- function(N, mu, sigma, mu_d,sigma_d, beta,gam){
  circ <- rand_norm_mat(N, mu, sigma)
  mat_1 <- change_diag_norm(N, circ, mu_d, sigma_d)
  mat_2 <- mat_1 - (beta + gam)*diag(N)
  return(mat_2)
}  

# Function which plot the eigenvalues:
eigen_mat <- function(mat){
  eigen_m <- as.complex(eigen(mat, only.values = TRUE)$values)
  df <- data.frame(re = Re(eigen_m), im = Im(eigen_m))
  return(df)
}

############COMPUTATIONS#############
## NGM Matrix:
mu_bi = 20          # Bivariate mean
sigma_bi = 10        # Bivariate variance
N = 100             # Size matrix
beta = 1.2          # gamma
mu_d = 1            # Media diagonal
sigma_d = 0.0000001 # Variance diagonal
gam = 1.4          # beta
##--------------------------------------------------------------#
# Check the eigenvalues with same distribution in the diagonal.
mat_2 <- NGM_matrix(N, mu_bi, sigma_bi, mu_d,sigma_d, beta, gam)
eig <- eigen_mat(mat_2)
plot_J <- ggplot(eig) + geom_point(aes(re,im), size = 0.05) 

plot_J +
   coord_fixed()
# plot_J +
#   xlim(c(-50,50))

center = (1/gam)*beta*(1-mu_bi)
radius = (1/gam)*beta*sigma_bi*sqrt(N)
outlier = (mu_bi*(N-1)+1)*(1/gam)*beta
# plot_J + 
#   geom_circle(aes(x0 = center,
#                    y0 = 0,
#                    r = radius, colour = "red")) +
#   geom_point(aes(outlier,0), colour =  "red") +
#   xlim(c(-50,50))

plot_J + 
  geom_circle(aes(x0 = center,
                  y0 = 0,
                  r = radius, colour = "red")) +
  geom_point(aes(outlier,0), colour =  "red") 
