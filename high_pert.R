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
  mat_2 <- beta*mat_1 -  gam*diag(N)
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
mu_bi = -2          # Bivariate mean
sigma_bi = 3        # Bivariate variance
N = 300             # Size matrix
beta = 1.2          # gamma
mu_d = 1            # Media diagonal
sigma_d = 0.00001 # Variance diagonal
gam = 13        # beta

# Case where there are two set totally disconected
mat_2 <- J_matrix(N, mu_bi, sigma_bi, mu_d,sigma_d, beta, gam)

mat_2[lower.tri(mat_2)] <- 0
eig <- eigen_mat(mat_2)
plot_J <- ggplot(eig) + geom_point(aes(re,im), size = 0.05) 

plot_J 

df <- data.frame(x1 = 0, x2 = 0, y1 =0, y2 = 5.0)
plot_J <- plot_J +
  geom_segment(aes(x = 0, y = -7, xend = 0, yend = 7,
                   colour = "segment"), data = df) 


plot_J
plot_J +
  coord_fixed()

center = (1/gam)*beta*(1-mu_bi)
radius = (1/gam)*beta*sigma_bi*sqrt(N)
outlier = (mu_bi*(N-1)+1)*(1/gam)*beta
