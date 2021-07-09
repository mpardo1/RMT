#sLibraries
library(matlib)
##########CHECKING CONSISTENCY RMT WITH J AND NGM###############

#####FUNCTIONS#######

# Function which returns the stability of the Jacobian matrix.
J_stab <- function(mat){
  eigen_max = max(eigen(mat)$values)
  stab = "unstable"
  if( eigen_max < 0){
    stab = "stable"
  }
  return(stab)
}
# Function which returns the stability of the NGM.
R0_stab <- function(mat){
  rho = max(Mod(eigen(mat)$values))
  stab = "unstable"
  if(rho < 1){
    stab = "stable"
  }
  return(stab)
}
# Conditions for stability:
# 1º condition J and NGM
cond_1 <- function(beta, gamma){
  cond = "false"
  if(beta < gamma){
    cond = "true"
  }
  return(cond)
}
# 2º condition J and NGM
cond_2 <- function(beta, gamma, mu, rho, N,sigma){
  cond = "false"
  if((sigma/mu) < (1/(sqrt(N)*(1+rho)))*(((gamma-beta)/mu)+N)){
    cond = "true"
  }
  return(cond)
}
# 3º condition NGM
cond_3 <- function(beta, gamma, mu, rho, N,sigma){
  cond = "false"
  if((sigma/mu) < (1/(sqrt(N)*(1+rho)))*(((gamma+beta)/mu)+(N-2))){
    cond = "true"
  }
  return(cond)
}
# 4º condition NGM
cond_4 <- function(beta, gamma, mu, rho, N,sigma){
  cond = "false"
  if((sigma/mu) < (1/(sqrt(N)*(1-rho)))*(((gamma)/mu)+(N-1))){
    cond = "true"
  }
  return(cond)
}

# Conditions J
cond_J <- function(cond1,cond2){
  cond = "Conditions not satisfied"
  if(cond1 = "true" && cond2 = "true"){
    cond = "Conditions satisfied"
  }
  return(cond)
}

# Conditions NGM
cond_NGM <- function(cond1,cond2,cond3,cond4){
  cond = "Conditions not satisfied"
  if(cond1 = "true" && cond2 = "true" && cond3 = "true" && cond4 = "true" ){
    cond = "Conditions satisfied"
  }
  return(cond)
}

# Plot eigenvalues:
plot_eigen <- function(mat){
  plot(eigen(mat)$values, xlab = "real", ylab = "imaginary")
}

# Function which create a normal random matrix.
rand_norm_mat <- function(N, mu, sigma){
  rmatrix <- matrix(rnorm(N^2, mu, sigma),nrow = N)
  return(rmatrix)
}

# 1º Function sample two numbers from a bivariate normal distribution with correlation rho.
biv_norm <- function(mu, sigma, rho){
    x = rnorm(2, mu, sigma) 
    x3 = rho*x[1] + sqrt(1-rho^2)*x[2]
    y1 = mu + sigma*x[1]
    y2 = mu + sigma*x3
    return(c(y1,y2))
}

# Function that creates a elliptic random matrix.
elip_mat <- function(N,mu, sigma, rho){
  mat = matrix(0,N,N)
  for(i in c(1:N)){
    for(j in c(i:N)){
      x = biv_norm(mu, sigma, rho)
      mat[i,j] = x[1]
      mat[j,i] = x[2]
    }
  }
  return(mat)
}

J_matrx <- function(N,mu,sigma,rho, beta, sigma2, gamma){
  mat <- elip_mat(N,mu,sigma,rho)
  diag(mat) <- rnorm(N, beta,sigma2)
  C = mu*(N-1)
  mat <- mat - (gamma+C)*diag(N)
  return(mat)
}

# Jacobian matrix for SIR.
J_matrix <- function(N,mu,sigma,rho, beta, sigma2, gamma,sigma3){
  mat <- elip_mat(N,mu,sigma,rho)
  mat_aux <- mat
  diag(mat_aux) = 0
  diag(mat) <- rnorm(N, beta,sigma2) - rowSums(mat_aux)
  V = diag(N)
  diag(V) = rnorm(N,gamma,sigma3)
  mat <- mat - V
  return(mat)
}
# NGM for SIR.
NGM_matrix <- function(N,mu,sigma,rho, beta, sigma2, gamma,sigma3){
  mat <- elip_mat(N,mu,sigma,rho)
  mat_aux <- mat
  diag(mat_aux) = 0
  diag(mat) <- rnorm(N, beta,sigma2) - rowSums(mat_aux)
  V <- diag(N)
  diag(V) = rnorm(N,gamma,sigma3)
  V = inv(V)
  mat <- mat%*%V
  return(mat)
}
# Example
#function(N,mu,sigma,rho, beta, sigma2, gamma)
J <- J_matrix(300, 0.5,0.3,-0.3,2,0.3,0.2,0.00002)
plot_eigen(J)
NGM <- NGM_matrix(300, 0.5,0.3,-0.3,2,0.3,0.2,0.00002)
plot_eigen(NGM)



