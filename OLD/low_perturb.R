rm(list =ls())
library(easypackages)
libraries("gdata", "ggExtra","ggplot2", "numbers",
          "tidyverse", "MASS", "bivariate", "barsurf",
          "ggforce", "MethylCapSig", "simstudy")
library("ggpubr")

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
NGM_matrix <- function(N, mu, sigma, mu_d,sigma_d, beta,gamma){
  circ <- rand_norm_mat(N, mu, sigma)
  mat1 <- change_diag_norm(N, circ, mu_d, sigma_d)
  mat2 <- (1/gamma)*beta*mat1
  return(mat2)
}

# Create the Jacobian:
J_matrix <- function(N, mu, sigma, mu_d,sigma_d, beta,gamma){
  circ <- rand_norm_mat(N, mu, sigma)
  mat_1 <- change_diag_norm(N, circ, mu_d, sigma_d)
  mat_2 <- beta*mat_1 -  gamma*diag(N)
  return(mat_2)
}  
# Function which plot the eigenvalues:
eigen_mat <- function(mat){
  eigen_m <- as.complex(eigen(mat, only.values = TRUE)$values)
  df <- data.frame(re = Re(eigen_m), im = Im(eigen_m))
  return(df)
}

# Create a random matrix
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

############COMPUTATIONS#############
# Normal distribution ######
muw = -2          # Bivariate mean
sw = 3        # Bivariate variance
N = 300             # Size matrix
beta = 1.2          # gamma
mu_d = 1            # Media diagonal
sigma_d = 0.00001 # Variance diagonal
gamma = 13        # beta
##--------------------------------------------------------------#
# NGM:.
mat_2 <- NGM_matrix(N, muw, sw, mu_d,sigma_d, beta, gamma)
eig <- eigen_mat(mat_2)
plot_J <- ggplot(eig) + geom_point(aes(re,im), size = 0.05) 

plot_J +
  coord_fixed()

center = (1/gamma)*beta*(1-muw)
radius = (1/gamma)*beta*sw*sqrt(N)
outlier = (muw*(N-1)+1)*(1/gamma)*beta

plot_J + 
  geom_circle(aes(x0 = center,
                  y0 = 0,
                  r = radius, colour = "red")) +
  geom_point(aes(outlier,0), colour =  "red") +
  coord_fixed()

# Jacobian:
mat_2 <- J_matrix(N, muw, sw, mu_d,sigma_d, beta, gamma)
eig <- eigen_mat(mat_2)
plot_J <- ggplot(eig) + geom_point(aes(re,im), size = 0.0005) 

plot_J +
  coord_fixed()

center = beta*(1-muw)-gamma
radius = beta*sw*sqrt(N)
outlier = (muw*(N-1)+1)*beta -gamma

# Gershgorin Theorem:
center_gresh <- beta-gamma
radius_gresh <- beta*muw*(N-1)
radius_gresh1 <- sum(mat_2[1,2:ncol(mat_2)])

plot_J + 
  geom_circle(aes(x0 = center,
                  y0 = 0,
                  r = radius), color = "red",
                  size = 0.2 ,
                  show.legend = NA) + 
  geom_circle(aes(x0 = center_gresh,
                  y0 = 0,
                  r = radius_gresh), colour = "lightgreen",
                  size = 0.2 ,
                  show.legend = NA) +
  geom_point(aes(outlier,0), color =  "blue",
             show.legend = NA) +
  coord_fixed() +
  theme_bw() + scale_fill_discrete(name = "", labels = c("Circular law", "Gershgorin theorem"))

#------------------------------------------------------------------------------------------------#
# Low perturbations of the matrix:
# 1. Perturbed flow for one single flow and a single direction.
muw = 2          # Bivariate mean
sw = 3        # Bivariate variance
N = 300             # Size matrix
beta = 1.2          # gamma
gamma = 13        # beta
muws = 10

betas <- matrix(rep(0,N^2), nrow = N)
diag(betas) <- rep(beta,N)
BIGT <- rand_mat(N, muw, sw, distrib = "gamma")
diag(BIGT) <- rep(1,N)
BIGT[1:2,1] <- BIGT[1:2,1]+muws-muw 
BIGT <- BIGT%*%betas
BIGS <- matrix(rep(0,N^2), nrow = N)
diag(BIGS) <- rep(-gamma, N)
jacobian <- BIGT+BIGS

eig <- eigen_mat(jacobian)

center = beta*(1-muw)-gamma
radius = beta*sw*sqrt(N)
outlier1 <- (muw*(N/2)+sqrt(muw^2*((N^2/4-1))+muw*muws))
outlier2 <- (muw*(N/2)-sqrt(muw^2*((N^2/4-1))+muw*muws))
outlier1 <- beta*(outlier1 + (1-muw)) - gamma
outlier2 <- beta*(outlier2 + (1-muw)) - gamma

plot_J <- ggplot(eig) + geom_point(aes(re,im), size = 0.05) 

plot_J + 
  geom_circle(aes(x0 = center,
                  y0 = 0,
                  r = radius), colour = "red",
              show.legend = NA,size = 0.2) +
  geom_point(aes(outlier1,0), colour =  "red",
             show.legend = NA) +
  geom_point(aes(outlier2,0), colour =  "red",
             show.legend = NA) +
 coord_fixed() +
  theme_bw() 

# 2. Perturbed flow a single patch and for all (income/outcome) directions.
muws <- 20 # New mean
betas <- matrix(rep(0,N^2), nrow = N)
diag(betas) <- rep(beta,N)
BIGT <- rand_mat(N, muw, sw, distrib = "gamma")
diag(BIGT) <- rep(1,N)
BIGT[2:N,1] <- BIGT[2:N,1]+muws-muw 
BIGT <- BIGT%*%betas
BIGS <- matrix(rep(0,N^2), nrow = N)
diag(BIGS) <- rep(-gamma, N)
jacobian <- BIGT+BIGS

eig <- eigen_mat(jacobian)

center = beta*(1-muw)-gamma
radius = beta*sw*sqrt(N)
outlier1 <- (muw*(N/2)+sqrt(muw^2*((N^2-2)/4)+muw*muws*(N-1)))
outlier2 <- (muw*(N/2)-sqrt(muw^2*((N^2-2)/4)+muw*muws*(N-1)))
outlier1 <- beta*(outlier1 + (1-muws)) - gamma
outlier2 <- beta*(outlier2 + (1-muws)) - gamma

plot_J <- ggplot(eig) + geom_point(aes(re,im), size = 0.05) 

plot_J + 
  geom_circle(aes(x0 = center,
                  y0 = 0,
                  r = radius),size = 0.2 , colour = "red",
              show.legend = NA) +
  geom_point(aes(outlier1,0), colour =  "blue",
             show.legend = NA) +
  geom_point(aes(outlier2,0), colour =  "green",
             show.legend = NA) +
  coord_fixed() +
  theme_bw() 

# 3. Perturbed flow for one single flow in both directions.
N <- 300 
muw <- 3
beta <- 2
gamma <- 4
mues <- 20 # New mean
sw = 3    # Bivariate variance

betas <- matrix(rep(0,N^2), nrow = N)
diag(betas) <- rep(beta,N)
BIGT <- rand_mat(N, muw, sw, distrib = "gamma")
diag(BIGT) <- rep(1,N)
BIGT[2,1] <- BIGT[2,1]+muws-muw 
BIGT[1,2] <- BIGT[1,2]+muws-muw 
BIGT <- BIGT%*%betas
BIGS <- matrix(rep(0,N^2), nrow = N)
diag(BIGS) <- rep(-gamma, N)
jacobian <- BIGT+BIGS

eig <- eigen_mat(jacobian)

center = beta*(1-muw)-gamma
radius = beta*sw*sqrt(N)

outlier1 <- (muw*((N-1)/2)+muws/2+sqrt((muw^2)*((N^2+2*N-7)/4)-muw*muws*((N-3)/2)+(muws^2/4)))
outlier2 <- (muw*((N-1)/2)+muws/2-sqrt((muw^2)*((N^2+2*N-7)/4)-muw*muws*((N-3)/2)+(muws^2/4)))
outlier3 <- muw - muws
outlier1 <- beta*(outlier1 + (1-muw)) - gamma
outlier2 <- beta*(outlier2 + (1-muw)) - gamma
outlier3 <- beta*(outlier3 + (1-muw)) - gamma

plot_J <- ggplot(eig) + geom_point(aes(re,im), size = 0.05) 

plot_J <- plot_J + 
  geom_circle(aes(x0 = center,
                  y0 = 0,
                  r = radius), size = 0.2 ,colour = "red",
              show.legend = NA) +
  geom_point(aes(outlier1,0), colour =  "blue",
             show.legend = NA) +
  geom_point(aes(outlier2,0), colour =  "forestgreen",
             show.legend = NA) +
  geom_point(aes(outlier3,0), colour =  "purple",
             show.legend = NA)

plot_J  +
  coord_fixed() + theme_bw()

# 4. Perturbed for two patches and for all (income/outcome) directions.
N <- 200 
muw <- 3
beta <- 2
gamma <- 4
muws <- 20 # New mean
sw = 3     # Bivariate variance

center = beta*(1-muw)-gamma
radius = beta*sw*sqrt(N)

betas <- matrix(rep(0,N^2), nrow = N)
diag(betas) <- rep(beta,N)
BIGT <- rand_mat(N, muw, sw, distrib = "gamma")
diag(BIGT) <- rep(1,N)
BIGT[2:N,1] <- BIGT[2:N,1]+muws-muw 
BIGT[1,2] <- BIGT[1,2]+muws-muw 
BIGT[3:N,2] <- BIGT[3:N,2]+muws-muw 
BIGT <- BIGT%*%betas
BIGS <- matrix(rep(0,N^2), nrow = N)
diag(BIGS) <- rep(-gamma, N)
jacobian <- BIGT+BIGS

eig <- eigen_mat(jacobian)

outlier1 <- (1/2)*((N-1)*muw+muws+sqrt(((N-3)^2)*muw^2 + (6*N-10)*muw*muws+muws^2))
outlier2 <- (1/2)*((N-1)*muw+muws-sqrt(((N-3)^2)*muw^2 + (6*N-10)*muw*muws+muws^2))
outlier1 <- beta*(outlier1 + (1-muw)) - gamma
outlier2 <- beta*(outlier2 + (1-muw)) - gamma
outlier3 <- beta*((muw - muws) + (1-muw)) - gamma

plot_J <- ggplot(eig) + geom_point(aes(re,im), size = 0.05) 

c(outlier1,outlier2,outlier3)

plot_J + 
  geom_circle(aes(x0 = center,
                  y0 = 0,
                  r = radius),size = 0.2 , colour = "red",
              show.legend = NA) +
  geom_point(aes(outlier1,0), colour =  "blue",
             show.legend = NA) +
  geom_point(aes(outlier2,0), colour =  "forestgreen",
             show.legend = NA) +
  geom_point(aes(outlier3,0), colour =  "purple",
             show.legend = NA) +
  coord_fixed() + theme_bw()

# 5. Perturbed for three patches and for all (income/outcome) directions.
N <- 150 
muw <- 0.3
beta <- 0.03
gamma <- 1.5
muws <- 3 # New mean
sw <- 1

center = beta*(1-muw)-gamma
radius = beta*sw*sqrt(N)

betas <- matrix(rep(0,N^2), nrow = N)
diag(betas) <- rep(beta,N)
BIGT <- rand_mat(N, muw, sw, distrib = "gamma")
diag(BIGT) <- rep(1,N)
BIGT <- BIGT%*%betas
BIGS <- matrix(rep(0,N^2), nrow = N)
diag(BIGS) <- rep(-gamma, N)
jacobian <- BIGT+BIGS
eig <- eigen_mat(jacobian)

outlier = (muw*(N-1)+1)*beta -gamma

plot_1 <- ggplot(eig) + geom_point(aes(re,im), size = 0.05)  +
  geom_segment(aes(x = 0, y = -0.5, xend = 0, yend = 0.5,
                   colour = "segment"), data = df)  + 
  geom_circle(aes(x0 = center,
                  y0 = 0,
                  r = radius),size = 0.2 , colour = "red",
              show.legend = NA) +
  geom_point(aes(outlier,0), colour =  "blue",
             show.legend = NA) +
  coord_fixed() +
  theme_bw()
plot_1

betas <- matrix(rep(0,N^2), nrow = N)
diag(betas) <- rep(beta,N)
BIGT <- rand_mat(N, muw, sw, distrib = "gamma")
diag(BIGT) <- rep(1,N)
BIGT[2:N,1] <- BIGT[2:N,1]+muws-muw 
BIGT[1,2] <- BIGT[1,2]+muws-muw 
BIGT[3:N,2] <- BIGT[3:N,2]+muws-muw 
BIGT[1:2,3] <- BIGT[1:2,3]+muws-muw 
BIGT[4:N,3] <- BIGT[4:N,3]+muws-muw 
BIGT <- BIGT%*%betas
BIGS <- matrix(rep(0,N^2), nrow = N)
diag(BIGS) <- rep(-gamma, N)
jacobian <- BIGT+BIGS

eig_pert <- eigen_mat(jacobian)

outlier1 <- (muw*((N-2)/2)+muws+sqrt((muw^2)*((N-4)^2/4)+muw*muws*(2*N-5)+(muws^2)))
outlier2 <- (muw*((N-2)/2)+muws-sqrt((muw^2)*((N-4)^2/4)+muw*muws*(2*N-5)+(muws^2)))
outlier3 <- muw - muws
outlier1 <- beta*(outlier1 + (1-muw)) - gamma
outlier2 <- beta*(outlier2 + (1-muw)) - gamma
outlier3 <- beta*(outlier3 + (1-muw)) - gamma
outlier = (muw*(N-1)+1)*beta -gamma


plot_J_pert <- ggplot(eig_pert) + geom_point(aes(re,im), size = 0.05) 

plot_J_pert

plot_J_pert <- plot_J_pert + 
  geom_circle(aes(x0 = center,
                  y0 = 0,
                  r = radius),size = 0.2 , colour = "red",
              show.legend = NA) +
  geom_point(aes(outlier1,0), colour =  "blue",
             show.legend = NA) +
  geom_point(aes(outlier2,0), colour =  "forestgreen",
             show.legend = NA) +
  geom_point(aes(outlier3,0), colour =  "purple",
             show.legend = NA) +
  coord_fixed() + theme_bw()

plot_J_pert

df <- data.frame(x1 = 0, x2 = 0, y1 =0, y2 = 5.0)
plot_2 <- plot_J_pert +
  geom_segment(aes(x = 0, y = -0.5, xend = 0, yend = 0.5,
                   colour = "segment"), data = df)

plot_2

ggarrange(plot_1, plot_2,  legend = "none", ncol = 1, nrow = 2)

# 6. Perturbed for four patches and for all (income/outcome) directions.
N <- 200 
muw <- 3
beta <- 2
gamma <- 2
muws <- 2 # New mean

center = beta*(1-muw)-gamma
radius = beta*sw*sqrt(N)

mat_2 <- J_matrix(N, muw, sw, mu_d,sigma_d, beta, gamma)
mat_2[2:N,1] <- rnorm(N-1,muws,sw)*beta
mat_2[1,2] <- rnorm(1,muws,sw)*beta
mat_2[3:N,2] <- rnorm(N-2,muws,sw)*beta
mat_2[1:2,3] <- rnorm(2,muws,sw)*beta
mat_2[4:N,3] <- rnorm(N-3,muws,sw)*beta
mat_2[1:3,4] <- rnorm(3,muws,sw)*beta
mat_2[5:N,4] <- rnorm(N-4,muws,sw)*beta
eig <- eigen_mat(mat_2)

outlier1 <- (muw*((N-3)/2)+(3/2)*muws+sqrt(muw^2*((N-5)^2/4)+muw*muws*(5*N-17))+((9/4)*muws^2))
outlier2 <- (muw*((N-3)/2)+(3/2)*muws-sqrt(muw^2*((N-5)^2/4)+muw*muws*(5*N-17))+((9/4)*muws^2))
outlier3 <- muw - muws
outlier1 <- beta*(outlier1 + (1-muw)) - gamma
outlier2 <- beta*(outlier2 + (1-muw)) - gamma
outlier3 <- beta*(outlier3 + (1-muw)) - gamma

plot_J <- ggplot(eig) + geom_point(aes(re,im), size = 0.05) 

plot_J + 
  geom_circle(aes(x0 = center,
                  y0 = 0,
                  r = radius),size = 0.2 , colour = "red",
              show.legend = NA) +
  geom_point(aes(outlier1,0), colour =  "blue",
             show.legend = NA) +
  geom_point(aes(outlier2,0), colour =  "forestgreen",
             show.legend = NA) +
  geom_point(aes(outlier3,0), colour =  "purple",
             show.legend = NA) +
  coord_fixed() + theme_bw()

outlier = (muw*(N-1)+1)*beta -gamma
plot_J + 
  geom_circle(aes(x0 = center,
                  y0 = 0,
                  r = radius),size = 0.2 , colour = "red",
              show.legend = NA) +
  geom_point(aes(outlier,0), colour =  "blue",
             show.legend = NA) +
  coord_fixed() + theme_bw()

# 7. Change trnasmision rates.
N <- 300 
muw <- 3
beta <- 2
gamma <- 2
muws <- 20 # New mean

center = beta*(1-muw)-gamma
radius = beta*sw*sqrt(N)

mat_2 <- J_matrix(N, muw, sw, mu_d,sigma_d, beta, gamma)
mat_2[2:N,1] <- rnorm(N-1,(muw*alpha/beta)+1,sw)*beta
mat_2[1,1] <- rnorm(1,1+(alpha/beta),sw)*beta
eig <- eigen_mat(mat_2)

outlier1 <- (muw*(N/2)+(1/2)*(alpha/beta)+sqrt(muw^2*((N^2)/4)+(N-1)*muw^2*(alpha/beta)-(N/2-1)*muw*(alpha/beta)+(1/4)*((alpha^2)/(beta^2))))
outlier2 <- (muw*(N/2)+(1/2)*(alpha/beta)-sqrt(muw^2*((N^2)/4)+(N-1)*muw^2*(alpha/beta)-(N/2-1)*muw*(alpha/beta)+(1/4)*((alpha^2)/(beta^2))))
outlier1 <- beta*(outlier1 + (1-muw)) - gamma
outlier2 <- beta*(outlier2 + (1-muw)) - gamma

plot_J <- ggplot(eig) + geom_point(aes(re,im), size = 0.05) 

plot_J + 
  geom_circle(aes(x0 = center,
                  y0 = 0,
                  r = radius, colour = "red")) +
  geom_point(aes(outlier1,0), colour =  "red") +
  geom_point(aes(outlier2,0), colour =  "red") +
  coord_fixed()

# Comparison of the outliers.

vec <- seq(50,1000,1)

# Functions of the ouliers:
out_simp <- function(N){
  out <- (muw*(N-1)+1)*beta -gamma
  return(out)
}


out1 <- function(N){
  out <- (muw*(N/2)+sqrt(muw^2*((N^2/4-1))+muw*muws))*beta-gamma
  return(out)
}

out2 <- function(N){
  out <- (muw*(N/2)+sqrt(muw^2*((N^2-2)/4)+muw*muws*(N-1)))*beta-gamma
  return(out)
}

out3 <- function(N){
  out <-  (muw*((N-1)/2)+muws/2+sqrt(muw^2*((N^2+2*N-7)/4)-muw*muws*(N-3)/2)+(muws^2/4))*beta-gamma
  return(out)
}

out4 <- function(N){
  out <- (muw*((N-1)/2)+muws/2+sqrt(muw^2*((N-3)^2/4)-muw*muws*(3*N-5)/2)+(muws^2/4))*beta-gamma
  return(out)
}

out5 <- function(N){
  out <- (muw*((N-2)/2)+muws+sqrt(muw^2*((N-4)^2/4)-muw*muws*(3*N-5))+(muws^2))*beta-gamma
  return(out)
}

out6 <- function(N){
  out <- (muw*((N-3)/2)+(3/2)*muws+sqrt(muw^2*((N-5)^2/4)-muw*muws*(5*N-17))+((9/4)*muws^2))*beta-gamma
  return(out)
}

# for(i in c(1:N)){
#   a <- (muw*(i/2)+sqrt(muw^2*((i^2/4-1))+muw*muws))*beta-gamma
#   if(is.na(a)){
#     print(paste0("Iteration: ",i))
#     print("na number")
#   }
# }
out_1 <- as.numeric(unlist(lapply(vec,out1)))
out_2 <- as.numeric(unlist(lapply(vec,out2)))
out_3 <- as.numeric(unlist(lapply(vec,out3)))
out_4 <- as.numeric(unlist(lapply(vec,out4)))
out_5 <- as.numeric(unlist(lapply(vec,out5)))
out_6 <- as.numeric(unlist(lapply(vec,out6)))
out_simp <- as.numeric(unlist(lapply(vec,out_simp)))
df_out <- data.frame(x =vec, out_simp = out_simp, outlier1 = out_1,
                     outlier2 = out_2, outlier3 = out_3,
                     outlier4 = out_4, outlier5 = out_5,
                     outlier6 = out_6)

df_plot <- reshape2::melt(df_out, id.vars = c("x"))
ggplot(df_plot,aes(x, value)) +
  geom_line(aes( colour = variable)) 
