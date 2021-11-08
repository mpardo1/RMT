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
mu_bi = 2          # Bivariate mean
sigma_bi = 3        # Bivariate variance
N = 300             # Size matrix
beta = 1.2          # gamma
mu_d = 1            # Media diagonal
sigma_d = 0.00001 # Variance diagonal
gam = 3        # beta
##--------------------------------------------------------------#
# NGM:.
mat_2 <- NGM_matrix(N, mu_bi, sigma_bi, mu_d,sigma_d, beta, gam)
eig <- eigen_mat(mat_2)
plot_J <- ggplot(eig) + geom_point(aes(re,im), size = 0.05) 

plot_J +
  coord_fixed()

center = (1/gam)*beta*(1-mu_bi)
radius = (1/gam)*beta*sigma_bi*sqrt(N)
outlier = (mu_bi*(N-1)+1)*(1/gam)*beta

plot_J + 
  geom_circle(aes(x0 = center,
                  y0 = 0,
                  r = radius, colour = "red")) +
  geom_point(aes(outlier,0), colour =  "red") +
  coord_fixed()

# Jacobian:
mat_2 <- J_matrix(N, mu_bi, sigma_bi, mu_d,sigma_d, beta, gam)
eig <- eigen_mat(mat_2)
plot_J <- ggplot(eig) + geom_point(aes(re,im), size = 0.0005) 

plot_J +
  coord_fixed()

center = beta*(1-mu_bi)-gam
radius = beta*sigma_bi*sqrt(N)
outlier = (mu_bi*(N-1)+1)*beta -gam

# Gershgorin Theorem:
center_gresh <- beta-gam
radius_gresh <- beta*mu_bi*(N-1)
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

# Low perturbations of the matrix:
# 1. Perturbed flow for one single flow and a single direction.
mu_new <- 0.5 # New mean
mat_2 <- J_matrix(N, mu_bi, sigma_bi, mu_d,sigma_d, beta, gam)
mat_2[1,2] <- rnorm(1,mu_new,sigma_bi)
eig <- eigen_mat(mat_2)

outlier1 <- (mu_bi*(N/2)+sqrt(mu_bi^2*((N^2/4-1))+mu_bi*mu_new))
outlier2 <- (mu_bi*(N/2)-sqrt(mu_bi^2*((N^2/4-1))+mu_bi*mu_new))
outlier1 <- beta*(outlier1 + (1-mu_bi)) - gam
outlier2 <- beta*(outlier2 + (1-mu_bi)) - gam

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
mu_new <- 2 # New mean
mat_2 <- J_matrix(N, mu_bi, sigma_bi, mu_d,sigma_d, beta, gam)
mat_2[2:N,1] <- rnorm(N-1,mu_new,sigma_bi)
eig <- eigen_mat(mat_2)

outlier1 <- (mu_bi*(N/2)+sqrt(mu_bi^2*((N^2-2)/4)+mu_bi*mu_new*(N-1)))
outlier2 <- (mu_bi*(N/2)-sqrt(mu_bi^2*((N^2-2)/4)+mu_bi*mu_new*(N-1)))
outlier1 <- beta*(outlier1 + (1-mu_bi)) - gam
outlier2 <- beta*(outlier2 + (1-mu_bi)) - gam

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
# N <- 700
N <- 300 
mu_bi <- 3
beta <- 2
gam <- 4
mu_new <- 3 # New mean

center = beta*(1-mu_bi)-gam
radius = beta*sigma_bi*sqrt(N)
mat_2 <- J_matrix(N, mu_bi, sigma_bi, mu_d,sigma_d, beta, gam)
mat_2[2,1] <- rnorm(1,mu_new,sigma_bi)
mat_2[1,2] <- rnorm(1,mu_new,sigma_bi)
eig <- eigen_mat(mat_2)

outlier1 <- (mu_bi*((N-1)/2)+mu_new/2+sqrt(mu_bi^2*((N^2+2*N-7)/4)+mu_bi*mu_new*(N-3)/2)+(mu_new^2/4))
outlier2 <- (mu_bi*((N-1)/2)+mu_new/2-sqrt(mu_bi^2*((N^2+2*N-7)/4)+mu_bi*mu_new*(N-3)/2)+(mu_new^2/4))
outlier3 <- mu_bi - mu_new
outlier1 <- beta*(outlier1 + (1-mu_bi)) - gam
outlier2 <- beta*(outlier2 + (1-mu_bi)) - gam
outlier3 <- beta*(outlier3 + (1-mu_bi)) - gam

plot_J <- ggplot(eig) + geom_point(aes(re,im), size = 0.05) 

plot_J + 
  geom_circle(aes(x0 = center,
                  y0 = 0,
                  r = radius), size = 0.2 ,colour = "red",
              show.legend = NA) +
  geom_point(aes(outlier1,0), colour =  "blue",
             show.legend = NA) +
  geom_point(aes(outlier2,0), colour =  "forestgreen",
             show.legend = NA) +
  geom_point(aes(outlier3,0), colour =  "purple",
             show.legend = NA) +
  coord_fixed() + theme_bw()

# 4. Perturbed for two patches and for all (income/outcome) directions.
N <- 300 
mu_bi <- 3
beta <- 2
gam <- 2
mu_new <- 4 # New mean

center = beta*(1-mu_bi)-gam
radius = beta*sigma_bi*sqrt(N)

mat_2 <- J_matrix(N, mu_bi, sigma_bi, mu_d,sigma_d, beta, gam)
mat_2[2:N,1] <- rnorm(N-1,mu_new,sigma_bi)
mat_2[1,2] <- rnorm(1,mu_new,sigma_bi)
mat_2[3:N,2] <- rnorm(N-2,mu_new,sigma_bi)
eig <- eigen_mat(mat_2)

outlier1 <- (mu_bi*((N-1)/2)+mu_new/2+sqrt(mu_bi^2*((N-3)^2/4)+mu_bi*mu_new*(3*N-5)/2)+(mu_new^2/4))
outlier2 <- (mu_bi*((N-1)/2)+mu_new/2-sqrt(mu_bi^2*((N-3)^2/4)+mu_bi*mu_new*(3*N-5)/2)+(mu_new^2/4))
outlier3 <- (mu_bi - mu_new)
outlier1 <- beta*(outlier1 + (1-mu_bi)) - gam
outlier2 <- beta*(outlier2 + (1-mu_bi)) - gam
outlier3 <- beta*(outlier3 + (1-mu_bi)) - gam


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

# 5. Perturbed for three patches and for all (income/outcome) directions.
N <- 400 
mu_bi <- 3
beta <- 2
gam <- 2
mu_new <- 4 # New mean

center = beta*(1-mu_bi)-gam
radius = beta*sigma_bi*sqrt(N)

mat_2 <- J_matrix(N, mu_bi, sigma_bi, mu_d,sigma_d, beta, gam)
mat_2[2:N,1] <- rnorm(N-1,mu_new,sigma_bi)
mat_2[1,2] <- rnorm(1,mu_new,sigma_bi)
mat_2[3:N,2] <- rnorm(N-2,mu_new,sigma_bi)
mat_2[1:2,3] <- rnorm(2,mu_new,sigma_bi)
mat_2[4:N,3] <- rnorm(N-3,mu_new,sigma_bi)
eig <- eigen_mat(mat_2)

outlier1 <- (mu_bi*((N-2)/2)+mu_new+sqrt(mu_bi^2*((N-4)^2/4)+mu_bi*mu_new*(3*N-5))+(mu_new^2))
outlier2 <- (mu_bi*((N-2)/2)+mu_new-sqrt(mu_bi^2*((N-4)^2/4)+mu_bi*mu_new*(3*N-5))+(mu_new^2))
outlier3 <- mu_bi - mu_new
outlier1 <- beta*(outlier1 + (1-mu_bi)) - gam
outlier2 <- beta*(outlier2 + (1-mu_bi)) - gam
outlier3 <- beta*(outlier3 + (1-mu_bi)) - gam


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


# 6. Perturbed for four patches and for all (income/outcome) directions.
N <- 200 
mu_bi <- 3
beta <- 2
gam <- 2
mu_new <- 4 # New mean

center = beta*(1-mu_bi)-gam
radius = beta*sigma_bi*sqrt(N)

mat_2 <- J_matrix(N, mu_bi, sigma_bi, mu_d,sigma_d, beta, gam)
mat_2[2:N,1] <- rnorm(N-1,mu_new,sigma_bi)
mat_2[1,2] <- rnorm(1,mu_new,sigma_bi)
mat_2[3:N,2] <- rnorm(N-2,mu_new,sigma_bi)
mat_2[1:2,3] <- rnorm(2,mu_new,sigma_bi)
mat_2[4:N,3] <- rnorm(N-3,mu_new,sigma_bi)
mat_2[1:3,4] <- rnorm(3,mu_new,sigma_bi)
mat_2[5:N,4] <- rnorm(N-4,mu_new,sigma_bi)
eig <- eigen_mat(mat_2)

outlier1 <- (mu_bi*((N-3)/2)+(3/2)*mu_new+sqrt(mu_bi^2*((N-5)^2/4)+mu_bi*mu_new*(5*N-17))+((9/4)*mu_new^2))
outlier2 <- (mu_bi*((N-3)/2)+(3/2)*mu_new-sqrt(mu_bi^2*((N-5)^2/4)+mu_bi*mu_new*(5*N-17))+((9/4)*mu_new^2))
outlier3 <- mu_bi - mu_new
outlier1 <- beta*(outlier1 + (1-mu_bi)) - gam
outlier2 <- beta*(outlier2 + (1-mu_bi)) - gam
outlier3 <- beta*(outlier3 + (1-mu_bi)) - gam

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

# 7. Change trnasmision rates.
mu_new <- 2 # New mean
alpha <- 2
mat_2 <- J_matrix(N, mu_bi, sigma_bi, mu_d,sigma_d, beta, gam)
mat_2[2:N,1] <- rnorm(N-1,(mu_bi*alpha/beta)+1,sigma_bi)
mat_2[1,1] <- rnorm(1,1+(alpha/beta),sigma_bi)
eig <- eigen_mat(mat_2)

outlier1 <- (mu_bi*(N/2)+(1/2)*(alpha/beta)+sqrt(mu_bi^2*((N^2)/4)+(N-1)*mu_bi^2*(alpha/beta)-(N/2-1)*mu_bi*(alpha/beta)+(1/4)*((alpha^2)/(beta^2))))
outlier2 <- (mu_bi*(N/2)+(1/2)*(alpha/beta)-sqrt(mu_bi^2*((N^2)/4)+(N-1)*mu_bi^2*(alpha/beta)-(N/2-1)*mu_bi*(alpha/beta)+(1/4)*((alpha^2)/(beta^2))))

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
  out <- (mu_bi*(N-1)+1)*beta -gam
  return(out)
}


out1 <- function(N){
  out <- (mu_bi*(N/2)+sqrt(mu_bi^2*((N^2/4-1))+mu_bi*mu_new))*beta-gam
  return(out)
}

out2 <- function(N){
  out <- (mu_bi*(N/2)+sqrt(mu_bi^2*((N^2-2)/4)+mu_bi*mu_new*(N-1)))*beta-gam
  return(out)
}

out3 <- function(N){
  out <-  (mu_bi*((N-1)/2)+mu_new/2+sqrt(mu_bi^2*((N^2+2*N-7)/4)-mu_bi*mu_new*(N-3)/2)+(mu_new^2/4))*beta-gam
  return(out)
}

out4 <- function(N){
  out <- (mu_bi*((N-1)/2)+mu_new/2+sqrt(mu_bi^2*((N-3)^2/4)-mu_bi*mu_new*(3*N-5)/2)+(mu_new^2/4))*beta-gam
  return(out)
}

out5 <- function(N){
  out <- (mu_bi*((N-2)/2)+mu_new+sqrt(mu_bi^2*((N-4)^2/4)-mu_bi*mu_new*(3*N-5))+(mu_new^2))*beta-gam
  return(out)
}

out6 <- function(N){
  out <- (mu_bi*((N-3)/2)+(3/2)*mu_new+sqrt(mu_bi^2*((N-5)^2/4)-mu_bi*mu_new*(5*N-17))+((9/4)*mu_new^2))*beta-gam
  return(out)
}

# for(i in c(1:N)){
#   a <- (mu_bi*(i/2)+sqrt(mu_bi^2*((i^2/4-1))+mu_bi*mu_new))*beta-gam
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
