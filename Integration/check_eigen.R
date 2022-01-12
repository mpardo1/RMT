rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")
library("ggpubr")
library("ggsci")
library("ggplot2")

# Parameters: 
dim <- 43
bet <- 1
mu_w <- 2
mu_c <- 2.4
alp <- 2
mat <- matrix(0, nrow = dim, ncol = dim)

a <- bet*mu_w +mu_c
b <- alp
c <- alp*mu_w

mat[,] <- a
mat[1,1] <- a+b
mat[2:dim,1] <- a+c

eigen_m <- as.complex(eigen(mat, only.values = TRUE)$values)
df <- data.frame(re = Re(eigen_m), im = Im(eigen_m))
plot_eig <- ggplot(df) + geom_point(aes(re,im), size = 0.05) 

bet_cte <- bet
mu_m <- mu_c
N <- dim
bet_new <- bet_cte + alp


outl <- (1/2)*(N*a + b + sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
outl2 <- (1/2)*(N*a + b - sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))

plot_eig + 
  geom_point(aes(outl,0), colour =  "blue",
             show.legend = NA) +
  geom_point(aes(outl2,0), colour =  "blue",
           show.legend = NA) +
  geom_point(aes(0,0), colour =  "blue",
             show.legend = NA)


