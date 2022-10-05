####### RANDOM MATRICES FOR METAPOPULATION MODELS #######
### 
### article plots
###
### plots for the article
### works over "RMT_parent.R"
###
rm(list = ls())
source("~/RMT/David/RMT_genrandom.R")
source("~/RMT/David/RMT_plotmobility.R")
source("~/RMT/David/d_functions_eigen_int.R")
library("viridis")

####### COMMUTING: PANEL ################################
library("copula")
library("copula")
library("viridis")
library("ggsci")
library("grid")
library("reshape")
library("ggpubr")
library("cowplot")
library("ggforce")

a <- 1
b <- 2
c <- 3

mat_func <- function(k,N){
  mat <- matrix(a,N,N)
  mat[,1:k] <- a + c
  diag(mat) <- a
  if(k > 1){
    diag(mat[1:k,1:k]) <- a + b
  }else{
    mat[1,1] <- a + b
  }
  mat
}

maxeig_mat <- function(k,N){
  (1/2)*(N*a + b + (k-1)*c + sqrt((N*a)^2-(2*N-4*k)*a*b + 
                                               ((2*k+2)*N-(4*k))*a*c +
                                       2*(k-1)*b*c + b^2 + ((k-1)*c)^2))
}

mat_f <- mat_func(2,3)

eig <- eigen_mat(mat_f)
max(eig$re) == maxeig_mat(2,3)
