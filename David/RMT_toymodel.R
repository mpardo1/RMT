rm(list = ls())
source("~/RMT/David/RMT_genrandom.R")
source("~/RMT/David/RMT_plotmobility.R")
source("~/RMT/David/d_functions_eigen_int.R")
library("viridis")
library("latex2exp")
library("ggforce")
####

N = 100
mat <- matrix(rnorm(N^2,0,1/10), nrow = N)
eig_m <- eigen_mat(mat)
ggplot(eig_m) +
  geom_point(aes(re, im), size = 0.5) +
  geom_ellipse(aes(x0 = 0, y0 = 0, 
                   a = 1,
                   b = 1, 
                   angle = 0), color = "blue") +
  theme_bw() + coord_fixed() + 
  theme(text = element_text(size = 20)) +
  xlab("") + 
  ylab("")

ggplot(eig_m) +
  geom_point(aes(re, im), size = 0.5) +
  theme_bw() + coord_fixed() + 
  theme(text = element_text(size = 20)) +
  xlab("") + 
  ylab("")
