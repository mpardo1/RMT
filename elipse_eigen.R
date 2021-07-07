rm(list = ls())
library(easypackages)
libraries("gdata", "ggExtra","ggplot2", "numbers","tidyverse",
          "data.table","multiplex","reshape","viridis","stats",
          "ggpubr","ggstatsplot","e1071","mlr3misc","deSolve",
          "gganimate", "matlib")

elipse <- function(t){
  x <- (mu_d - mu)/gamma + (sigma*sqrt(N)*(1+rho)/gamma)*cos(t)
  y <- (sigma*sqrt(N)*(1-rho)/gamma)*sin(t)
  vec <- c(x,y)
  return(vec)
}


rho = 0.5
sigma = 0.4
N = 100
mu_d = 3
mu = 2
gamma = 0.4
min_seq <- (sigma*sqrt(N)*(1+rho) - (mu_d-mu))*gamma^(-1)
max_seq <- (sigma*sqrt(N)*(1+rho) + (mu_d-mu))*gamma^(-1)
min_seq <- 0
max_seq <- 5
x <- seq(min_seq, max_seq, 0.1)
vec <- unlist(lapply(x,elipse))

df <- data.frame( x , vec)
df_1 <- data.frame(x,vec = -vec)
df_join <- rbind(df,df_1)
ggplot(df_join) +
  geom_point(aes(x,vec))

