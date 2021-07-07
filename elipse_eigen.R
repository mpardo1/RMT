rm(list = ls())
library(easypackages)
libraries("gdata", "ggExtra","ggplot2", "numbers","tidyverse",
          "data.table","multiplex","reshape","viridis","stats",
          "ggpubr","ggstatsplot","e1071","mlr3misc","deSolve",
          "gganimate", "matlib")

elipse <- function(x){
  y <- (1-rho)*sqrt((1-(x-((1/gamma)*(mu_d-mu)))^2/(1+rho)))
  return(y)
}

min_seq <- 1
max_seq <- 5
x <- seq(min_seq, max_seq, 0.1)
rho = 0.5
sigma = 0.4
N = 100
mu_d = 3
mu = 2
gamma = 0.4
vec <- unlist(lapply(x,elipse))

df <- data.frame( x , vec)
ggplot(df) +
  geom_line(aes(x,vec))

