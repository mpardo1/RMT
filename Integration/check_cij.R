rm(list = ls())
library("tidyverse")
library("deSolve")
library("ggplot2")
library("gifski")
library("ggforce")
library("ggpubr")

#----------------------------------------------------------------------------#
source("~/RMT/Integration/functions_eigen_int.R")
#-------------------LOGIC----------------------
# Mobility parameter, 0: Just commuting, 1: just migration 2: migration & commuting.
MOB <- 0
# Integration parameter, 0: No integration, 1: integration.
INT <- 1
# Parameter for initial population. 0: No cte, 1: cte.
CTE_POP <- 1
# Parameter for transmission rate. 0: No cte, 1: cte.
BETA_CTE <- 0
# Parameter for initial infected ind. 0: No cte , 1: cte.
CTE_INF <- 1

#-----------------------------------------------------------------------------#

sig <- seq(0,1/4,0.01)
l <- length(sig)
df = data.frame(mean = 0, sigma = 0, cond = FALSE)
for(i in c(1:l)){
  mu_w <- seq(0,1,0.01)
  list <- sapply(mu_w, validate_mu_s, sigma = sig[i])
  df1 <- data.frame(mean = mu_w, sigma = sig[i], cond = list)
  df <- rbind(df,df1)
}
 
# # Color map TRUE posible values for mu and sigma in a beta distribution.
ggplot(df) +
  geom_point(aes(mean,sigma, colour = cond))

# Filter the values where sigma = 0 and the beta dist properties holds.
df_filt <- df %>% filter(mean != 0 & sigma != 0 & cond == TRUE)

beta_a_b_mat <- function(mat){
  beta_a_b(mat[1], mat[2])
}

mat_ab <- as.matrix(df_filt[,c(1,2)])
vec_ab <- apply(mat_ab,1, beta_a_b_mat)
# df_filt a df with the possible values for the mean and variance of the beta distribution
# with the respectives values for a and b (parameters of the beta distribution). 
df_filt$a <- vec_ab[1,]
df_filt$b <- vec_ab[2,]

# --------------------------BETA DIST WITH A & B-----------------------------#
# Parameters:
N = 200
beta_ct = 0.5
beta_ct <- rep(beta_ct,N)
gamma_ct = 1.4
gamma_ct <- rep(gamma_ct,N)
MOB <- 2
mu_w <- 0.1 
s_w <- 0.01
vec<-  beta_a_b(mu_w,s_w)
a_w <- vec[1]
b_w <- vec[2]
mu_c <- 0.01 
s_c <- 0.000001
vec<-  beta_a_b(mu_c,s_c)
a_c <- vec[1]
b_c <- vec[2]

df_filt$s <- sqrt(df_filt$sigma)

mu_c <- a_c/(a_c + b_c)
s_c <-  sqrt((a_c*b_c)/(((a_c + b_c)^2)*(1+a_c+b_c)))

mu_w <- a_w/(a_w + b_w)
s_w <-  sqrt((a_w*b_w)/(((a_w + b_w)^2)*(1+a_w+b_w)))

DIST <-  "beta"
mig_mat <- mat_conect(N,a_c,b_c,DIST)
com_mat <- mat_conect(N,a_w,b_w,DIST)

jac <- jacobian(N,beta_ct,gamma_ct, com_mat, mig_mat,0, MOB)
eig <- eigen_mat(jac)

# the 2 is to use the model with commuting and migration:
outl <- pred_outlier(N, beta_ct, gamma_ct, mu_w, MOB)
# Compute the predicted radius and center
rad <- pred_radius(N, beta_ct[1], gamma_ct[1], 0, mu_c, s_c, mu_w, s_w, MOB)
cent <- pred_center(N, beta_ct[1], gamma_ct[1], 0, mu_c, s_c, mu_w, s_w, MOB) 

plot_eigen(eig, cent, rad, outl,MOB) +
  ggtitle(paste0("Parameters:\n \n", 
                 "N: ",N,",", " gamma: ", gamma_ct[1],",", " beta: ", beta_ct[1], "\n",
                 "Migration:    ", expression(mu_c), ": ", round(mu_c, 4),","," ",
                 expression(Sigma_c), ": ",round(s_c, 4),"\n",
                 "Commuting:    ",expression(mu_w), ": ", round(mu_w, 2)," ",
                 expression(sigma_w), ": ", round(s_w, 2)))

plot_eig <- ggplot(eig) + geom_point(aes(re,im), size = 0.05) 

colSums(mig_mat)
mean(colSums(mig_mat))
var(colSums(mig_mat))
(N-1)*mu_c


ind <-  which(eig$re == max(eig$re))
eig_filt <-  eig[-ind,]
plot_eig <- ggplot(eig_filt) + geom_point(aes(re,im), size = 0.05)  + 
  geom_circle(aes(x0 = cent,
                  y0 = 0,
                  r = rad), colour = "blue",
              show.legend = NA,size = 0.2) +
  coord_fixed() +
  ggtitle(paste0("Parameters:\n \n", 
                 "N: ",N,",", " gamma: ", gamma_ct[1],",", " beta: ", beta_ct[1], "\n",
                 "Migration:    ", expression(mu_c), ": ", round(mu_c, 4),","," ",
                 expression(Sigma_c), ": ",round(s_c, 4),"\n",
                 "Commuting:    ",expression(mu_w), ": ", round(mu_w, 2)," ",
                 expression(sigma_w), ": ", round(s_w, 2)))
plot_eig

plot_no_circ <- plot_no_circ + labs(title="a")
plot_circ <- plot_circ + labs(title="b")
ggarrange(plot_no_circ,plot_circ, ncol = 1, nrow = 2)

gamma_ct_w <- format(round(gamma_ct,2), decimal.mark = ',')
beta_ct_w <- format(round(beta_ct,2), decimal.mark = ',')
mu_w_w <- format(round(mu_w,2), decimal.mark = ',')
s_w_w <- format(round(s_w,2), decimal.mark = ',')
mu_c_w <- format(round(mu_c,2), decimal.mark = ',')
s_c_w <- format(round(s_c,2), decimal.mark = ',')
Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/cij_issue/"
path <- paste0(Path,"N",N,"g",gamma_ct_w,"b",beta_ct_w,"mw",
               mu_w_w,"sw",s_w_w,"mc",mu_c_w,"sc","0.0316","sc_b","0,001",".png")
# path <- paste0(Path,"mod1patchtrans","N",N,"g",gamma_ct,"b",
# beta_ct,"mw","bnew","14",
# mu_w,"sw",s_w,"mm",mu_m,"sm",s_m,".png")
png(file = path, width = 8000, height = 6000, res = 1100)
ggarrange(plot_no_circ,plot_circ, ncol = 1, nrow = 2)
dev.off()
