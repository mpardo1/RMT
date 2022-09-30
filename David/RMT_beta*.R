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

# Eigenvalues:

stragBET <- function(muc){
  alp <- alp_bet
  a <- mub*muw +muc
  b <- alp
  c <- alp*muw
  outl <- (1/2)*(N*a + b + sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
  outl + (mub*(1-muw) - N*muc - mug)
} 
#----------------------------------------------------------------------------#
####### Parameter values #####

# N <- 100
N <- 50
# para los plots
sus_init <- rep(100000, N) # initial susceptibles
inf_init <- runif(N, min = 50,100)  # initial infecteds
end_time <- 200
end_time_rand <- 20

Deltas <- rep(0.1, N) # birth rate
mub <- 0.1
sb <- 0.001
betas <- rep(mub, N) # transmission rates
#betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.1, N) # loss of immunity rates
mud <- 0.5
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.45
alphas <- rep(mua, N) # recovery rates
mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates
gammas = deaths + alphas + deltas
mug <- gammas[1]
alp_bet <- 1
  
muw <- 0.1
sw <- 0.07/3
rhow <- 0             # original rho (Gamma of baron et al)

muc <- 0.0001
sc <- 0.00001
rhoc <- .001

mub*muw*(N-1)+mub-mua-mud-mudel

col_unstab <- "#0D0C18"
colA <- "#F69A79"
colB <- "#F26430"
colC <- "#D2430F"
colD <- "#85DCFF"
colE <- "#0ABAFF"
colF <- "#0084B8"

COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
MIGRATION <- rand_mat_ell(N, muc, sc, rhoc, distrib = "beta")
# MIGRATION <- matrix(0,N,N)
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)

################################################################
##### Change K:
stragBET <- function(k){
  (N/2)*(mub*muw+muc) + (alp_bet/2)*(1+(k-1)*muw) + mub*(1-muw) -
    mug - N*muc + (1/2)*sqrt((N*(mub*muw+muc))^2 +
                               alp_bet^2*(1+(k-1)*muw)^2 +
                 2*alp_bet*(mub*muw+muc)*(N*(1+(k-1)*muw) +
                                            2*(N-k)*(muw-1)))
} 

betas[1] <- betas[1] + alp_bet
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))
plot_eigen(jacobian) + 
  geom_point(aes(stragBET(1),0), color = "red")

####################################################################
# Patches and new mean:
k <- 1
nu <- 0.5

eig_stab <- eigen_mat(jacobianC)
plot_eigen(jacobianC) +
  geom_point(aes(stragC(muc),0), color = "red")

### Outliers:
vec <- seq(0,1,0.01)
out_BET <- sapply(vec, stragBET)

df_out <- data.frame(muc = vec, 
                     outBET = out_BET)

mig_plot <- ggplot(df_out) + 
  geom_line(aes(muc,outBET)) +
  ylab("s(J)") + xlab(TeX("$\\mu_m$")) + 
  theme_bw()

######### CHANGING the PERTURBATION #############
# Eigenvalues:
stragBET <- function(muw){
  alp <- alp_bet
  a <- mub*muw +muc
  b <- alp
  c <- alp*muw
  outl <- (1/2)*(N*a + b + sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
  outl + (mub*(1-muw) - N*muc - mug)
} 

vec <- seq(0,1,0.01)
out_BET <- sapply(vec, stragBET)

df_out <- data.frame(muw = vec, 
                     outBET = out_BET)

com_plot <- ggplot(df_out) + 
  geom_line(aes(muw,outBET)) +
  ylab("s(J)") + xlab(TeX("$\\mu_c$")) + 
  theme_bw()

######### CHANGING the PERTURBATION #############
# Eigenvalues:
stragBET <- function(mub){
  alp <- alp_bet
  a <- mub*muw +muc
  b <- alp
  c <- alp*muw
  outl <- (1/2)*(N*a + b + sqrt((N*a)^2 - (2*N-4)*a*b + (4*N-4)*a*c + b^2))
  outl + (mub*(1-muw) - N*muc - mug)
} 

vec <- seq(0,1,0.01)
out_BET <- sapply(vec, stragBET)

df_out <- data.frame(mub = vec, 
                     outBET = out_BET)

bet_plot <- ggplot(df_out) + 
  geom_line(aes(mub,outBET)) +
  ylab("s(J)") + xlab(TeX("$\\beta$")) + 
  theme_bw()

##### Change K:
stragBET <- function(k){
  (N/2)*(mub*muw+muc) + (alp_bet/2)*(1+(k-1)*muw) + mub*(1-muw) - mug - N*muc +
    (1/2)*sqrt((N*(mub*muw+muc))^2 + alp_bet^2*(1+(k-1)*muw)^2 +
                 2*alp_bet*(mub*muw+muc)*(N*(1+(k-1)*muw) +
                                            2*(N-k)*(muw-1)))
} 

vec <- seq(0,N/2,1)
out_BET <- sapply(vec, stragBET)

df_out <- data.frame(k = vec, 
                     outBET = out_BET)

k_plot <- ggplot(df_out) + 
  geom_line(aes(k,outBET)) +
  ylab("s(J)") + xlab("k") + 
  theme_bw()
### panel
col <- "#5C5D8D"
plot_grid(mig_plot + 
            scale_colour_manual(values = c(col))  ,
          com_plot + rremove("ylab") + 
            scale_colour_manual(values = c(col)) ,
          bet_plot + 
            scale_colour_manual(values = c(col)) ,
          k_plot + rremove("ylab") + 
            scale_colour_manual(values = c(col)) )
Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/SM/"
muw
mub
N
mug
muc
ggsave(file=paste0(Path,"StrategiesSJ50N0_05mum0_1muc0_1mub0_95mug.pdf"))
