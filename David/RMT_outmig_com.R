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

stragAB <- function(muc){
  (N/2)*(mub*muw+muc) + ((k-1)/2)*(mub*nu) + mub*(1-muw) - mug - N*muc + 
    (1/2)*sqrt(N^2*(mub*muw + muc)^2 + ((k-1)*mub*nu)^2 + 2*(N+N*k-4*k)*(mub*muw + muc)*mub*nu)
  
} 

stragC <- function(muc){
  N*(mub*muw+muc)/2 + ((k-1)/2)*mub*nu +
    sqrt((N*(mub*muw+muc))^2+
           2*mub*nu*(mub*muw+muc)*(N+k*(3*N-(2*k)-2)) +
           (mub*nu)^2*(4*k*(N-k) + (k-1)^2))/2 +
    mub*(1-muw) - mug - N*muc
} 

stragDEF <- function(k, muwstar){
  mub*(muw*(N-k)+muwstar*k/N)*(N-1)/N + mub - mug
} 

#----------------------------------------------------------------------------#

####### Parameter values #####

# N <- 100
N <- 40
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
mud <- 0.1
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.45
alphas <- rep(mua, N) # recovery rates
mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates
gammas = deaths + alphas + deltas
mug <- gammas[1]
muw <- 0.12
sw <- 0.07/3
rhow <- 0             # original rho (Gamma of baron et al)

muc <- 0.0001
sc <- 0.00001
rhoc <- .001

mub*muw*(N-1)+mub-mua-mud-mudel

col_unstab = "#0D0C18"   # negro
col_stabA = "#06D622"  # verde
col_stabB = "#E96F1D"  # naranja
col_stabC = "#F8DC22"  # amarillo
col_stabD = "#B62F2F"  # rojo
col_stabE = "#2C5CE1"  # azul
col_stabF = "#AA2FB5"  # violeta
col_circ = "#010C0C" # negro

COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
MIGRATION <- rand_mat_ell(N, muc, sc, rhoc, distrib = "beta")
# MIGRATION <- matrix(0,N,N)
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)

# Unstable scenario
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

vec <- seq(0,1,0.01)
k <- 2
nu <- 0.5
out_AB <- sapply(vec, stragAB)
out_C <- sapply(vec, stragC)

df_out <- data.frame(muc = vec, 
                     AB = out_AB,
                     "C" = out_C)

ggplot(df_out) +
  geom_line(aes(muc,outC))
df_plot <- reshape2::melt(df_out, id.vars = c("muc"))
library("latex2exp")
colnames(df_plot) <- c("muc", "Scenario","value")
ggplot(df_plot) + 
  geom_line(aes(muc,value, 
                colour = Scenario,
                linetype = Scenario ), size = 1) + 
  ylab("s(J)") + xlab(TeX("$\\mu_m$")) + 
  theme_bw()

Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/SM/"
muw
mub
N
mug
ggsave(file=paste0(Path,"MUCvsSJ40N0_12muc0_1mub0_55mug.pdf"))
