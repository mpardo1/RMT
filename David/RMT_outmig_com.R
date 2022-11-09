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
  c <- (mub*muw+muc)
  (N/2)*c + ((k-1)/2)*(mub*nu) + mub*(1-muw) - mug - N*muc + 
    (1/2)*sqrt(N^2*c^2 + ((k-1)*mub*nu)^2 + 2*(N+N*k-4*k)*c*mub*nu)
  
} 

stragC <- function(muc){
  c <- (mub*muw+muc)
  (N/2)*c + ((k-1)/2)*(mub*nu) + mub*(1-muw) - mug - N*muc + 
    (1/2)*sqrt((N*c)^2 + 2*mub*nu*c*(N+k*(3*N-2*k-2)) +
                 ((mub*nu)^2)*(4*k*(N-k) + (k-1)^2))
} 


stragDEF <- function(muc){
  k0 <- k/N
  (N-1)*mub*(muw + k0*nu) + mub - mug
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
# Patches and new mean:
k <- 1
nu <- 0.5
### Compute the eigenvalues:
fs <- sample(1:N,2)
COMMUTINGA <- COMMUTINGC <- COMMUTING 
# STRATEGY A:
COMMUTINGA[fs,] <- COMMUTINGA[fs,] + nu
jacobianA <- (COMMUTINGA + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

eig_stab <- eigen_mat(jacobianA)
plot_eigen(jacobianA) +
  geom_point(aes(stragAB(muc),0), color = "red")

# STRATEGY C:
COMMUTINGC[fs[1],] <- COMMUTINGC[fs[1],] + nu
COMMUTINGC[,fs[1]] <- COMMUTINGC[,fs[1]] + nu
jacobianC <- (COMMUTINGC + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

eig_stab <- eigen_mat(jacobianC)
plot_eigen(jacobianC) +
  geom_point(aes(stragC(muc),0), color = "red")

### Outliers:
vec <- seq(0,1,0.01)
out_AB <- sapply(vec, stragAB)
out_C <- sapply(vec, stragC)
out_DEF <- sapply(vec, stragDEF)

df_out <- data.frame(muc = vec, 
                     AB = out_AB,
                     C = out_C,
                     DEF = out_DEF)

ggplot(df_out) + geom_line(aes(muc,C))
df_plot <- reshape2::melt(df_out, id.vars = c("muc"))
library("latex2exp")
text_size = 16
colnames(df_plot) <- c("muc", "Strategy","value")
mig_plot <- ggplot(df_plot) + 
  geom_line(aes(muc,value, 
                colour = Strategy ), size = 1) + 
  ylab("s(J)") + xlab(TeX("$\\mu_m$")) + 
  scale_colour_manual(values = c(colA,colC,colD)) +
  theme_bw() + theme(text = element_text(size = text_size))

Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/SM/"
muw
mub
N
mug
ggsave(file=paste0(Path,"MUCvsSJ40N0_12muc0_1mub0_55mug.pdf"))

######### CHANGING the PERTURBATION #############
# Eigenvalues:
# Changing nu:
stragAB <- function(nu){
  c <- (mub*muw+muc)
  (N/2)*c + ((k-1)/2)*(mub*nu) + mub*(1-muw) - mug - N*muc + 
    (1/2)*sqrt(N^2*c^2 + ((k-1)*mub*nu)^2 + 2*(N+N*k-4*k)*c*mub*nu)
  
} 

stragC <- function(nu){
  c <- (mub*muw+muc)
  (N/2)*c + ((k-1)/2)*(mub*nu) + mub*(1-muw) - mug - N*muc + 
    (1/2)*sqrt((N*c)^2 + 2*mub*nu*c*(N+k*(3*N-2*k-2)) +
                 ((mub*nu)^2)*(4*k*(N-k) + (k-1)^2))
} 

stragDEF <- function(nu){
  k0 <- k/N
  (N-1)*mub*(muw + k0*nu) + mub - mug
} 

vec <- seq(-muw,1,0.01)
# k <- 2
# muc <- 0.05
out_AB <- sapply(vec, stragAB)
out_C <- sapply(vec, stragC)
out_DEF <- sapply(vec, stragDEF)

df_out <- data.frame(nu = vec, 
                     AB = out_AB,
                     C = out_C,
                     DEF = out_DEF)

ggplot(df_out) +
  geom_line(aes(nu,C))
df_plot <- reshape2::melt(df_out, id.vars = c("nu"))
library("latex2exp")
colnames(df_plot) <- c("nu", "Strategy","value")
nu_plot <- ggplot(df_plot) + 
  geom_line(aes(nu,value, 
                colour = Strategy ), size = 1) + 
  scale_colour_manual(values = c(colA,colC,colD)) +
  ylab("s(J)") + xlab(TeX("$\\nu$")) + 
  theme_bw() + theme(text = element_text(size = text_size))

nu <- 0.5
# Changing k:
stragAB <- function(k){
  c <- (mub*muw+muc)
  (N/2)*c + ((k-1)/2)*(mub*nu) + mub*(1-muw) - mug - N*muc + 
    (1/2)*sqrt(N^2*c^2 + ((k-1)*mub*nu)^2 + 2*(N+N*k-4*k)*c*mub*nu)
  
} 

stragC <- function(k){
  c <- (mub*muw+muc)
  (N/2)*c + ((k-1)/2)*(mub*nu) + mub*(1-muw) - mug - N*muc + 
    (1/2)*sqrt((N*c)^2 + 2*mub*nu*c*(N+k*(3*N-2*k-2)) +
                 ((mub*nu)^2)*(4*k*(N-k) + (k-1)^2))
}

stragDEF <- function(k){
  k0 <- k/N
  (N-1)*mub*(muw + k0*nu) + mub - mug
} 

vec <- seq(0,N/2,1)
# nu <- 0.5
out_AB <- sapply(vec, stragAB)
out_C <- sapply(vec, stragC)
out_DEF <- sapply(vec, stragDEF)

df_out <- data.frame(k = vec, 
                     AB = out_AB,
                     C = out_C,
                     DEF = out_DEF)

ggplot(df_out) +
  geom_line(aes(k,C))
df_plot <- reshape2::melt(df_out, id.vars = c("k"))
library("latex2exp")
colnames(df_plot) <- c("k", "Strategy","value")
k_plot <- ggplot(df_plot) + 
  geom_line(aes(k,value, 
                colour = Strategy ), size = 1) + 
  scale_colour_manual(values = c(colA,colC,colD)) +
  ylab("s(J)") + xlab("k") + 
  theme_bw() + theme(text = element_text(size = text_size))

k <- 1
### Changing muw:
stragAB <- function(muw){
  c <- (mub*muw+muc)
  (N/2)*c + ((k-1)/2)*(mub*nu) + mub*(1-muw) - mug - N*muc + 
    (1/2)*sqrt(N^2*c^2 + ((k-1)*mub*nu)^2 + 2*(N+N*k-4*k)*c*mub*nu)
  
} 

stragC <- function(muw){
  c <- (mub*muw+muc)
  (N/2)*c + ((k-1)/2)*(mub*nu) + mub*(1-muw) - mug - N*muc + 
    (1/2)*sqrt((N*c)^2 + 2*mub*nu*c*(N+k*(3*N-2*k-2)) +
                 ((mub*nu)^2)*(4*k*(N-k) + (k-1)^2))
} 

stragDEF <- function(muw){
  k0 <-  k/N
  (N-1)*mub*(muw + k0*nu) + mub - mug
} 

vec <- seq(0,1,0.01)
# nu <- 0.5
# muc <- 0.05
out_AB <- sapply(vec, stragAB)
out_C <- sapply(vec, stragC)
out_DEF <- sapply(vec, stragDEF)

df_out <- data.frame(muw = vec, 
                     AB = out_AB,
                     C = out_C,
                     DEF = out_DEF)

ggplot(df_out) +
  geom_line(aes(muw,C))
df_plot <- reshape2::melt(df_out, id.vars = c("muw"))
library("latex2exp")
colnames(df_plot) <- c("muw", "Strategy","value")
muw_plot <- ggplot(df_plot) + 
  geom_line(aes(muw,value, 
                colour = Strategy ), size = 1) + 
  ylab("s(J)") + xlab(TeX("$\\mu_c$")) + 
  scale_colour_manual(values = c(colA,colC,colD)) +
  theme_bw() + theme(text = element_text(size = text_size))
leg <- get_legend(ggplot(df_plot) + 
                    geom_line(aes(muw,value, 
                                  colour = Strategy ), size = 1) + 
                    ylab("s(J)") + xlab(TeX("$\\mu_c$")) + 
                    scale_colour_manual(values = c(colA,colC,colD)) +
                    theme_bw() + theme(text = element_text(size = 15)))

muw <- 0.1
### panel
ggarrange(nu_plot ,
          mig_plot + rremove("ylab"),
          k_plot,
          common.legend = TRUE)

### Changing beta:
stragAB <- function(mub){
  c <- (mub*muw+muc)
  (N/2)*c + ((k-1)/2)*(mub*nu) + mub*(1-muw) - mug - N*muc + 
    (1/2)*sqrt(N^2*c^2 + ((k-1)*mub*nu)^2 + 2*(N+N*k-4*k)*c*mub*nu)
  
} 

stragC <- function(mub){
  c <- (mub*muw+muc)
  (N/2)*c + ((k-1)/2)*(mub*nu) + mub*(1-muw) - mug - N*muc + 
    (1/2)*sqrt((N*c)^2 + 2*mub*nu*c*(N+k*(3*N-2*k-2)) +
                 ((mub*nu)^2)*(4*k*(N-k) + (k-1)^2))
} 

stragDEF <- function(mub){
  k0 <- k/N
  (N-1)*mub*(muw + k0*nu) + mub - mug
}

vec <- seq(0,1,0.01)
# nu <- 0.5
# muc <- 0.05
out_AB <- sapply(vec, stragAB)
out_C <- sapply(vec, stragC)
out_DEF <- sapply(vec, stragDEF)

df_out <- data.frame(mub = vec, 
                     AB = out_AB,
                     C = out_C,
                     DEF = out_DEF)

ggplot(df_out) +
  geom_line(aes(mub,C))
df_plot <- reshape2::melt(df_out, id.vars = c("mub"))
library("latex2exp")
colnames(df_plot) <- c("mub", "Strategy","value")
beta_plot <- ggplot(df_plot) + 
  geom_line(aes(mub,value, 
                colour = Strategy ), size = 1) + 
  ylab("s(J)") + xlab(TeX("$\\beta$")) + 
  scale_colour_manual(values = c(colA,colC,colD)) +
  theme_bw() + theme(text = element_text(size = text_size))

mub <- 0.1
### panel
library("cowplot")

text_size = 16
plot_sJ <- plot_grid(nu_plot + geom_hline(yintercept = 0, color = "red",
                                          linetype = "dashed") +
                       theme(legend.position = "none", 
                             text = element_text(size = text_size)),
          mig_plot  +  geom_hline(yintercept = 0, color = "red", 
                                  linetype = "dashed") +
            theme(legend.position = "none",
                  text = element_text(size = text_size)) + rremove("ylab"),
          muw_plot + geom_hline(yintercept = 0, color = "red", 
                                linetype = "dashed") +
            theme(legend.position = "none",
                  text = element_text(size = text_size)),
          beta_plot + geom_hline(yintercept = 0, color = "red", 
                                 linetype = "dashed") +
            theme(legend.position = "none", 
                  text = element_text(size = text_size)) + 
            rremove("ylab"),
          ncol = 2, nrow = 2, align = "v")

plot_SJ <- ggdraw() +
  draw_plot(plot_sJ) +
  draw_plot(leg, x = 0.8, y = .75, width = .2, height = .2)

plot_SJ

# plot_grid( leg, plot_sJ, ncol = 1, rel_heighs = c(0.1,2.2))
Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/SM/"
muw
mub
N
mug
muc
ggsave(file=paste0(Path,"Strategies1SJ50N0_05mum0_1muc0_1mub0_95mug.pdf"))

